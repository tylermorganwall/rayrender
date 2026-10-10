#include "wavefront.h"
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include "sobol_directions.h"
namespace wavefront_tables {
#include "samplerBlueNoise.h"
}
#include <chrono>
#include <cstring>
#include <stdexcept>
#include <thread>
#include <array>

namespace {
const char *source =
#include "wavefront_kernels.inc"
    ;

std::runtime_error MetalError(const char *context, NSError *error) {
  return std::runtime_error(
      std::string(context) + ": " +
      (error ? [[error localizedDescription] UTF8String] : "Metal resource allocation failed"));
}

// Pipeline and immutable sampler tables are shared between renders. Per-scene
// buffers, acceleration structures and queues remain owned by one session.
struct Pipelines {
  id<MTLDevice> device;
  id<MTLComputePipelineState> primary, intersect, shade, shadow, dispatch;
  id<MTLBuffer> directions, blue;
  Pipelines() {
    device = MTLCreateSystemDefaultDevice();
    if (!device || !device.supportsRaytracing)
      throw std::runtime_error("this Metal device does not support ray tracing");
    NSError *error = nil;
    MTLCompileOptions *options = [MTLCompileOptions new];
#if __MAC_OS_X_VERSION_MAX_ALLOWED >= 150000
    if (@available(macOS 15.0, *))
      options.mathMode = MTLMathModeSafe;
    else
#endif
    {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
      options.fastMathEnabled = NO;
#pragma clang diagnostic pop
    }
    options.languageVersion = MTLLanguageVersion2_3;
    id<MTLLibrary> library = [device newLibraryWithSource:[NSString stringWithUTF8String:source]
                                                  options:options
                                                    error:&error];
    if (!library)
      throw MetalError("wavefront shader compilation", error);
    auto pipeline = [&](NSString *name) {
      id<MTLFunction> function = [library newFunctionWithName:name];
      id<MTLComputePipelineState> result = [device newComputePipelineStateWithFunction:function
                                                                                 error:&error];
      if (!result)
        throw MetalError("wavefront pipeline compilation", error);
      return result;
    };
    primary = pipeline(@"primary");
    intersect = pipeline(@"intersect_paths");
    shade = pipeline(@"shade");
    shadow = pipeline(@"shadows");
    dispatch = pipeline(@"dispatch_size");
    directions = [device newBufferWithBytes:SPACEFILLR_SOBOL_DIRECTIONS
                                     length:21200 * 32 * sizeof(uint32_t)
                                    options:MTLResourceStorageModeShared];
    std::vector<uint32_t> values;
    values.insert(values.end(), std::begin(wavefront_tables::sobol_256spp_256d),
                  std::end(wavefront_tables::sobol_256spp_256d));
    values.insert(values.end(), std::begin(wavefront_tables::rankingTile),
                  std::end(wavefront_tables::rankingTile));
    values.insert(values.end(), std::begin(wavefront_tables::scramblingTile),
                  std::end(wavefront_tables::scramblingTile));
    blue = [device newBufferWithBytes:values.data()
                               length:values.size() * sizeof(uint32_t)
                              options:MTLResourceStorageModeShared];
    if (!directions || !blue)
      throw MetalError("wavefront sampler tables", nil);
  }
};

class MetalWavefront final : public WavefrontSession {
  Pipelines &kernels;
  id<MTLCommandQueue> queue;
  id<MTLAccelerationStructure> acceleration;
  // Slots match SCENE_ARGUMENTS. Parameters use setBytes at 12; the acceleration
  // structure uses setAccelerationStructure at 16, so neither occupies a buffer.
  std::array<id<MTLBuffer>, 16> buffers;
  id<MTLBuffer> indirect;
  size_t capacity;

  id<MTLBuffer> Buffer(size_t size, const void *data = nullptr) {
    size_t allocation_size = std::max(size, size_t(16));
    id<MTLBuffer> result = [kernels.device newBufferWithLength:allocation_size
                                                       options:MTLResourceStorageModeShared];
    if (!result)
      throw MetalError("wavefront buffer allocation", nil);
    if (data)
      std::memcpy(result.contents, data, size);
    return result;
  }
  template <class T> id<MTLBuffer> Upload(const std::vector<T> &data) {
    if (data.empty())
      return Buffer(16);
    return Buffer(data.size() * sizeof(T), data.data());
  }

  bool Wait(id<MTLCommandBuffer> command, const std::function<bool()> &cancelled) {
    bool stop = false;
    std::exception_ptr interrupted;
    [command commit];
    while (command.status < MTLCommandBufferStatusCompleted) {
      // Keep resources alive even after cancellation until this bounded GPU
      // batch completes. Never run R APIs or preview callbacks on a Metal thread.
      if (!interrupted) {
        try {
          stop = cancelled() || stop;
        } catch (...) {
          interrupted = std::current_exception();
        }
      }
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    if (interrupted)
      std::rethrow_exception(interrupted);
    if (command.status == MTLCommandBufferStatusError)
      throw MetalError("wavefront command", command.error);
    return !(stop || cancelled());
  }

  void DispatchSize(id<MTLCommandBuffer> command, uint32_t stage) {
    auto encoder = [command computeCommandEncoder];
    encoder.label = @"Wavefront queue counts";
    [encoder setComputePipelineState:kernels.dispatch];
    [encoder setBuffer:buffers[11] offset:0 atIndex:0];
    [encoder setBuffer:indirect offset:0 atIndex:1];
    [encoder setBytes:&stage length:sizeof(stage) atIndex:2];
    [encoder dispatchThreadgroups:MTLSizeMake(1, 1, 1) threadsPerThreadgroup:MTLSizeMake(1, 1, 1)];
    [encoder endEncoding];
  }

  void Dispatch(id<MTLCommandBuffer> command, id<MTLComputePipelineState> pipeline,
                const WFParameters &parameters, bool primary = false) {
    auto encoder = [command computeCommandEncoder];
    encoder.label = primary ? @"Wavefront primary rays"
                    : pipeline == kernels.intersect ? @"Wavefront intersections"
                    : pipeline == kernels.shade ? @"Wavefront surface shading"
                                                : @"Wavefront shadow rays";
    [encoder setComputePipelineState:pipeline];
    for (size_t i = 0; i < buffers.size(); ++i)
      if (buffers[i])
        [encoder setBuffer:buffers[i] offset:0 atIndex:i];
    [encoder setBytes:&parameters length:sizeof(parameters) atIndex:12];
    [encoder setAccelerationStructure:acceleration atBufferIndex:16];
    if (primary)
      [encoder dispatchThreadgroups:MTLSizeMake(std::max(1u, (parameters.active + 63) / 64), 1, 1)
              threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
    else
      [encoder dispatchThreadgroupsWithIndirectBuffer:indirect
                                 indirectBufferOffset:0
                                threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
    [encoder endEncoding];
  }

public:
  MetalWavefront(Pipelines &kernels, const WavefrontScene &scene, size_t capacity)
      : kernels(kernels), capacity(capacity) {
    queue = [kernels.device newCommandQueue];
    if (!queue)
      throw MetalError("wavefront command queue", nil);
    buffers[0] = Upload(scene.triangles);
    buffers[1] = Upload(scene.materials);
    buffers[2] = Upload(scene.textures);
    buffers[3] = Upload(scene.texels);
    buffers[4] = Upload(scene.lights);
    buffers[5] = Buffer(capacity * 80);
    buffers[6] = Buffer(capacity * sizeof(WFPixel));
    buffers[7] = Buffer(capacity * 4);
    buffers[8] = Buffer(capacity * 4);
    buffers[9] = Buffer(capacity * 16);
    buffers[10] = Buffer(capacity * 64);
    buffers[11] = Buffer(5 * 4);
    buffers[13] = Buffer(capacity * 4);
    buffers[14] = kernels.directions;
    buffers[15] = kernels.blue;
    indirect = Buffer(16);
    // Reuse packed triangle vertices rather than allocating a second positions
    // array. Ten float4 records per triangle; the first three are its vertices.
    if (scene.triangles.size() > UINT32_MAX / 10)
      throw std::runtime_error("Metal vertex index capacity exceeded");
    std::vector<uint32_t> indices(scene.triangles.size() * 3);
    for (size_t i = 0; i < scene.triangles.size(); ++i)
      for (size_t k = 0; k < 3; ++k)
        indices[3 * i + k] = 10 * i + k;
    id<MTLBuffer> index_buffer = Upload(indices);
    auto geometry = [MTLAccelerationStructureTriangleGeometryDescriptor descriptor];
    geometry.vertexBuffer = buffers[0];
    geometry.vertexStride = sizeof(WFVector);
    geometry.indexBuffer = index_buffer;
    geometry.indexType = MTLIndexTypeUInt32;
    geometry.triangleCount = scene.triangles.size();
    geometry.opaque = YES;
    auto descriptor = [MTLPrimitiveAccelerationStructureDescriptor descriptor];
    descriptor.geometryDescriptors = @[ geometry ];
    MTLAccelerationStructureSizes sizes =
        [kernels.device accelerationStructureSizesWithDescriptor:descriptor];
    acceleration =
        [kernels.device newAccelerationStructureWithSize:sizes.accelerationStructureSize];
    id<MTLBuffer> scratch = Buffer(sizes.buildScratchBufferSize);
    if (!acceleration)
      throw MetalError("wavefront acceleration structure", nil);
    auto command = [queue commandBuffer];
    auto encoder = [command accelerationStructureCommandEncoder];
    [encoder buildAccelerationStructure:acceleration
                             descriptor:descriptor
                          scratchBuffer:scratch
                    scratchBufferOffset:0];
    [encoder endEncoding];
    Wait(command, [] { return false; });
  }

  bool Render(const WFParameters &input, const std::vector<uint32_t> &active,
              const std::vector<WFLight> &lights, const std::function<bool()> &cancelled,
              const WFPixel *&pixels) override {
    @autoreleasepool {
      if (active.size() > capacity || input.width * size_t(input.height) > capacity)
        throw std::runtime_error("wavefront image exceeds queue capacity");
      if (cancelled())
        return false;
      if (!active.empty())
        std::memcpy(buffers[13].contents, active.data(), active.size() * sizeof(uint32_t));
      if (!lights.empty())
        std::memcpy(buffers[4].contents, lights.data(), lights.size() * sizeof(WFLight));
      WFParameters parameters = input;
      // Four bounces per command buffer bound cancellation latency. Queue
      // transitions and indirect dispatches are entirely on the GPU.
      for (uint32_t first = 0; first <= input.max_depth; first += 4) {
        auto command = [queue commandBuffer];
        command.label = @"Wavefront transport batch";
        if (first == 0)
          Dispatch(command, kernels.primary, parameters, true);
        for (uint32_t depth = first; depth < std::min(first + 4, input.max_depth + 1); ++depth) {
          parameters.depth = depth;
          DispatchSize(command, 0);
          Dispatch(command, kernels.intersect, parameters);
          DispatchSize(command, 1);
          Dispatch(command, kernels.shade, parameters);
          DispatchSize(command, 2);
          Dispatch(command, kernels.shadow, parameters);
          DispatchSize(command, 3);
          std::swap(buffers[7], buffers[8]);
        }
        if (!Wait(command, cancelled))
          return false;
      }
      auto *counts = static_cast<uint32_t *>(buffers[11].contents);
      if (counts[4])
        throw std::runtime_error("Metal wavefront queue overflow; sample discarded");
      pixels = static_cast<const WFPixel *>(buffers[6].contents);
      return true;
    }
  }
};
} // namespace

std::unique_ptr<WavefrontSession> MakeMetalWavefront(const WavefrontScene &scene, size_t capacity) {
  @autoreleasepool {
    if (@available(macOS 11.0, *)) {
      static Pipelines pipelines;
      return std::make_unique<MetalWavefront>(pipelines, scene, capacity);
    }
    throw std::runtime_error("Metal wavefront requires macOS 11 or later");
  }
}
