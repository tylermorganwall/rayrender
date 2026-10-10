#include "sobol_directions.h"
#include "wavefront.h"
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
namespace wavefront_tables {
#include "samplerBlueNoise.h"
}
#include <array>
#include <chrono>
#include <cstring>
#include <mutex>
#include <stdexcept>
#include <thread>

namespace {
const char *source = "#include <metal_stdlib>\nusing namespace metal;\nnamespace wf_openpbr {\n"
#include "wavefront_openpbr_generated.inc"
                     "\n}\n"
#include "wavefront_kernels.inc"
#ifdef NOT_CRAN
#include "wavefront_probes.inc"
#endif
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
  id<MTLLibrary> library;
  id<MTLComputePipelineState> primary, intersect, shade, shade_materials, shade_pbr,
      shade_pbr_lights, shade_diffusion, shadow, dispatch;
  id<MTLBuffer> directions, blue;
  bool boundary_support;
  std::mutex material_mutex;

  id<MTLComputePipelineState> Compile(NSString *name) {
    NSError *error = nil;
    auto constants = [MTLFunctionConstantValues new];
    [constants setConstantValue:&boundary_support type:MTLDataTypeBool atIndex:0];
    auto function = [library newFunctionWithName:name constantValues:constants error:&error];
    if (!function)
      throw MetalError("wavefront function specialization", error);
    auto descriptor = [MTLComputePipelineDescriptor new];
    descriptor.computeFunction = function;
    if (boundary_support) {
      auto linked = [MTLLinkedFunctions linkedFunctions];
      linked.functions = @[ [library newFunctionWithName:@"containment_face"] ];
      descriptor.linkedFunctions = linked;
    }
    auto result = [device newComputePipelineStateWithDescriptor:descriptor
                                                        options:MTLPipelineOptionNone
                                                     reflection:nil
                                                          error:&error];
    if (!result)
      throw MetalError("wavefront pipeline compilation", error);
    return result;
  }
  void PrepareMaterials(bool basic, bool pbr, bool emission, bool diffusion) {
    std::lock_guard<std::mutex> lock(material_mutex);
    if (basic && !shade_materials)
      shade_materials = Compile(@"shade_materials");
    if (pbr && !shade_pbr)
      shade_pbr = Compile(@"shade_pbr");
    if (emission && !shade_pbr_lights)
      shade_pbr_lights = Compile(@"shade_pbr_lights");
    if (diffusion && !shade_diffusion)
      shade_diffusion = Compile(@"shade_diffusion");
  }
  explicit Pipelines(bool boundary_support) : boundary_support(boundary_support) {
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
    library = [device newLibraryWithSource:[NSString stringWithUTF8String:source]
                                   options:options
                                     error:&error];
    if (!library)
      throw MetalError("wavefront shader compilation", error);
    primary = Compile(@"primary");
    intersect = Compile(@"intersect_paths");
    shade = Compile(@"shade");
    shadow = Compile(@"shadows");
    dispatch = Compile(@"advance_queues");
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

Pipelines &SharedPipelines(bool boundaries = false) {
  if (boundaries) {
    static Pipelines pipelines(true);
    return pipelines;
  }
  static Pipelines pipelines(false);
  return pipelines;
}

class MetalWavefront final : public WavefrontSession {
  Pipelines &kernels;
  id<MTLCommandQueue> queue;
  id<MTLAccelerationStructure> acceleration;
  // Slots match SCENE_ARGUMENTS. Parameters use setBytes at 12; the acceleration
  // structure uses setAccelerationStructure at 16, so neither occupies a buffer.
  std::array<id<MTLBuffer>, 25> buffers;
  id<MTLBuffer> indirect;
  std::vector<std::pair<id<MTLComputePipelineState>, id<MTLIntersectionFunctionTable>>>
      containment_tables;
  size_t capacity;
  bool has_materials = false, has_pbr = false, has_pbr_emission = false, has_diffusion = false;

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

  void AdvanceQueues(id<MTLCommandBuffer> command) {
    auto encoder = [command computeCommandEncoder];
    encoder.label = @"Wavefront queue counts";
    [encoder setComputePipelineState:kernels.dispatch];
    [encoder setBuffer:buffers[11] offset:0 atIndex:0];
    [encoder setBuffer:indirect offset:0 atIndex:1];
    [encoder dispatchThreadgroups:MTLSizeMake(1, 1, 1) threadsPerThreadgroup:MTLSizeMake(1, 1, 1)];
    [encoder endEncoding];
  }

  void Dispatch(id<MTLCommandBuffer> command, id<MTLComputePipelineState> pipeline,
                const WFParameters &parameters, bool primary = false) {
    auto encoder = [command computeCommandEncoder];
    encoder.label = primary                                ? @"Wavefront primary rays"
                    : pipeline == kernels.intersect        ? @"Wavefront intersections"
                    : pipeline == kernels.shade            ? @"Wavefront surface shading"
                    : pipeline == kernels.shade_materials  ? @"Wavefront material shading"
                    : pipeline == kernels.shade_pbr        ? @"Wavefront OpenPBR shading"
                    : pipeline == kernels.shade_pbr_lights ? @"Wavefront OpenPBR light shading"
                    : pipeline == kernels.shade_diffusion  ? @"Wavefront normalized diffusion"
                                                           : @"Wavefront shadow rays";
    [encoder setComputePipelineState:pipeline];
    for (size_t i = 0; i < buffers.size(); ++i)
      if (buffers[i])
        [encoder setBuffer:buffers[i] offset:0 atIndex:i];
    [encoder setBytes:&parameters length:sizeof(parameters) atIndex:12];
    [encoder setAccelerationStructure:acceleration atBufferIndex:16];
    for (const auto &entry : containment_tables)
      if (entry.first == pipeline)
        [encoder setIntersectionFunctionTable:entry.second atBufferIndex:25];
    if (primary)
      [encoder dispatchThreadgroups:MTLSizeMake(std::max(1u, (parameters.active + 63) / 64), 1, 1)
              threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
    else
      [encoder dispatchThreadgroupsWithIndirectBuffer:indirect
                                 indirectBufferOffset:(pipeline == kernels.shade ||
                                                       pipeline == kernels.shade_materials ||
                                                       pipeline == kernels.shade_pbr ||
                                                       pipeline == kernels.shade_pbr_lights ||
                                                       pipeline == kernels.shade_diffusion)
                                                          ? 12
                                                          : (pipeline == kernels.shadow ? 24 : 0)
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
    buffers[9] = Buffer(capacity * 32);
    buffers[10] = Buffer(capacity * 144);
    buffers[11] = Buffer(7 * 4);
    buffers[13] = Buffer(capacity * 4);
    buffers[14] = kernels.directions;
    buffers[15] = kernels.blue;
    buffers[17] = Upload(scene.boundaries);
    buffers[18] = Buffer(scene.boundaries.size() > 1 ? capacity * 80 : 80);
    indirect = Buffer(9 * sizeof(uint32_t));
    buffers[19] = indirect;
    buffers[20] = Upload(scene.material_data);
    has_diffusion = !scene.diffusion_exits.empty();
    buffers[21] = Upload(scene.diffusion_exits);
    buffers[22] = Upload(scene.diffusion_samples);
    buffers[23] = Upload(scene.diffusion_pairs);
    buffers[24] = Buffer(has_diffusion ? capacity * 32 : 32);
    for (const auto &material : scene.materials) {
      has_materials |= material.type > 0 && material.type < 7;
      has_pbr |= material.type == 7;
      has_pbr_emission |= material.type == 7 && material.emissive;
    }
    kernels.PrepareMaterials(has_materials, has_pbr, has_pbr_emission, has_diffusion);
    for (id<MTLComputePipelineState> pipeline :
         {kernels.primary, kernels.intersect, kernels.shade, kernels.shade_materials,
          kernels.shade_pbr, kernels.shade_pbr_lights, kernels.shade_diffusion, kernels.shadow}) {
      if (!pipeline || !kernels.boundary_support)
        continue;
      auto table_descriptor =
          [MTLIntersectionFunctionTableDescriptor intersectionFunctionTableDescriptor];
      table_descriptor.functionCount = 1;
      auto table = [pipeline newIntersectionFunctionTableWithDescriptor:table_descriptor];
      auto handle = [pipeline
          functionHandleWithFunction:[kernels.library newFunctionWithName:@"containment_face"]];
      if (!table || !handle)
        throw MetalError("wavefront containment table", nil);
      [table setFunction:handle atIndex:0];
      [table setBuffer:buffers[0] offset:0 atIndex:0];
      containment_tables.emplace_back(pipeline, table);
    }
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
              const WFPixel *&pixels, WavefrontReport &report) override {
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
      // Media carry depth per path: internal flights and hidden crossings must
      // not spend ordinary bounce depth. Run bounded event batches until their
      // queue is empty, with no finite walk cap. Only batch completion crosses
      // to the CPU, so cancellation stays responsive during very long walks.
      const uint32_t batch = input.boundaries ? 16 : 4;
      for (uint64_t first = 0; input.boundaries || first <= input.max_depth; first += batch) {
        auto command = [queue commandBuffer];
        command.label = @"Wavefront transport batch";
        if (first == 0)
          Dispatch(command, kernels.primary, parameters, true);
        const uint32_t events =
            input.boundaries ? batch : std::min<uint64_t>(batch, input.max_depth + 1ULL - first);
        for (uint32_t event = 0; event < events; ++event) {
          parameters.depth = first + event;
          Dispatch(command, kernels.intersect, parameters);
          Dispatch(command,
                   has_pbr_emission ? kernels.shade_pbr_lights
                   : has_materials  ? kernels.shade_materials
                                    : kernels.shade,
                   parameters);
          if (has_pbr)
            Dispatch(command, kernels.shade_pbr, parameters);
          if (has_diffusion)
            Dispatch(command, kernels.shade_diffusion, parameters);
          Dispatch(command, kernels.shadow, parameters);
          AdvanceQueues(command);
          std::swap(buffers[7], buffers[8]);
        }
        if (!Wait(command, cancelled))
          return false;
        auto *counts = static_cast<uint32_t *>(buffers[11].contents);
        if (counts[4])
          throw std::runtime_error("Metal wavefront queue overflow; sample discarded");
        if (input.boundaries && counts[0] == 0)
          break;
      }
      auto *counts = static_cast<uint32_t *>(buffers[11].contents);
      if (counts[4])
        throw std::runtime_error("Metal wavefront queue overflow; sample discarded");
      report.discarded_paths += counts[5];
      report.scattering_events += counts[6];
      pixels = static_cast<const WFPixel *>(buffers[6].contents);
      return true;
    }
  }
};
} // namespace

std::unique_ptr<WavefrontSession> MakeMetalWavefront(const WavefrontScene &scene, size_t capacity) {
  @autoreleasepool {
    if (@available(macOS 11.0, *)) {
      return std::make_unique<MetalWavefront>(SharedPipelines(scene.boundaries.size() > 1), scene,
                                              capacity);
    }
    throw std::runtime_error("Metal wavefront requires macOS 11 or later");
  }
}

#ifdef NOT_CRAN
std::vector<WFBsdfResult> ProbeMetalMaterials(const std::vector<WFBsdfProbe> &inputs,
                                              const std::vector<WFVector> &data) {
  static_assert(sizeof(WFBsdfProbe) == 240 && sizeof(WFBsdfResult) == 48);
  std::vector<WFBsdfResult> output(inputs.size());
  if (inputs.empty())
    return output;
  @autoreleasepool {
    auto &shared = SharedPipelines();
    NSError *error = nil;
    auto function = [shared.library newFunctionWithName:@"probe_materials"];
    auto pipeline = [shared.device newComputePipelineStateWithFunction:function error:&error];
    if (!pipeline)
      throw MetalError("material probe pipeline", error);
    auto input = [shared.device newBufferWithBytes:inputs.data()
                                            length:inputs.size() * sizeof(WFBsdfProbe)
                                           options:MTLResourceStorageModeShared];
    auto result = [shared.device newBufferWithLength:output.size() * sizeof(WFBsdfResult)
                                             options:MTLResourceStorageModeShared];
    WFVector empty;
    auto parameters =
        [shared.device newBufferWithBytes:data.empty() ? &empty : data.data()
                                   length:std::max(size_t(1), data.size()) * sizeof(WFVector)
                                  options:MTLResourceStorageModeShared];
    auto queue = [shared.device newCommandQueue];
    auto command = [queue commandBuffer];
    auto encoder = [command computeCommandEncoder];
    [encoder setComputePipelineState:pipeline];
    [encoder setBuffer:input offset:0 atIndex:0];
    [encoder setBuffer:result offset:0 atIndex:1];
    [encoder setBuffer:parameters offset:0 atIndex:2];
    [encoder dispatchThreads:MTLSizeMake(inputs.size(), 1, 1)
        threadsPerThreadgroup:MTLSizeMake(32, 1, 1)];
    [encoder endEncoding];
    [command commit];
    [command waitUntilCompleted];
    if (command.status != MTLCommandBufferStatusCompleted)
      throw MetalError("material probe", command.error);
    std::memcpy(output.data(), result.contents, output.size() * sizeof(WFBsdfResult));
  }
  return output;
}
#endif
