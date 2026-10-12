#include "wavefront_medium_probe.h"
#include "sobol_directions.h"
#include "wavefront.h"
#include "wavefront_motion.h"
#include <skymodelr/prague.hpp>
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
namespace wavefront_tables {
#include "samplerBlueNoise.h"
}
#include <array>
#include <chrono>
#include <cstring>
#include <map>
#include <mutex>
#include <stdexcept>
#include <thread>

namespace {
const char *source = "#include <metal_stdlib>\nusing namespace metal;\n"
#include <skymodelr/prague_metal.inc>
                     "\nnamespace wf_openpbr {\n"
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
  id<MTLFunction> analytic_function;
  id<MTLComputePipelineState> primary, intersect, shade, shade_materials, shade_pbr,
      shade_pbr_lights, shade_diffusion, shadow, shadow_transport, shadow_lights,
      dispatch, atmosphere_paths, atmosphere_dispatch, coverage;
  id<MTLBuffer> directions, blue;
  std::array<id<MTLComputePipelineState>, 4> motion{};
  bool boundary_support, analytic_support, instance_support, mixed_support, atmosphere_support;
  uint32_t texture_features;
  std::mutex material_mutex;

  id<MTLComputePipelineState> Compile(NSString *name, bool defer_shadow_lighting = false,
                                    bool primary_coverage = false) {
    NSError *error = nil;
    auto constants = [MTLFunctionConstantValues new];
    [constants setConstantValue:&boundary_support type:MTLDataTypeBool atIndex:0];
    [constants setConstantValue:&texture_features type:MTLDataTypeUInt atIndex:1];
    [constants setConstantValue:&analytic_support type:MTLDataTypeBool atIndex:2];
    [constants setConstantValue:&instance_support type:MTLDataTypeBool atIndex:3];
    [constants setConstantValue:&atmosphere_support type:MTLDataTypeBool atIndex:4];
    [constants setConstantValue:&defer_shadow_lighting type:MTLDataTypeBool atIndex:5];
    [constants setConstantValue:&primary_coverage type:MTLDataTypeBool atIndex:6];
    if (analytic_support && !analytic_function) {
      analytic_function=[library newFunctionWithName:@"intersect_quadric"
                                     constantValues:constants error:&error];
      if (!analytic_function) throw MetalError("analytic intersection specialization",error);
    }
    auto function = [library newFunctionWithName:name constantValues:constants error:&error];
    if (!function)
      throw MetalError("wavefront function specialization", error);
    auto descriptor = [MTLComputePipelineDescriptor new];
    descriptor.computeFunction = function;
    if (boundary_support || mixed_support) {
      auto linked = [MTLLinkedFunctions linkedFunctions];
      NSMutableArray<id<MTLFunction>> *functions = [NSMutableArray array];
      [functions
          addObject:[library newFunctionWithName:mixed_support ? @"containment_face_instanced"
                                                               : @"containment_face"]];
      if (analytic_support)
        [functions addObject:analytic_function];
      linked.functions = functions;
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
  void PrepareMaterials(bool basic, bool pbr, bool emission, bool diffusion, bool transparent) {
    std::lock_guard<std::mutex> lock(material_mutex);
    if (basic && !shade_materials)
      shade_materials = Compile(@"shade_materials");
    if (pbr && !shade_pbr)
      shade_pbr = Compile(@"shade_pbr");
    if (emission && !shade_pbr_lights)
      shade_pbr_lights = Compile(@"shade_pbr_lights");
    if (diffusion && !shade_diffusion)
      shade_diffusion = Compile(@"shade_diffusion");
    if (transparent && !coverage)
      coverage = Compile(@"shadows", false, true);
  }
  void PrepareMotion() {
    std::lock_guard<std::mutex> lock(material_mutex);
    if (motion[0] && motion[1] && motion[2] && motion[3]) return;
    NSArray<NSString *> *names = @[@"motion_triangles", @"motion_quadrics",
                                  @"motion_surfaces", @"motion_instances"];
    // A failed pipeline allocation must not leave a partially populated cache
    // looking ready when a later render retries the Metal backend.
    for (size_t i=0; i<4; ++i)
      if (!motion[i]) motion[i]=Compile(names[i]);
  }
  explicit Pipelines(bool boundary_support, uint32_t textures, bool analytics, bool instances, bool atmosphere)
      : boundary_support(boundary_support), analytic_support(analytics),
        instance_support(instances), mixed_support(analytics || instances),
        atmosphere_support(atmosphere), texture_features(textures) {
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
    if (atmosphere) {
      atmosphere_paths = Compile(@"atmosphere_paths");
      atmosphere_dispatch = Compile(@"atmosphere_dispatch");
      // Opaque visibility already has a short any-hit path. An additional pass
      // pays off when boundary/alpha traversal also computes RGB attenuation.
      if (boundary_support || (texture_features & (8 | 64))) {
        shadow_lights = Compile(@"shadow_lights");
        // Specialization removes unused spectral/haze state from scalar work.
        shadow_transport = Compile(@"shadows", true);
      }
    }
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

Pipelines &SharedPipelines(bool boundaries = false, uint32_t textures = 0, bool analytics = false,
                           bool instances = false, bool atmosphere = false) {
  // Specialization removes graph registers, alpha traversal and bump footprint
  // storage from ordinary opaque scenes. Sessions retain their own resources.
  static std::mutex mutex;
  static std::map<uint32_t, std::unique_ptr<Pipelines>> pipelines;
  std::lock_guard<std::mutex> lock(mutex);
  uint32_t key = (textures << 4) | uint32_t(boundaries) | (uint32_t(analytics) << 1) |
                 (uint32_t(instances) << 2) | (uint32_t(atmosphere) << 3);
  auto &entry = pipelines[key];
  if (!entry)
    entry = std::make_unique<Pipelines>(boundaries, textures, analytics, instances, atmosphere);
  return *entry;
}

class MetalWavefront final : public WavefrontSession {
  Pipelines &kernels;
  id<MTLCommandQueue> queue;
  id<MTLAccelerationStructure> acceleration, mixed_acceleration;
  NSArray<id<MTLAccelerationStructure>> *bottom_levels;
  // Slots match SCENE_ARGUMENTS. Parameters use setBytes at 12; the acceleration
  // structure uses setAccelerationStructure at 16, so neither occupies a buffer.
  std::array<id<MTLBuffer>, 29> buffers;
  id<MTLBuffer> indirect;
  std::vector<std::pair<id<MTLComputePipelineState>, id<MTLIntersectionFunctionTable>>>
      intersection_tables;
  size_t capacity;
  WFGeometryStorage geometry_storage;
  bool has_materials = false, has_pbr = false, has_pbr_emission = false, has_diffusion = false;
  bool has_haze = false;
  std::unique_ptr<WavefrontAtmosphere> atmosphere;
  const WavefrontScene &scene_data;
  bool gpu_atmosphere = false;
  std::vector<id<MTLAccelerationStructure>> refit_structures;
  std::vector<MTLAccelerationStructureDescriptor *> refit_descriptors;
  id<MTLBuffer> refit_scratch = nil, motion_boxes = nil, motion_instances = nil;
  bool moving_triangles = false, moving_quadrics = false, moving_instances = false;
  bool whole_scene_frame = false;
  id<MTLBuffer> reference_vertices = nil;
  std::array<id<MTLBuffer>,4> motion_references{};
  std::array<uint32_t,4> motion_counts{};
  id<MTLBuffer> motion_poses = nil;

  void EncodeMotion(id<MTLCommandBuffer> command) {
    if (!scene_data.motion || !scene_data.motion_changed)
      return;
    const auto &s = scene_data;
    // Only placement matrices and small transport metadata cross the host
    // boundary. All primitive records and refit inputs are written on device.
    std::memcpy(motion_poses.contents, s.motion->poses.data(),
                s.motion->poses.size()*sizeof(WavefrontMotion::Pose));
    std::memcpy(buffers[17].contents, s.boundaries.data(),
                s.boundaries.size() * sizeof(WFBoundary));
    auto *bytes = static_cast<char *>(buffers[28].contents);
    if (whole_scene_frame)
      std::memcpy(bytes + geometry_storage.motion_offset, s.motion->inverse_frame,
                  3 * sizeof(WFVector));
    std::memcpy(bytes + geometry_storage.medium_offset, s.media.data(), s.media.size()*sizeof(WFMedium));
    for (size_t kind=0; kind<4; ++kind) {
      if (!motion_counts[kind]) continue;
      auto encoder=[command computeCommandEncoder];
      encoder.label=@"Wavefront geometry motion";
      [encoder setComputePipelineState:kernels.motion[kind]];
      id<MTLBuffer> target=kind==0 ? buffers[0] : kind==2 ? buffers[26] : buffers[28];
      size_t offset=kind==1 ? geometry_storage.quadric_offset : kind==3 ? geometry_storage.instance_offset : 0;
      [encoder setBuffer:target offset:offset atIndex:0];
      [encoder setBuffer:motion_references[kind] offset:0 atIndex:1];
      [encoder setBuffer:motion_poses offset:0 atIndex:2];
      id<MTLBuffer> descriptors=kind==1 ? motion_boxes : motion_instances;
      [encoder setBuffer:descriptors ? descriptors : motion_poses offset:0 atIndex:3];
      [encoder setBytes:&motion_counts[kind] length:sizeof(uint32_t) atIndex:4];
      bool refit=descriptors != nil;
      [encoder setBytes:&refit length:sizeof(bool) atIndex:5];
      [encoder dispatchThreads:MTLSizeMake(motion_counts[kind],1,1)
          threadsPerThreadgroup:MTLSizeMake(64,1,1)];
      [encoder endEncoding];
    }
    // Separate encoders order child bounds before the parent TLAS refit.
    for (size_t i = 0; i < refit_structures.size(); ++i) {
      auto encoder = [command accelerationStructureCommandEncoder];
      [encoder refitAccelerationStructure:refit_structures[i]
                               descriptor:refit_descriptors[i]
                              destination:nil
                            scratchBuffer:refit_scratch
                      scratchBufferOffset:0];
      [encoder endEncoding];
    }
  }

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

  id<MTLAccelerationStructure> BuildAcceleration(MTLAccelerationStructureDescriptor *descriptor,
                                                bool changes = false) {
    if (changes) descriptor.usage |= MTLAccelerationStructureUsageRefit;
    auto sizes = [kernels.device accelerationStructureSizesWithDescriptor:descriptor];
    auto result = [kernels.device newAccelerationStructureWithSize:sizes.accelerationStructureSize];
    auto scratch = Buffer(sizes.buildScratchBufferSize);
    if (!result)
      throw MetalError("wavefront acceleration structure", nil);
    auto command = [queue commandBuffer];
    auto encoder = [command accelerationStructureCommandEncoder];
    [encoder buildAccelerationStructure:result
                             descriptor:descriptor
                          scratchBuffer:scratch
                    scratchBufferOffset:0];
    [encoder endEncoding];
    Wait(command, [] { return false; });
    if (changes) {
      refit_structures.push_back(result);
      refit_descriptors.push_back(descriptor);
      if (!refit_scratch || refit_scratch.length < sizes.refitScratchBufferSize)
        refit_scratch = Buffer(sizes.refitScratchBufferSize);
    }
    return result;
  }

  void RegisterIntersections(id<MTLComputePipelineState> pipeline) {
    auto table_descriptor =
        [MTLIntersectionFunctionTableDescriptor intersectionFunctionTableDescriptor];
    table_descriptor.functionCount = kernels.analytic_support ? 2 : 1;
    auto table = [pipeline newIntersectionFunctionTableWithDescriptor:table_descriptor];
    auto handle = [pipeline
        functionHandleWithFunction:[kernels.library
                                       newFunctionWithName:kernels.mixed_support
                                                               ? @"containment_face_instanced"
                                                               : @"containment_face"]];
    if (!table || !handle)
      throw MetalError("wavefront containment table", nil);
    [table setFunction:handle atIndex:0];
    [table setBuffer:buffers[0] offset:0 atIndex:0];
    if (kernels.analytic_support) {
      auto function = kernels.analytic_function;
      auto analytic = [pipeline functionHandleWithFunction:function];
      if (!analytic)
        throw MetalError("wavefront analytic intersection table", nil);
      [table setFunction:analytic atIndex:1];
      [table setBuffer:buffers[28]
                offset:whole_scene_frame ? geometry_storage.motion_offset + 3 * sizeof(WFVector)
                                         : geometry_storage.quadric_offset
               atIndex:1];
    }
    [table setBuffer:buffers[28] offset:geometry_storage.instance_offset atIndex:2];
    [table setBuffer:buffers[11] offset:7*sizeof(uint32_t) atIndex:3];
    intersection_tables.emplace_back(pipeline, table);
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
    bool spectral = gpu_atmosphere &&
        (pipeline == kernels.atmosphere_paths || pipeline == kernels.shadow_lights ||
         pipeline == kernels.shadow || pipeline == kernels.coverage);
    if (spectral) {
      auto encoder = [command computeCommandEncoder];
      [encoder setComputePipelineState:kernels.atmosphere_dispatch];
      [encoder setBuffer:buffers[11] offset:0 atIndex:0];
      [encoder setBuffer:indirect offset:0 atIndex:1];
      [encoder dispatchThreadgroups:MTLSizeMake(1,1,1) threadsPerThreadgroup:MTLSizeMake(1,1,1)];
      [encoder endEncoding];
    }
    auto encoder = [command computeCommandEncoder];
    encoder.label = primary                                ? @"Wavefront primary rays"
                    : pipeline == kernels.intersect        ? @"Wavefront intersections"
                    : pipeline == kernels.shade            ? @"Wavefront surface shading"
                    : pipeline == kernels.shade_materials  ? @"Wavefront material shading"
                    : pipeline == kernels.shade_pbr        ? @"Wavefront OpenPBR shading"
                    : pipeline == kernels.shade_pbr_lights ? @"Wavefront OpenPBR light shading"
                    : pipeline == kernels.shade_diffusion  ? @"Wavefront normalized diffusion"
                    : pipeline == kernels.atmosphere_paths ? @"Wavefront spectral atmosphere"
                    : pipeline == kernels.shadow_lights    ? @"Wavefront shadow lighting"
                    : pipeline == kernels.coverage         ? @"Wavefront camera coverage"
                                                           : @"Wavefront shadow rays";
    [encoder setComputePipelineState:pipeline];
    for (size_t i = 0; i < buffers.size(); ++i)
      if (buffers[i])
        [encoder setBuffer:buffers[i] offset:0 atIndex:i];
    [encoder setBytes:&parameters length:sizeof(parameters) atIndex:12];
    // TLAS references do not make their BLAS resources resident automatically.
    // Declare every indirect child for each encoder, including after preview
    // edits or denoising has submitted work on another command queue.
    if (kernels.mixed_support)
      for (id<MTLAccelerationStructure> child in bottom_levels)
        [encoder useResource:child usage:MTLResourceUsageRead];
    [encoder setAccelerationStructure:acceleration atBufferIndex:16];
    if (kernels.mixed_support)
      [encoder setAccelerationStructure:mixed_acceleration atBufferIndex:29];
    for (const auto &entry : intersection_tables)
      if (entry.first == pipeline)
        [encoder setIntersectionFunctionTable:entry.second
                                atBufferIndex:kernels.mixed_support ? 30 : 25];
    if (primary)
      [encoder dispatchThreadgroups:MTLSizeMake(std::max(1u, (parameters.active + 63) / 64), 1, 1)
              threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
    else {
      // Each indirect dispatch record contains three uint32 dimensions.
      // Spectral queues expand each connection into one 32-lane SIMD group.
      NSUInteger offset = 0;
      if (spectral)
        offset = pipeline == kernels.atmosphere_paths || pipeline == kernels.coverage ? 36 : 48;
      else if (pipeline == kernels.shade || pipeline == kernels.shade_materials ||
               pipeline == kernels.shade_pbr || pipeline == kernels.shade_pbr_lights ||
               pipeline == kernels.shade_diffusion)
        offset = 12;
      else if (pipeline == kernels.shadow || pipeline == kernels.shadow_transport)
        offset = 24;
      [encoder dispatchThreadgroupsWithIndirectBuffer:indirect
                                 indirectBufferOffset:offset
                                threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
    }
    [encoder endEncoding];
  }

public:
  MetalWavefront(Pipelines &kernels, const WavefrontScene &scene, size_t capacity)
      : kernels(kernels), capacity(capacity), scene_data(scene) {
    if (scene.motion) {
      whole_scene_frame = scene.motion->whole_scene_frame;
      auto changes = [&](size_t kind) {
        const auto &owners = scene.motion->owners[kind];
        return std::any_of(owners.begin(), owners.end(), [](uint32_t owner) { return owner != 0; });
      };
      moving_triangles = changes(0);
      moving_quadrics = changes(1);
      moving_instances = changes(3);
    }
    for (const auto &material : scene.materials) {
      has_materials |= material.type > 0 && material.type < 7;
      has_pbr |= material.type == 7;
      has_pbr_emission |= material.type == 7 && material.emissive;
    }
    has_diffusion = !scene.diffusion_exits.empty();
    kernels.PrepareMaterials(has_materials, has_pbr, has_pbr_emission, has_diffusion,
                            scene.transparent_background);
    gpu_atmosphere = scene.gpu_atmosphere;
    if (gpu_atmosphere && (kernels.atmosphere_paths.threadExecutionWidth != 32 ||
                           kernels.shadow.threadExecutionWidth != 32 ||
                           (scene.transparent_background && kernels.coverage.threadExecutionWidth != 32) ||
                           (kernels.shadow_lights && kernels.shadow_lights.threadExecutionWidth != 32)))
      gpu_atmosphere = false;
    has_haze = scene.transport_atmosphere != nullptr;
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
    buffers[11] = Buffer(8 * 4);
    buffers[13] = Buffer(capacity * 4);
    buffers[14] = kernels.directions;
    buffers[15] = kernels.blue;
    buffers[17] = Upload(scene.boundaries);
    buffers[18] = Buffer(scene.boundaries.size() > 1 ? capacity * 80 : 80);
    indirect = Buffer(15 * sizeof(uint32_t));
    buffers[19] = indirect;
    buffers[20] = Upload(scene.material_data);
    buffers[21] = Upload(scene.diffusion_exits);
    buffers[22] = Upload(scene.diffusion_samples);
    buffers[23] = Upload(scene.diffusion_pairs);
    // Keep diffusion-only scenes at 32 bytes per pixel; append atmosphere requests only when needed.
    buffers[24] = Buffer((has_diffusion || scene.atmosphere || (scene.texture_features & 256) ? capacity : 1) * 32 +
                         (scene.atmosphere ? capacity * (sizeof(WFAtmosphereRequest) + 64) : 0));
    if (scene.atmosphere && !gpu_atmosphere)
      atmosphere = std::make_unique<WavefrontAtmosphere>(scene, capacity);
    buffers[26] = Upload(scene.surfaces);
    // Metal exposes 31 buffer slots, all already in use. Pack geometry, media,
    // camera and atmosphere metadata with explicit offsets, not a pointer ABI.
    // Moving scenes also retain immutable quadrics for reference-frame tracing.
    auto instances = scene.instances;
    if (instances.empty())
      instances.resize((!scene.triangles.empty()) + (!scene.quadrics.empty()));
    size_t csg_bytes=scene.csg_instructions.size()*sizeof(WFCsgInstruction);
    size_t curve_bytes = scene.curves.size() * sizeof(WFCurve);
    size_t instance_offset = sizeof(WFGeometryStorage) + scene.quadrics.size() * sizeof(WFQuadric)+csg_bytes+curve_bytes;
    if (instance_offset > UINT32_MAX)
      throw std::runtime_error("Metal analytic metadata exceeds the 32-bit buffer offset capacity");
    geometry_storage.instance_offset = uint32_t(instance_offset);
    size_t storage_size = geometry_storage.instance_offset +
                          std::max(size_t(1), instances.size()) * sizeof(WFMeshInstance);
    if (storage_size + scene.media.size() * sizeof(WFMedium) +
            scene.fields.size() * sizeof(float) > UINT32_MAX)
      throw std::runtime_error("Metal medium metadata exceeds the 32-bit buffer capacity");
    geometry_storage.medium_offset = uint32_t(storage_size);
    storage_size += scene.media.size() * sizeof(WFMedium);
    geometry_storage.field_offset = uint32_t(storage_size);
    storage_size += scene.fields.size() * sizeof(float);
    storage_size = (storage_size + 15) & ~size_t(15);
    if (whole_scene_frame) {
      size_t motion_bytes = 3 * sizeof(WFVector) + scene.quadrics.size() * sizeof(WFQuadric)+csg_bytes+curve_bytes;
      if (storage_size + motion_bytes > UINT32_MAX)
        throw std::runtime_error("Metal motion metadata exceeds the 32-bit scene buffer capacity");
      geometry_storage.motion_offset = uint32_t(storage_size);
      storage_size += motion_bytes;
    }
    geometry_storage.camera_offset = uint32_t(storage_size);
    storage_size += scene.camera_data.size() * sizeof(WFVector);
    if (storage_size + scene.light_nodes.size()*sizeof(WFLightNode) > UINT32_MAX)
      throw std::runtime_error("Metal light metadata exceeds the 32-bit buffer capacity");
    geometry_storage.light_offset = uint32_t(storage_size);
    storage_size += scene.light_nodes.size()*sizeof(WFLightNode);
    geometry_storage.atmosphere_offset = uint32_t(storage_size);
    size_t coefficient_count = 0;
    const float *coefficients = nullptr;
    if (gpu_atmosphere) {
      coefficients = scene.atmosphere_snapshot->data(coefficient_count);
      storage_size += sizeof(WFAtmosphereData);
    }
    geometry_storage.atmosphere_data_offset = uint32_t(storage_size);
    if (coefficient_count > (UINT32_MAX - storage_size) / sizeof(float))
      throw std::runtime_error("Metal atmosphere exceeds the 32-bit scene buffer capacity");
    storage_size += coefficient_count * sizeof(float);
    buffers[28] = Buffer(storage_size);
    auto *bytes = static_cast<unsigned char *>(buffers[28].contents);
    std::memcpy(bytes, &geometry_storage, sizeof(geometry_storage));
    if (!scene.light_nodes.empty())
      std::memcpy(bytes + geometry_storage.light_offset, scene.light_nodes.data(),
                  scene.light_nodes.size()*sizeof(WFLightNode));
    if (whole_scene_frame) {
      std::memcpy(bytes + geometry_storage.motion_offset, scene.motion->inverse_frame,
                  3 * sizeof(WFVector));
      if (!scene.quadrics.empty())
        std::memcpy(bytes + geometry_storage.motion_offset + 3 * sizeof(WFVector),
                    scene.quadrics.data(), scene.quadrics.size() * sizeof(WFQuadric));
      if(csg_bytes) std::memcpy(bytes+geometry_storage.motion_offset+3*sizeof(WFVector)+
          scene.quadrics.size()*sizeof(WFQuadric),scene.csg_instructions.data(),csg_bytes);
      if (curve_bytes) std::memcpy(bytes + geometry_storage.motion_offset + 3*sizeof(WFVector) +
          scene.quadrics.size()*sizeof(WFQuadric) + csg_bytes, scene.curves.data(), curve_bytes);
    }
    if (gpu_atmosphere) {
      std::memcpy(bytes + geometry_storage.atmosphere_offset, &scene.atmosphere_data, sizeof(WFAtmosphereData));
      std::memcpy(bytes + geometry_storage.atmosphere_data_offset, coefficients, coefficient_count*sizeof(float));
    }
    if (!scene.quadrics.empty())
      std::memcpy(bytes + geometry_storage.quadric_offset, scene.quadrics.data(),
                  scene.quadrics.size() * sizeof(WFQuadric));
    if(csg_bytes) std::memcpy(bytes+geometry_storage.quadric_offset+
        scene.quadrics.size()*sizeof(WFQuadric),scene.csg_instructions.data(),csg_bytes);
    if (curve_bytes) std::memcpy(bytes + geometry_storage.quadric_offset +
        scene.quadrics.size()*sizeof(WFQuadric) + csg_bytes, scene.curves.data(), curve_bytes);
    if (!instances.empty())
      std::memcpy(bytes + geometry_storage.instance_offset, instances.data(),
                  instances.size() * sizeof(WFMeshInstance));
    std::memcpy(bytes + geometry_storage.medium_offset, scene.media.data(),
                scene.media.size() * sizeof(WFMedium));
    if (!scene.fields.empty())
      std::memcpy(bytes + geometry_storage.field_offset, scene.fields.data(),
                  scene.fields.size() * sizeof(float));

    buffers[27] = Buffer(scene.texture_features & 4 ? capacity * 64 : 64);
    for (id<MTLComputePipelineState> pipeline :
         {kernels.primary, kernels.intersect, kernels.shade, kernels.shade_materials,
          kernels.shade_pbr, kernels.shade_pbr_lights, kernels.shade_diffusion,
          kernels.shadow, kernels.shadow_transport, kernels.coverage}) {
      if (!pipeline || !(kernels.boundary_support || kernels.mixed_support))
        continue;
      RegisterIntersections(pipeline);
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
    NSMutableArray<MTLAccelerationStructureGeometryDescriptor *> *geometries =
        [NSMutableArray array];
    size_t loose = scene.meshes.empty() ? scene.triangles.size() : scene.loose_triangles;
    if (whole_scene_frame && moving_triangles)
      reference_vertices = Buffer(loose * sizeof(WFTriangle), scene.triangles.data());
    if (loose) {
      auto geometry = [MTLAccelerationStructureTriangleGeometryDescriptor descriptor];
      geometry.vertexBuffer = reference_vertices ? reference_vertices : buffers[0];
      geometry.vertexStride = sizeof(WFVector);
      geometry.indexBuffer = index_buffer;
      geometry.indexType = MTLIndexTypeUInt32;
      geometry.triangleCount = loose;
      geometry.opaque = YES;
      geometry.intersectionFunctionTableOffset = 0;
      [geometries addObject:geometry];
    }
    id<MTLBuffer> box_buffer = nil;
    if (!scene.quadrics.empty()) {
      std::vector<MTLAxisAlignedBoundingBox> boxes(scene.quadrics.size());
      for (size_t i = 0; i < boxes.size(); ++i) {
        const auto &lo = scene.quadric_bounds[2 * i];
        const auto &hi = scene.quadric_bounds[2 * i + 1];
        boxes[i].min = MTLPackedFloat3{lo.x, lo.y, lo.z};
        boxes[i].max = MTLPackedFloat3{hi.x, hi.y, hi.z};
      }
      box_buffer = Upload(boxes);
      if (moving_quadrics && !whole_scene_frame)
        motion_boxes = box_buffer;
      auto geometry = [MTLAccelerationStructureBoundingBoxGeometryDescriptor descriptor];
      geometry.boundingBoxBuffer = box_buffer;
      geometry.boundingBoxCount = boxes.size();
      geometry.boundingBoxStride = sizeof(MTLAxisAlignedBoundingBox);
      geometry.intersectionFunctionTableOffset = 1;
      [geometries addObject:geometry];
    }
    for (const auto &mesh : scene.meshes) {
      auto geometry = [MTLAccelerationStructureTriangleGeometryDescriptor descriptor];
      geometry.vertexBuffer = buffers[0];
      geometry.vertexStride = sizeof(WFVector);
      geometry.indexBuffer = index_buffer;
      geometry.indexBufferOffset = mesh.start * 3 * sizeof(uint32_t);
      geometry.indexType = MTLIndexTypeUInt32;
      geometry.triangleCount = mesh.count;
      geometry.opaque = YES;
      geometry.intersectionFunctionTableOffset = 0;
      [geometries addObject:geometry];
    }
    NSMutableArray<id<MTLAccelerationStructure>> *children = [NSMutableArray array];
    size_t child_index = 0;
    for (MTLAccelerationStructureGeometryDescriptor *geometry in geometries) {
      auto descriptor = [MTLPrimitiveAccelerationStructureDescriptor descriptor];
      descriptor.geometryDescriptors = @[ geometry ];
      bool changes =
          (loose && child_index == 0 && moving_triangles) ||
          (!scene.quadrics.empty() && child_index == size_t(loose != 0) && moving_quadrics);
      changes = changes && !whole_scene_frame;
      [children addObject:BuildAcceleration(descriptor, changes)];
      ++child_index;
    }
    bottom_levels = children;
    acceleration = children[0];
    if (kernels.mixed_support) {
      // Mesh prototypes each have one BLAS. Placements reference it through a
      // TLAS transform; loose triangles and world-space quadrics use identities.
      size_t prefix = (loose != 0) + !scene.quadrics.empty();
      size_t count = scene.instances.empty() ? prefix : scene.instances.size();
      std::vector<MTLAccelerationStructureInstanceDescriptor> instances(count);
      for (size_t i = 0; i < count; ++i) {
        auto &instance = instances[i];
        instance.mask = 0xff;
        if (i < prefix) {
          instance.transformationMatrix.columns[0].x = 1;
          instance.transformationMatrix.columns[1].y = 1;
          instance.transformationMatrix.columns[2].z = 1;
          instance.accelerationStructureIndex = i;
        } else {
          const auto &placement = scene.instances[i];
          const auto *a = placement.forward;
          instance.transformationMatrix.columns[0] = MTLPackedFloat3{a[0].x, a[1].x, a[2].x};
          instance.transformationMatrix.columns[1] = MTLPackedFloat3{a[0].y, a[1].y, a[2].y};
          instance.transformationMatrix.columns[2] = MTLPackedFloat3{a[0].z, a[1].z, a[2].z};
          instance.transformationMatrix.columns[3] = MTLPackedFloat3{a[0].w, a[1].w, a[2].w};
          instance.accelerationStructureIndex = prefix + placement.mesh;
        }
      }
      auto descriptor = [MTLInstanceAccelerationStructureDescriptor descriptor];
      descriptor.instancedAccelerationStructures = children;
      descriptor.instanceCount = instances.size();
      descriptor.instanceDescriptorBuffer = Upload(instances);
      if (moving_instances && !whole_scene_frame)
        motion_instances = descriptor.instanceDescriptorBuffer;
      mixed_acceleration = BuildAcceleration(
          descriptor,
          !whole_scene_frame && (moving_triangles || moving_quadrics || moving_instances));
    }
    if (scene.motion) {
      static_assert(sizeof(WavefrontMotion::Record<WFTriangle>) == 176 &&
                    sizeof(WavefrontMotion::Record<WFQuadric>) == 192 &&
                    sizeof(WavefrontMotion::Record<WFSurface>) == 256 &&
                    sizeof(WavefrontMotion::Record<WFMeshInstance>) == 144 &&
                    sizeof(MTLAccelerationStructureInstanceDescriptor) == 64,
                    "motion reference ABI mismatch");
      kernels.PrepareMotion();
      auto upload=[&](size_t kind,const auto &references) {
        if (references.size()>UINT32_MAX) throw std::runtime_error("Too many moving Metal primitives");
        motion_counts[kind]=references.size(); motion_references[kind]=Upload(references);
      };
      upload(0,scene.motion->TriangleReferences()); upload(1,scene.motion->QuadricReferences());
      upload(2,scene.motion->SurfaceReferences()); upload(3,scene.motion->InstanceReferences());
      motion_poses=Buffer((scene.motion->nodes.size()+1)*sizeof(WavefrontMotion::Pose));
      scene.motion->EnableGPU();
    }
  }

#ifdef NOT_CRAN
  std::vector<WFGeometryResult> Probe(const std::vector<WFGeometryProbe> &input) {
    auto pipeline = kernels.Compile(@"probe_geometry");
    RegisterIntersections(pipeline);
    auto probes = Upload(input);
    auto results = Buffer(input.size() * sizeof(WFGeometryResult));
    auto command = [queue commandBuffer];
    EncodeMotion(command);
    auto encoder = [command computeCommandEncoder];
    [encoder setComputePipelineState:pipeline];
    [encoder setBuffer:buffers[0] offset:0 atIndex:0];
    [encoder setBuffer:probes offset:0 atIndex:1];
    [encoder setBuffer:results offset:0 atIndex:2];
    // TLAS references do not make their BLAS resources resident automatically.
    // Declare every indirect child for each encoder, including after preview
    // edits or denoising has submitted work on another command queue.
    if (kernels.mixed_support)
      for (id<MTLAccelerationStructure> child in bottom_levels)
        [encoder useResource:child usage:MTLResourceUsageRead];
    [encoder setAccelerationStructure:acceleration atBufferIndex:16];
    if (kernels.mixed_support)
      [encoder setAccelerationStructure:mixed_acceleration atBufferIndex:29];
    [encoder setIntersectionFunctionTable:intersection_tables.back().second
                            atBufferIndex:kernels.mixed_support ? 30 : 25];
    [encoder setBuffer:buffers[26] offset:0 atIndex:26];
    [encoder setBuffer:buffers[28] offset:0 atIndex:28];
    [encoder dispatchThreads:MTLSizeMake(input.size(), 1, 1)
        threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];
    [encoder endEncoding];
    Wait(command, [] { return false; });
    const auto *data = static_cast<const WFGeometryResult *>(results.contents);
    return std::vector<WFGeometryResult>(data, data + input.size());
  }
#endif

  bool Render(const WFParameters &input, const std::vector<uint32_t> &active,
              const std::vector<WFLight> &lights, const std::function<bool()> &cancelled,
              const WFPixel *&pixels, WavefrontReport &report) override {
    @autoreleasepool {
      if (active.size() > capacity || input.width * size_t(input.height) > capacity)
        throw std::runtime_error("wavefront image exceeds queue capacity");
      if (cancelled())
        return false;
      // Primary rays can query CSG for initial membership. Reset on the host
      // before dispatch, rather than racing the first lane against those queries.
      static_cast<uint32_t *>(buffers[11].contents)[7]=0;
      if (scene_data.motion_changed && !scene_data.light_nodes.empty())
        std::memcpy(static_cast<char *>(buffers[28].contents) + geometry_storage.light_offset,
                    scene_data.light_nodes.data(),scene_data.light_nodes.size()*sizeof(WFLightNode));
      if (!active.empty())
        std::memcpy(buffers[13].contents, active.data(), active.size() * sizeof(uint32_t));
      if (!lights.empty())
        std::memcpy(buffers[4].contents, lights.data(), lights.size() * sizeof(WFLight));
      WFParameters parameters = input;
      if (!scene_data.camera_data.empty())
        std::memcpy(static_cast<char *>(buffers[28].contents) + geometry_storage.camera_offset,
                    scene_data.camera_data.data(), scene_data.camera_data.size()*sizeof(WFVector));
      parameters.camera.padding[0] = static_cast<uint32_t>(capacity);
      parameters.camera.padding[1] |= uint32_t(has_haze) | (uint32_t(gpu_atmosphere) << 1);
      if (gpu_atmosphere)
        std::memcpy(static_cast<char *>(buffers[28].contents) + geometry_storage.atmosphere_offset,
                    &scene_data.atmosphere_data, sizeof(WFAtmosphereData));
      if (atmosphere)
        atmosphere->Begin(active);
      // Media carry depth per path: internal flights and hidden crossings must
      // not spend ordinary bounce depth. Run bounded event batches until their
      // queue is empty, with no finite walk cap. Only batch completion crosses
      // to the CPU, so cancellation stays responsive during very long walks.
      const uint32_t batch = atmosphere ? 1 : (input.boundaries ? 16 : 4);
      for (uint64_t first = 0; input.boundaries || first <= input.max_depth; first += batch) {
        auto command = [queue commandBuffer];
        command.label = @"Wavefront transport batch";
        if (first == 0) {
          // Motion, child/parent refits, and traversal share one ordered command
          // buffer. A separate submit/wait costs more than small poses need.
          EncodeMotion(command);
          Dispatch(command, kernels.primary, parameters, true);
          if (parameters.camera.padding[1] & 4) {
            if (atmosphere && has_haze) {
              // The existing CPU spectral service also owns coverage for a
              // hybrid sky. Read the generated rays only after GPU completion.
              if (!Wait(command, cancelled)) return false;
              auto atmosphere_started = std::chrono::steady_clock::now();
              if (!atmosphere->Coverage(parameters, active,
                      static_cast<WFPath *>(buffers[5].contents),
                      static_cast<WFPixel *>(buffers[6].contents), cancelled)) return false;
              report.atmosphere_seconds += std::chrono::duration<double>(
                  std::chrono::steady_clock::now() - atmosphere_started).count();
              command = [queue commandBuffer];
            } else
              Dispatch(command, kernels.coverage, parameters);
          }
        }
        const uint32_t events =
            input.boundaries ? batch : std::min<uint64_t>(batch, input.max_depth + 1ULL - first);
        for (uint32_t event = 0; event < events; ++event) {
          parameters.depth = first + event;
          Dispatch(command, kernels.intersect, parameters);
          if (gpu_atmosphere) Dispatch(command, kernels.atmosphere_paths, parameters);
          if (atmosphere) {
            if (!Wait(command, cancelled))
              return false;
            auto *counts = static_cast<uint32_t *>(buffers[11].contents);
            auto atmosphere_started = std::chrono::steady_clock::now();
            if (!atmosphere->Paths(parameters, static_cast<uint32_t *>(buffers[7].contents), counts[0],
                              reinterpret_cast<WFAtmosphereRequest *>(
                                  static_cast<char *>(buffers[24].contents) + capacity * 32),
                              static_cast<WFPath *>(buffers[5].contents),
                              static_cast<WFPixel *>(buffers[6].contents),cancelled)) return false;
            report.atmosphere_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - atmosphere_started).count();
            if (cancelled()) return false;
            command = [queue commandBuffer];
          }
          Dispatch(command,
                   has_pbr_emission ? kernels.shade_pbr_lights
                   : has_materials  ? kernels.shade_materials
                                    : kernels.shade,
                   parameters);
          if (has_pbr)
            Dispatch(command, kernels.shade_pbr, parameters);
          if (has_diffusion)
            Dispatch(command, kernels.shade_diffusion, parameters);
          bool split_shadows = gpu_atmosphere && !has_haze && kernels.shadow_transport;
          Dispatch(command, split_shadows ? kernels.shadow_transport : kernels.shadow,
                   parameters);
          if (split_shadows)
            Dispatch(command, kernels.shadow_lights, parameters);
          if (atmosphere) {
            if (!Wait(command, cancelled))
              return false;
            auto *counts = static_cast<uint32_t *>(buffers[11].contents);
            auto atmosphere_started = std::chrono::steady_clock::now();
            if (!atmosphere->Shadows(parameters, static_cast<WFShadow *>(buffers[10].contents), counts[2],
                                static_cast<WFPixel *>(buffers[6].contents),cancelled)) return false;
            const auto *shadow_results = static_cast<const WFShadow *>(buffers[10].contents);
            for (size_t j=0; j<counts[2]; ++j)
              report.discarded_paths += shadow_results[j].padding[2] != 0;
            report.atmosphere_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - atmosphere_started).count();
            if (cancelled()) return false;
            command = [queue commandBuffer];
          }
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
      report.csg_failed_queries += counts[7];
      pixels = static_cast<const WFPixel *>(buffers[6].contents);
      return true;
    }
  }
};
} // namespace

#if defined(NOT_CRAN)
std::vector<WFGeometryResult> ProbeMetalGeometry(const WavefrontScene &scene,
                                                 const std::vector<WFGeometryProbe> &input) {
  @autoreleasepool {
    auto &pipelines = SharedPipelines(scene.boundaries.size() > 1, scene.texture_features,
                                      !scene.quadrics.empty(), !scene.meshes.empty());
    MetalWavefront session(pipelines, scene, 1);
    return session.Probe(input);
  }
}
#endif

std::unique_ptr<WavefrontSession> MakeMetalWavefront(const WavefrontScene &scene, size_t capacity) {
  @autoreleasepool {
    if (@available(macOS 11.0, *)) {
      return std::make_unique<MetalWavefront>(
          SharedPipelines(scene.boundaries.size() > 1, scene.texture_features,
                          !scene.quadrics.empty(), !scene.meshes.empty(), scene.atmosphere != nullptr),
          scene, capacity);
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
std::vector<WFVector> ProbeMetalTextures(const WavefrontScene &scene,
                                         const std::vector<WFTextureProbe> &inputs) {
  std::vector<WFVector> output(inputs.size());
  if (inputs.empty())
    return output;
  @autoreleasepool {
    auto &shared = SharedPipelines(false, scene.texture_features);
    auto pipeline = shared.Compile(@"probe_textures");
    auto upload = [&](const void *data, size_t bytes) {
      WFVector empty;
      return [shared.device newBufferWithBytes:bytes ? data : &empty
                                        length:std::max(bytes, sizeof(empty))
                                       options:MTLResourceStorageModeShared];
    };
    auto input = upload(inputs.data(), inputs.size() * sizeof(WFTextureProbe));
    auto result = upload(output.data(), output.size() * sizeof(WFVector));
    auto textures = upload(scene.textures.data(), scene.textures.size() * sizeof(WFTexture));
    auto texels = upload(scene.texels.data(), scene.texels.size() * sizeof(WFVector));
    auto queue = [shared.device newCommandQueue];
    auto command = [queue commandBuffer];
    auto encoder = [command computeCommandEncoder];
    [encoder setComputePipelineState:pipeline];
    [encoder setBuffer:input offset:0 atIndex:0];
    [encoder setBuffer:result offset:0 atIndex:1];
    [encoder setBuffer:textures offset:0 atIndex:2];
    [encoder setBuffer:texels offset:0 atIndex:3];
    [encoder dispatchThreads:MTLSizeMake(inputs.size(), 1, 1)
        threadsPerThreadgroup:MTLSizeMake(32, 1, 1)];
    [encoder endEncoding];
    [command commit];
    [command waitUntilCompleted];
    if (command.status != MTLCommandBufferStatusCompleted)
      throw MetalError("texture probe", command.error);
    std::memcpy(output.data(), result.contents, output.size() * sizeof(WFVector));
  }
  return output;
}
std::vector<WFVector> ProbeMetalEnvironment(const WavefrontScene &scene, const WFLight &light,
                                            const std::vector<WFEnvironmentProbe> &inputs) {
  std::vector<WFVector> output(inputs.size());
  if (inputs.empty())
    return output;
  @autoreleasepool {
    auto &shared = SharedPipelines(false, scene.texture_features);
    auto pipeline = shared.Compile(@"probe_environment");
    auto upload = [&](const void *data, size_t bytes) {
      WFVector empty;
      return [shared.device newBufferWithBytes:bytes ? data : &empty
                                        length:std::max(bytes, sizeof(empty))
                                       options:MTLResourceStorageModeShared];
    };
    auto input = upload(inputs.data(), inputs.size() * sizeof(WFEnvironmentProbe));
    auto result = upload(output.data(), output.size() * sizeof(WFVector));
    auto texels = upload(scene.texels.data(), scene.texels.size() * sizeof(WFVector));
    auto queue = [shared.device newCommandQueue];
    auto command = [queue commandBuffer];
    auto encoder = [command computeCommandEncoder];
    [encoder setComputePipelineState:pipeline];
    [encoder setBuffer:input offset:0 atIndex:0];
    [encoder setBuffer:result offset:0 atIndex:1];
    [encoder setBytes:&light length:sizeof(light) atIndex:2];
    [encoder setBuffer:texels offset:0 atIndex:3];
    [encoder dispatchThreads:MTLSizeMake(inputs.size(), 1, 1)
        threadsPerThreadgroup:MTLSizeMake(32, 1, 1)];
    [encoder endEncoding];
    [command commit];
    [command waitUntilCompleted];
    if (command.status != MTLCommandBufferStatusCompleted)
      throw MetalError("environment probe", command.error);
    std::memcpy(output.data(), result.contents, output.size() * sizeof(WFVector));
  }
  return output;
}

std::vector<WFVector> ProbeMetalLightSelection(const WavefrontScene &scene,
                                             const std::vector<WFVector> &inputs) {
  std::vector<WFVector> output(inputs.size());
  if (inputs.empty()) return output;
  @autoreleasepool {
    auto &shared=SharedPipelines(false,scene.texture_features);
    auto pipeline=shared.Compile(@"probe_light_selection");
    auto input=[shared.device newBufferWithBytes:inputs.data() length:inputs.size()*sizeof(WFVector)
                       options:MTLResourceStorageModeShared];
    auto result=[shared.device newBufferWithLength:output.size()*sizeof(WFVector)
                        options:MTLResourceStorageModeShared];
    WFLightNode empty;
    auto nodes=[shared.device newBufferWithBytes:scene.light_nodes.empty() ? &empty : scene.light_nodes.data()
                      length:std::max(size_t(1),scene.light_nodes.size())*sizeof(WFLightNode)
                      options:MTLResourceStorageModeShared];
    WFParameters p; p.lights=scene.lights.size(); p.environments=scene.environments;
    auto queue=[shared.device newCommandQueue];
    auto command=[queue commandBuffer];
    auto encoder=[command computeCommandEncoder];
    [encoder setComputePipelineState:pipeline];
    [encoder setBuffer:input offset:0 atIndex:0];
    [encoder setBuffer:result offset:0 atIndex:1];
    [encoder setBytes:&p length:sizeof(p) atIndex:2];
    [encoder setBuffer:nodes offset:0 atIndex:3];
    [encoder dispatchThreads:MTLSizeMake(inputs.size(),1,1) threadsPerThreadgroup:MTLSizeMake(32,1,1)];
    [encoder endEncoding]; [command commit]; [command waitUntilCompleted];
    if (command.status != MTLCommandBufferStatusCompleted)
      throw MetalError("light selection probe",command.error);
    std::memcpy(output.data(),result.contents,output.size()*sizeof(WFVector));
  }
  return output;
}
#endif

#ifdef NOT_CRAN
namespace {
template <class Input>
std::vector<WFCameraResult> camera_haze_probe(NSString *name, const std::vector<Input> &input,
                                              const WFCamera &camera,
                                              const std::vector<WFVector> &data) {
  std::vector<WFCameraResult> result(input.size());
  if (input.empty())
    return result;
  @autoreleasepool {
    auto &shared = SharedPipelines();
    auto pipeline = shared.Compile(name);
    auto source = [shared.device newBufferWithBytes:input.data()
                                             length:input.size() * sizeof(Input)
                                            options:MTLResourceStorageModeShared];
    auto output = [shared.device newBufferWithLength:result.size() * sizeof(WFCameraResult)
                                             options:MTLResourceStorageModeShared];
    WFVector empty;
    auto values =
        [shared.device newBufferWithBytes:data.empty() ? &empty : data.data()
                                   length:std::max(size_t(1), data.size()) * sizeof(WFVector)
                                  options:MTLResourceStorageModeShared];
    auto queue = [shared.device newCommandQueue];
    auto command = [queue commandBuffer];
    auto encoder = [command computeCommandEncoder];
    [encoder setComputePipelineState:pipeline];
    [encoder setBuffer:source offset:0 atIndex:0];
    [encoder setBuffer:output offset:0 atIndex:1];
    [encoder setBuffer:values offset:0 atIndex:2];
    [encoder setBytes:&camera length:sizeof(camera) atIndex:3];
    [encoder dispatchThreads:MTLSizeMake(input.size(), 1, 1)
        threadsPerThreadgroup:MTLSizeMake(32, 1, 1)];
    [encoder endEncoding];
    [command commit];
    [command waitUntilCompleted];
    if (command.status != MTLCommandBufferStatusCompleted)
      throw MetalError("camera/haze probe", command.error);
    std::memcpy(result.data(), output.contents, result.size() * sizeof(WFCameraResult));
  }
  return result;
}
}
std::vector<WFCameraResult> ProbeMetalCamera(const WFCamera &camera,
                                             const std::vector<WFVector> &data,
                                             const std::vector<WFCameraProbe> &input) {
  return camera_haze_probe(@"probe_camera", input, camera, data);
}
std::vector<WFCameraResult> ProbeMetalHaze(const std::vector<WFHazeProbe> &input) {
  return camera_haze_probe(@"probe_haze", input, WFCamera(), {});
}
std::vector<WFCameraResult> ProbeMetalSubsurface(const std::vector<WFSubsurfaceProbe> &input) {
  return camera_haze_probe(@"probe_subsurface", input, WFCamera(), {});
}
#endif
#ifdef NOT_CRAN
std::vector<WFMediumResult> ProbeMetalMedium(const std::vector<WFMediumProbe> &input,
                                           const std::vector<float> &values) {
  std::vector<WFMediumResult> result(input.size());
  if (input.empty()) return result;
  @autoreleasepool {
    auto &shared = SharedPipelines();
    auto pipeline = shared.Compile(@"probe_medium");
    auto source = [shared.device newBufferWithBytes:input.data()
        length:input.size() * sizeof(WFMediumProbe) options:MTLResourceStorageModeShared];
    auto output = [shared.device newBufferWithLength:result.size() * sizeof(WFMediumResult)
        options:MTLResourceStorageModeShared];
    auto field = [shared.device newBufferWithBytes:values.data()
        length:values.size() * sizeof(float) options:MTLResourceStorageModeShared];
    auto queue = [shared.device newCommandQueue];
    auto command = [queue commandBuffer];
    auto encoder = [command computeCommandEncoder];
    [encoder setComputePipelineState:pipeline];
    [encoder setBuffer:source offset:0 atIndex:0];
    [encoder setBuffer:output offset:0 atIndex:1];
    [encoder setBuffer:field offset:0 atIndex:2];
    [encoder dispatchThreads:MTLSizeMake(input.size(), 1, 1)
        threadsPerThreadgroup:MTLSizeMake(32, 1, 1)];
    [encoder endEncoding];
    [command commit];
    [command waitUntilCompleted];
    if (command.status != MTLCommandBufferStatusCompleted)
      throw MetalError("medium flight probe", command.error);
    std::memcpy(result.data(), output.contents, result.size() * sizeof(WFMediumResult));
  }
  return result;
}
#endif
