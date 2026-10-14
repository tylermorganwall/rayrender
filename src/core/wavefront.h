#ifndef RAYRENDER_WAVEFRONT_H
#define RAYRENDER_WAVEFRONT_H

#include <cmath>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <skymodelr/prague_metal.h>
namespace skymodelr { class PragueMetalSnapshot; }

class hitable_list;
class RayCamera;
class adaptive_sampler;
class InfiniteLight;
class hitable;
class Transform;
class Atmosphere;
class Medium;
class VolumeScene;
class InfiniteAreaLight;
class WavefrontMotion;

// Explicit float4 storage keeps the C++/MSL ABI independent of native Float,
// SIMD configuration and Objective-C. No GPU record owns a native pointer.
struct alignas(16) WFVector {
  float x = 0, y = 0, z = 0, w = 0;
};
struct WFTriangle {
  // Vertex-colored native meshes retain source UVs in p[k].w/n[k].w.
  // Native curve contacts synthesize derivative vertices only for shading.
  WFVector p[3], n[3], uv01, uv2, geometric;
  uint32_t material = 0, light = UINT32_MAX, boundary = 0, flags = 0;
};
struct WFTexture {
  // Types: constant, checker, u-gradient, v-gradient, vertex color, image,
  // lat-long image. Image coordinates apply mapping.xy*uv + mapping.zw.
  WFVector a, b, c, mapping{1, 1, 0, 0};
  uint32_t type = 0, offset = 0, width = 0, height = 0;
};
struct WFSurface {
  // Affine world-to-object map. One record per placement/material mapping,
  // shared by all its faces; triangle uv2.w holds its integer index as bits.
  WFVector inverse[3];
  WFVector forward[3];
  // Bump heights are evaluated before native instance transforms. Keep that
  // frame separately from the object's texture coordinates.
  WFVector bump_inverse[3], bump_forward[3];
  WFVector axes; // Analytic type in w: mesh=0, sphere=1, ellipsoid=2.
  uint32_t alpha = UINT32_MAX, bump = UINT32_MAX, bump_transformed = 0, padding = 0;
  WFVector bump_parameters; // Height intensity; remaining lanes reserved.
};
struct WFQuadric {
  // Unit sphere <-> world affine maps, including ellipsoid axes and instances.
  // Rays retain their world-distance parameter when transformed; do not normalize
  // the object-space direction. One record replaces an entire tessellated sphere.
  WFVector inverse[3], forward[3];
  // |det(linear map)|, shading orientation, uniform world radius (or zero),
  // boundary orientation. The determinant supplies the area-sampling Jacobian.
  WFVector measure;
  uint32_t material = 0, light = UINT32_MAX, boundary = 0, surface = 0;
  // Nonzero count selects a distance-field program, stored immediately after
  // the analytic records. Offset is in bytes from the start of that array.
  uint32_t csg_offset = 0, csg_count = 0, curve_offset = 0, reserved = 0;
  WFVector csg_lo, csg_hi;
};
// Native cubic segment in object space. Width endpoints and ribbon normals
// belong to the original curve; width.zw retains this segment's global u range.
// Appended to the analytic buffer, addressed relative to the Quadric array.
struct WFCurve {
  WFVector cp[4], normal[2], width, fallback;
  uint32_t type = 0, padding[3]{};
};
class curve;
WFCurve ExportWavefrontCurve(const curve &);
struct WFCsgInstruction {
  WFVector a, b, c, d;
  uint32_t operation = 0, padding[3]{};
};
class ImplicitShape;
void ExportWavefrontCSG(const ImplicitShape &, std::vector<WFCsgInstruction> &);
void WavefrontCSGBounds(const WFQuadric &, WFVector &, WFVector &);
// One affine placement references a shared triangle range and its own mapping,
// emitter and medium metadata. IDs are TLAS descriptor indices plus one (0=flat).
struct WFMeshInstance {
  WFVector inverse[3], forward[3];
  uint32_t triangle_start = 0, triangle_count = 0, surface_base = 0, boundary = 0;
  uint32_t light_base = 0, mesh = 0;
  float orientation = 1;
  uint32_t padding = 0;
};
struct WFMeshRange {
  uint32_t start = 0, count = 0;
};
struct WFGeometryStorage {
  uint32_t quadric_offset = 48, instance_offset = 48, medium_offset = 48, field_offset = 48;
  uint32_t atmosphere_offset = 48, atmosphere_data_offset = 48, camera_offset = 48, motion_offset = 0;
  uint32_t light_offset = 48, padding[3]{};
};
// Leaves precede internal nodes; parent links allow emitter-hit MIS to evaluate
// the same position-dependent PMF as downward sampling, without a search.
struct WFLightNode {
  WFVector center_radius2;
  float power = 0;
  uint32_t left = UINT32_MAX, right = UINT32_MAX, parent = UINT32_MAX;
};
struct WFAtmosphereData {
  skymodelr_prague_gpu_layout model{};
  WFVector transform[3], origin, gain, settings, options;
  WFVector weights[16], transmission[16], solar[16], lunar[16], filter[129];
  uint32_t channels = 0, sky = 0, sun = 0, haze = 0;
  uint32_t query_altitude = 0, haze_filter = 0, haze_in_volumes = 0, padding = 0;
};
static_assert(sizeof(skymodelr_prague_gpu_layout) == 144 && sizeof(WFAtmosphereData) == 3376,
              "Prague Metal layout mismatch");
struct WFField {
  uint32_t offset = 0, nx = 0, ny = 0, nz = 0;
};
// Fields retain their native sample lattice. Interpolate temperature before
// converting it to RGB, and interpolate emission and its scale independently.
struct WFMedium {
  WFVector inverse[3], lo, hi, emission, controls;
  WFField density, temperature, red, green, blue, scale, blackbody, majorant;
};
struct WFMaterial {

  WFVector diffuse{1, 0, 0, 0}; // EON A, B, average loss; emission intensity in w.
  WFVector optical{1, 0, 0, 0}, attenuation;
  uint32_t texture = 0, emissive = 0, dielectric = 0, priority = 0;
  // Tagged BSDF parameters. Diffuse retains its small, separate shader path.
  WFVector eta, k, parameters, extra;
  uint32_t type = 0, distribution = 0, second_texture = UINT32_MAX, data = 0;
};
struct WFBoundary {
  WFVector sigma_a, sigma_s; // roughness in a.w, HG anisotropy in s.w.
  uint32_t material = 0, subsurface = 0, diffusion_offset = 0, medium = 0;
  // Diffusion uses RGB color/radius. Guided random walks use color.xyz for the
  // proposal poles and radius.w as eligibility; explicit media use radius.w
  // for keep_surface. Tags keep these layouts disjoint without enlarging them.
  WFVector color, radius, lo{INFINITY, INFINITY, INFINITY, 0},
      hi{-INFINITY, -INFINITY, -INFINITY, 0};
};
struct WFDiffusionExit {
  WFVector shape; // interior/exterior IOR, critical cosine squared, support, mass.
  WFVector normalization;
  uint32_t offset = 0, padding[3]{};
};
struct WFLight {
  // Types: image environment, infinite disk, point/spot, emissive primitive.
  WFVector forward, right, up, color;
  uint32_t type = 0, texture = 0, primitive = 0, flags = 0;
  uint32_t instance = 0, alias_offset = UINT32_MAX, alias_width = 0, alias_height = 0;
};
struct WFCamera {
  WFVector origin, lower_left, horizontal, vertical, right, up, forward;
  uint32_t orthographic = 0, flip_x = 0, padding[2]{};
  // Camera kind: perspective=0, orthographic=1, panorama=2, physical lens=3.
  WFVector shutter, motion_origin[2], motion_rotation[2], lens;
  uint32_t moving = 0, elements = 0, pupils = 0, simple_weighting = 0;
};
void ExportWavefrontCamera(RayCamera &, WFCamera &, std::vector<WFVector> &);
struct WFParameters {
  WFCamera camera;
  uint32_t width = 0, height = 0, active = 0, sample = 0;
  uint32_t lights = 0, depth = 0, roulette = 0, sampler = 0;
  uint32_t seed = 0, max_depth = 0, boundaries = 0, environments = 0;
  WFVector boundary_lo, boundary_hi;
};
struct WFPixel {
  // radiance.w is camera-ray transparency (1 = escaped, 0 = opaque), not
  // opacity or RGB path throughput. normal.w carries camera importance weight.
  WFVector radiance, normal, albedo;
};
struct WFPath {
  WFVector origin, direction, beta, previous;
  // Flags: delta predecessor=1, guide written=2, invalid membership=4,
  // diffusion entry=8, any earlier non-specular scattering=16.
  uint32_t dimension, depth, internal, flags;
};
struct WFAtmosphereRequest {
  WFVector origin, direction, weight, radiance; // atmospheric segment request
  uint32_t body = 0, glass = 0, event = 0, padding = 0;
};
struct WFShadow {
  WFVector origin, direction, contribution;
  uint32_t pixel, padding[3], ids[16], count, state_padding[3];
};
static_assert(sizeof(WFPath) == 80 && sizeof(WFAtmosphereRequest) == 80 && sizeof(WFShadow) == 144);
static_assert(sizeof(WFMedium) == 240 && sizeof(WFField) == 16 &&
              sizeof(WFMeshInstance) == 128 && sizeof(WFGeometryStorage) == 48 &&
              sizeof(WFLightNode) == 32 &&
              sizeof(WFQuadric) == 176 && sizeof(WFCurve) == 144 && sizeof(WFCsgInstruction) == 80 && sizeof(WFSurface) == 240 && sizeof(WFTriangle) == 160 &&
              sizeof(WFTexture) == 80 && sizeof(WFMaterial) == 144 && sizeof(WFBoundary) == 112 &&
              sizeof(WFDiffusionExit) == 48 && sizeof(WFLight) == 96 && sizeof(WFCamera) == 240 &&
              sizeof(WFParameters) == 320 && sizeof(WFPixel) == 48);

struct WavefrontScene {
  bool transparent_background = false;
  std::shared_ptr<WavefrontMotion> motion;
  std::vector<WFVector> camera_data;
  std::shared_ptr<const skymodelr::PragueMetalSnapshot> atmosphere_snapshot;
  WFAtmosphereData atmosphere_data;
  bool gpu_atmosphere = false;
  std::vector<WFTriangle> triangles;
  std::vector<WFQuadric> quadrics;
  std::vector<WFCsgInstruction> csg_instructions;
  std::vector<WFCurve> curves;
  std::vector<WFMeshRange> meshes;
  std::vector<WFMeshInstance> instances;
  size_t loose_triangles = 0;
  std::vector<WFVector> quadric_bounds; // Two entries (lo, hi) per primitive.
  std::vector<WFMaterial> materials;
  std::vector<WFTexture> textures;
  std::vector<WFVector> texels;
  std::vector<WFSurface> surfaces{WFSurface()};
  // Graphs=1, surface context=2, bump=4, alpha=8, general media=16,
  // invisible area emitters=64, finite light BVH=128, guided random walks=256,
  // native curves=512, independent geometry shutter times=1024.
  // Shared pipeline specialization key.
  uint32_t texture_features = 0;
  std::vector<WFLight> lights;
  std::vector<WFLightNode> light_nodes;
  std::vector<double> emission_estimates;
  std::vector<WFVector> material_data;
  std::vector<WFMedium> media{WFMedium()};
  std::vector<float> fields;
  std::vector<const Medium *> native_media{nullptr};
  const Atmosphere *atmosphere = nullptr;
  const Atmosphere *transport_atmosphere = nullptr;
  const VolumeScene *volume_scene = nullptr;
  const hitable_list *world = nullptr; // Native atmospheric shadows and camera coverage.
  size_t workers = 1;
  std::vector<WFDiffusionExit> diffusion_exits;
  std::vector<WFVector> diffusion_samples;
  std::vector<uint32_t> diffusion_pairs;
  uint32_t environments = 0; // Infinite emitters form the prefix of lights.
  // Index zero means vacuum/no boundary. Each placement gets its own ID, even
  // when instances share material and geometry pointers.
  std::vector<WFBoundary> boundaries{WFBoundary()};
  WFVector boundary_lo, boundary_hi;
  // Borrowed only on the host: refresh environment rotations after UI edits.
  std::vector<const InfiniteLight *> light_sources;
  std::vector<std::shared_ptr<InfiniteLight>> owned_lights;
  std::vector<std::pair<const InfiniteAreaLight *, const InfiniteLight *>> environment_versions;
};

struct WavefrontReport {
  bool requested = false, used = false, scene_reused = false;
  std::string fallback;
  double upload_seconds = 0, sample_seconds = 0, atmosphere_seconds = 0;
  size_t triangles = 0, completed_samples = 0;
  size_t tessellated_primitives = 0, discarded_paths = 0, scattering_events = 0;
  size_t analytic_primitives = 0, unique_meshes = 0, mesh_instances = 0;
  size_t csg_failed_queries = 0;
};

// Device work is transactional: cancellation never commits a partial sample.
// Queues have at most one continuation and one shadow per active pixel.
class WavefrontSession {
public:
  virtual ~WavefrontSession() = default;
  virtual bool Render(const WFParameters &, const std::vector<uint32_t> &active,
                      const std::vector<WFLight> &lights, const std::function<bool()> &cancelled,
                      const WFPixel *&pixels, WavefrontReport &report) = 0;
};

// Own the host descriptors and their dependent device session together. An
// animation keeps this object alive across frames; ordinary renders use a local
// one. The scene must not move while a session holds references into it.
class WavefrontContext {
public:
  WavefrontScene scene;
  std::unique_ptr<WavefrontSession> session;
  WavefrontContext() = default;
  WavefrontContext(const WavefrontContext &) = delete;
  WavefrontContext &operator=(const WavefrontContext &) = delete;
  void BeginFrame(const hitable_list &, const hitable_list &, RayCamera &,
                  size_t capacity, int sampler, WavefrontReport &, size_t workers);
private:
  const hitable_list *world = nullptr, *lights = nullptr;
  size_t pixel_capacity = 0, camera_capacity = 0;
  int sampling_method = -1;
  WavefrontReport setup_report;
};
// Empty output means this primitive has no faithful static triangulation.
std::vector<WFTriangle> TriangulateWavefrontPrimitive(const hitable &, const Transform &);
std::unique_ptr<WavefrontSession> MakeMetalWavefront(const WavefrontScene &, size_t capacity);
// Capability probe only: do not compile pipelines or allocate scene resources.
bool MetalWavefrontAvailable();
std::unique_ptr<WavefrontSession> PrepareWavefront(const hitable_list &, const hitable_list &,
                                                   RayCamera &, size_t capacity, int sampler,
                                                   WavefrontScene &, WavefrontReport &, size_t workers = 1);
bool RenderWavefrontSample(WavefrontSession &, WavefrontScene &, RayCamera &, adaptive_sampler &,
                           size_t sample, int sampler, uint32_t seed, size_t depth, size_t roulette,
                           float clamp, const std::function<bool()> &cancelled, WavefrontReport &);
bool WavefrontLightsChanged(const WavefrontScene &);
// Exact native Prague queries are an explicit wavefront service, not a fallback
// integrator. Shared Metal buffers are accessed only after command completion.
class WavefrontAtmosphere {
public:
  WavefrontAtmosphere(const WavefrontScene &, size_t capacity);
  ~WavefrontAtmosphere();
  void Begin(const std::vector<uint32_t> &active, const float *times = nullptr);
  bool Coverage(const WFParameters &, const std::vector<uint32_t> &active,
                WFPath *, WFPixel *, const std::function<bool()> &cancelled);
  bool Paths(const WFParameters &, const uint32_t *queue, size_t count,
             WFAtmosphereRequest *, WFPath *, WFPixel *, const std::function<bool()> &cancelled);
  bool Shadows(const WFParameters &, WFShadow *, size_t count, WFPixel *,
               const std::function<bool()> &cancelled);
private:
  struct Data;
  std::unique_ptr<Data> data;
};

#if defined(NOT_CRAN) && defined(RAY_HAS_METAL_BVH)
struct WFCameraProbe { WFVector film_lens, time; };
struct WFCameraResult { WFVector origin, direction; };
std::vector<WFCameraResult> ProbeMetalCamera(const WFCamera &,const std::vector<WFVector> &,
                                            const std::vector<WFCameraProbe> &);
struct WFHazeProbe { WFVector values[2], start, delta; };
std::vector<WFCameraResult> ProbeMetalHaze(const std::vector<WFHazeProbe> &);
struct WFSubsurfaceProbe { WFVector pole, axis, random, wi, rates; };
std::vector<WFCameraResult> ProbeMetalSubsurface(const std::vector<WFSubsurfaceProbe> &);
// Development-only numerical cross-backend probes. No R API or CPU hot-path hook.
struct WFGeometryProbe {
  WFVector origin, direction; // w holds minimum/maximum ray distance.
  uint32_t body = 0, padding[3]{}; // padding[1]: absolute shutter time as float bits.
};
struct WFGeometryResult {
  WFVector point, normal, uv; // point.w is ray distance; normal.w is hit/miss.
  uint32_t primitive = UINT32_MAX, material = 0, boundary = 0, padding = 0;
};
std::vector<WFGeometryResult> ProbeMetalGeometry(const WavefrontScene &,
                                                 const std::vector<WFGeometryProbe> &);
struct WFBsdfProbe {
  WFMaterial material;
  WFVector rho, transmission, view, direction, random, context;
};
struct WFBsdfResult {
  WFVector evaluated, direction, weight;
};
struct WFEnvironmentProbe {
  uint32_t bits = 0;
  float latitude = 0;
  uint32_t padding[2]{};
};
std::vector<WFVector> ProbeMetalLightSelection(const WavefrontScene &, const std::vector<WFVector> &);
std::vector<WFVector> ProbeMetalEnvironment(const WavefrontScene &, const WFLight &,
                                            const std::vector<WFEnvironmentProbe> &);
struct WFTextureProbe {
  WFVector uv, position, object_position, normal, object_normal, footprint;
};
std::vector<WFVector> ProbeMetalTextures(const WavefrontScene &,
                                         const std::vector<WFTextureProbe> &);
std::vector<WFBsdfResult> ProbeMetalMaterials(const std::vector<WFBsdfProbe> &,
                                              const std::vector<WFVector> &data);
#endif

#endif
