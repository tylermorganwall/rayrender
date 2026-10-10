#ifndef RAYRENDER_WAVEFRONT_H
#define RAYRENDER_WAVEFRONT_H

#include <cmath>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

class hitable_list;
class RayCamera;
class adaptive_sampler;
class InfiniteLight;
class hitable;
class Transform;

// Explicit float4 storage keeps the C++/MSL ABI independent of native Float,
// SIMD configuration and Objective-C. No GPU record owns a native pointer.
struct alignas(16) WFVector {
  float x = 0, y = 0, z = 0, w = 0;
};
struct WFTriangle {
  // Curve-only metadata in unused float4 lanes: p[].w transforms a unit error
  // cube; geometric.w/uv2.z hold endpoint widths, n[].w the world radius.
  // Positions remain xyz float3s.
  // Vertex-colored native meshes instead retain source UVs in p[k].w/n[k].w.
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
  uint32_t material = 0, subsurface = 0, diffusion_offset = 0, padding = 0;
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
};
struct WFCamera {
  WFVector origin, lower_left, horizontal, vertical, right, up, forward;
  uint32_t orthographic = 0, flip_x = 0, padding[2]{};
};
struct WFParameters {
  WFCamera camera;
  uint32_t width = 0, height = 0, active = 0, sample = 0;
  uint32_t lights = 0, depth = 0, roulette = 0, sampler = 0;
  uint32_t seed = 0, max_depth = 0, boundaries = 0, environments = 0;
  WFVector boundary_lo, boundary_hi;
};
struct WFPixel {
  WFVector radiance, normal, albedo;
};
static_assert(sizeof(WFQuadric) == 128 && sizeof(WFSurface) == 240 && sizeof(WFTriangle) == 160 &&
              sizeof(WFTexture) == 80 && sizeof(WFMaterial) == 144 && sizeof(WFBoundary) == 112 &&
              sizeof(WFDiffusionExit) == 48 && sizeof(WFLight) == 80 && sizeof(WFCamera) == 128 &&
              sizeof(WFParameters) == 208 && sizeof(WFPixel) == 48);

struct WavefrontScene {
  std::vector<WFTriangle> triangles;
  std::vector<WFQuadric> quadrics;
  std::vector<WFVector> quadric_bounds; // Two entries (lo, hi) per primitive.
  std::vector<WFMaterial> materials;
  std::vector<WFTexture> textures;
  std::vector<WFVector> texels;
  std::vector<WFSurface> surfaces{WFSurface()};
  uint32_t texture_features = 0; // Graphs=1, surface context=2, bump=4, alpha=8.
  std::vector<WFLight> lights;
  std::vector<WFVector> material_data;
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
};

struct WavefrontReport {
  bool requested = false, used = false;
  std::string fallback;
  double upload_seconds = 0, sample_seconds = 0;
  size_t triangles = 0, completed_samples = 0;
  size_t tessellated_primitives = 0, discarded_paths = 0, scattering_events = 0;
  size_t analytic_primitives = 0;
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
// Empty output means this primitive has no faithful static triangulation.
std::vector<WFTriangle> TriangulateWavefrontPrimitive(const hitable &, const Transform &);
std::unique_ptr<WavefrontSession> MakeMetalWavefront(const WavefrontScene &, size_t capacity);
std::unique_ptr<WavefrontSession> PrepareWavefront(const hitable_list &, const hitable_list &,
                                                   RayCamera &, size_t capacity, int sampler,
                                                   WavefrontScene &, WavefrontReport &);
bool RenderWavefrontSample(WavefrontSession &, WavefrontScene &, RayCamera &, adaptive_sampler &,
                           size_t sample, int sampler, uint32_t seed, size_t depth, size_t roulette,
                           float clamp, const std::function<bool()> &cancelled, WavefrontReport &);

#if defined(NOT_CRAN) && defined(RAY_HAS_METAL_BVH)
// Development-only numerical cross-backend probes. No R API or CPU hot-path hook.
struct WFGeometryProbe {
  WFVector origin, direction; // w holds minimum/maximum ray distance.
  uint32_t body = 0, padding[3]{};
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
struct WFTextureProbe {
  WFVector uv, position, object_position, normal, object_normal, footprint;
};
std::vector<WFVector> ProbeMetalTextures(const WavefrontScene &,
                                         const std::vector<WFTextureProbe> &);
std::vector<WFBsdfResult> ProbeMetalMaterials(const std::vector<WFBsdfProbe> &,
                                              const std::vector<WFVector> &data);
#endif

#endif
