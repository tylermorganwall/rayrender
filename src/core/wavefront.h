#ifndef RAYRENDER_WAVEFRONT_H
#define RAYRENDER_WAVEFRONT_H

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

class hitable_list;
class RayCamera;
class adaptive_sampler;
class InfiniteLight;

// Explicit float4 storage keeps the C++/MSL ABI independent of native Float,
// SIMD configuration and Objective-C. No GPU record owns a native pointer.
struct alignas(16) WFVector {
  float x = 0, y = 0, z = 0, w = 0;
};
struct WFTriangle {
  WFVector p[3], n[3], uv01, uv2, geometric;
  uint32_t material = 0, light = UINT32_MAX, padding[2]{};
};
struct WFTexture {
  // Types: constant, checker, u-gradient, v-gradient, vertex color, image,
  // lat-long image. Image coordinates apply mapping.xy*uv + mapping.zw.
  WFVector a, b, c, mapping{1, 1, 0, 0};
  uint32_t type = 0, offset = 0, width = 0, height = 0;
};
struct WFMaterial {
  WFVector diffuse{1, 0, 0, 0}; // EON A, B, average loss; emission intensity in w.
  uint32_t texture = 0, emissive = 0, padding[2]{};
};
struct WFLight {
  // Types: image environment, infinite disk, point/spot, emissive triangle.
  WFVector forward, right, up, color;
  uint32_t type = 0, texture = 0, triangle = 0, flags = 0;
};
struct WFCamera {
  WFVector origin, lower_left, horizontal, vertical, right, up, forward;
  uint32_t orthographic = 0, flip_x = 0, padding[2]{};
};
struct WFParameters {
  WFCamera camera;
  uint32_t width = 0, height = 0, active = 0, sample = 0;
  uint32_t lights = 0, depth = 0, roulette = 0, sampler = 0;
  uint32_t seed = 0, max_depth = 0, padding[2]{};
};
struct WFPixel {
  WFVector radiance, normal, albedo;
};
static_assert(sizeof(WFTriangle) == 160 && sizeof(WFTexture) == 80 && sizeof(WFMaterial) == 32 &&
              sizeof(WFLight) == 80 && sizeof(WFCamera) == 128 && sizeof(WFParameters) == 176 &&
              sizeof(WFPixel) == 48);

struct WavefrontScene {
  std::vector<WFTriangle> triangles;
  std::vector<WFMaterial> materials;
  std::vector<WFTexture> textures;
  std::vector<WFVector> texels;
  std::vector<WFLight> lights;
  // Borrowed only on the host: refresh environment rotations after UI edits.
  std::vector<const InfiniteLight *> light_sources;
  std::vector<std::shared_ptr<InfiniteLight>> owned_lights;
};

struct WavefrontReport {
  bool requested = false, used = false;
  std::string fallback;
  double upload_seconds = 0, sample_seconds = 0;
  size_t triangles = 0, completed_samples = 0;
};

// Device work is transactional: cancellation never commits a partial sample.
// Queues have at most one continuation and one shadow per active pixel.
class WavefrontSession {
public:
  virtual ~WavefrontSession() = default;
  virtual bool Render(const WFParameters &, const std::vector<uint32_t> &active,
                      const std::vector<WFLight> &lights, const std::function<bool()> &cancelled,
                      const WFPixel *&pixels) = 0;
};
std::unique_ptr<WavefrontSession> MakeMetalWavefront(const WavefrontScene &, size_t capacity);
std::unique_ptr<WavefrontSession> PrepareWavefront(const hitable_list &, const hitable_list &,
                                                   RayCamera &, size_t capacity, int sampler,
                                                   WavefrontScene &, WavefrontReport &);
bool RenderWavefrontSample(WavefrontSession &, WavefrontScene &, RayCamera &, adaptive_sampler &,
                           size_t sample, int sampler, uint32_t seed, size_t depth, size_t roulette,
                           float clamp, const std::function<bool()> &cancelled, WavefrontReport &);

#endif
