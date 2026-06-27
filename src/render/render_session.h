#ifndef RAYRENDER_RENDER_SESSION_H
#define RAYRENDER_RENDER_SESSION_H

#include "../core/buildscene.h"
#include "../core/camera.h"
#include "../core/color.h"
#include "../hitables/hitable.h"
#include "../hitables/hitablelist.h"
#include "../materials/material.h"
#include "../materials/texture.h"
#include "../materials/texturecache.h"
#include "../math/RayMatrix.h"
#include "../math/rng.h"
#include "../math/transformcache.h"

#include <memory>
#include <string>
#include <vector>

namespace rayrender {
namespace render {

struct RenderOptions {
  int nx = 0;
  int ny = 0;
  int samples = 0;
  int numbercores = 1;
  bool progress_bar = false;
  bool verbose = false;
  int debug_channel = 0;
  Float min_variance = 0;
  int min_adaptive_size = 0;
  IntegratorType integrator_type = IntegratorType::ShadowRays;
};

struct CameraFrameState {
  point3f lookfrom;
  point3f lookat;
  vec3f camera_up;
  Float fov = 40;
  Float aperture = 0;
  Float focus_distance = 1;
  Float orthox = 1;
  Float orthoy = 1;
};

struct LegacyCameraStatic {
  int nx = 0;
  int ny = 0;
  Float shutteropen = 0;
  Float shutterclose = 1;
  Float film_size = 1;
  Float camera_scale = 1;
  Float iso = 100;
  std::vector<Float> lens_data;
  std::string missing_lens_error = "No lens data passed in lens descriptor file.";
};

struct LegacyEnvironmentInput {
  bool ambient_light = false;
  bool hasbackground = false;
  point3f backgroundhigh;
  point3f backgroundlow;
  std::string background;
  Float rotate_env = 0;
  Float intensity_env = 1;
  Float world_radius = 1;
  vec3f world_center;
};

struct LegacyEnvironment {
  bool ambient_light = false;
  bool hasbackground = false;
  point3f backgroundhigh;
  point3f backgroundlow;
  std::shared_ptr<texture> background_texture;
  std::shared_ptr<material> background_material;
  std::shared_ptr<hitable> background_sphere;
};

struct LegacyOutputBuffers {
  LegacyOutputBuffers(int nx, int ny);

  RayMatrix rgb_output;
  RayMatrix draw_rgb_output;
  RayMatrix alpha_output;
  RayMatrix normal_output;
  RayMatrix albedo_output;
};

struct RenderFrameSetup {
  std::unique_ptr<RayCamera> camera;
  LegacyOutputBuffers outputs;
};

class CompiledScene {
public:
  TransformCache transform_cache;
  TextureCache texture_cache;
  std::vector<Float*> textures;
  std::vector<unsigned char*> alpha_textures;
  std::vector<unsigned char*> bump_textures;
  std::vector<unsigned char*> roughness_textures;
  std::vector<std::shared_ptr<material> > shared_materials;
  hitable_list importance_sample_objects;
  std::vector<std::shared_ptr<hitable> > instanced_objects;
  std::vector<std::shared_ptr<hitable_list> > instance_importance_sampled;
  std::vector<std::shared_ptr<alpha_texture> > alpha;
  std::vector<std::shared_ptr<bump_texture> > bump;
  std::vector<std::shared_ptr<roughness_texture> > roughness;
  std::vector<int> texture_idx;
  std::shared_ptr<hitable> world_bvh;
  hitable_list world;

  void Compile(Rcpp::List& scene,
               Rcpp::IntegerVector& shape,
               Float shutteropen,
               Float shutterclose,
               int bvh_type,
               bool verbose,
               random_gen& rng);

  void AssembleWorld(const std::shared_ptr<hitable>& background_sphere);
  void AddEnvironmentToImportanceSampling(
    const std::shared_ptr<hitable>& background_sphere,
    bool importance_only_background,
    bool hasbackground);
};

class RenderSession {
public:
  RenderSession(RenderOptions options, CompiledScene& scene);

  RenderFrameSetup RenderFrame(const CameraFrameState& frame,
                               const LegacyCameraStatic& camera_static,
                               TransformCache& transform_cache) const;
  std::unique_ptr<RayCamera> CreateCamera(
    const CameraFrameState& frame,
    const LegacyCameraStatic& camera_static,
    TransformCache& transform_cache) const;
  LegacyOutputBuffers CreateOutputBuffers() const;

private:
  RenderOptions options_;
  CompiledScene* scene_;
};

std::vector<Float> FlattenLensData(const Rcpp::NumericMatrix& real_camera_info);

std::unique_ptr<RayCamera> MakeLegacyCamera(
  const CameraFrameState& frame,
  const LegacyCameraStatic& camera_static,
  TransformCache& transform_cache);

LegacyEnvironment MakeLegacyEnvironment(const LegacyEnvironmentInput& input,
                                        TextureCache& texture_cache,
                                        TransformCache& transform_cache);

} // namespace render
} // namespace rayrender

#endif
