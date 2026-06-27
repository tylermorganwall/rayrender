#include "render_session.h"

#include "../hitables/infinite_area_light.h"
#include "../materials/constant.h"
#include "../math/transform.h"

#include <cmath>
#include <stdexcept>

namespace rayrender {
namespace render {

LegacyOutputBuffers::LegacyOutputBuffers(int nx, int ny)
  : rgb_output(nx, ny, 3),
    draw_rgb_output(nx, ny, 3),
    alpha_output(nx, ny, 1),
    normal_output(nx, ny, 3),
    albedo_output(nx, ny, 3) {}

void CompiledScene::Compile(Rcpp::List& scene,
                            Rcpp::IntegerVector& shape,
                            Float shutteropen,
                            Float shutterclose,
                            int bvh_type,
                            bool verbose,
                            random_gen& rng) {
  world_bvh = build_scene(scene,
                          shape,
                          shutteropen,
                          shutterclose,
                          textures,
                          alpha_textures,
                          bump_textures,
                          roughness_textures,
                          &shared_materials,
                          alpha,
                          bump,
                          roughness,
                          bvh_type,
                          transform_cache,
                          texture_cache,
                          importance_sample_objects,
                          instanced_objects,
                          instance_importance_sampled,
                          texture_idx,
                          verbose,
                          rng);
}

void CompiledScene::AssembleWorld(
  const std::shared_ptr<hitable>& background_sphere) {
  world.objects.clear();
  world.add(world_bvh);
  world.add(background_sphere);
}

void CompiledScene::AddEnvironmentToImportanceSampling(
  const std::shared_ptr<hitable>& background_sphere,
  bool importance_only_background,
  bool hasbackground) {
  if(importance_only_background || hasbackground) {
    importance_sample_objects.add(background_sphere);
  }
}

RenderSession::RenderSession(RenderOptions options, CompiledScene& scene)
  : options_(options),
    scene_(&scene) {}

RenderFrameSetup RenderSession::RenderFrame(
  const CameraFrameState& frame,
  const LegacyCameraStatic& camera_static,
  TransformCache& transform_cache) const {
  return RenderFrameSetup{
    CreateCamera(frame, camera_static, transform_cache),
    CreateOutputBuffers()
  };
}

std::unique_ptr<RayCamera> RenderSession::CreateCamera(
  const CameraFrameState& frame,
  const LegacyCameraStatic& camera_static,
  TransformCache& transform_cache) const {
  (void)scene_;
  return MakeLegacyCamera(frame, camera_static, transform_cache);
}

LegacyOutputBuffers RenderSession::CreateOutputBuffers() const {
  return LegacyOutputBuffers(options_.nx, options_.ny);
}

std::vector<Float> FlattenLensData(const Rcpp::NumericMatrix& real_camera_info) {
  std::vector<Float> lens_data;
  lens_data.reserve(static_cast<std::size_t>(real_camera_info.rows()) *
                    static_cast<std::size_t>(real_camera_info.cols()));
  for(int i = 0; i < real_camera_info.rows(); i++) {
    for(int j = 0; j < real_camera_info.cols(); j++) {
      lens_data.push_back(real_camera_info.at(i, j));
    }
  }
  return lens_data;
}

std::unique_ptr<RayCamera> MakeLegacyCamera(
  const CameraFrameState& frame,
  const LegacyCameraStatic& camera_static,
  TransformCache& transform_cache) {
  if(frame.fov < 0) {
    Transform cam_transform = LookAt(frame.lookfrom,
                                     frame.lookat,
                                     frame.camera_up).GetInverseMatrix();
    Transform* camera_transform = transform_cache.Lookup(cam_transform);
    AnimatedTransform animated_camera_transform(camera_transform,
                                                0,
                                                camera_transform,
                                                0);
    std::vector<Float> lens_data = camera_static.lens_data;
    if(lens_data.empty()) {
      throw std::runtime_error(camera_static.missing_lens_error);
    }
    return std::unique_ptr<RayCamera>(new RealisticCamera(
      animated_camera_transform,
      camera_static.shutteropen,
      camera_static.shutterclose,
      frame.aperture,
      camera_static.nx,
      camera_static.ny,
      frame.focus_distance,
      false,
      lens_data,
      camera_static.film_size,
      camera_static.camera_scale,
      camera_static.iso,
      frame.camera_up,
      cam_transform,
      frame.lookat));
  }
  if(frame.fov == 0) {
    return std::unique_ptr<RayCamera>(new ortho_camera(frame.lookfrom,
                                                       frame.lookat,
                                                       frame.camera_up,
                                                       frame.orthox,
                                                       frame.orthoy,
                                                       camera_static.shutteropen,
                                                       camera_static.shutterclose,
                                                       camera_static.iso));
  }
  if(frame.fov == 360) {
    return std::unique_ptr<RayCamera>(new environment_camera(
      frame.lookfrom,
      frame.lookat,
      frame.camera_up,
      camera_static.shutteropen,
      camera_static.shutterclose,
      camera_static.iso));
  }
  return std::unique_ptr<RayCamera>(new camera(frame.lookfrom,
                                               frame.lookat,
                                               frame.camera_up,
                                               frame.fov,
                                               Float(camera_static.nx) /
                                                 Float(camera_static.ny),
                                               frame.aperture,
                                               frame.focus_distance,
                                               camera_static.shutteropen,
                                               camera_static.shutterclose,
                                               camera_static.iso));
}

LegacyEnvironment MakeLegacyEnvironment(const LegacyEnvironmentInput& input,
                                        TextureCache& texture_cache,
                                        TransformCache& transform_cache) {
  LegacyEnvironment environment;
  environment.ambient_light = input.ambient_light;
  environment.hasbackground = input.hasbackground;
  environment.backgroundhigh = input.backgroundhigh;
  environment.backgroundlow = input.backgroundlow;

  Matrix4x4 identity;
  Transform background_angle(identity);
  if(input.rotate_env != 0) {
    background_angle = Translate(input.world_center) * RotateY(input.rotate_env);
  } else {
    background_angle = Translate(input.world_center);
  }
  Transform* background_transform = transform_cache.Lookup(background_angle);
  Transform* background_transform_inv =
    transform_cache.Lookup(background_angle.GetInverseMatrix());

  int nx1 = 0;
  int ny1 = 0;
  int nn1 = 0;
  Float* background_texture_data = nullptr;
  if(environment.hasbackground) {
    background_texture_data = texture_cache.LookupFloat(input.background,
                                                        nx1,
                                                        ny1,
                                                        nn1,
                                                        3);
    if(background_texture_data) {
      bool has_env_light = false;
      const std::size_t env_size = static_cast<std::size_t>(nx1) *
        static_cast<std::size_t>(ny1) *
        static_cast<std::size_t>(nn1);
      for(std::size_t i = 0; i < env_size; i++) {
        if(background_texture_data[i] > 0) {
          has_env_light = true;
          break;
        }
      }
      if(has_env_light) {
        environment.background_texture =
          std::make_shared<image_texture_float>(background_texture_data,
                                                nx1,
                                                ny1,
                                                nn1,
                                                1,
                                                1,
                                                input.intensity_env);
        environment.background_material =
          std::make_shared<diffuse_light>(environment.background_texture,
                                          1.0,
                                          false);
        environment.background_sphere =
          std::make_shared<InfiniteAreaLight>(
            nx1,
            ny1,
            input.world_radius * 2,
            convert_to_point3(input.world_center),
            environment.background_texture,
            environment.background_material,
            background_transform,
            background_transform_inv,
            false);
        return environment;
      }
      environment.hasbackground = false;
      environment.ambient_light = true;
    } else {
      Rcpp::Rcout << "Failed to load background image at " <<
        input.background << "\n";
      environment.hasbackground = false;
      environment.ambient_light = true;
    }
  }

  if(environment.hasbackground || environment.ambient_light) {
    if(!environment.hasbackground &&
       environment.backgroundhigh.length() == 0 &&
       environment.backgroundlow.length() == 0) {
      environment.backgroundhigh = point3f(FLT_MIN, FLT_MIN, FLT_MIN);
      environment.backgroundlow = point3f(FLT_MIN, FLT_MIN, FLT_MIN);
    }
    if(!environment.hasbackground && input.hasbackground) {
      environment.backgroundhigh = point3f(FLT_MIN, FLT_MIN, FLT_MIN);
      environment.backgroundlow = point3f(FLT_MIN, FLT_MIN, FLT_MIN);
    }
    environment.background_texture =
      std::make_shared<gradient_texture>(environment.backgroundlow,
                                         environment.backgroundhigh,
                                         false,
                                         false);
    environment.background_material =
      std::make_shared<diffuse_light>(environment.background_texture,
                                      1.0,
                                      false);
    environment.background_sphere =
      std::make_shared<InfiniteAreaLight>(
        100,
        100,
        input.world_radius * 2,
        convert_to_point3(input.world_center),
        environment.background_texture,
        environment.background_material,
        background_transform,
        background_transform_inv,
        false);
    return environment;
  }

  environment.background_texture =
    std::make_shared<constant_texture>(point3f(FLT_MIN, FLT_MIN, FLT_MIN));
  environment.background_material =
    std::make_shared<diffuse_light>(environment.background_texture,
                                    1.0,
                                    false);
  environment.background_sphere =
    std::make_shared<InfiniteAreaLight>(
      100,
      100,
      input.world_radius * 2,
      convert_to_point3(input.world_center),
      environment.background_texture,
      environment.background_material,
      background_transform,
      background_transform_inv,
      false);
  return environment;
}

} // namespace render
} // namespace rayrender
