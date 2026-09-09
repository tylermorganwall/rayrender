#include "volumes/boundary.h"
#define RCPP_USE_UNWIND_PROTECT

#include "math/float.h"
#include "math/vectypes.h"
#include "math/mathinline.h"
#include "core/camera.h"
#include "math/float.h"
#include "core/buildscene.h"
#include "math/rng.h"
#include "hitables/infinite_area_light.h"
#include "core/adaptivesampler.h"
#include "math/sampler.h"
#include "core/color.h"
#include "core/integrator.h"
#include "core/oidn_aux.h"
#include "core/oidn_denoiser.h"
#include "math/matrix.h"
#include "math/transform.h"
#include "math/transformcache.h"
#include "materials/texturecache.h"
#include "utils/debug.h"
#include "core/bvh.h"
#include "core/PreviewDisplay.h"

#ifdef HAS_OIDN
#undef None
#include <OpenImageDenoise/oidn.hpp>
#endif

#include "Rcpp.h"
#include "math/RayMatrix.h"
using namespace Rcpp;
#include "RcppThread.h"
#include "RProgress.h"


using namespace std;

namespace {

int animation_camera_type(Float fov) {
  if(fov < 0) {
    return(-1);
  }
  if(fov == 0) {
    return(0);
  }
  if(fov == 360) {
    return(360);
  }
  return(1);
}

bool can_update_animation_camera(int camera_type) {
  return(camera_type == 0 || camera_type == 1 || camera_type == 360);
}

bool same_animation_up(const vec3f& lhs, const vec3f& rhs) {
  return(lhs.xyz.x == rhs.xyz.x &&
         lhs.xyz.y == rhs.xyz.y &&
         lhs.xyz.z == rhs.xyz.z);
}

std::unique_ptr<RayCamera> make_animation_camera(
    const point3f& lookfrom,
    const point3f& lookat,
    const vec3f& camera_up,
    Float fov,
    Float aperture,
    Float focus_distance,
    Float orthox,
    Float orthoy,
    int nx,
    int ny,
    Float shutteropen,
    Float shutterclose,
    NumericMatrix realCameraInfo,
    Float film_size,
    Float camera_scale,
    Float iso,
    Float shutter_speed,
    TransformCache& transformCache) {
  std::unique_ptr<RayCamera> cam;
  if(fov < 0) {
    Transform CamTransform = LookAt(lookfrom,
                                    lookat,
                                    camera_up).GetInverseMatrix();
    Transform* CameraTransform = transformCache.Lookup(CamTransform);
    AnimatedTransform CamTr(CameraTransform,0,CameraTransform,0);
    
    std::vector<Float> lensData;
    for(int lens_row = 0; lens_row < realCameraInfo.rows(); lens_row++) {
      for(int lens_col = 0; lens_col < realCameraInfo.cols(); lens_col++) {
        lensData.push_back(realCameraInfo.at(lens_row,lens_col));
      }
    }
    
    if(lensData.size() == 0) {
      throw std::runtime_error("No lense data passed in lens descriptor file.");
    }
    
    cam = std::unique_ptr<RayCamera>(new RealisticCamera(CamTr, shutteropen, shutterclose,
                                                         aperture, nx, ny, focus_distance, false,
                                                         lensData, film_size, camera_scale, iso,
                                                         camera_up, CamTransform, lookat));
  } else if(fov == 0) {
    cam = std::unique_ptr<RayCamera>(new ortho_camera(lookfrom, lookat, camera_up,
                                                      orthox, orthoy,
                                                      shutteropen, shutterclose, iso));
  } else if(fov == 360) {
    cam = std::unique_ptr<RayCamera>(new environment_camera(lookfrom, lookat, camera_up,
                                                            shutteropen, shutterclose, iso));
  } else {
    cam = std::unique_ptr<RayCamera>(new camera(lookfrom, lookat, camera_up, fov,
                                                Float(nx)/Float(ny), aperture, focus_distance,
                                                shutteropen, shutterclose, iso));
  }
  cam->set_shutter_speed(shutter_speed);
  return cam;
}

void update_animation_camera(
    RayCamera* cam,
    int camera_type,
    const point3f& lookfrom,
    const point3f& lookat,
    Float aperture,
    Float fov,
    Float focus_distance,
    Float orthox,
    Float orthoy) {
  cam->update_position_absolute(lookfrom);
  cam->update_lookat(lookat);
  if(camera_type == 1) {
    cam->update_aperture_absolute(aperture);
    cam->update_focal_absolute(focus_distance);
    cam->update_fov_absolute(fov);
  } else if(camera_type == 0) {
    cam->update_ortho_absolute(vec2f(orthox, orthoy));
  }
}

bool frame_camera_motion_blur_enabled(
    int frame,
    bool default_camera_motion_blur,
    bool has_camera_motion_blur_column,
    const LogicalVector& camera_motion_blur_values) {
  if(!has_camera_motion_blur_column ||
     frame >= camera_motion_blur_values.size()) {
    return default_camera_motion_blur;
  }
  int value = camera_motion_blur_values[frame];
  if(value == NA_LOGICAL) {
    return default_camera_motion_blur;
  }
  return value == TRUE;
}

int camera_motion_blur_end_frame(
    int frame,
    int frame_count,
    bool has_camera_motion_blur_group,
    const IntegerVector& camera_motion_blur_group,
    const NumericVector& cam_fov) {
  int end_frame = frame + 1;
  if(end_frame >= frame_count) {
    return frame;
  }
  if(animation_camera_type(cam_fov(frame)) !=
     animation_camera_type(cam_fov(end_frame))) {
    return frame;
  }
  if(has_camera_motion_blur_group &&
     frame < camera_motion_blur_group.size() &&
     end_frame < camera_motion_blur_group.size() &&
     camera_motion_blur_group(frame) != camera_motion_blur_group(end_frame)) {
    return frame;
  }
  return end_frame;
}

void set_animation_camera_motion_blur(
    RayCamera* cam,
    int frame,
    int end_frame,
    bool enabled,
    const NumericVector& cam_x,
    const NumericVector& cam_y,
    const NumericVector& cam_z,
    const NumericVector& cam_dx,
    const NumericVector& cam_dy,
    const NumericVector& cam_dz,
    const NumericVector& cam_upx,
    const NumericVector& cam_upy,
    const NumericVector& cam_upz,
    const NumericVector& cam_focal) {
  point3f start_origin(cam_x(frame), cam_y(frame), cam_z(frame));
  point3f start_lookat(cam_dx(frame), cam_dy(frame), cam_dz(frame));
  vec3f start_up(cam_upx(frame), cam_upy(frame), cam_upz(frame));
  point3f end_origin(cam_x(end_frame), cam_y(end_frame), cam_z(end_frame));
  point3f end_lookat(cam_dx(end_frame), cam_dy(end_frame), cam_dz(end_frame));
  vec3f end_up(cam_upx(end_frame), cam_upy(end_frame), cam_upz(end_frame));

  cam->set_camera_motion_blur(enabled);
  cam->set_camera_motion_blur_range(start_origin, start_lookat, start_up,
                                    cam_focal(frame),
                                    end_origin, end_lookat, end_up,
                                    cam_focal(end_frame));
}

} // namespace

// [[Rcpp::export]]
List render_animation_rcpp(List scene, List camera_info, List scene_info, List render_info,
                           List camera_movement, 
                           int start_frame, int end_frame,
                           CharacterVector filenames, Function post_process_frame, 
						   CharacterVector tonemap,
                           bool bloom, bool write_image, bool transparent_background) {
  //Unpack scene info
  IntegerVector shape = as<IntegerVector>(scene_info["shape"]);
  
  //Unpack render info
  bool ambient_light = as<bool>(render_info["ambient_light"]);
  NumericVector bghigh  = as<NumericVector>(render_info["bghigh"]);
  NumericVector bglow = as<NumericVector>(render_info["bglow"]);
  Float clampval = as<Float>(render_info["clampval"]);
  bool progress_bar = as<bool>(render_info["progress_bar"]);
  int numbercores = as<int>(render_info["numbercores"]);
  bool hasbackground = as<bool>(render_info["hasbackground"]);
  Float rotate_env = as<Float>(render_info["rotate_env"]);
  bool verbose = as<bool>(render_info["verbose"]);
  int debug_channel = as<int>(render_info["debug_channel"]);
  bool plot_scene = as<bool>(render_info["plot_scene"]);
  Float min_variance = as<Float>(render_info["min_variance"]);
  int min_adaptive_size = as<int>(render_info["min_adaptive_size"]);
  IntegratorType integrator_type = static_cast<IntegratorType>(as<int>(render_info["integrator_type"]));
#ifdef HAS_OIDN
  bool denoise = as<bool>(render_info["denoise"]);
#endif
  bool has_frame_seed = render_info.containsElementNamed("frame_seed");
  unsigned int frame_seed = 0;
  if(has_frame_seed) {
    frame_seed = static_cast<unsigned int>(as<int>(render_info["frame_seed"]));
  }

  Environment pkg = Environment::namespace_env("rayrender");
  Function print_time = pkg["print_time"];


  //Unpack Camera Info
  int nx = as<int>(camera_info["nx"]);
  int ny = as<int>(camera_info["ny"]);
  int ns = as<int>(camera_info["ns"]);
  Float shutteropen = as<Float>(camera_info["shutteropen"]);
  Float shutterclose = as<Float>(camera_info["shutterclose"]);
  Float shutter_speed = camera_info.containsElementNamed("shutter_speed") ?
    as<Float>(camera_info["shutter_speed"]) :
    static_cast<Float>(2);
  std::size_t max_depth = as<std::size_t>(camera_info["max_depth"]);
  std::size_t roulette_active = as<std::size_t>(camera_info["roulette_active_depth"]);
  int sample_method = as<int>(camera_info["sample_method"]);
  NumericVector stratified_dim = as<NumericVector>(camera_info["stratified_dim"]);
  NumericVector light_direction = as<NumericVector>(camera_info["light_direction"]);
  int stratified_x = static_cast<int>(stratified_dim(0));
  int stratified_y = static_cast<int>(stratified_dim(1));
  vec3f preview_light_direction(light_direction(0), light_direction(1), light_direction(2));
  Float preview_exponent = light_direction.size() > 3 ? static_cast<Float>(light_direction(3)) : 0;
  int bvh_type = as<int>(camera_info["bvh"]);
  NumericMatrix realCameraInfo = as<NumericMatrix>(camera_info["real_camera_info"]);
  Float film_size = as<Float>(camera_info["film_size"]);
  Float camera_scale = as<Float>(camera_info["camera_scale"]);
  Float sample_dist = as<Float>(camera_info["sample_dist"]);
  bool keep_colors = as<bool>(camera_info["keep_colors"]);
  bool preview = as<bool>(camera_info["preview"]);
  bool camera_motion_blur = camera_info.containsElementNamed("camera_motion_blur") ?
    as<bool>(camera_info["camera_motion_blur"]) :
    false;
  Float iso = as<Float>(camera_info["iso"]);
  
  std::unique_ptr<RayCamera> cam;

  //unpack motion info
  NumericVector cam_x        = as<NumericVector>(camera_movement["x"]);
  NumericVector cam_y        = as<NumericVector>(camera_movement["y"]);
  NumericVector cam_z        = as<NumericVector>(camera_movement["z"]);
  NumericVector cam_dx       = as<NumericVector>(camera_movement["dx"]);
  NumericVector cam_dy       = as<NumericVector>(camera_movement["dy"]);
  NumericVector cam_dz       = as<NumericVector>(camera_movement["dz"]);
  NumericVector cam_upx      = as<NumericVector>(camera_movement["upx"]);
  NumericVector cam_upy      = as<NumericVector>(camera_movement["upy"]);
  NumericVector cam_upz      = as<NumericVector>(camera_movement["upz"]);
  NumericVector cam_aperture = as<NumericVector>(camera_movement["aperture"]);
  NumericVector cam_fov      = as<NumericVector>(camera_movement["fov"]);
  NumericVector cam_focal    = as<NumericVector>(camera_movement["focal"]);
  NumericVector cam_orthox   = as<NumericVector>(camera_movement["orthox"]);
  NumericVector cam_orthoy   = as<NumericVector>(camera_movement["orthoy"]);
  bool has_camera_motion_blur_column = camera_movement.containsElementNamed("camera_motion_blur");
  LogicalVector cam_motion_blur =
    has_camera_motion_blur_column ?
    as<LogicalVector>(camera_movement["camera_motion_blur"]) :
    LogicalVector();
  bool has_camera_motion_blur_group = camera_movement.containsElementNamed("camera_motion_blur_group");
  IntegerVector cam_motion_blur_group =
    has_camera_motion_blur_group ?
    as<IntegerVector>(camera_movement["camera_motion_blur_group"]) :
    IntegerVector();
  int motion_frame_count = cam_x.size();
  int n_frames = end_frame;
  List output_frames(n_frames - start_frame);
  int output_frame_index = 0;

  point3f backgroundhigh(bghigh[0],bghigh[1],bghigh[2]);
  point3f backgroundlow(bglow[0],bglow[1],bglow[2]);

  RcppThread::ThreadPool pool(numbercores);
  GetRNGstate();
  random_gen rng(unif_rand() * std::pow(2,32));
  print_time(verbose, "Loaded Data");

  std::vector<Float* > textures;
  std::vector<unsigned char * > alpha_textures;
  std::vector<unsigned char * > bump_textures;
  std::vector<unsigned char * > roughness_textures;
  //Shared material vector
  std::vector<std::shared_ptr<material> >* shared_materials = new std::vector<std::shared_ptr<material> >;


  //Initialize transformation cache
  TransformCache transformCache;
  TransformCache transformCacheBg;
  
  //Initialize texture cache
  TextureCache texCache;
  
  hitable_list imp_sample_objects;
  if(integrator_type == IntegratorType::ShadowRays) {
    imp_sample_objects.volume_scene = std::make_shared<VolumeScene>();
    imp_sample_objects.volume_scene->transparent_background = render_info.containsElementNamed("transparent_background") && as<bool>(render_info["transparent_background"]);
  }
  std::vector<std::shared_ptr<hitable> > instanced_objects;
  std::vector<std::shared_ptr<hitable_list> > instance_importance_sampled;
  std::vector<std::shared_ptr<alpha_texture> > alpha;
  std::vector<std::shared_ptr<bump_texture> > bump;
  std::vector<std::shared_ptr<roughness_texture> > roughness;
  std::vector<int> texture_idx;

  std::shared_ptr<hitable> worldbvh = build_scene(scene, shape, 
                                                  static_cast<Float>(0), static_cast<Float>(1),
                                                  textures, 
                                                  alpha_textures,
                                                  bump_textures,
                                                  roughness_textures, 
                                                  shared_materials, 
                                                  alpha, bump, roughness,
                                                  bvh_type,
                                                  transformCache, 
                                                  texCache,
                                                  imp_sample_objects,
                                                  instanced_objects,
                                                  instance_importance_sampled,
                                                  texture_idx,
                                                  verbose, rng);
  bool has_atmosphere = render_info.containsElementNamed("has_atmosphere") &&
                        Rcpp::as<bool>(render_info["has_atmosphere"]);
  if (has_atmosphere && imp_sample_objects.volume_scene)
    imp_sample_objects.volume_scene->has_media = true;
  bool has_media = imp_sample_objects.volume_scene && imp_sample_objects.volume_scene->has_media;
  print_time(verbose, "Built Scene BVH" );
  

  //Calculate world bounds
  aabb bounding_box_world;
  worldbvh->bounding_box(0,0,bounding_box_world);
  Float world_radius = bounding_box_world.Diag().length() ;
  vec3f world_center  = convert_to_vec3(bounding_box_world.Centroid());
  for(int i = 0; i < cam_x.length(); i++) {
    vec3f lf(cam_x(i),cam_y(i),cam_z(i));
    world_radius = world_radius > (lf - world_center).length() ? world_radius : 
      1.1*(lf - world_center).length();
  }

  
  std::shared_ptr<texture> background_texture = nullptr;
  std::shared_ptr<material> background_material = nullptr;
  std::shared_ptr<hitable> background_sphere = nullptr;
  Matrix4x4 Identity;
  Transform BackgroundAngle(Identity);
  if(rotate_env != 0) {
    BackgroundAngle = Translate(world_center) * RotateY(rotate_env);
  } else {
    BackgroundAngle = Translate(world_center);
  }
  Transform* BackgroundTransform = transformCacheBg.Lookup(BackgroundAngle);
  Transform* BackgroundTransformInv = transformCacheBg.Lookup(BackgroundAngle.GetInverseMatrix());

  if(hasbackground) {
    auto infinite_lights = BuildInfiniteLights(
        Rcpp::as<Rcpp::List>(render_info["infinite_lights"]), texCache);
    if (imp_sample_objects.volume_scene)
      imp_sample_objects.volume_scene->atmosphere = infinite_lights->GetAtmosphere();
    background_sphere = std::make_shared<InfiniteAreaLight>(
        infinite_lights, world_radius * 2, convert_to_point3(world_center),
        BackgroundTransform, BackgroundTransformInv);
  } else if(ambient_light) {
    //Check if both high and low are black, and set to FLT_MIN
    if(backgroundhigh.length() == 0 && backgroundlow.length() == 0) {
      backgroundhigh = point3f(FLT_MIN,FLT_MIN,FLT_MIN);
      backgroundlow = point3f(FLT_MIN,FLT_MIN,FLT_MIN);
    }
    background_texture = std::make_shared<gradient_texture>(backgroundlow, backgroundhigh, false, false);
    background_material = std::make_shared<diffuse_light>(background_texture, 1.0, false);
    background_sphere = std::make_shared<InfiniteAreaLight>(100, 100, world_radius*2, convert_to_point3(world_center),
                                                            background_texture, background_material,
                                                            BackgroundTransform, BackgroundTransformInv, false);
    
  } else {
    //Minimum intensity FLT_MIN so the CDF isn't NAN
    background_texture = std::make_shared<constant_texture>(point3f(FLT_MIN,FLT_MIN,FLT_MIN));
    background_material = std::make_shared<diffuse_light>(background_texture, 1.0, false);
    background_sphere = std::make_shared<InfiniteAreaLight>(100, 100, world_radius*2, convert_to_point3(world_center),
                                                            background_texture, background_material,
                                                            BackgroundTransform,
                                                            BackgroundTransformInv, false);
  }
  //Initialize background
  print_time(verbose, "Loaded Background" );
  hitable_list world;
  world.add(worldbvh);

  bool impl_only_bg = false;
  world.add(background_sphere);
  if(((imp_sample_objects.size() == 0 && !(imp_sample_objects.volume_scene && imp_sample_objects.volume_scene->has_emission)) || hasbackground || ambient_light) && debug_channel != 18) {
    impl_only_bg = true;
  }
  
  if(impl_only_bg || hasbackground) {
    imp_sample_objects.add(background_sphere);
  }

  if(verbose && !progress_bar) {
    Rcpp::message(CharacterVector("Starting Raytracing"));
  }
  RProgress::RProgress pb_sampler("Generating Samples [:bar] :percent%");
  pb_sampler.set_width(70);
  RProgress::RProgress pb("Adaptive Raytracing [:bar] :percent%");
  pb.set_width(70);
  RProgress::RProgress pb_frames("Frame :current/:total [:bar] :percent%");

  pb_frames.set_width(70);

  if(progress_bar) {
    pb_sampler.set_total(ny);
    pb.set_total(ns);
    pb_frames.set_total(n_frames - start_frame);
  }

  if(debug_channel != 0) {
    for(int i = start_frame; i < n_frames; i++ ) {
      random_gen frame_rng(frame_seed);
      random_gen* rng_for_frame = has_frame_seed ? &frame_rng : &rng;
      if(progress_bar) {
        pb_frames.tick();
      }
      point3f lookfrom(cam_x(i),cam_y(i),cam_z(i));
      point3f lookat(cam_dx(i),cam_dy(i),cam_dz(i));
      Float fov = cam_fov(i);
      Float aperture = cam_aperture(i);
      Float focus_distance = cam_focal(i);
      Float orthox = cam_orthox(i);
      Float orthoy = cam_orthoy(i);
      vec3f camera_up = vec3f(cam_upx(i),cam_upy(i),cam_upz(i));


      cam = make_animation_camera(lookfrom, lookat, camera_up, fov, aperture,
                                  focus_distance, orthox, orthoy, nx, ny,
                                  shutteropen, shutterclose, realCameraInfo,
                                  film_size, camera_scale, iso, shutter_speed,
                                  transformCache);
      bool blur_enabled = frame_camera_motion_blur_enabled(
        i,
        camera_motion_blur,
        has_camera_motion_blur_column,
        cam_motion_blur
      );
      int blur_end_frame = camera_motion_blur_end_frame(
        i,
        motion_frame_count,
        has_camera_motion_blur_group,
        cam_motion_blur_group,
        cam_fov
      );
      set_animation_camera_motion_blur(cam.get(), i, blur_end_frame, blur_enabled,
                                       cam_x, cam_y, cam_z,
                                       cam_dx, cam_dy, cam_dz,
                                       cam_upx, cam_upy, cam_upz,
                                       cam_focal);

      // world_radius = world_radius > (lookfrom - world_center).length() ? world_radius : (lookfrom - world_center).length()*2;

      if(fov == 0) {
        Float ortho_diag = sqrt(pow(orthox,2) + pow(orthoy,2));
        world_radius += ortho_diag;
      }

      //Initialize output matrices
      RayMatrix rgb_output(nx,ny, 3);
      RayMatrix normalOutput(nx,ny, 3);
      RayMatrix albedoOutput(nx,ny, 3);

      debug_scene(numbercores, nx, ny, ns, debug_channel,
                  min_variance, min_adaptive_size,
                  rgb_output, normalOutput, albedoOutput,
                  progress_bar, sample_method, stratified_x, stratified_y,
                  verbose, cam.get(), fov,
                  world, imp_sample_objects, 
                  clampval, max_depth, roulette_active,
                  preview_light_direction, preview_exponent, *rng_for_frame, sample_dist, keep_colors, backgroundhigh);
      List temp = List::create(_["r"] = rgb_output.ConvertRcpp(0), 
                               _["g"] = rgb_output.ConvertRcpp(1), 
                               _["b"] = rgb_output.ConvertRcpp(2));
      if(imp_sample_objects.volume_scene && imp_sample_objects.volume_scene->collect_statistics)
        temp.attr("volume_statistics")=imp_sample_objects.volume_scene->Statistics();
      std::string frame_filename = as<std::string>(filenames(i));
      bool write_current_image = write_image && !frame_filename.empty();
      RObject frame_output = post_process_frame(temp, debug_channel, frame_filename, 
                         as<std::string>(tonemap(0)), bloom,
                         transparent_background, write_current_image, plot_scene);
      output_frames[output_frame_index++] = frame_output;
    }
  } else {
    bool persistent_preview = preview;
    std::unique_ptr<RayCamera> preview_cam;
    bool preview_cam_initialized = false;
    int preview_camera_type = 0;
    vec3f preview_camera_up;
    std::unique_ptr<PreviewDisplay> preview_display;

    for(int i = start_frame; i < n_frames; i++ ) {
      random_gen frame_rng(frame_seed);
      random_gen* rng_for_frame = has_frame_seed ? &frame_rng : nullptr;
      if(progress_bar) {
        pb_frames.tick();
      }
      point3f lookfrom(cam_x(i),cam_y(i),cam_z(i));
      point3f lookat(cam_dx(i),cam_dy(i),cam_dz(i));
      vec3f camera_up = vec3f(cam_upx(i),cam_upy(i),cam_upz(i));

      Float fov = cam_fov(i);
      Float aperture = cam_aperture(i);
      Float focus_distance = cam_focal(i);
      Float orthox = cam_orthox(i);
      Float orthoy = cam_orthoy(i);

      int current_camera_type = animation_camera_type(fov);
      RayCamera* frame_cam = nullptr;
      std::unique_ptr<RayCamera> frame_cam_storage;
      if(persistent_preview) {
        if(preview_cam_initialized &&
           current_camera_type == preview_camera_type &&
           can_update_animation_camera(current_camera_type) &&
           same_animation_up(camera_up, preview_camera_up)) {
          update_animation_camera(preview_cam.get(), current_camera_type,
                                  lookfrom, lookat, aperture, fov,
                                  focus_distance, orthox, orthoy);
        } else {
          preview_cam = make_animation_camera(lookfrom, lookat, camera_up, fov, aperture,
                                              focus_distance, orthox, orthoy, nx, ny,
                                              shutteropen, shutterclose, realCameraInfo,
                                              film_size, camera_scale, iso, shutter_speed,
                                              transformCache);
          preview_cam_initialized = true;
          preview_camera_type = current_camera_type;
          preview_camera_up = camera_up;
        }
        frame_cam = preview_cam.get();
      } else {
        frame_cam_storage = make_animation_camera(lookfrom, lookat, camera_up, fov, aperture,
                                                 focus_distance, orthox, orthoy, nx, ny,
                                                 shutteropen, shutterclose, realCameraInfo,
                                                 film_size, camera_scale, iso, shutter_speed,
                                                 transformCache);
        frame_cam = frame_cam_storage.get();
      }
      bool blur_enabled = frame_camera_motion_blur_enabled(
        i,
        camera_motion_blur,
        has_camera_motion_blur_column,
        cam_motion_blur
      );
      int blur_end_frame = camera_motion_blur_end_frame(
        i,
        motion_frame_count,
        has_camera_motion_blur_group,
        cam_motion_blur_group,
        cam_fov
      );
      set_animation_camera_motion_blur(frame_cam, i, blur_end_frame, blur_enabled,
                                       cam_x, cam_y, cam_z,
                                       cam_dx, cam_dy, cam_dz,
                                       cam_upx, cam_upy, cam_upz,
                                       cam_focal);

      world_radius = world_radius > (lookfrom - world_center).length() ? world_radius : (lookfrom - world_center).length()*2;

      if(fov == 0) {
        Float ortho_diag = sqrt(pow(orthox,2) + pow(orthoy,2));
        world_radius += ortho_diag;
      }

      //Initialize output matrices
      RayMatrix rgb_output(nx,ny, 3);
      RayMatrix normalOutput(nx,ny, 3);
      RayMatrix albedoOutput(nx,ny, 3);
      RayMatrix alpha_output(nx,ny, 1);
      RayMatrix draw_rgb_output(nx,ny, 3);

#ifdef HAS_OIDN
      bool denoise_frame = denoise;
      RayMatrix oidn_normal_output(nx, ny, 3);
      RayMatrix oidn_albedo_output(nx, ny, 3);
      RayOidnDenoiser oidn_denoiser;
      if(denoise_frame) {
        oidn_denoiser.Setup(rgb_output,
                            oidn_albedo_output,
                            oidn_normal_output,
                            draw_rgb_output,
                            nx,
                            ny,
                            RayOidnQuality::Balanced,
                            false,
                            false, !has_media);
      }
#endif

      bool terminated = false;
#ifdef HAS_OIDN
      if(persistent_preview) {
        if(!preview_display) {
          preview_display = std::unique_ptr<PreviewDisplay>(new PreviewDisplay(
            nx, ny, preview, false, false, 20.0f, frame_cam,
            background_sphere->ObjectToWorld,
            background_sphere->WorldToObject,
            &oidn_denoiser,
            &oidn_albedo_output,
            &oidn_normal_output,
            denoise_frame, false));
        } else {
          preview_display->SetCamera(frame_cam);
          preview_display->SetDenoiser(&oidn_denoiser,
                                       &oidn_albedo_output,
                                       &oidn_normal_output,
                                       denoise_frame);
        }
        pathtracer(numbercores, nx, ny, ns, debug_channel,
                   min_variance, min_adaptive_size,
                   rgb_output, normalOutput, albedoOutput,
                   alpha_output, draw_rgb_output,
                   progress_bar, sample_method, stratified_x, stratified_y,
                   verbose, frame_cam, fov,
                   world, imp_sample_objects,
                   clampval, max_depth, roulette_active, *preview_display,
                   integrator_type, rng_for_frame);
        terminated = preview_display->terminate;
      } else {
        PreviewDisplay d(nx, ny, preview, false, 
                     false,
                     20.0f, frame_cam, 
                     background_sphere->ObjectToWorld,
                     background_sphere->WorldToObject,
                     &oidn_denoiser,
                     &oidn_albedo_output,
                     &oidn_normal_output,
                     denoise_frame, false);
        pathtracer(numbercores, nx, ny, ns, debug_channel,
                   min_variance, min_adaptive_size,
                   rgb_output, normalOutput, albedoOutput,
                   alpha_output, draw_rgb_output,
                   progress_bar, sample_method, stratified_x, stratified_y,
                   verbose, frame_cam, fov,
                   world, imp_sample_objects,
                   clampval, max_depth, roulette_active, d, integrator_type, rng_for_frame);
        terminated = d.terminate;
      }
#else
      if(persistent_preview) {
        if(!preview_display) {
          preview_display = std::unique_ptr<PreviewDisplay>(new PreviewDisplay(
            nx, ny, preview, false, false, 20.0f, frame_cam,
            background_sphere->ObjectToWorld,
            background_sphere->WorldToObject,
            false));
        } else {
          preview_display->SetCamera(frame_cam);
        }
        pathtracer(numbercores, nx, ny, ns, debug_channel,
                   min_variance, min_adaptive_size,
                   rgb_output, normalOutput, albedoOutput,
                   alpha_output, draw_rgb_output,
                   progress_bar, sample_method, stratified_x, stratified_y,
                   verbose, frame_cam, fov,
                   world, imp_sample_objects,
                   clampval, max_depth, roulette_active, *preview_display,
                   integrator_type, rng_for_frame);
        terminated = preview_display->terminate;
      } else {
        PreviewDisplay d(nx,ny, preview, false, 
                         false,
                         20.0f, frame_cam,
                         background_sphere->ObjectToWorld,
                         background_sphere->WorldToObject,
                         false);
        pathtracer(numbercores, nx, ny, ns, debug_channel,
                   min_variance, min_adaptive_size,
                   rgb_output, normalOutput, albedoOutput,
                   alpha_output, draw_rgb_output,
                   progress_bar, sample_method, stratified_x, stratified_y,
                   verbose, frame_cam,  fov,
                   world, imp_sample_objects,
                   clampval, max_depth, roulette_active, d, integrator_type, rng_for_frame);
        terminated = d.terminate;
      }
#endif
      if(terminated) {
        break;
      }
#ifdef HAS_OIDN
      if(denoise_frame) {
        OidnAuxRenderOptions oidn_aux_options;
        oidn_aux_options.samples = static_cast<std::size_t>(std::max(ns, 1));
        oidn_aux_options.max_depth = max_depth;
        oidn_aux_options.max_dielectric_splits = 1;
        oidn_aux_options.sample_method = sample_method;
        oidn_aux_options.stratified_x = stratified_x;
        oidn_aux_options.stratified_y = stratified_y;
        if(!has_media) render_oidn_aux_features(numbercores,
                                 nx,
                                 ny,
                                 frame_cam,
                                 fov,
                                 &world,
                                 oidn_aux_options,
                                 oidn_normal_output,
                                 oidn_albedo_output);
        oidn_denoiser.Setup(rgb_output,
                            oidn_albedo_output,
                            oidn_normal_output,
                            draw_rgb_output,
                            nx,
                            ny,
                            RayOidnQuality::High,
                            true,
                            true, !has_media);
        oidn_denoiser.Execute();
        oidn_denoiser.ReportError();
      }
      RayMatrix final_output = denoise_frame ? draw_rgb_output : rgb_output;
      List temp = List::create(_["r"] = final_output.ConvertRcpp(0), 
                               _["g"] = final_output.ConvertRcpp(1), 
                               _["b"] = final_output.ConvertRcpp(2),
                               _["a"] = alpha_output.ConvertRcpp(),
                               _["premultiplied"] = integrator_type == IntegratorType::ShadowRays);
      if(imp_sample_objects.volume_scene && imp_sample_objects.volume_scene->collect_statistics)
        temp.attr("volume_statistics")=imp_sample_objects.volume_scene->Statistics();
      std::string frame_filename = as<std::string>(filenames(i));
      bool write_current_image = write_image && !frame_filename.empty();
      RObject frame_output = post_process_frame(temp, debug_channel, frame_filename, as<std::string>(tonemap(0)), bloom,
                       transparent_background, write_current_image, plot_scene);
      output_frames[output_frame_index++] = frame_output;
#else
      List temp = List::create(_["r"] = rgb_output.ConvertRcpp(0), 
                               _["g"] = rgb_output.ConvertRcpp(1), 
                               _["b"] = rgb_output.ConvertRcpp(2),
                               _["a"] = alpha_output.ConvertRcpp(),
                               _["premultiplied"] = integrator_type == IntegratorType::ShadowRays);
      if(imp_sample_objects.volume_scene && imp_sample_objects.volume_scene->collect_statistics)
        temp.attr("volume_statistics")=imp_sample_objects.volume_scene->Statistics();
      std::string frame_filename = as<std::string>(filenames(i));
      bool write_current_image = write_image && !frame_filename.empty();
      RObject frame_output = post_process_frame(temp, debug_channel, frame_filename, as<std::string>(tonemap(0)), bloom,
                       transparent_background, write_current_image, plot_scene);
      output_frames[output_frame_index++] = frame_output;
#endif
    }
  }

  delete shared_materials;
  PutRNGstate();
  print_time(verbose, "Finished rendering" );
  if(output_frame_index < output_frames.size()) {
    List trimmed_output_frames(output_frame_index);
    for(int i = 0; i < output_frame_index; i++) {
      trimmed_output_frames[i] = output_frames[i];
    }
    return trimmed_output_frames;
  }
  return output_frames;
}
