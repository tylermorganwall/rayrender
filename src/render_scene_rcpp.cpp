#include "core/wavefront.h"
#include "animation_rcpp.h"
#include "core/wavefront_report.h"
#include "core/animation_report.h"
#include <chrono>
#include "volumes/boundary.h"

#include "math/float.h"
#include "math/vectypes.h"
#include "math/mathinline.h"
#include "math/transform.h"
#include "math/transformcache.h"
#include "core/camera.h"
#include "math/float.h"
#include "core/buildscene.h"
#include "core/bvh_timing.h"
#include "rng.h"
#include "hitables/infinite_area_light.h"
#include "core/adaptivesampler.h"
#include "math/sampler.h"
#include "core/color.h"
#include "core/integrator.h"
#include "core/oidn_aux.h"
#include "core/oidn_denoiser.h"
#include "utils/debug.h"

#include "materials/texturecache.h"
#include "hitables/box.h"
#include "hitables/sphere.h"
#include "core/PreviewDisplay.h"
#include "utils/raylog.h"
#include <cfenv>
#include <array>

// #define DEBUG

#ifdef DEBUG
#include <iostream>
#endif

#ifdef HAS_OIDN
#undef None
#include <OpenImageDenoise/oidn.hpp>
#endif


#include "RProgress.h"
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
#include "math/RayMatrix.h"
#include "RcppThread.h"

using namespace std;

#include "screen_overlays.h"
using namespace rayrender_screen;

// [[Rcpp::export]]
List render_scene_rcpp(List scene, List camera_info, List scene_info, List render_info) {
  const bool benchmark_timing = render_info.containsElementNamed("benchmark_timing") &&
    as<bool>(render_info["benchmark_timing"]);
  BVHBuildTiming bvh_timing(benchmark_timing);
  RESET_RAYLOG();
  START_TIMER("Overall Time");
  feclearexcept(FE_ALL_EXCEPT);

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
  Float min_variance = as<Float>(render_info["min_variance"]);
  int min_adaptive_size = as<int>(render_info["min_adaptive_size"]);
  IntegratorType integrator_type = static_cast<IntegratorType>(as<int>(render_info["integrator_type"]));
  WavefrontReport wavefront_report;
  wavefront_report.requested = render_info.containsElementNamed("metal_wavefront") &&
      as<bool>(render_info["metal_wavefront"]);
  bool print_debug_info = as<bool>(render_info["print_debug_info"]);
#ifdef HAS_OIDN
  bool denoise = as<bool>(render_info["denoise"]);
#endif

  Environment pkg = Environment::namespace_env("rayrender");
  Function print_time = pkg["print_time"];


  //Unpack Camera Info
  int nx = as<int>(camera_info["nx"]);
  int ny = as<int>(camera_info["ny"]);
  int ns = as<int>(camera_info["ns"]);
  Float fov = as<Float>(camera_info["fov"]);
  NumericVector lookfromvec = as<NumericVector>(camera_info["lookfrom"]);
  NumericVector lookatvec = as<NumericVector>(camera_info["lookat"]);
  Float aperture = as<Float>(camera_info["aperture"]);
  NumericVector camera_up = as<NumericVector>(camera_info["camera_up"]);
  Float shutteropen = as<Float>(camera_info["shutteropen"]);
  Float shutterclose = as<Float>(camera_info["shutterclose"]);
  Float shutter_speed = camera_info.containsElementNamed("shutter_speed") ?
    as<Float>(camera_info["shutter_speed"]) :
    static_cast<Float>(2);
  Float focus_distance = as<Float>(camera_info["focal_distance"]);
  NumericVector ortho_dimensions = as<NumericVector>(camera_info["ortho_dimensions"]);
  std::size_t max_depth = as<std::size_t>(camera_info["max_depth"]);
  std::size_t roulette_active = as<std::size_t>(camera_info["roulette_active_depth"]);
  int sample_method = as<int>(camera_info["sample_method"]);
  NumericVector stratified_dim = as<NumericVector>(camera_info["stratified_dim"]);
  NumericVector light_direction = as<NumericVector>(camera_info["light_direction"]);
  int stratified_x = static_cast<int>(stratified_dim(0));
  int stratified_y = static_cast<int>(stratified_dim(1));
  vec3f preview_light_direction(light_direction(0), light_direction(1), light_direction(2));
  Float preview_exponent = light_direction.size() > 3 ? static_cast<Float>(light_direction(3)) : 0;
  NumericMatrix realCameraInfo = as<NumericMatrix>(camera_info["real_camera_info"]);
  Float film_size = as<Float>(camera_info["film_size"]);
  Float camera_scale = as<Float>(camera_info["camera_scale"]);
  Float sample_dist = as<Float>(camera_info["sample_dist"]);
  bool keep_colors = as<bool>(camera_info["keep_colors"]);
  bool preview     = as<bool>(camera_info["preview"]);
  bool interactive = as<bool>(camera_info["interactive"]);
  bool deferred_render = as<bool>(camera_info["deferred_render"]);
  bool auto_exposure = as<bool>(camera_info["auto_exposure"]);
  bool camera_motion_blur = camera_info.containsElementNamed("camera_motion_blur") ?
    as<bool>(camera_info["camera_motion_blur"]) :
    false;
  List keyframe_motion_args = camera_info.containsElementNamed("keyframe_motion_args") ?
    as<List>(camera_info["keyframe_motion_args"]) :
    List::create(_["type"] = "spline",
                 _["smooth_orientation"] = true,
                 _["damp_motion"] = true);
  std::string snapshot_filename;
  if(camera_info.containsElementNamed("snapshot_filename")) {
    CharacterVector snapshot_filename_value =
      as<CharacterVector>(camera_info["snapshot_filename"]);
    if(snapshot_filename_value.size() > 0 &&
       !CharacterVector::is_na(snapshot_filename_value[0])) {
      snapshot_filename = as<std::string>(snapshot_filename_value[0]);
    }
  }
  Float iso = as<Float>(camera_info["iso"]);
  BVHBuildOptions bvh_type{static_cast<BVHBuildMethod>(as<int>(camera_info["bvh"])),
                           static_cast<unsigned>(std::max(1, numbercores))};

  
  //Initialize transformation cache
  TransformCache transformCache;
  TransformCache transformCacheBg;
  
  //Initialize texture cache
  TextureCache texCache;

  //Initialize output matrices
  RayMatrix rgb_output(nx,ny, 3);
  RayMatrix draw_rgb_output(nx,ny, 3);
  RayMatrix alpha_output(nx,ny, 1);
  RayMatrix normalOutput(nx,ny, 3);
  RayMatrix albedoOutput(nx,ny, 3);

#ifdef HAS_OIDN
  RayMatrix oidn_normal_output(nx, ny, 3);
  RayMatrix oidn_albedo_output(nx, ny, 3);
  RayOidnDenoiser oidn_denoiser;

#endif
  
  point3f lookfrom(lookfromvec[0],lookfromvec[1],lookfromvec[2]);
  point3f lookat(lookatvec[0],lookatvec[1],lookatvec[2]);
  point3f backgroundhigh(bghigh[0],bghigh[1],bghigh[2]);
  point3f backgroundlow(bglow[0],bglow[1],bglow[2]);
  Float dist_to_focus = focus_distance;
  
  std::vector<bool> has_image;
  std::vector<bool> has_alpha;
  std::vector<bool> has_bump;
  std::vector<bool> has_roughness;
  
  std::unique_ptr<RayCamera> cam;


  RcppThread::ThreadPool pool(numbercores);
  GetRNGstate();
  random_gen rng(unif_rand() * std::pow(2,32));
  if(fov < 0) {
    Transform CamTransform = LookAt(lookfrom,
                                    lookat,
                                    vec3f(camera_up(0),camera_up(1),camera_up(2))).GetInverseMatrix();
    Transform* CameraTransform = transformCache.Lookup(CamTransform);
    
    AnimatedTransform CamTr(CameraTransform,0,CameraTransform,0);
    
    std::vector<Float> lensData;
    for(int i = 0; i < realCameraInfo.rows(); i++) {
      for(int j = 0; j < realCameraInfo.cols(); j++) {
        lensData.push_back(realCameraInfo.at(i,j));
      }
    }
    
    if(fov < 0 && lensData.size() == 0) {
      throw std::runtime_error("No lens data passed in lens descriptor file.");
    }
    
    cam = std::unique_ptr<RayCamera>(new RealisticCamera(CamTr,shutteropen, shutterclose,
                         aperture, nx,ny, focus_distance, false, lensData,
                         film_size, camera_scale, iso, vec3f(camera_up(0),camera_up(1),camera_up(2)),
                         CamTransform, lookat));
  } else if(fov == 0) {
    cam = std::unique_ptr<RayCamera>(new ortho_camera(lookfrom, lookat, vec3f(camera_up(0),camera_up(1),camera_up(2)),
                      ortho_dimensions(0), ortho_dimensions(1),
                      shutteropen, shutterclose, iso));
  } else if (fov == 360) {
    cam = std::unique_ptr<RayCamera>(new environment_camera(lookfrom, lookat, vec3f(camera_up(0),camera_up(1),camera_up(2)),
                            shutteropen, shutterclose, iso));
  } else {
    cam = std::unique_ptr<RayCamera>(new camera(lookfrom, lookat, vec3f(camera_up(0),camera_up(1),camera_up(2)), fov, Float(nx)/Float(ny),
                                     aperture, dist_to_focus,
                                     shutteropen, shutterclose, iso));
  }
  cam->set_camera_motion_blur(camera_motion_blur);
  cam->set_shutter_speed(shutter_speed);
  cam->set_camera_flip_x(camera_info.containsElementNamed("camera_flip_x") &&
                         as<bool>(camera_info["camera_flip_x"]));
  if (camera_motion_blur && camera_info.containsElementNamed("motion_endpoint")) {
    List endpoint = camera_info["motion_endpoint"];
    NumericVector origin = endpoint["lookfrom"], target = endpoint["lookat"], up = endpoint["up"];
    cam->set_camera_motion_blur_range(lookfrom, lookat,
        vec3f(camera_up(0), camera_up(1), camera_up(2)), dist_to_focus,
        point3f(origin[0], origin[1], origin[2]), point3f(target[0], target[1], target[2]),
        vec3f(up[0], up[1], up[2]), as<Float>(endpoint["focal"]));
  }
  print_time(verbose, "Generated Camera" );



  std::vector<Float* > textures;
  std::vector<unsigned char * > alpha_textures;
  std::vector<std::shared_ptr<const HeightImage>> bump_textures;
  std::vector<unsigned char * > roughness_textures;
  
  //Shared material vector
  std::vector<std::shared_ptr<material> >* shared_materials = new std::vector<std::shared_ptr<material> >;

  
  hitable_list imp_sample_objects;
  if(integrator_type == IntegratorType::ShadowRays) {
    imp_sample_objects.volume_scene = std::make_shared<VolumeScene>();
    if (render_info.containsElementNamed("point_lights"))
      imp_sample_objects.volume_scene->point_lights = PointLightSet(as<List>(render_info["point_lights"]));
    imp_sample_objects.volume_scene->transparent_background = render_info.containsElementNamed("transparent_background") && as<bool>(render_info["transparent_background"]);
  }
  std::vector<std::shared_ptr<hitable> > instanced_objects;
  std::vector<std::shared_ptr<hitable_list> > instance_importance_sampled;
  std::vector<std::shared_ptr<alpha_texture> > alpha;
  std::vector<std::shared_ptr<bump_texture> > bump;
  std::vector<std::shared_ptr<roughness_texture> > roughness;
  std::vector<int> texture_idx;

  const bool has_timeline = camera_info.containsElementNamed("timeline") && !Rf_isNull(camera_info["timeline"]);
  const double time_origin = camera_info.containsElementNamed("scene_time") && !Rf_isNull(camera_info["scene_time"]) ?
    as<double>(camera_info["scene_time"]) : 0;
  const bool explicit_time = camera_info.containsElementNamed("scene_time") && !Rf_isNull(camera_info["scene_time"]);
  auto scene_time = explicit_time ? std::make_shared<SceneTimeContext>() : nullptr;
  if (scene_time) {
    scene_time->epoch = time_origin;
    scene_time->benchmark_timing = benchmark_timing;
    scene_time->full_pose_blur = render_info.containsElementNamed("animation_motion_blur") &&
      as<std::string>(render_info["animation_motion_blur"]) == "full";
    scene_time->frozen_pose_notice = [] {
      Rcpp::message(Rcpp::CharacterVector(
        "Internal model poses are sampled at frame time; placement and other scene motion retain motion blur."));
    };
  }
  cam->scene_camera_time = scene_time;
  std::shared_ptr<const animation::CameraTrack> camera_track;
  if (camera_info.containsElementNamed("timeline_camera") && !Rf_isNull(camera_info["timeline_camera"])) {
    DataFrame motion = as<DataFrame>(camera_info["timeline_camera"]);
    if (motion.nrows() > 1 && (has_timeline || shutteropen != 0 || shutterclose != 0)) {
      camera_track = read_camera_track(motion, 1, 0);
      camera_track->ValidateExposure(shutteropen,shutterclose);
      cam->scene_camera_track = camera_track;
      cam->scene_camera_time = scene_time;
      cam->scene_camera_aspect = Float(nx)/Float(ny);
    }
  }
  NumericVector timeline_range = has_timeline ? as<NumericVector>(camera_info["timeline"]) : NumericVector();
  Float bounds_open = shutteropen, bounds_close = shutterclose;
  if (has_timeline) {
    scene_time->epoch = time_origin;
    // Include the zero-width dragging exposure as well as the authored shutter.
    bounds_open = std::nextafter(Float(timeline_range[0] - time_origin + std::min(Float(0), shutteropen)), -std::numeric_limits<Float>::infinity());
    bounds_close = std::nextafter(Float(timeline_range[1] - time_origin + std::max(Float(0), shutterclose)), std::numeric_limits<Float>::infinity());
  }
  const auto compile_started = benchmark_timing ? std::chrono::steady_clock::now() :
    std::chrono::steady_clock::time_point{};
  std::shared_ptr<hitable> worldbvh = build_scene(scene, 
                                                   shape, 
                                                   bounds_open,
                                                   bounds_close,
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
                                                   verbose, 
                                                   rng, nullptr,
                                                   time_origin, scene_time);
  if (scene_time && benchmark_timing) scene_time->compile_seconds =
    std::chrono::duration<double>(std::chrono::steady_clock::now()-compile_started).count();
  if (scene_time) scene_time->UpdateGeometry(time_origin, 0, shutteropen, shutterclose);
  bool has_atmosphere = render_info.containsElementNamed("has_atmosphere") &&
                        Rcpp::as<bool>(render_info["has_atmosphere"]);
  if (has_atmosphere && imp_sample_objects.volume_scene)
    imp_sample_objects.volume_scene->has_media = true;
#ifdef HAS_OIDN
  bool has_media = imp_sample_objects.volume_scene && imp_sample_objects.volume_scene->has_media;
  if(denoise) {
    oidn_denoiser.Setup(rgb_output,
                        oidn_albedo_output,
                        oidn_normal_output,
                        draw_rgb_output,
                        nx,
                        ny,
                        RayOidnQuality::Balanced,
                        false,
                        false, true);
  }
#endif
  print_time(verbose, "Built Scene BVH" );
  if(print_debug_info) {
    worldbvh->hitable_info_bounds(static_cast<Float>(0), static_cast<Float>(1));
  }
  //Calculate world bounds and ensure camera is inside infinite area light
  aabb bounding_box_world;
  worldbvh->bounding_box(0,0,bounding_box_world);
  Float world_radius = bounding_box_world.Diag().length() ;
  vec3f world_center  = convert_to_vec3(bounding_box_world.Centroid());
  world_radius = world_radius > (lookfrom - world_center).length() ? world_radius : 
     1.1*(lookfrom - world_center).length();
  world_radius *= interactive ? 100 : 1;
  
  if(fov == 0 || (interactive && fov > 0 && fov < 180)) {
    // The preview may switch to orthographic rays without rebuilding the sky.
    Float ortho_diag = sqrt(pow(ortho_dimensions(0),2) + pow(ortho_dimensions(1),2));
    world_radius += ortho_diag;
  }

  std::shared_ptr<texture> background_texture = nullptr;
  std::shared_ptr<material> background_material = nullptr;
  std::shared_ptr<hitable> background_sphere = nullptr;

  //Background rotation
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
      imp_sample_objects.volume_scene->atmosphere = infinite_lights->GetTransportAtmosphere();
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
  print_time(verbose, "Loaded background" );
  hitable_list world;
  world.add(worldbvh);

  bool impl_only_bg = false;
  world.add(background_sphere);
  if(((imp_sample_objects.size() == 0 && !(imp_sample_objects.volume_scene && imp_sample_objects.volume_scene->has_emission)) || hasbackground || ambient_light || interactive) && debug_channel != 18) {
    impl_only_bg = true;
  }
  LogicalVector screen_text_visible = compute_screen_text_visibility(
    render_info,
    cam.get(),
    worldbvh.get(),
    rng
  );
  LogicalVector screen_line_visible = compute_screen_line_visibility(
    render_info,
    cam.get(),
    worldbvh.get(),
    rng
  );
  std::vector<PreviewTextOverlay> text_overlays =
    parse_preview_text_overlays(render_info);
  std::vector<PreviewLineOverlay> line_overlays =
    parse_preview_line_overlays(render_info);
  preview = preview && debug_channel == 0;
#ifdef HAS_OIDN
  PreviewDisplay Display(nx,ny, preview, interactive, 
                         deferred_render, cam.get(),
                         background_sphere->ObjectToWorld,
                         background_sphere->WorldToObject,
                         &oidn_denoiser,
                         &oidn_albedo_output,
                         &oidn_normal_output,
                         denoise, auto_exposure);
#else
  PreviewDisplay Display(nx,ny, preview, interactive, 
                         deferred_render, cam.get(),
                         background_sphere->ObjectToWorld,
                         background_sphere->WorldToObject,
                         auto_exposure);
#endif
  Display.navigation_bounds = bounding_box_world;
  Display.ConfigureInteractiveCameras(&cam, vec2f(ortho_dimensions(0), ortho_dimensions(1)));
  if (has_timeline) {
    Display.ConfigureTimeline(timeline_range[0], timeline_range[1], time_origin,
      as<double>(camera_info["timeline_fps"]), scene_time, shutteropen, shutterclose);
    if (camera_track) Display.ConfigureTimelineCamera(camera_track, &cam);
  }
  // Keep complete, immutable light variants for interactive atmosphere changes.
  // Rebuilding the mixture also updates altitude-dependent sampling weights and
  // celestial filtering. Coefficients and image textures use their existing caches.
  if(interactive && has_atmosphere && hasbackground) {
    auto environment_light = std::static_pointer_cast<InfiniteAreaLight>(background_sphere);
    Rcpp::List descriptions = Rcpp::as<Rcpp::List>(render_info["infinite_lights"]);
    for(R_xlen_t i = 0; i < descriptions.size(); ++i) {
      Rcpp::List description = descriptions[i];
      if(Rcpp::as<std::string>(description["type"]) != "prague") continue;
      const bool haze = !description.containsElementNamed("haze") || Rcpp::as<bool>(description["haze"]);
      const bool query_altitude = !description.containsElementNamed("query_altitude") ||
                                  Rcpp::as<bool>(description["query_altitude"]);
      auto variants = std::make_shared<std::array<std::shared_ptr<InfiniteLight>, 4>>();
      (*variants)[2 * query_altitude + haze] = environment_light->light;
      auto volume_scene = imp_sample_objects.volume_scene;
      Display.SetAtmosphereControls(haze, query_altitude,
        [descriptions, i, variants, environment_light, volume_scene, &texCache](bool haze, bool altitude) {
          auto &light = (*variants)[2 * altitude + haze];
          if(!light) {
            Rcpp::List updated = Rcpp::clone(descriptions);
            Rcpp::List atmosphere = updated[i];
            atmosphere["haze"] = haze;
            atmosphere["query_altitude"] = altitude;
            light = BuildInfiniteLights(updated, texCache);
          }
          environment_light->SetLight(light);
          volume_scene->atmosphere = light->GetTransportAtmosphere();
        });
      break;
    }
  }
  Display.SetSnapshotFilename(snapshot_filename);
  Display.SetKeyframeMotionArgs(keyframe_motion_args);
  // Restore this camera's authored keys, without consulting R's last-render
  // keyframe buffer or moving away from the selected initial camera.
  if(camera_info.containsElementNamed("keyframes")) {
    Rcpp::List saved_keys = camera_info["keyframes"];
    for(R_xlen_t i = 0; i < saved_keys.size(); ++i) {
      Display.Keyframes.push_back(Rcpp::clone(Rcpp::as<Rcpp::List>(saved_keys[i])));
    }
  }
  Display.SetTextOverlays(text_overlays);
  Display.SetLineOverlays(line_overlays);
  
  if(impl_only_bg || hasbackground) {
    imp_sample_objects.add(background_sphere);
  }

  QUERY_MEMORY_USAGE();
  PRINT_CURRENT_MEMORY("Before raytracing");
  
  const auto trace_start = benchmark_timing ? std::chrono::steady_clock::now() :
                                            std::chrono::steady_clock::time_point{};
  if(debug_channel != 0) {
    if (scene_time && scene_time->full_pose_blur && scene_time->has_model_poses &&
        shutteropen != shutterclose)
      Rcpp::message(Rcpp::CharacterVector(
        "Debug channels retain frame-frozen model poses; full-pose shutter integration is available for color renders."));
    debug_scene(numbercores, nx, ny, ns, debug_channel,
                min_variance, min_adaptive_size,
                rgb_output, normalOutput, albedoOutput,
                progress_bar, sample_method, stratified_x, stratified_y,
                verbose, cam.get(), fov,
                world, imp_sample_objects, 
                clampval, max_depth, roulette_active,
                preview_light_direction, preview_exponent, rng, sample_dist, keep_colors,
                backgroundhigh);
  } else {
    pathtracer(numbercores, nx, ny, ns, debug_channel,
               min_variance, min_adaptive_size,
               rgb_output, normalOutput, albedoOutput,
               alpha_output,
               draw_rgb_output,
               progress_bar, sample_method, stratified_x, stratified_y,
               verbose, cam.get(),  fov,
               world, imp_sample_objects,
               clampval, max_depth, roulette_active, Display, integrator_type, nullptr, &wavefront_report);
    fov = cam->get_fov();
  }
  const double trace_seconds = benchmark_timing ?
    std::chrono::duration<double>(std::chrono::steady_clock::now() - trace_start).count() : 0;
  PRINT_CURRENT_MEMORY("After raytracing");
#ifdef HAS_OIDN
  Display.PollCloseEvent();
  bool denoise_output = denoise && !Display.terminate;
  bool use_cached_denoised_preview = denoise &&
    Display.terminate &&
    Display.HasDenoisedPreview();
  if(denoise_output) {
    OidnAuxRenderOptions oidn_aux_options;
    oidn_aux_options.samples = static_cast<std::size_t>(std::max(ns, 1));
    oidn_aux_options.max_depth = max_depth;
    oidn_aux_options.max_dielectric_splits = 1;
    oidn_aux_options.sample_method = sample_method;
    oidn_aux_options.stratified_x = stratified_x;
    oidn_aux_options.stratified_y = stratified_y;
    if (has_media || wavefront_report.used || (scene_time && scene_time->shutter_samples > 0)) {
      std::copy(normalOutput.begin(), normalOutput.end(), oidn_normal_output.begin());
      std::copy(albedoOutput.begin(), albedoOutput.end(), oidn_albedo_output.begin());
    }
    if(!has_media && !wavefront_report.used && !(scene_time && scene_time->shutter_samples > 0)) render_oidn_aux_features(numbercores,
                             nx,
                             ny,
                             cam.get(),
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
                        true, true);
    oidn_denoiser.Execute();
    oidn_denoiser.ReportError();
  }
#endif
  delete shared_materials;
  PutRNGstate();
  print_time(verbose, "Finished rendering" );
  RayMatrix final_output = rgb_output;
#ifdef HAS_OIDN
  if(denoise_output || use_cached_denoised_preview) {
    final_output = draw_rgb_output;
  }
#endif
  #ifdef RAY_COLOR_DEBUG
  final_output.print();
  #endif
  List final_image = List::create(_["r"] = final_output.ConvertRcpp(0), 
                                  _["g"] = final_output.ConvertRcpp(1), 
                                  _["b"] = final_output.ConvertRcpp(2),

                                  _["nx"] = normalOutput.ConvertRcpp(0), 
                                  _["ny"] = normalOutput.ConvertRcpp(1), 
                                  _["nz"] = normalOutput.ConvertRcpp(2),

                                  _["cx"] = albedoOutput.ConvertRcpp(0), 
                                  _["cy"] = albedoOutput.ConvertRcpp(1), 
                                  _["cz"] = albedoOutput.ConvertRcpp(2),

                                  _["a"] = alpha_output.ConvertRcpp(),
                                  _["premultiplied"] = integrator_type == IntegratorType::ShadowRays);
  if(imp_sample_objects.volume_scene && imp_sample_objects.volume_scene->collect_statistics)
    final_image.attr("volume_statistics")=imp_sample_objects.volume_scene->Statistics();
  if (imp_sample_objects.volume_scene)
    final_image.attr("path_warnings") = imp_sample_objects.volume_scene->path_diagnostics.Take();
  if(Display.Keyframes.size() > 0) {
    List keyframes(Display.Keyframes.size());
    for(unsigned int i = 0; i < Display.Keyframes.size(); i++ ) {
      keyframes(i) = Display.Keyframes[i];
    }
    final_image.attr("keyframes") = keyframes;
  }
  if (benchmark_timing) {
    final_image.attr("bvh_build_seconds") = bvh_timing.count ? bvh_timing.seconds : NA_REAL;
    final_image.attr("bvh_build_count") = bvh_timing.count;
    final_image.attr("trace_seconds") = trace_seconds;
    if (scene_time) final_image.attr("animation_timing") = AnimationTimingToR(*scene_time);
  }
  if(wavefront_report.requested) {
    final_image.attr("wavefront") = WavefrontReportToR(wavefront_report);
  }
  final_image.attr("render_cancelled") = Display.terminate;
  final_image.attr("preview_exposure") = Display.preview_exposure_adjustment;
  final_image.attr("screen_camera_info") = get_screen_camera_info(cam.get());
  if (Display.timeline) {
    const auto &published = Display.timeline->Published();
    const auto &displayed = Display.timeline->Displayed();
    final_image.attr("timeline") = List::create(Named("published_time") = published.time,
      Named("displayed_time") = displayed.time, Named("revision") = double(published.revision),
      Named("displayed_revision") = double(displayed.revision), Named("rate") = Display.timeline->Rate());
  }
  if(screen_text_visible.size() > 0) {
    final_image.attr("screen_text_visible") = screen_text_visible;
  }
  if(screen_line_visible.size() > 0) {
    final_image.attr("screen_line_visible") = screen_line_visible;
  }
  if(should_return_screen_text_overlay(render_info)) {
    final_image.attr("screen_text_overlay") = composite_screen_text_overlay(
      text_overlays,
      static_cast<unsigned int>(nx),
      static_cast<unsigned int>(ny),
      cam.get(),
      worldbvh.get(),
      rng
    );
  }
  if(should_return_screen_line_overlay(render_info)) {
    final_image.attr("screen_line_overlay") = composite_screen_line_overlay(
      line_overlays,
      static_cast<unsigned int>(nx),
      static_cast<unsigned int>(ny),
      cam.get(),
      worldbvh.get(),
      rng
    );
  }
  STOP_TIMER("Overall Time");

  PRINT_LOG_REPORT(numbercores);

  PRINT_CURRENT_MEMORY("After cleanup");
  
  texCache.ReportPtex(verbose);
  return(final_image);
}
