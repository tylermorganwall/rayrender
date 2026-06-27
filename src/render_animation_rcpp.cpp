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
#include "math/matrix.h"
#include "math/transform.h"
#include "math/transformcache.h"
#include "materials/texturecache.h"
#include "utils/debug.h"
#include "core/bvh.h"
#include "core/PreviewDisplay.h"
#include "render/render_session.h"

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

// [[Rcpp::export]]
void render_animation_rcpp(List scene, List camera_info, List scene_info, List render_info,
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
  std::string background = as<std::string>(render_info["background"]);
  Float rotate_env = as<Float>(render_info["rotate_env"]);
  Float intensity_env = as<Float>(render_info["intensity_env"]);
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
  Float iso = as<Float>(camera_info["iso"]);
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
  int n_frames = end_frame;

  point3f backgroundhigh(bghigh[0],bghigh[1],bghigh[2]);
  point3f backgroundlow(bglow[0],bglow[1],bglow[2]);

  rayrender::render::RenderOptions render_options;
  render_options.nx = nx;
  render_options.ny = ny;
  render_options.samples = ns;
  render_options.numbercores = numbercores;
  render_options.progress_bar = progress_bar;
  render_options.verbose = verbose;
  render_options.debug_channel = debug_channel;
  render_options.min_variance = min_variance;
  render_options.min_adaptive_size = min_adaptive_size;
  render_options.integrator_type = integrator_type;

  rayrender::render::LegacyCameraStatic camera_static;
  camera_static.nx = nx;
  camera_static.ny = ny;
  camera_static.shutteropen = shutteropen;
  camera_static.shutterclose = shutterclose;
  camera_static.film_size = film_size;
  camera_static.camera_scale = camera_scale;
  camera_static.iso = iso;
  camera_static.lens_data =
    rayrender::render::FlattenLensData(realCameraInfo);
  camera_static.missing_lens_error =
    "No lense data passed in lens descriptor file.";

  rayrender::render::CompiledScene compiled_scene;
  rayrender::render::RenderSession render_session(render_options,
                                                  compiled_scene);

  RcppThread::ThreadPool pool(numbercores);
  GetRNGstate();
  random_gen rng(unif_rand() * std::pow(2,32));
  print_time(verbose, "Loaded Data");

  compiled_scene.Compile(scene,
                         shape,
                         shutteropen,
                         shutterclose,
                         bvh_type,
                         verbose,
                         rng);
  print_time(verbose, "Built Scene BVH" );
  

  //Calculate world bounds
  aabb bounding_box_world;
  compiled_scene.world_bvh->bounding_box(0,0,bounding_box_world);
  Float world_radius = bounding_box_world.Diag().length() ;
  vec3f world_center  = convert_to_vec3(bounding_box_world.Centroid());
  for(int i = 0; i < cam_x.length(); i++) {
    vec3f lf(cam_x(i),cam_y(i),cam_z(i));
    world_radius = world_radius > (lf - world_center).length() ? world_radius : 
      1.1*(lf - world_center).length();
  }

  rayrender::render::LegacyEnvironmentInput environment_input;
  environment_input.ambient_light = ambient_light;
  environment_input.hasbackground = hasbackground;
  environment_input.backgroundhigh = backgroundhigh;
  environment_input.backgroundlow = backgroundlow;
  environment_input.background = background;
  environment_input.rotate_env = rotate_env;
  environment_input.intensity_env = intensity_env;
  environment_input.world_radius = world_radius;
  environment_input.world_center = world_center;
  rayrender::render::LegacyEnvironment environment =
    rayrender::render::MakeLegacyEnvironment(
      environment_input,
      compiled_scene.texture_cache,
      compiled_scene.transform_cache);
  ambient_light = environment.ambient_light;
  hasbackground = environment.hasbackground;
  backgroundhigh = environment.backgroundhigh;
  backgroundlow = environment.backgroundlow;
  std::shared_ptr<hitable> background_sphere = environment.background_sphere;
  //Initialize background
  print_time(verbose, "Loaded Background" );
  compiled_scene.AssembleWorld(background_sphere);
  hitable_list& world = compiled_scene.world;
  hitable_list& imp_sample_objects =
    compiled_scene.importance_sample_objects;

  bool impl_only_bg = false;
  if((imp_sample_objects.size() == 0 || hasbackground || ambient_light) && debug_channel != 18) {
    impl_only_bg = true;
  }
  
  compiled_scene.AddEnvironmentToImportanceSampling(background_sphere,
                                                    impl_only_bg,
                                                    hasbackground);

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

      rayrender::render::CameraFrameState frame_state;
      frame_state.lookfrom = lookfrom;
      frame_state.lookat = lookat;
      frame_state.camera_up = camera_up;
      frame_state.fov = fov;
      frame_state.aperture = aperture;
      frame_state.focus_distance = focus_distance;
      frame_state.orthox = orthox;
      frame_state.orthoy = orthoy;
      rayrender::render::RenderFrameSetup frame_setup =
        render_session.RenderFrame(frame_state,
                                   camera_static,
                                   compiled_scene.transform_cache);
      std::unique_ptr<RayCamera>& cam = frame_setup.camera;

      // world_radius = world_radius > (lookfrom - world_center).length() ? world_radius : (lookfrom - world_center).length()*2;

      if(fov == 0) {
        Float ortho_diag = sqrt(pow(orthox,2) + pow(orthoy,2));
        world_radius += ortho_diag;
      }

      RayMatrix& rgb_output = frame_setup.outputs.rgb_output;
      RayMatrix& normalOutput = frame_setup.outputs.normal_output;
      RayMatrix& albedoOutput = frame_setup.outputs.albedo_output;

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
      post_process_frame(temp, debug_channel, as<std::string>(filenames(i)), 
                         as<std::string>(tonemap(0)), bloom,
                         transparent_background, write_image, plot_scene);
    }
  } else {
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

      rayrender::render::CameraFrameState frame_state;
      frame_state.lookfrom = lookfrom;
      frame_state.lookat = lookat;
      frame_state.camera_up = camera_up;
      frame_state.fov = fov;
      frame_state.aperture = aperture;
      frame_state.focus_distance = focus_distance;
      frame_state.orthox = orthox;
      frame_state.orthoy = orthoy;
      rayrender::render::RenderFrameSetup frame_setup =
        render_session.RenderFrame(frame_state,
                                   camera_static,
                                   compiled_scene.transform_cache);
      std::unique_ptr<RayCamera>& cam = frame_setup.camera;

      world_radius = world_radius > (lookfrom - world_center).length() ? world_radius : (lookfrom - world_center).length()*2;

      if(fov == 0) {
        Float ortho_diag = sqrt(pow(orthox,2) + pow(orthoy,2));
        world_radius += ortho_diag;
      }

      RayMatrix& rgb_output = frame_setup.outputs.rgb_output;
      RayMatrix& normalOutput = frame_setup.outputs.normal_output;
      RayMatrix& albedoOutput = frame_setup.outputs.albedo_output;
      RayMatrix& alpha_output = frame_setup.outputs.alpha_output;
      RayMatrix& draw_rgb_output = frame_setup.outputs.draw_rgb_output;

#ifdef HAS_OIDN
      // Create an Open Image Denoise device
      oidn::DeviceRef device = oidn::newDevice(); // CPU or GPU if available
      // oidn::DeviceRef device = oidn::newDevice(oidn::DeviceType::CPU);
      device.commit();
      // Create buffers for input/output images accessible by both host (CPU) and device (CPU/GPU)
      oidn::BufferRef colorBuf  = device.newBuffer(rgb_output.begin(), nx * ny * 3 * sizeof(Float));
      oidn::BufferRef albedoBuf = device.newBuffer(albedoOutput.begin(), nx * ny * 3 * sizeof(Float));
      oidn::BufferRef normalBuf = device.newBuffer(normalOutput.begin(), nx * ny * 3 * sizeof(Float));
      oidn::BufferRef colorBuf2 = device.newBuffer(draw_rgb_output.begin(), nx * ny * 3 * sizeof(Float));

      // Create a filter for denoising a beauty (color) image using optional auxiliary images too
      // This can be an expensive operation, so try no to create a new filter for every image!
      oidn::FilterRef filter = device.newFilter("RT"); // generic ray tracing filter
      filter.setImage("color",  colorBuf,  oidn::Format::Float3, nx, ny); // beauty
      filter.setImage("albedo", albedoBuf, oidn::Format::Float3, nx, ny); // auxiliary
      filter.setImage("normal", normalBuf, oidn::Format::Float3, nx, ny); // auxiliary
      filter.setImage("output", colorBuf2,  oidn::Format::Float3, nx, ny); // denoised beauty
      filter.set("hdr", true); // beauty image is HDR
      filter.commit();
#endif

#ifdef HAS_OIDN
      PreviewDisplay d(nx, ny, preview, false, 
                   false,
                   20.0f, cam.get(), 
                   background_sphere->ObjectToWorld,
                   background_sphere->WorldToObject,
                   filter, denoise, false);
#else
  PreviewDisplay d(nx,ny, preview, false, 
                         false,
                         20.0f, cam.get(),
                         background_sphere->ObjectToWorld,
                         background_sphere->WorldToObject,
                         false);
#endif
      pathtracer(numbercores, nx, ny, ns, debug_channel,
                 min_variance, min_adaptive_size,
                 rgb_output, normalOutput, albedoOutput,
                 alpha_output, draw_rgb_output,
                 progress_bar, sample_method, stratified_x, stratified_y,
                 verbose, cam.get(),  fov,
                 world, imp_sample_objects,
                 clampval, max_depth, roulette_active, d, integrator_type, rng_for_frame);
      if(d.terminate) {
        break;
      }
#ifdef HAS_OIDN
      filter.execute();
      const char* errorMessage;
      if (device.getError(errorMessage) != oidn::Error::None) {
        Rcpp::Rcout << "Error: " << errorMessage << std::endl;
      }
      List temp = List::create(_["r"] = draw_rgb_output.ConvertRcpp(0), 
                               _["g"] = draw_rgb_output.ConvertRcpp(1), 
                               _["b"] = draw_rgb_output.ConvertRcpp(2),
                               _["a"] = alpha_output.ConvertRcpp());
      post_process_frame(temp, debug_channel, as<std::string>(filenames(i)), as<std::string>(tonemap(0)), bloom,
                       transparent_background, write_image, plot_scene);
#else
      List temp = List::create(_["r"] = rgb_output.ConvertRcpp(0), 
                               _["g"] = rgb_output.ConvertRcpp(1), 
                               _["b"] = rgb_output.ConvertRcpp(2),
                               _["a"] = alpha_output.ConvertRcpp());
      post_process_frame(temp, debug_channel, as<std::string>(filenames(i)), as<std::string>(tonemap(0)), bloom,
                       transparent_background, write_image, plot_scene);
#endif
    }
  }

  PutRNGstate();
  print_time(verbose, "Finished rendering" );
  
}
