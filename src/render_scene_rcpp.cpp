#ifndef STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION 
#endif

#include "float.h"
#include "vectypes.h"
#include "vec2.h"
#include "point2.h"
#include "mathinline.h"
#include "transform.h"
#include "transformcache.h"
#include "camera.h"
#include "float.h"
#include "buildscene.h"
#include "rng.h"
#include "tonemap.h"
#include "infinite_area_light.h"
#include "adaptivesampler.h"
#include "sampler.h"
#include "color.h"
#include "integrator.h"
#include "debug.h"

#include "texturecache.h"
#include "box.h"
#include "sphere.h"
#include "PreviewDisplay.h"
#include "raylog.h"
#include <cfenv>

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
#include "RayMatrix.h"
#include "RcppThread.h"

using namespace std;


// [[Rcpp::export]]
List render_scene_rcpp(List scene, List camera_info, List scene_info, List render_info) {
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
  std::string background = as<std::string>(render_info["background"]);
  Float rotate_env = as<Float>(render_info["rotate_env"]);
  Float intensity_env = as<Float>(render_info["intensity_env"]);
  bool verbose = as<bool>(render_info["verbose"]);
  int debug_channel = as<int>(render_info["debug_channel"]);
  Float min_variance = as<Float>(render_info["min_variance"]);
  int min_adaptive_size = as<int>(render_info["min_adaptive_size"]);
  IntegratorType integrator_type = static_cast<IntegratorType>(as<int>(render_info["integrator_type"]));
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
  Float focus_distance = as<Float>(camera_info["focal_distance"]);
  NumericVector ortho_dimensions = as<NumericVector>(camera_info["ortho_dimensions"]);
  std::size_t max_depth = as<std::size_t>(camera_info["max_depth"]);
  std::size_t roulette_active = as<std::size_t>(camera_info["roulette_active_depth"]);
  int sample_method = as<int>(camera_info["sample_method"]);
  NumericVector stratified_dim = as<NumericVector>(camera_info["stratified_dim"]);
  NumericVector light_direction = as<NumericVector>(camera_info["light_direction"]);
  NumericMatrix realCameraInfo = as<NumericMatrix>(camera_info["real_camera_info"]);
  Float film_size = as<Float>(camera_info["film_size"]);
  Float camera_scale = as<Float>(camera_info["camera_scale"]);
  Float sample_dist = as<Float>(camera_info["sample_dist"]);
  bool keep_colors = as<bool>(camera_info["keep_colors"]);
  bool preview     = as<bool>(camera_info["preview"]);
  bool interactive = as<bool>(camera_info["interactive"]);
  Float iso = as<Float>(camera_info["iso"]);
  int bvh_type = as<int>(camera_info["bvh"]);

  
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
    std::shared_ptr<Transform> CameraTransform = transformCache.Lookup(CamTransform);
    
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
                      shutteropen, shutterclose));
  } else if (fov == 360) {
    cam = std::unique_ptr<RayCamera>(new environment_camera(lookfrom, lookat, vec3f(camera_up(0),camera_up(1),camera_up(2)),
                            shutteropen, shutterclose));
  } else {
    cam = std::unique_ptr<RayCamera>(new camera(lookfrom, lookat, vec3f(camera_up(0),camera_up(1),camera_up(2)), fov, Float(nx)/Float(ny),
                                     aperture, dist_to_focus,
                                     shutteropen, shutterclose));
  }
  print_time(verbose, "Generated Camera" );


  int nx1, ny1, nn1;

  std::vector<Float* > textures;
  std::vector<unsigned char * > alpha_textures;
  std::vector<unsigned char * > bump_textures;
  std::vector<unsigned char * > roughness_textures;
  
  //Shared material vector
  std::vector<std::shared_ptr<material> >* shared_materials = new std::vector<std::shared_ptr<material> >;

  
  hitable_list imp_sample_objects;
  std::vector<std::shared_ptr<hitable> > instanced_objects;
  std::vector<std::shared_ptr<hitable_list> > instance_importance_sampled;
  std::vector<std::shared_ptr<alpha_texture> > alpha;
  std::vector<std::shared_ptr<bump_texture> > bump;
  std::vector<std::shared_ptr<roughness_texture> > roughness;
  std::vector<int> texture_idx;

  std::shared_ptr<hitable> worldbvh = build_scene(scene, 
                                                   shape, 
                                                   shutteropen,
                                                   shutterclose,
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
                                                   rng);
  print_time(verbose, "Built Scene BVH" );
  if(print_debug_info) {
    worldbvh->hitable_info_bounds(shutteropen,shutterclose);
  }
  //Calculate world bounds and ensure camera is inside infinite area light
  aabb bounding_box_world;
  worldbvh->bounding_box(0,0,bounding_box_world);
  Float world_radius = bounding_box_world.Diag().length() ;
  vec3f world_center  = convert_to_vec3(bounding_box_world.Centroid());
  world_radius = world_radius > (lookfrom - world_center).length() ? world_radius : 
     1.1*(lookfrom - world_center).length();
  world_radius *= interactive ? 100 : 1;
  
  if(fov == 0) {
    Float ortho_diag = sqrt(pow(ortho_dimensions(0),2) + pow(ortho_dimensions(1),2));
    world_radius += ortho_diag;
  }

  std::shared_ptr<texture> background_texture = nullptr;
  std::shared_ptr<material> background_material = nullptr;
  std::shared_ptr<hitable> background_sphere = nullptr;
  Float *background_texture_data = nullptr;

  //Background rotation
  Matrix4x4 Identity;
  Transform BackgroundAngle(Identity);
  if(rotate_env != 0) {
    BackgroundAngle = Translate(world_center) * RotateY(rotate_env);
  } else {
    BackgroundAngle = Translate(world_center);
  }

  std::shared_ptr<Transform> BackgroundTransform = transformCacheBg.Lookup(BackgroundAngle);
  std::shared_ptr<Transform> BackgroundTransformInv = transformCacheBg.Lookup(BackgroundAngle.GetInverseMatrix());
  if(hasbackground) {
    background_texture_data = texCache.LookupFloat(background, nx1, ny1, nn1, 3);
    // nn1 = 3;
    // texture_bytes += nx1 * ny1 * nn1;
    
    if(background_texture_data) {
      background_texture = std::make_shared<image_texture_float>(background_texture_data, nx1, ny1, nn1,
                                                                 1, 1, intensity_env);
      background_material = std::make_shared<diffuse_light>(background_texture, 1.0, false);
      background_sphere = std::make_shared<InfiniteAreaLight>(nx1, ny1, world_radius*2, convert_to_point3(world_center),
                                                              background_texture, background_material,
                                                              BackgroundTransform,
                                                              BackgroundTransformInv, false);
    } else {
      Rcpp::Rcout << "Failed to load background image at " << background << "\n";
      hasbackground = false;
      ambient_light = true;
      backgroundhigh = point3f(FLT_MIN,FLT_MIN,FLT_MIN);
      backgroundlow = point3f(FLT_MIN,FLT_MIN,FLT_MIN);
      background_texture = std::make_shared<gradient_texture>(backgroundlow, backgroundhigh, false, false);
      background_material = std::make_shared<diffuse_light>(background_texture, 1.0, false);
      background_sphere = std::make_shared<InfiniteAreaLight>(100, 100, world_radius*2, convert_to_point3(world_center),
                                                              background_texture, background_material,
                                                              BackgroundTransform,BackgroundTransformInv,false);
    }
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
  if((imp_sample_objects.size() == 0 || hasbackground || ambient_light || interactive) && debug_channel != 18) {
    impl_only_bg = true;
  }
  preview = preview && debug_channel == 0;
#ifdef HAS_OIDN
  PreviewDisplay Display(nx,ny, preview, interactive, 
                         (lookat-lookfrom).length(), cam.get(),
                         background_sphere->ObjectToWorld.get(),
                         background_sphere->WorldToObject.get(),
                         filter, denoise);
#else
  PreviewDisplay Display(nx,ny, preview, interactive, 
                         (lookat-lookfrom).length(), cam.get(),
                         background_sphere->ObjectToWorld.get(),
                         background_sphere->WorldToObject.get());
#endif
  
  if(impl_only_bg || hasbackground) {
    imp_sample_objects.add(background_sphere);
  }

  QUERY_MEMORY_USAGE();
  PRINT_CURRENT_MEMORY("Before raytracing");
  
  // Rcpp::Rcout << "Total world size: " << world.GetSize() + texture_bytes << " (Textures: " << texture_bytes << ") \n";
  if(debug_channel != 0) {
    debug_scene(numbercores, nx, ny, ns, debug_channel,
                min_variance, min_adaptive_size,
                rgb_output, normalOutput, albedoOutput,
                progress_bar, sample_method, stratified_dim,
                verbose, cam.get(), fov,
                world, imp_sample_objects, 
                clampval, max_depth, roulette_active,
                light_direction, rng, sample_dist, keep_colors,
                backgroundhigh);
  } else {
    pathtracer(numbercores, nx, ny, ns, debug_channel,
               min_variance, min_adaptive_size,
               rgb_output, normalOutput, albedoOutput,
               alpha_output,
               draw_rgb_output,
               progress_bar, sample_method, stratified_dim,
               verbose, cam.get(),  fov,
               world, imp_sample_objects,
               clampval, max_depth, roulette_active, Display, integrator_type);
  }
  PRINT_CURRENT_MEMORY("After raytracing");
#ifdef HAS_OIDN
  filter.execute();
  const char* errorMessage;
  if (device.getError(errorMessage) != oidn::Error::None) {
    Rcpp::Rcout << "Error: " << errorMessage << std::endl;
  }
#endif
  delete shared_materials;
  PutRNGstate();
  print_time(verbose, "Finished rendering" );
  RayMatrix &final_output = rgb_output;
#ifdef HAS_OIDN
  if(denoise) {
    final_output = draw_rgb_output;
  }
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

                                  _["a"] = alpha_output.ConvertRcpp());
  if(Display.Keyframes.size() > 0) {
    List keyframes(Display.Keyframes.size());
    for(unsigned int i = 0; i < Display.Keyframes.size(); i++ ) {
      keyframes(i) = Display.Keyframes[i];
    }
    final_image.attr("keyframes") = keyframes;
  }
  STOP_TIMER("Overall Time");

  PRINT_LOG_REPORT(numbercores);

  PRINT_CURRENT_MEMORY("After cleanup");
  
  return(final_image);
}

#include "testthat.h"

context("Test overall functionality") {
  #ifdef HAS_OIDN
  unsigned char denoised_final_render[432] =
{
	0x53, 0x5e, 0x4b, 0x3f, 0x1a, 0xeb, 0x60, 0x3f, 
	0x00, 0x00, 0x80, 0x3f, 0xbe, 0xf8, 0x48, 0x3f, 
	0xc1, 0x80, 0x5f, 0x3f, 0x00, 0x00, 0x80, 0x3f, 
	0xad, 0x8d, 0x4a, 0x3f, 0xe3, 0x6f, 0x60, 0x3f, 
	0x00, 0x00, 0x80, 0x3f, 0xa4, 0xa8, 0x4a, 0x3f, 
	0xd0, 0x7f, 0x60, 0x3f, 0x00, 0x00, 0x80, 0x3f, 
	0x0c, 0x4c, 0x4b, 0x3f, 0x50, 0xe0, 0x60, 0x3f, 
	0x00, 0x00, 0x80, 0x3f, 0x8a, 0xed, 0x49, 0x3f, 
	0x52, 0x11, 0x60, 0x3f, 0x00, 0x00, 0x80, 0x3f, 
	0x8e, 0x54, 0x47, 0x3f, 0x9c, 0x88, 0x5e, 0x3f, 
	0x00, 0x00, 0x80, 0x3f, 0x07, 0x99, 0x48, 0x3f, 
	0x3a, 0x48, 0x5f, 0x3f, 0x00, 0x00, 0x80, 0x3f, 
	0xe1, 0x9d, 0x81, 0x3f, 0x61, 0x50, 0x8d, 0x3f, 
	0x60, 0x2f, 0x9e, 0x3f, 0x60, 0xcf, 0xfd, 0x3e, 
	0x0b, 0x6f, 0x1c, 0x3f, 0xf6, 0x05, 0x47, 0x3f, 
	0xae, 0x15, 0x47, 0x3f, 0x7a, 0x63, 0x5e, 0x3f, 
	0x00, 0x00, 0x80, 0x3f, 0xb0, 0x55, 0x47, 0x3f, 
	0x48, 0x89, 0x5e, 0x3f, 0x00, 0x00, 0x80, 0x3f, 
	0x27, 0xc7, 0x44, 0x3f, 0xbe, 0x06, 0x5d, 0x3f, 
	0x00, 0x00, 0x80, 0x3f, 0x51, 0xb5, 0x45, 0x3f, 
	0x64, 0x93, 0x5d, 0x3f, 0x00, 0x00, 0x80, 0x3f, 
	0xa3, 0x82, 0xa1, 0x3f, 0x87, 0x1c, 0xad, 0x3f, 
	0x07, 0xd8, 0xbd, 0x3f, 0xb0, 0xfd, 0x78, 0x3f, 
	0x94, 0x23, 0x96, 0x3f, 0xf0, 0x1f, 0xbb, 0x3f, 
	0x4c, 0xe9, 0x7a, 0x3f, 0x40, 0x27, 0x8b, 0x3f, 
	0xd9, 0xe8, 0x9e, 0x3f, 0x0f, 0x98, 0x45, 0x3f, 
	0x1e, 0x82, 0x5d, 0x3f, 0x00, 0x00, 0x80, 0x3f, 
	0xcc, 0xda, 0x43, 0x3f, 0x2a, 0x7b, 0x5c, 0x3f, 
	0x00, 0x00, 0x80, 0x3f, 0xb6, 0x59, 0x42, 0x3f, 
	0xc0, 0x97, 0x5b, 0x3f, 0x00, 0x00, 0x80, 0x3f, 
	0x8c, 0x23, 0x86, 0x3f, 0x77, 0x47, 0xa1, 0x3f, 
	0x7f, 0x6c, 0xc8, 0x3f, 0x00, 0x00, 0x00, 0x00, 
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
	0x00, 0x00, 0x00, 0x00, 0xb6, 0x2c, 0x42, 0x3f, 
	0x2c, 0x7d, 0x5b, 0x3f, 0x00, 0x00, 0x80, 0x3f, 
	0xe8, 0xb2, 0x3f, 0x3f, 0xe2, 0x06, 0x5a, 0x3f, 
	0x00, 0x00, 0x80, 0x3f, 0x54, 0xe7, 0x40, 0x3f, 
	0x06, 0xbd, 0x5a, 0x3f, 0x00, 0x00, 0x80, 0x3f, 
	0x22, 0x63, 0x80, 0x3f, 0x89, 0x2d, 0xa1, 0x3f, 
	0xd9, 0x78, 0xd0, 0x3f, 0xb6, 0x92, 0x71, 0x3f, 
	0x63, 0x58, 0x9d, 0x3f, 0xff, 0x12, 0xd2, 0x3f, 
	0xec, 0x7b, 0x41, 0x3f, 0xc6, 0x14, 0x5b, 0x3f, 
	0x00, 0x00, 0x80, 0x3f, 0xc0, 0x6c, 0x40, 0x3f, 
	0xa2, 0x74, 0x5a, 0x3f, 0x00, 0x00, 0x80, 0x3f, 
	0x19, 0x16, 0x3f, 0x3f, 0x47, 0xaa, 0x59, 0x3f, 
	0x00, 0x00, 0x80, 0x3f, 0x28, 0x7b, 0x3d, 0x3f, 
	0x99, 0xb7, 0x58, 0x3f, 0x00, 0x00, 0x80, 0x3f, 
	0x3b, 0xd4, 0x3e, 0x3f, 0x62, 0x83, 0x59, 0x3f, 
	0x00, 0x00, 0x80, 0x3f, 0x20, 0xb0, 0x3d, 0x3f, 
	0xe0, 0xd6, 0x58, 0x3f, 0x00, 0x00, 0x80, 0x3f, 
	0xc6, 0xaf, 0x3e, 0x3f, 0xda, 0x6d, 0x59, 0x3f, 
	0x00, 0x00, 0x80, 0x3f, 0x18, 0xe7, 0x3d, 0x3f, 
	0x57, 0xf7, 0x58, 0x3f, 0x00, 0x00, 0x80, 0x3f, 
};

  

  test_that("OIDN Test") {
    vec3f v; v.e[0] = 1.0f; v.e[1] = 2.0f; v.e[2] = 3.0f;
    point3f p = convert_to_point3(v);

    expect_true(p.x() == Approx(1.0f));
    expect_true(p.y() == Approx(2.0f));
    expect_true(p.z() == Approx(3.0f));
  } 
  #endif
}