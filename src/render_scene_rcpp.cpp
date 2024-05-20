#ifndef STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION 
#endif

#include <Rcpp.h>

// #define DEBUG_MEMORY
#ifdef DEBUG_MEMORY

#include <cstdlib>
#include <iostream>
#include <new>


static std::size_t alloc{0};
static std::size_t dealloc{0};
static std::size_t memory_used{0};
static std::size_t mem_counter(0);

void* operator new(std::size_t sz){
  alloc += 1;
  memory_used += sz;
  Rcpp::Rcout << alloc << " : " << sz << " : " << memory_used << "\n";
  return std::malloc(sz);
}

void operator delete(void* ptr) noexcept{
  dealloc += 1;
  std::free(ptr);
}

void* operator new[](std::size_t sz){
  alloc += 1;
  memory_used += sz;
  // Rcpp::Rcout << alloc << " : " << sz << " : " << memory_used << "\n";
  
  return std::malloc(sz);
}

void operator delete[](void* ptr) noexcept {
  dealloc += 1;
  std::free(ptr);
}

void getInfo(){
  Rcpp::Rcout << "Number of allocations  : " << alloc << std::endl;
  // Rcpp::Rcout << "Number of deallocations: " << dealloc << std::endl;
  Rcpp::Rcout << "Total memory used      : " << memory_used <<      " bytes (";

  if(memory_used < 1024) {
  } else if (memory_used <= 1024*1024) {
    Rcpp::Rcout << (float)memory_used/(float)1024 << " KB)" << std::endl;
  } else if (memory_used <= 1024*1024*1024) {
    Rcpp::Rcout << (float)memory_used/(float)(1024*1024) <<      " MB)" << std::endl;
  } else {
    Rcpp::Rcout << (float)memory_used/float(1024*1024*1024) <<      " GB)" << std::endl;
  }
}
#endif

#include "float.h"
#include "vec3.h"
#include "vec2.h"
#include "point3.h"
#include "point2.h"
#include "normal.h" 
#include "RayMatrix.h"
#include "mathinline.h"
#include "transform.h"
#include "transformcache.h"
#include "camera.h"
#include "float.h"
#include "buildscene.h"
#include "RProgress.h"
#include "rng.h"
#include "tonemap.h"
#include "infinite_area_light.h"
#include "adaptivesampler.h"
#include "sampler.h"
#include "color.h"
#include "integrator.h"
#include "debug.h"
#include "bvh_node.h"
#include "texturecache.h"
#include "box.h"
#include "sphere.h"
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
#include "RcppThread.h"
#include "PreviewDisplay.h"

// #define DEBUG

#ifdef DEBUG
#include <iostream>
#endif

using namespace std;



// [[Rcpp::export]]
List render_scene_rcpp(List scene, List camera_info, List scene_info, List render_info) {
#ifdef DEBUG_MEMORY
  alloc = 0;
  // dealloc = 0;
  memory_used = 0;
  mem_counter = 0;
#endif

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
  RayMatrix routput(nx,ny);
  RayMatrix goutput(nx,ny);
  RayMatrix boutput(nx,ny);
  RayMatrix alpha_output(nx,ny);
  
  vec3f lookfrom(lookfromvec[0],lookfromvec[1],lookfromvec[2]);
  vec3f lookat(lookatvec[0],lookatvec[1],lookatvec[2]);
  vec3f backgroundhigh(bghigh[0],bghigh[1],bghigh[2]);
  vec3f backgroundlow(bglow[0],bglow[1],bglow[2]);
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
  
  std::shared_ptr<bvh_node> worldbvh = build_scene(scene, 
                                                   shape, 
                                                   shutteropen,
                                                   shutterclose,
                                                   textures, 
                                                   alpha_textures,
                                                   bump_textures,
                                                   roughness_textures, 
                                                   shared_materials, 
                                                   bvh_type,
                                                   transformCache, 
                                                   texCache,
                                                   imp_sample_objects,
                                                   instanced_objects,
                                                   instance_importance_sampled,
                                                   verbose, 
                                                   rng);
  print_time(verbose, "Built Scene BVH" );

  //Calculate world bounds and ensure camera is inside infinite area light
  aabb bounding_box_world;
  worldbvh->bounding_box(0,0,bounding_box_world);
  Float world_radius = bounding_box_world.Diag().length() ;
  vec3f world_center  = bounding_box_world.Centroid();
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
      background_sphere = std::make_shared<InfiniteAreaLight>(nx1, ny1, world_radius*2, world_center,
                                                              background_texture, background_material,
                                                              BackgroundTransform,
                                                              BackgroundTransformInv, false);
    } else {
      Rcpp::Rcout << "Failed to load background image at " << background << "\n";
      hasbackground = false;
      ambient_light = true;
      backgroundhigh = vec3f(FLT_MIN,FLT_MIN,FLT_MIN);
      backgroundlow = vec3f(FLT_MIN,FLT_MIN,FLT_MIN);
      background_texture = std::make_shared<gradient_texture>(backgroundlow, backgroundhigh, false, false);
      background_material = std::make_shared<diffuse_light>(background_texture, 1.0, false);
      background_sphere = std::make_shared<InfiniteAreaLight>(100, 100, world_radius*2, world_center,
                                                              background_texture, background_material,
                                                              BackgroundTransform,BackgroundTransformInv,false);
    }
  } else if(ambient_light) {
    //Check if both high and low are black, and set to FLT_MIN
    if(backgroundhigh.length() == 0 && backgroundlow.length() == 0) {
      backgroundhigh = vec3f(FLT_MIN,FLT_MIN,FLT_MIN);
      backgroundlow = vec3f(FLT_MIN,FLT_MIN,FLT_MIN);
    }
    background_texture = std::make_shared<gradient_texture>(backgroundlow, backgroundhigh, false, false);
    background_material = std::make_shared<diffuse_light>(background_texture, 1.0, false);
    background_sphere = std::make_shared<InfiniteAreaLight>(100, 100, world_radius*2, world_center,
                                              background_texture, background_material,
                                              BackgroundTransform, BackgroundTransformInv, false);

  } else {
    //Minimum intensity FLT_MIN so the CDF isn't NAN
    background_texture = std::make_shared<constant_texture>(vec3f(FLT_MIN,FLT_MIN,FLT_MIN));
    background_material = std::make_shared<diffuse_light>(background_texture, 1.0, false);
    background_sphere = std::make_shared<InfiniteAreaLight>(100, 100, world_radius*2, world_center,
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
  
  PreviewDisplay Display(nx,ny, preview, interactive, (lookat-lookfrom).length(), cam.get(),
                         background_sphere->ObjectToWorld.get(),
                         background_sphere->WorldToObject.get());
  
  if(impl_only_bg || hasbackground) {
    imp_sample_objects.add(background_sphere);
  }

  if(min_variance == 0) {
    min_adaptive_size = 1;
    min_variance = 10E-8;
  }
  // Rcpp::Rcout << "Total world size: " << world.GetSize() + texture_bytes << " (Textures: " << texture_bytes << ") \n";
  if(debug_channel != 0) {
    debug_scene(numbercores, nx, ny, ns, debug_channel,
                min_variance, min_adaptive_size,
                routput, goutput,boutput,
                progress_bar, sample_method, stratified_dim,
                verbose, cam.get(), fov,
                world, imp_sample_objects, 
                clampval, max_depth, roulette_active,
                light_direction, rng, sample_dist, keep_colors,
                backgroundhigh);
  } else {
    pathtracer(numbercores, nx, ny, ns, debug_channel,
               min_variance, min_adaptive_size,
               routput, goutput,boutput, alpha_output,
               progress_bar, sample_method, stratified_dim,
               verbose, cam.get(),  fov,
               world, imp_sample_objects,
               clampval, max_depth, roulette_active, Display);
  }
  
  delete shared_materials;
  PutRNGstate();
  print_time(verbose, "Finished rendering" );
  List final_image = List::create(_["r"] = routput.ConvertRcpp(), 
                                  _["g"] = goutput.ConvertRcpp(), 
                                  _["b"] = boutput.ConvertRcpp(),
                                  _["a"] = alpha_output.ConvertRcpp());
  if(Display.Keyframes.size() > 0) {
    List keyframes(Display.Keyframes.size());
    for(unsigned int i = 0; i < Display.Keyframes.size(); i++ ) {
      keyframes(i) = Display.Keyframes[i];
    }
    final_image.attr("keyframes") = keyframes;
  }
#ifdef DEBUG_MEMORY
  Rcpp::Rcout << "Test alloc #" << mem_counter << "\n";
  mem_counter++;
  getInfo();
#endif
  return(final_image);
}

