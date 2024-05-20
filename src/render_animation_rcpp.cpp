#define RCPP_USE_UNWIND_PROTECT

#include "float.h"
#include "vec3.h"
#include "vec2.h"
#include "RayMatrix.h"
#include "mathinline.h"
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
#include "matrix.h"
#include "transform.h"
#include "transformcache.h"
#include "texturecache.h"
#include "debug.h"
#include "bvh_node.h"
using namespace Rcpp;
#include "RcppThread.h"
#include "PreviewDisplay.h"

using namespace std;

// [[Rcpp::export]]
void render_animation_rcpp(List scene, List camera_info, List scene_info, List render_info,
                           List camera_movement, 
                           int start_frame, int end_frame,
                           CharacterVector filenames, Function post_process_frame, int toneval,
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
  Float min_variance = as<Float>(render_info["min_variance"]);
  int min_adaptive_size = as<int>(render_info["min_adaptive_size"]);
  
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
  int bvh_type = as<int>(camera_info["bvh"]);
  NumericMatrix realCameraInfo = as<NumericMatrix>(camera_info["real_camera_info"]);
  Float film_size = as<Float>(camera_info["film_size"]);
  Float camera_scale = as<Float>(camera_info["camera_scale"]);
  Float sample_dist = as<Float>(camera_info["sample_dist"]);
  bool keep_colors = as<bool>(camera_info["keep_colors"]);
  bool preview = as<bool>(camera_info["preview"]);
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
  int n_frames = end_frame;

  vec3f backgroundhigh(bghigh[0],bghigh[1],bghigh[2]);
  vec3f backgroundlow(bglow[0],bglow[1],bglow[2]);

  RcppThread::ThreadPool pool(numbercores);
  GetRNGstate();
  random_gen rng(unif_rand() * std::pow(2,32));
  int nx1, ny1, nn1;
  print_time(verbose, "Loaded Data");

  std::vector<Float* > textures;
  std::vector<unsigned char * > alpha_textures;
  std::vector<unsigned char * > bump_textures;
  std::vector<unsigned char * > roughness_textures;
  //Shared material vector
  std::vector<std::shared_ptr<material> >* shared_materials = new std::vector<std::shared_ptr<material> >;


  //Initialize transformation cache
  TransformCache transformCache;
  
  //Initialize texture cache
  TextureCache texCache;
  
  hitable_list imp_sample_objects;
  std::vector<std::shared_ptr<hitable> > instanced_objects;
  std::vector<std::shared_ptr<hitable_list> > instance_importance_sampled;
  
  std::shared_ptr<bvh_node> worldbvh = build_scene(scene, shape, 
                                                  shutteropen,shutterclose,
                                                  textures, 
                                                  alpha_textures,
                                                  bump_textures,
                                                  roughness_textures, 
                                                  shared_materials, bvh_type,
                                                  transformCache, 
                                                  texCache,
                                                  imp_sample_objects,
                                                  instanced_objects,
                                                  instance_importance_sampled,
                                                  verbose, rng);
  print_time(verbose, "Built Scene BVH" );
  

  //Calculate world bounds
  aabb bounding_box_world;
  worldbvh->bounding_box(0,0,bounding_box_world);
  Float world_radius = bounding_box_world.Diag().length() ;
  vec3f world_center  = bounding_box_world.Centroid();
  for(int i = 0; i < cam_x.length(); i++) {
    vec3f lf(cam_x(i),cam_y(i),cam_z(i));
    world_radius = world_radius > (lf - world_center).length() ? world_radius : 
      1.1*(lf - world_center).length();
  }

  
  std::shared_ptr<texture> background_texture = nullptr;
  std::shared_ptr<material> background_material = nullptr;
  std::shared_ptr<hitable> background_sphere = nullptr;
  Float *background_texture_data = nullptr;
  Matrix4x4 Identity;
  Transform BackgroundAngle(Identity);
  if(rotate_env != 0) {
    BackgroundAngle = Translate(world_center) * RotateY(rotate_env);
  } else {
    BackgroundAngle = Translate(world_center);
  }
  std::shared_ptr<Transform> BackgroundTransform = transformCache.Lookup(BackgroundAngle);
  std::shared_ptr<Transform> BackgroundTransformInv = transformCache.Lookup(BackgroundAngle.GetInverseMatrix());

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
  //Initialize background
  print_time(verbose, "Loaded Background" );
  hitable_list world;
  world.add(worldbvh);

  bool impl_only_bg = false;
  world.add(background_sphere);
  if((imp_sample_objects.size() == 0 || hasbackground || ambient_light) && debug_channel != 18) {
    impl_only_bg = true;
  }
  
  PreviewDisplay d(nx, ny, preview, false, 20.0f, cam.get(), 
                   background_sphere->ObjectToWorld.get(),
                   background_sphere->WorldToObject.get());
  
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
  if(min_variance == 0) {
    min_adaptive_size = 1;
    min_variance = 10E-8;
  }

  if(debug_channel != 0) {
    for(int i = start_frame; i < n_frames; i++ ) {
      if(progress_bar) {
        pb_frames.tick();
      }
      vec3f lookfrom = vec3f(cam_x(i),cam_y(i),cam_z(i));
      vec3f lookat = vec3f(cam_dx(i),cam_dy(i),cam_dz(i));
      Float fov = cam_fov(i);
      Float aperture = cam_aperture(i);
      Float focus_distance = cam_focal(i);
      Float orthox = cam_orthox(i);
      Float orthoy = cam_orthoy(i);
      vec3f camera_up = vec3f(cam_upx(i),cam_upy(i),cam_upz(i));


      if(fov < 0) {
        Transform CamTransform = LookAt(lookfrom,
                                        lookat,
                                        camera_up).GetInverseMatrix();
        std::shared_ptr<Transform> CameraTransform = transformCache.Lookup(CamTransform);
        AnimatedTransform CamTr(CameraTransform,0,CameraTransform,0);
        
        std::vector<Float> lensData;
        for(int i = 0; i < realCameraInfo.rows(); i++) {
          for(int j = 0; j < realCameraInfo.cols(); j++) {
            lensData.push_back(realCameraInfo.at(i,j));
          }
        }
        
        if(fov < 0 && lensData.size() == 0) {
          throw std::runtime_error("No lense data passed in lens descriptor file.");
        }
        
        cam = std::unique_ptr<RayCamera>(new RealisticCamera(CamTr,shutteropen, shutterclose,
                                                             aperture, nx,ny, focus_distance, false, lensData,
                                                             film_size, camera_scale, iso, camera_up, CamTransform,
                                                             lookat));
      } else if(fov == 0) {
        cam = std::unique_ptr<RayCamera>(new ortho_camera(lookfrom, lookat, camera_up,
                                                          orthox, orthoy,
                                                          shutteropen, shutterclose));
      } else if (fov == 360) {
        cam = std::unique_ptr<RayCamera>(new environment_camera(lookfrom, lookat, camera_up,
                                                                shutteropen, shutterclose));
      } else {
        cam = std::unique_ptr<RayCamera>(new camera(lookfrom, lookat, camera_up, fov, Float(nx)/Float(ny),
                                                    aperture, focus_distance,
                                                    shutteropen, shutterclose));
      }

      // world_radius = world_radius > (lookfrom - world_center).length() ? world_radius : (lookfrom - world_center).length()*2;

      if(fov == 0) {
        Float ortho_diag = sqrt(pow(orthox,2) + pow(orthoy,2));
        world_radius += ortho_diag;
      }

      //Initialize output matrices
      RayMatrix routput(nx,ny);
      RayMatrix goutput(nx,ny);
      RayMatrix boutput(nx,ny);
      debug_scene(numbercores, nx, ny, ns, debug_channel,
                  min_variance, min_adaptive_size,
                  routput, goutput,boutput,
                  progress_bar, sample_method, stratified_dim,
                  verbose, cam.get(), fov,
                  world, imp_sample_objects,
                  clampval, max_depth, roulette_active,
                  light_direction, rng, sample_dist, keep_colors, backgroundhigh);
      List temp = List::create(_["r"] = routput.ConvertRcpp(), 
                               _["g"] = goutput.ConvertRcpp(), 
                               _["b"] = boutput.ConvertRcpp());
      post_process_frame(temp, debug_channel, as<std::string>(filenames(i)), toneval);
    }
  } else {
    for(int i = start_frame; i < n_frames; i++ ) {
      if(progress_bar) {
        pb_frames.tick();
      }
      vec3f lookfrom = vec3f(cam_x(i),cam_y(i),cam_z(i));
      vec3f lookat = vec3f(cam_dx(i),cam_dy(i),cam_dz(i));
      vec3f camera_up = vec3f(cam_upx(i),cam_upy(i),cam_upz(i));

      Float fov = cam_fov(i);
      Float aperture = cam_aperture(i);
      Float focus_distance = cam_focal(i);
      Float orthox = cam_orthox(i);
      Float orthoy = cam_orthoy(i);

      std::unique_ptr<RayCamera> cam;
      if(fov < 0) {
        Transform CamTransform = LookAt(lookfrom,
                                        lookat,
                                        camera_up).GetInverseMatrix();
        std::shared_ptr<Transform> CameraTransform = transformCache.Lookup(CamTransform);
        
        AnimatedTransform CamTr(CameraTransform,0,CameraTransform,0);
        
        std::vector<Float> lensData;
        for(int i = 0; i < realCameraInfo.rows(); i++) {
          for(int j = 0; j < realCameraInfo.cols(); j++) {
            lensData.push_back(realCameraInfo.at(i,j));
          }
        }
        
        if(fov < 0 && lensData.size() == 0) {
          throw std::runtime_error("No lense data passed in lens descriptor file.");
        }
        
        cam = std::unique_ptr<RayCamera>(new RealisticCamera(CamTr,shutteropen, shutterclose,
                                                             aperture, nx,ny, focus_distance, false, lensData,
                                                             film_size, camera_scale, iso,camera_up,CamTransform,
                                                             lookat));
      } else if(fov == 0) {
        cam = std::unique_ptr<RayCamera>(new ortho_camera(lookfrom, lookat, camera_up,
                                                          orthox, orthoy,
                                                          shutteropen, shutterclose));
      } else if (fov == 360) {
        cam = std::unique_ptr<RayCamera>(new environment_camera(lookfrom, lookat, camera_up,
                                                                shutteropen, shutterclose));
      } else {
        cam = std::unique_ptr<RayCamera>(new camera(lookfrom, lookat, camera_up, fov, Float(nx)/Float(ny),
                                                    aperture, focus_distance,
                                                    shutteropen, shutterclose));
      }

      world_radius = world_radius > (lookfrom - world_center).length() ? world_radius : (lookfrom - world_center).length()*2;

      if(fov == 0) {
        Float ortho_diag = sqrt(pow(orthox,2) + pow(orthoy,2));
        world_radius += ortho_diag;
      }

      //Initialize output matrices
      RayMatrix routput(nx,ny);
      RayMatrix goutput(nx,ny);
      RayMatrix boutput(nx,ny);
      RayMatrix alpha_output(nx,ny);
      
      pathtracer(numbercores, nx, ny, ns, debug_channel,
                 min_variance, min_adaptive_size,
                 routput, goutput, boutput, alpha_output,
                 progress_bar, sample_method, stratified_dim,
                 verbose, cam.get(),  fov,
                 world, imp_sample_objects,
                 clampval, max_depth, roulette_active, d);
      if(d.terminate) {
        break;
      }
      List temp = List::create(_["r"] = routput.ConvertRcpp(), 
                               _["g"] = goutput.ConvertRcpp(), 
                               _["b"] = boutput.ConvertRcpp(),
                               _["a"] = alpha_output.ConvertRcpp());
      post_process_frame(temp, debug_channel, as<std::string>(filenames(i)), toneval, bloom,
                       transparent_background, write_image);
    }
  }

  delete shared_materials;
  PutRNGstate();
  print_time(verbose, "Finished rendering" );
  
}
