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
#include "debug.h"
using namespace Rcpp;
#include "RcppThread.h"
#include "PreviewDisplay.h"

using namespace std;

// [[Rcpp::export]]
void render_animation_rcpp(List scene, List camera_info, List scene_info, List camera_movement, 
                           int start_frame, int end_frame,
                           CharacterVector filenames, Function post_process_frame, int toneval,
                           bool bloom, bool write_image, bool transparent_background) {
  
  size_t n = scene.length();
  //Unpack scene info
  bool ambient_light = as<bool>(scene_info["ambient_light"]);
  IntegerVector shape = as<IntegerVector>(scene_info["shape"]);
  List position_list = as<List>(scene_info["position_list"]);
  NumericVector bghigh  = as<NumericVector>(scene_info["bghigh"]);
  NumericVector bglow = as<NumericVector>(scene_info["bglow"]);
  LogicalVector isimage = as<LogicalVector>(scene_info["isimage"]);
  CharacterVector filelocation = as<CharacterVector>(scene_info["filelocation"]);
  List alphalist = as<List>(scene_info["alphalist"]);
  Float clampval = as<Float>(scene_info["clampval"]);
  bool progress_bar = as<bool>(scene_info["progress_bar"]);
  int numbercores = as<int>(scene_info["numbercores"]);
  bool hasbackground = as<bool>(scene_info["hasbackground"]);
  CharacterVector background = as<CharacterVector>(scene_info["background"]);
  Float rotate_env = as<Float>(scene_info["rotate_env"]);
  Float intensity_env = as<Float>(scene_info["intensity_env"]);
  bool verbose = as<bool>(scene_info["verbose"]);
  int debug_channel = as<int>(scene_info["debug_channel"]);
  Float min_variance = as<Float>(scene_info["min_variance"]);
  int min_adaptive_size = as<int>(scene_info["min_adaptive_size"]);
  List roughness_list = as<List>(scene_info["roughness_list"]);


  auto startfirst = std::chrono::high_resolution_clock::now();
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
  CharacterVector alpha_files = as<CharacterVector>(alphalist["alpha_temp_file_names"]);
  LogicalVector has_alpha = as<LogicalVector>(alphalist["alpha_tex_bool"]);

  CharacterVector bump_files = as<CharacterVector>(alphalist["bump_temp_file_names"]);
  LogicalVector has_bump = as<LogicalVector>(alphalist["bump_tex_bool"]);

  CharacterVector roughness_files = as<CharacterVector>(roughness_list["rough_temp_file_names"]);
  LogicalVector has_roughness = as<LogicalVector>(roughness_list["rough_tex_bool"]);

  RcppThread::ThreadPool pool(numbercores);
  GetRNGstate();
  random_gen rng(unif_rand() * std::pow(2,32));
  int nx1, ny1, nn1;
  auto start = std::chrono::high_resolution_clock::now();
  if(verbose) {
    Rcpp::Rcout << "Building BVH: ";
  }

  std::vector<Float* > textures;
  std::vector<int* > nx_ny_nn;

  std::vector<unsigned char * > alpha_textures;
  std::vector<int* > nx_ny_nn_alpha;

  std::vector<unsigned char * > bump_textures;
  std::vector<int* > nx_ny_nn_bump;

  std::vector<unsigned char * > roughness_textures;
  std::vector<int* > nx_ny_nn_roughness;
  //Shared material vector
  std::vector<std::shared_ptr<material> >* shared_materials = new std::vector<std::shared_ptr<material> >;

  for(size_t i = 0; i < n; i++) {
    if(isimage(i)) {
      int nx, ny, nn;
      Float* tex_data = stbi_loadf(filelocation(i), &nx, &ny, &nn, 4);
      nn = 4;
      // texture_bytes += nx * ny * nn;
      textures.push_back(tex_data);
      nx_ny_nn.push_back(new int[3]);
      nx_ny_nn[i][0] = nx;
      nx_ny_nn[i][1] = ny;
      nx_ny_nn[i][2] = nn;
    } else {
      textures.push_back(nullptr);
      nx_ny_nn.push_back(nullptr);
    }
    if(has_alpha(i)) {
      int nxa, nya, nna;
      unsigned char * tex_data_alpha = stbi_load(alpha_files(i), &nxa, &nya, &nna, 0);
      // texture_bytes += nxa * nya * nna;
      
      alpha_textures.push_back(tex_data_alpha);
      nx_ny_nn_alpha.push_back(new int[3]);
      nx_ny_nn_alpha[i][0] = nxa;
      nx_ny_nn_alpha[i][1] = nya;
      nx_ny_nn_alpha[i][2] = nna;
    } else {
      alpha_textures.push_back(nullptr);
      nx_ny_nn_alpha.push_back(nullptr);
    }
    if(has_bump(i)) {
      int nxb, nyb, nnb;
      unsigned char * tex_data_bump = stbi_load(bump_files(i), &nxb, &nyb, &nnb, 0);
      bump_textures.push_back(tex_data_bump);
      nx_ny_nn_bump.push_back(new int[3]);
      nx_ny_nn_bump[i][0] = nxb;
      nx_ny_nn_bump[i][1] = nyb;
      nx_ny_nn_bump[i][2] = nnb;
    } else {
      bump_textures.push_back(nullptr);
      nx_ny_nn_bump.push_back(nullptr);
    }
    if(has_roughness(i)) {
      List material = as<List>(scene(i))["material"];
      NumericVector temp_glossy = as<NumericVector>(material["glossyinfo"]);
      int nxr, nyr, nnr;
      unsigned char * tex_data_roughness = stbi_load(roughness_files(i), &nxr, &nyr, &nnr, 0);
      // texture_bytes += nxr * nyr * nnr;
      
      Float min = temp_glossy(9), max = temp_glossy(10);
      Float rough_range = max-min;
      Float maxr = 0, minr = 1;
      for(int ii = 0; ii < nxr; ii++) {
        for(int jj = 0; jj < nyr; jj++) {
          Float temp_rough = tex_data_roughness[nnr*ii + nnr*nxr*jj];
          maxr = maxr < temp_rough ? temp_rough : maxr;
          minr = minr > temp_rough ? temp_rough : minr;
          if(nnr > 1) {
            temp_rough = tex_data_roughness[nnr*ii + nnr*nxr*jj+1];
            maxr = maxr < temp_rough ? temp_rough : maxr;
            minr = minr > temp_rough ? temp_rough : minr;
          }
        }
      }
      Float data_range = maxr-minr;
      for(int ii = 0; ii < nxr; ii++) {
        for(int jj = 0; jj < nyr; jj++) {
          if(!temp_glossy(11)) {
            tex_data_roughness[nnr*ii + nnr*nxr*jj] =
              (tex_data_roughness[nnr*ii + nnr*nxr*jj]-minr)/data_range * rough_range + min;
            if(nnr > 1) {
              tex_data_roughness[nnr*ii + nnr*nxr*jj+1] =
                (tex_data_roughness[nnr*ii + nnr*nxr*jj+1]-minr)/data_range * rough_range + min;
            }
          } else {
            tex_data_roughness[nnr*ii + nnr*nxr*jj] =
              (1.0-(tex_data_roughness[nnr*ii + nnr*nxr*jj]-minr)/data_range) * rough_range + min;
            if(nnr > 1) {
              tex_data_roughness[nnr*ii + nnr*nxr*jj+1] =
                (1.0-(tex_data_roughness[nnr*ii + nnr*nxr*jj+1]-minr)/data_range) * rough_range + min;
            }
          }
        }
      }
      roughness_textures.push_back(tex_data_roughness);
      nx_ny_nn_roughness.push_back(new int[3]);
      nx_ny_nn_roughness[i][0] = nxr;
      nx_ny_nn_roughness[i][1] = nyr;
      nx_ny_nn_roughness[i][2] = nnr;
    } else {
      roughness_textures.push_back(nullptr);
      nx_ny_nn_roughness.push_back(nullptr);
    }
  }

  //Initialize transformation cache
  TransformCache transformCache;
  hitable_list imp_sample_objects;
  
  std::shared_ptr<hitable> worldbvh = build_scene(scene, shape, position_list,
                                                  shutteropen,shutterclose,
                                                  textures, nx_ny_nn,
                                                  alpha_textures, nx_ny_nn_alpha,
                                                  bump_textures, nx_ny_nn_bump,
                                                  roughness_textures, nx_ny_nn_roughness,
                                                  shared_materials, bvh_type,
                                                  transformCache, imp_sample_objects,
                                                  verbose, rng);
  
  auto finish = std::chrono::high_resolution_clock::now();
  if(verbose) {
    std::chrono::duration<double> elapsed = finish - start;
    Rcpp::Rcout << elapsed.count() << " seconds" << "\n";
  }

  //Calculate world bounds
  aabb bounding_box_world;
  worldbvh->bounding_box(0,0,bounding_box_world);
  Float world_radius = bounding_box_world.Diag().length()/2 ;
  vec3f world_center  = bounding_box_world.Centroid();
  for(int i = 0; i < cam_x.length(); i++) {
    vec3f lf(cam_x(i),cam_y(i),cam_z(i));
    world_radius = world_radius > (lf - world_center).length() ? world_radius : 
      1.1*(lf - world_center).length();
  }

  //Initialize background
  if(verbose && hasbackground) {
    Rcpp::Rcout << "Loading Environment Image: ";
  }
  start = std::chrono::high_resolution_clock::now();
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
    background_texture_data = stbi_loadf(background[0], &nx1, &ny1, &nn1, 4);
    nn1 = 4;
    if(background_texture_data) {
      background_texture = std::make_shared<image_texture_float>(background_texture_data, nx1, ny1, nn1, 1, 1, intensity_env);
      background_material = std::make_shared<diffuse_light>(background_texture, 1.0, false);
      background_sphere = std::make_shared<InfiniteAreaLight>(nx1, ny1, world_radius*2, world_center,
                                                              background_texture, background_material,
                                                              BackgroundTransform,
                                                              BackgroundTransformInv, false);
    } else {
      Rcpp::Rcout << "Failed to load background image at " << background(0) << "\n";
      if(stbi_failure_reason()) {
        Rcpp::Rcout << stbi_failure_reason() << "\n";
      }
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
                                                            BackgroundTransform,BackgroundTransformInv,false);
  } else {
    //Minimum intensity FLT_MIN so the CDF isn't NAN
    background_texture = std::make_shared<constant_texture>(vec3f(FLT_MIN,FLT_MIN,FLT_MIN));
    background_material = std::make_shared<diffuse_light>(background_texture, 1.0, false);
    background_sphere = std::make_shared<InfiniteAreaLight>(100, 100, world_radius*2, world_center,
                                                            background_texture, background_material,
                                                            BackgroundTransform,BackgroundTransformInv,false);
  }
  finish = std::chrono::high_resolution_clock::now();
  if(verbose && hasbackground) {
    std::chrono::duration<double> elapsed = finish - start;
    Rcpp::Rcout << elapsed.count() << " seconds" << "\n";
  }
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
  
  if(verbose) {
    std::chrono::duration<double> elapsed = finish - start;
    Rcpp::Rcout << elapsed.count() << " seconds" << "\n";
  }
  if(impl_only_bg || hasbackground) {
    imp_sample_objects.add(background_sphere);
  }

  if(verbose && !progress_bar) {
    Rcpp::Rcout << "Starting Raytracing:\n ";
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

  if(verbose) {
    Rcpp::Rcout << "Cleaning up memory..." << "\n";
  }
  if(hasbackground) {
    stbi_image_free(background_texture_data);
  }
  for(int i = 0; i < n; i++) {
    if(isimage(i)) {
      stbi_image_free(textures[i]);
      delete[] nx_ny_nn[i];
    }
    if(has_alpha(i)) {
      stbi_image_free(alpha_textures[i]);
      delete[] nx_ny_nn_alpha[i];
    }
    if(has_bump(i)) {
      stbi_image_free(bump_textures[i]);
      delete[] nx_ny_nn_bump[i];
    }
    if(has_roughness(i)) {
      stbi_image_free(roughness_textures[i]);
      delete[] nx_ny_nn_roughness[i];
    }
  }
  delete shared_materials;
  PutRNGstate();
  finish = std::chrono::high_resolution_clock::now();
  if(verbose) {
    std::chrono::duration<double> elapsed = finish - startfirst;
    Rcpp::Rcout << "Total time elapsed: " << elapsed.count() << " seconds" << "\n";
  }
}
