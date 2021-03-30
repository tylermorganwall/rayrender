#define RCPP_USE_UNWIND_PROTECT

#include "vec3.h"
#include "vec2.h"
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
#include "debug.h"
using namespace Rcpp;
#include "RcppThread.h"

using namespace std;

// [[Rcpp::export]]
void render_animation_rcpp(List camera_info, List scene_info, List camera_movement, int start_frame,
                           CharacterVector filenames, Function post_process_frame, int toneval,
                           bool bloom) {
  
  //Unpack scene info
  bool ambient_light = as<bool>(scene_info["ambient_light"]);
  IntegerVector type = as<IntegerVector>(scene_info["type"]);
  NumericVector radius = as<NumericVector>(scene_info["radius"]);
  IntegerVector shape = as<IntegerVector>(scene_info["shape"]);
  List position_list = as<List>(scene_info["position_list"]);
  List properties = as<List>(scene_info["properties"]);
  List velocity = as<List>(scene_info["velocity"]);
  LogicalVector moving = as<LogicalVector>(scene_info["moving"]);
  int n = as<int>(scene_info["n"]);
  NumericVector bghigh  = as<NumericVector>(scene_info["bghigh"]);
  NumericVector bglow = as<NumericVector>(scene_info["bglow"]);
  LogicalVector ischeckered = as<LogicalVector>(scene_info["ischeckered"]);
  List checkercolors = as<List>(scene_info["checkercolors"]);
  List gradient_info = as<List>(scene_info["gradient_info"]);
  NumericVector noise = as<NumericVector>(scene_info["noise"]);
  LogicalVector isnoise = as<LogicalVector>(scene_info["isnoise"]);
  NumericVector noisephase = as<NumericVector>(scene_info["noisephase"]);
  NumericVector noiseintensity = as<NumericVector>(scene_info["noiseintensity"]);
  List noisecolorlist = as<List>(scene_info["noisecolorlist"]);
  List angle = as<List>(scene_info["angle"]);
  LogicalVector isimage = as<LogicalVector>(scene_info["isimage"]);
  CharacterVector filelocation = as<CharacterVector>(scene_info["filelocation"]);
  List alphalist = as<List>(scene_info["alphalist"]);
  NumericVector lightintensity = as<NumericVector>(scene_info["lightintensity"]);
  LogicalVector isflipped = as<LogicalVector>(scene_info["isflipped"]);
  LogicalVector isvolume = as<LogicalVector>(scene_info["isvolume"]);
  NumericVector voldensity = as<NumericVector>(scene_info["voldensity"]);
  LogicalVector implicit_sample = as<LogicalVector>(scene_info["implicit_sample"]);
  List order_rotation_list = as<List>(scene_info["order_rotation_list"]);
  float clampval = as<float>(scene_info["clampval"]);
  LogicalVector isgrouped = as<LogicalVector>(scene_info["isgrouped"]);
  List group_pivot = as<List>(scene_info["group_pivot"]);
  List group_translate = as<List>(scene_info["group_translate"]);
  List group_angle = as<List>(scene_info["group_angle"]);
  List group_order_rotation = as<List>(scene_info["group_order_rotation"]);
  List group_scale = as<List>(scene_info["group_scale"]);
  LogicalVector tri_normal_bools = as<LogicalVector>(scene_info["tri_normal_bools"]);
  LogicalVector is_tri_color = as<LogicalVector>(scene_info["is_tri_color"]);
  List tri_color_vert = as<List>(scene_info["tri_color_vert"]);
  CharacterVector fileinfo = as<CharacterVector>(scene_info["fileinfo"]);
  CharacterVector filebasedir = as<CharacterVector>(scene_info["filebasedir"]);
  bool progress_bar = as<bool>(scene_info["progress_bar"]);
  int numbercores = as<int>(scene_info["numbercores"]);
  bool hasbackground = as<bool>(scene_info["hasbackground"]);
  CharacterVector background = as<CharacterVector>(scene_info["background"]);
  List scale_list = as<List>(scene_info["scale_list"]);
  NumericVector sigmavec = as<NumericVector>(scene_info["sigmavec"]);
  float rotate_env = as<float>(scene_info["rotate_env"]);
  float intensity_env = as<float>(scene_info["intensity_env"]);
  bool verbose = as<bool>(scene_info["verbose"]);
  int debug_channel = as<int>(scene_info["debug_channel"]);
  IntegerVector shared_id_mat = as<IntegerVector>(scene_info["shared_id_mat"]);
  LogicalVector is_shared_mat = as<LogicalVector>(scene_info["is_shared_mat"]);
  float min_variance = as<float>(scene_info["min_variance"]);
  int min_adaptive_size = as<int>(scene_info["min_adaptive_size"]);
  List glossyinfo = as<List>(scene_info["glossyinfo"]);
  List image_repeat = as<List>(scene_info["image_repeat"]);
  List csg_info = as<List>(scene_info["csg_info"]);
  List mesh_list = as<List>(scene_info["mesh_list"]);
  
  
  auto startfirst = std::chrono::high_resolution_clock::now();
  //Unpack Camera Info
  int nx = as<int>(camera_info["nx"]);
  int ny = as<int>(camera_info["ny"]);
  int ns = as<int>(camera_info["ns"]);
  Float shutteropen = as<Float>(camera_info["shutteropen"]);
  Float shutterclose = as<Float>(camera_info["shutterclose"]);
  size_t max_depth = as<size_t>(camera_info["max_depth"]);
  size_t roulette_active = as<size_t>(camera_info["roulette_active_depth"]);
  int sample_method = as<int>(camera_info["sample_method"]);
  NumericVector stratified_dim = as<NumericVector>(camera_info["stratified_dim"]);
  NumericVector light_direction = as<NumericVector>(camera_info["light_direction"]);
  int bvh_type = as<int>(camera_info["bvh"]);
  
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
  int n_frames = cam_x.size();
  
  vec3 backgroundhigh(bghigh[0],bghigh[1],bghigh[2]);
  vec3 backgroundlow(bglow[0],bglow[1],bglow[2]);
  CharacterVector alpha_files = as<CharacterVector>(alphalist["alpha_temp_file_names"]);
  LogicalVector has_alpha = as<LogicalVector>(alphalist["alpha_tex_bool"]);
  
  CharacterVector bump_files = as<CharacterVector>(alphalist["bump_temp_file_names"]);
  LogicalVector has_bump = as<LogicalVector>(alphalist["bump_tex_bool"]);
  NumericVector bump_intensity = as<NumericVector>(alphalist["bump_intensity"]);
  
  
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
  
  std::vector<Float* > alpha_textures;
  std::vector<int* > nx_ny_nn_alpha;
  
  std::vector<Float* > bump_textures;
  std::vector<int* > nx_ny_nn_bump;
  //Shared material vector
  std::vector<std::shared_ptr<material> >* shared_materials = new std::vector<std::shared_ptr<material> >;
  
  for(int i = 0; i < n; i++) {
    if(isimage(i)) {
      int nx, ny, nn;
      Float* tex_data = stbi_loadf(filelocation(i), &nx, &ny, &nn, 4);
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
      Float* tex_data_alpha = stbi_loadf(alpha_files(i), &nxa, &nya, &nna, 4);
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
      Float* tex_data_bump = stbi_loadf(bump_files(i), &nxb, &nyb, &nnb, 4);
      bump_textures.push_back(tex_data_bump);
      nx_ny_nn_bump.push_back(new int[3]);
      nx_ny_nn_bump[i][0] = nxb;
      nx_ny_nn_bump[i][1] = nyb;
      nx_ny_nn_bump[i][2] = nnb;
    } else {
      bump_textures.push_back(nullptr);
      nx_ny_nn_bump.push_back(nullptr);
    }
  }
  
  
  std::shared_ptr<hitable> worldbvh = build_scene(type, radius, shape, position_list,
                                                  properties, velocity, moving,
                                                  n,shutteropen,shutterclose,
                                                  ischeckered, checkercolors, 
                                                  gradient_info,
                                                  noise, isnoise, noisephase, noiseintensity, noisecolorlist,
                                                  angle, 
                                                  isimage, has_alpha, alpha_textures, nx_ny_nn_alpha,
                                                  textures, nx_ny_nn, has_bump, bump_textures, nx_ny_nn_bump,
                                                  bump_intensity,
                                                  lightintensity, isflipped,
                                                  isvolume, voldensity, order_rotation_list, 
                                                  isgrouped, group_pivot, group_translate,
                                                  group_angle, group_order_rotation, group_scale,
                                                  tri_normal_bools, is_tri_color, tri_color_vert, 
                                                  fileinfo, filebasedir, 
                                                  scale_list, sigmavec, glossyinfo,
                                                  shared_id_mat, is_shared_mat, shared_materials,
                                                  image_repeat, csg_info, mesh_list, bvh_type, rng);
  auto finish = std::chrono::high_resolution_clock::now();
  if(verbose) {
    std::chrono::duration<double> elapsed = finish - start;
    Rcpp::Rcout << elapsed.count() << " seconds" << "\n";
  }
  
  //Calculate world bounds
  aabb bounding_box_world;
  worldbvh->bounding_box(0,0,bounding_box_world);
  Float world_radius = bounding_box_world.diag.length()/2 ;
  vec3 world_center  = bounding_box_world.centroid;
  
  //Initialize background
  if(verbose && hasbackground) {
    Rcpp::Rcout << "Loading Environment Image: ";
  }
  start = std::chrono::high_resolution_clock::now();
  std::shared_ptr<texture> background_texture = nullptr;
  std::shared_ptr<material> background_material = nullptr;
  std::shared_ptr<hitable> background_sphere = nullptr;
  Float *background_texture_data = nullptr;
  
  if(hasbackground) {
    background_texture_data = stbi_loadf(background[0], &nx1, &ny1, &nn1, 0);
    background_texture = std::make_shared<image_texture>(background_texture_data, nx1, ny1, nn1, 1, 1, intensity_env);
    background_material = std::make_shared<diffuse_light>(background_texture, 1.0, false);
    background_sphere = std::make_shared<InfiniteAreaLight>(nx1, ny1, world_radius*2, world_center,
                                                            background_texture, background_material);
    if(rotate_env != 0) {
      background_sphere = std::make_shared<rotate_y>(background_sphere, rotate_env);
    }
  } else if(ambient_light) {
    //Check if both high and low are black, and set to FLT_MIN
    if(backgroundhigh.length() == 0 && backgroundlow.length() == 0) {
      backgroundhigh = vec3(FLT_MIN,FLT_MIN,FLT_MIN);
      backgroundlow = vec3(FLT_MIN,FLT_MIN,FLT_MIN);
    }
    background_texture = std::make_shared<gradient_texture>(backgroundlow, backgroundhigh, false, false);
    background_material = std::make_shared<diffuse_light>(background_texture, 1.0, false);
    background_sphere = std::make_shared<InfiniteAreaLight>(100, 100, world_radius*2, world_center,
                                                            background_texture, background_material);
  } else {
    //Minimum intensity FLT_MIN so the CDF isn't NAN
    background_texture = std::make_shared<constant_texture>(vec3(FLT_MIN,FLT_MIN,FLT_MIN));
    background_material = std::make_shared<diffuse_light>(background_texture, 1.0, false);
    background_sphere = std::make_shared<InfiniteAreaLight>(100, 100, world_radius*2, world_center,
                                                            background_texture, background_material);
  }
  finish = std::chrono::high_resolution_clock::now();
  if(verbose && hasbackground) {
    std::chrono::duration<double> elapsed = finish - start;
    Rcpp::Rcout << elapsed.count() << " seconds" << "\n";
  }
  int numbertosample = 0;
  for(int i = 0; i < implicit_sample.size(); i++) {
    if(implicit_sample(i)) {
      numbertosample++;
    }
  }
  hitable_list world;
  
  world.add(worldbvh);
  
  bool impl_only_bg = false;
  if(numbertosample == 0 || hasbackground || ambient_light) {
    world.add(background_sphere);
    impl_only_bg = true;
  }
  
  hitable_list hlist;
  if(verbose) {
    Rcpp::Rcout << "Building Importance Sampling List: ";
  }
  start = std::chrono::high_resolution_clock::now();
  for(int i = 0; i < n; i++)  {
    if(implicit_sample(i)) {
      hlist.add(build_imp_sample(type, radius, shape, position_list,
                                 properties, velocity,
                                 n, shutteropen, shutterclose,
                                 angle, i, order_rotation_list,
                                 isgrouped, group_pivot, group_translate,
                                 group_angle, group_order_rotation, group_scale,
                                 fileinfo, filebasedir, scale_list, 
                                 mesh_list,bvh_type,  moving, rng));
    }
  }
  finish = std::chrono::high_resolution_clock::now();
  if(verbose) {
    std::chrono::duration<double> elapsed = finish - start;
    Rcpp::Rcout << elapsed.count() << " seconds" << "\n";
  }
  if(impl_only_bg || hasbackground) {
    hlist.add(background_sphere);
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
      vec3 lookfrom = vec3(cam_x(i),cam_y(i),cam_z(i));
      vec3 lookat = vec3(cam_dx(i),cam_dy(i),cam_dz(i));
      Float fov = cam_fov(i);
      Float aperture = cam_aperture(i);
      Float dist_to_focus = cam_focal(i);
      Float orthox = cam_orthox(i);
      Float orthoy = cam_orthoy(i);
      vec3 camera_up = vec3(cam_upx(i),cam_upy(i),cam_upz(i));
      
      
      camera cam(lookfrom, lookat, camera_up, fov, float(nx)/float(ny), 
                 aperture, dist_to_focus,
                 shutteropen, shutterclose);
      ortho_camera ocam(lookfrom, lookat, camera_up,
                        orthox, orthoy,
                        shutteropen, shutterclose);
      environment_camera ecam(lookfrom, lookat, camera_up,
                              shutteropen, shutterclose);
      
      world_radius = world_radius > (lookfrom - world_center).length() ? world_radius : (lookfrom - world_center).length()*2;
      
      if(fov == 0) {
        Float ortho_diag = sqrt(pow(orthox,2) + pow(orthoy,2));
        world_radius += ortho_diag;
      }
      
      //Initialize output matrices
      NumericMatrix routput(nx,ny);
      NumericMatrix goutput(nx,ny);
      NumericMatrix boutput(nx,ny);
      debug_scene(numbercores, nx, ny, ns, debug_channel,
                  min_variance, min_adaptive_size, 
                  routput, goutput,boutput,
                  progress_bar, sample_method, stratified_dim,
                  verbose, ocam, cam, ecam, fov,
                  world, hlist,
                  clampval, max_depth, roulette_active, 
                  light_direction, rng);
      List temp = List::create(_["r"] = routput, _["g"] = goutput, _["b"] = boutput);
      post_process_frame(temp, debug_channel, as<std::string>(filenames(i)), toneval);
    }
  } else {
    for(int i = start_frame; i < n_frames; i++ ) {
      if(progress_bar) {
        pb_frames.tick();
      }
      vec3 lookfrom = vec3(cam_x(i),cam_y(i),cam_z(i));
      vec3 lookat = vec3(cam_dx(i),cam_dy(i),cam_dz(i));
      vec3 camera_up = vec3(cam_upx(i),cam_upy(i),cam_upz(i));
      
      Float fov = cam_fov(i);
      Float aperture = cam_aperture(i);
      Float dist_to_focus = cam_focal(i);
      Float orthox = cam_orthox(i);
      Float orthoy = cam_orthoy(i);
      
      camera cam(lookfrom, lookat, camera_up, fov, float(nx)/float(ny), 
                 aperture, dist_to_focus,
                 shutteropen, shutterclose);
      ortho_camera ocam(lookfrom, lookat, camera_up,
                        orthox, orthoy,
                        shutteropen, shutterclose);
      environment_camera ecam(lookfrom, lookat, camera_up,
                              shutteropen, shutterclose);
      
      world_radius = world_radius > (lookfrom - world_center).length() ? world_radius : (lookfrom - world_center).length()*2;
      
      if(fov == 0) {
        Float ortho_diag = sqrt(pow(orthox,2) + pow(orthoy,2));
        world_radius += ortho_diag;
      }
      
      //Initialize output matrices
      NumericMatrix routput(nx,ny);
      NumericMatrix goutput(nx,ny);
      NumericMatrix boutput(nx,ny);
      pathtracer(numbercores, nx, ny, ns, debug_channel,
                 min_variance, min_adaptive_size, 
                 routput, goutput,boutput,
                 progress_bar, sample_method, stratified_dim,
                 verbose, ocam, cam, ecam, fov,
                 world, hlist,
                 clampval, max_depth, roulette_active);
      List temp = List::create(_["r"] = routput, _["g"] = goutput, _["b"] = boutput);
      post_process_frame(temp, debug_channel, as<std::string>(filenames(i)), toneval, bloom);
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
      delete nx_ny_nn[i];
    } 
    if(has_alpha(i)) {
      stbi_image_free(alpha_textures[i]);
      delete nx_ny_nn_alpha[i];
    }
    if(has_bump(i)) {
      stbi_image_free(bump_textures[i]);
      delete nx_ny_nn_bump[i];
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
