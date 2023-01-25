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
List render_scene_rcpp(List camera_info, List scene_info) {
#ifdef DEBUG_MEMORY
  alloc = 0;
  // dealloc = 0;
  memory_used = 0;
  mem_counter = 0;
#endif

  //Unpack scene info
  bool ambient_light = as<bool>(scene_info["ambient_light"]);
  IntegerVector type = as<IntegerVector>(scene_info["type"]);
  NumericVector radius = as<NumericVector>(scene_info["radius"]);
  IntegerVector shape = as<IntegerVector>(scene_info["shape"]);
  List position_list = as<List>(scene_info["position_list"]);
  List properties = as<List>(scene_info["properties"]);
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
  Float clampval = as<Float>(scene_info["clampval"]);
  LogicalVector isgrouped = as<LogicalVector>(scene_info["isgrouped"]);
  List group_transform = as<List>(scene_info["group_transform"]);
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
  Float rotate_env = as<Float>(scene_info["rotate_env"]);
  Float intensity_env = as<Float>(scene_info["intensity_env"]);
  bool verbose = as<bool>(scene_info["verbose"]);
  int debug_channel = as<int>(scene_info["debug_channel"]);
  IntegerVector shared_id_mat = as<IntegerVector>(scene_info["shared_id_mat"]);
  LogicalVector is_shared_mat = as<LogicalVector>(scene_info["is_shared_mat"]);
  Float min_variance = as<Float>(scene_info["min_variance"]);
  int min_adaptive_size = as<int>(scene_info["min_adaptive_size"]);
  List glossyinfo = as<List>(scene_info["glossyinfo"]);
  List image_repeat = as<List>(scene_info["image_repeat"]);
  List csg_info = as<List>(scene_info["csg_info"]);
  List mesh_list = as<List>(scene_info["mesh_list"]);
  List roughness_list = as<List>(scene_info["roughness_list"]);
  List animation_info = as<List>(scene_info["animation_info"]);

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
  

  //Initialize output matrices
  RayMatrix routput(nx,ny);
  RayMatrix goutput(nx,ny);
  RayMatrix boutput(nx,ny);

  vec3f lookfrom(lookfromvec[0],lookfromvec[1],lookfromvec[2]);
  vec3f lookat(lookatvec[0],lookatvec[1],lookatvec[2]);
  vec3f backgroundhigh(bghigh[0],bghigh[1],bghigh[2]);
  vec3f backgroundlow(bglow[0],bglow[1],bglow[2]);
  Float dist_to_focus = focus_distance;
  CharacterVector alpha_files = as<CharacterVector>(alphalist["alpha_temp_file_names"]);
  LogicalVector has_alpha = as<LogicalVector>(alphalist["alpha_tex_bool"]);

  CharacterVector bump_files = as<CharacterVector>(alphalist["bump_temp_file_names"]);
  LogicalVector has_bump = as<LogicalVector>(alphalist["bump_tex_bool"]);
  NumericVector bump_intensity = as<NumericVector>(alphalist["bump_intensity"]);

  CharacterVector roughness_files = as<CharacterVector>(roughness_list["rough_temp_file_names"]);
  LogicalVector has_roughness = as<LogicalVector>(roughness_list["rough_tex_bool"]);
  
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
  std::vector<int* > nx_ny_nn;

  std::vector<unsigned char * > alpha_textures;
  std::vector<int* > nx_ny_nn_alpha;

  std::vector<unsigned char * > bump_textures;
  std::vector<int* > nx_ny_nn_bump;

  std::vector<unsigned char * > roughness_textures;
  std::vector<int* > nx_ny_nn_roughness;
  //Shared material vector
  std::vector<std::shared_ptr<material> >* shared_materials = new std::vector<std::shared_ptr<material> >;

  
  // size_t texture_bytes = 0;
  for(int i = 0; i < n; i++) {
    if(isimage(i)) {
      int nx, ny, nn;
      Float* tex_data = stbi_loadf(filelocation(i), &nx, &ny, &nn, 0);
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
      NumericVector temp_glossy = as<NumericVector>(glossyinfo(i));
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
  print_time(verbose, "Loaded Textures" );
  
  hitable_list imp_sample_objects;
  std::shared_ptr<hitable> worldbvh = build_scene(type, radius, shape, position_list,
                                properties,
                                n,shutteropen,shutterclose,
                                ischeckered, checkercolors,
                                gradient_info,
                                noise, isnoise, noisephase, noiseintensity, noisecolorlist,
                                angle,
                                isimage, has_alpha, alpha_textures, nx_ny_nn_alpha,
                                textures, nx_ny_nn, has_bump, bump_textures, nx_ny_nn_bump,
                                bump_intensity,
                                roughness_textures, nx_ny_nn_roughness, has_roughness,
                                lightintensity, isflipped,
                                isvolume, voldensity, order_rotation_list,
                                isgrouped, group_transform,
                                tri_normal_bools, is_tri_color, tri_color_vert,
                                fileinfo, filebasedir,
                                scale_list, sigmavec, glossyinfo,
                                shared_id_mat, is_shared_mat, shared_materials,
                                image_repeat, csg_info, mesh_list, bvh_type, transformCache,
                                animation_info, implicit_sample, imp_sample_objects,
                                verbose, rng);
  print_time(verbose, "Built Scene BVH" );

  //Calculate world bounds and ensure camera is inside infinite area light
  aabb bounding_box_world;
  worldbvh->bounding_box(0,0,bounding_box_world);
  Float world_radius = bounding_box_world.Diag().length()/2 ;
  vec3f world_center  = bounding_box_world.Centroid();
  world_radius = world_radius > (lookfrom - world_center).length() ? world_radius : (lookfrom - world_center ).length();
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
    background_texture_data = stbi_loadf(background[0], &nx1, &ny1, &nn1, 0);
    // texture_bytes += nx1 * ny1 * nn1;
    
    if(background_texture_data) {
      background_texture = std::make_shared<image_texture_float>(background_texture_data, nx1, ny1, nn1, 1, 1, intensity_env);
      background_material = std::make_shared<diffuse_light>(background_texture, 1.0, false);
      background_sphere = std::make_shared<InfiniteAreaLight>(nx1, ny1, world_radius*2, vec3f(0.f),
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
      background_sphere = std::make_shared<InfiniteAreaLight>(100, 100, world_radius*2, vec3f(0.f),
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
    background_sphere = std::make_shared<InfiniteAreaLight>(100, 100, world_radius*2, vec3f(0.f),
                                              background_texture, background_material,
                                              BackgroundTransform,BackgroundTransformInv,false);

  } else {
    //Minimum intensity FLT_MIN so the CDF isn't NAN
    background_texture = std::make_shared<constant_texture>(vec3f(FLT_MIN,FLT_MIN,FLT_MIN));
    background_material = std::make_shared<diffuse_light>(background_texture, 1.0, false);
    background_sphere = std::make_shared<InfiniteAreaLight>(100, 100, world_radius*2, vec3f(0.f),
                                              background_texture, background_material,
                                              BackgroundTransform,
                                              BackgroundTransformInv,false);
  }
  print_time(verbose, "Loaded background" );
  hitable_list world;
  world.add(worldbvh);

  bool impl_only_bg = false;
  if((imp_sample_objects.size() == 0 || hasbackground || ambient_light || interactive) && debug_channel != 18) {
    world.add(background_sphere);
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
               routput, goutput,boutput,
               progress_bar, sample_method, stratified_dim,
               verbose, cam.get(),  fov,
               world, imp_sample_objects,
               clampval, max_depth, roulette_active, Display);
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
  print_time(verbose, "Finished rendering" );
  List final_image = List::create(_["r"] = routput.ConvertRcpp(), 
                                  _["g"] = goutput.ConvertRcpp(), 
                                  _["b"] = boutput.ConvertRcpp());
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


// [[Rcpp::export]] 
void PrintClassSizes() {
  Rcpp::Rcout << "hitable                  : " << sizeof(hitable) << "\n";
  Rcpp::Rcout << "---sphere                : " << sizeof(sphere) << "\n";
  Rcpp::Rcout << "---xy_rect               : " << sizeof(xy_rect) << "\n";
  Rcpp::Rcout << "---xz_rect               : " << sizeof(xz_rect) << "\n";
  Rcpp::Rcout << "---yz_rect               : " << sizeof(yz_rect) << "\n";
  Rcpp::Rcout << "---box                   : " << sizeof(box) << "\n";
  Rcpp::Rcout << "---AnimatedHitable       : " << sizeof(AnimatedHitable) << "\n";
  Rcpp::Rcout << "---triangle              : " << sizeof(triangle) << "\n";
  Rcpp::Rcout << "---trimesh               : " << sizeof(trimesh) << "\n";
  Rcpp::Rcout << "---disk                  : " << sizeof(disk) << "\n";
  Rcpp::Rcout << "---cylinder              : " << sizeof(cylinder) << "\n";
  Rcpp::Rcout << "---ellipsoid             : " << sizeof(ellipsoid) << "\n";
  Rcpp::Rcout << "---cone                  : " << sizeof(cone) << "\n";
  Rcpp::Rcout << "---curve                 : " << sizeof(curve) << "\n";
  Rcpp::Rcout << "---csg                   : " << sizeof(csg) << "\n";
  Rcpp::Rcout << "---plymesh               : " << sizeof(plymesh) << "\n";
  Rcpp::Rcout << "---mesh3d                : " << sizeof(mesh3d) << "\n\n";
  Rcpp::Rcout << "---InfiniteAreaLight     : " << sizeof(InfiniteAreaLight) << "\n\n";
  
  Rcpp::Rcout << "Float                    : " << sizeof(Float) << "\n";
  Rcpp::Rcout << "Matrix4x4                : " << sizeof(Matrix4x4) << "\n";
  Rcpp::Rcout << "Transform                : " << sizeof(Transform) << "\n";
  Rcpp::Rcout << "TransformCache           : " << sizeof(TransformCache) << "\n";
  Rcpp::Rcout << "RayMatrix                : " << sizeof(RayMatrix) << "\n";
  Rcpp::Rcout << "vec3f                    : " << sizeof(vec3f) << "\n";
  Rcpp::Rcout << "vec3i                    : " << sizeof(vec3i) << "\n";
  Rcpp::Rcout << "point3f                  : " << sizeof(point3f) << "\n";
  Rcpp::Rcout << "point3f                  : " << sizeof(point3f) << "\n";
  Rcpp::Rcout << "vec2f                    : " << sizeof(vec2f) << "\n";
  Rcpp::Rcout << "vec2i                    : " << sizeof(vec2i) << "\n";
  Rcpp::Rcout << "random_gen               : " << sizeof(random_gen) << "\n";
  Rcpp::Rcout << "aabb                     : " << sizeof(aabb) << "\n\n";
  Rcpp::Rcout << "bvh_node                 : " << sizeof(bvh_node) << "\n\n";
  Rcpp::Rcout << "hitable_list             : " << sizeof(hitable_list) << "\n\n";
  
  Rcpp::Rcout << "RayCamera                : " << sizeof(RayCamera) << "\n";
  Rcpp::Rcout << "---camera                : " << sizeof(camera) << "\n";
  Rcpp::Rcout << "---ortho_camera          : " << sizeof(ortho_camera) << "\n";
  Rcpp::Rcout << "---environment_camera    : " << sizeof(RayCamera) << "\n";
  Rcpp::Rcout << "---RealisticCamera       : " << sizeof(RealisticCamera) << "\n\n";
  Rcpp::Rcout << "material                 : " << sizeof(material) << "\n";
  Rcpp::Rcout << "---lambertian            : " << sizeof(lambertian) << "\n";
  Rcpp::Rcout << "---metal                 : " << sizeof(metal) << "\n";
  Rcpp::Rcout << "---dielectric            : " << sizeof(dielectric) << "\n";
  Rcpp::Rcout << "---diffuse_light         : " << sizeof(diffuse_light) << "\n";
  Rcpp::Rcout << "---spot_light            : " << sizeof(spot_light) << "\n";
  Rcpp::Rcout << "---isotropic             : " << sizeof(isotropic) << "\n";
  Rcpp::Rcout << "---orennayar             : " << sizeof(orennayar) << "\n";
  Rcpp::Rcout << "---MicrofacetReflection  : " << sizeof(MicrofacetReflection) << "\n";
  Rcpp::Rcout << "---MicrofacetTransmission: " << sizeof(MicrofacetTransmission) << "\n";
  Rcpp::Rcout << "---glossy                : " << sizeof(glossy) << "\n";
  Rcpp::Rcout << "---hair                  : " << sizeof(hair) << "\n";
}

