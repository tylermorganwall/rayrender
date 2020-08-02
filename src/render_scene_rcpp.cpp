#define STB_IMAGE_IMPLEMENTATION 


#ifndef FLOATDEF
#define FLOATDEF
#ifdef RAY_FLOAT_AS_DOUBLE
typedef double Float;
#else
typedef float Float;
#endif 
#endif

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
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
#include "RcppThread.h"

// #define DEBUG

#ifdef DEBUG
#include <iostream>
#include <fstream>
#endif

using namespace std;

inline vec3 de_nan(const vec3& c) {
  vec3 temp = c;
  if(std::isnan(c[0])) temp.e[0] = 0.0f;
  if(std::isnan(c[1])) temp.e[1] = 0.0f;
  if(std::isnan(c[2])) temp.e[2] = 0.0f;
  return(temp);
}


vec3 color(const ray& r, hitable *world, hitable_list *hlist,
           size_t max_depth, size_t roulette_activate, random_gen& rng, Sampler* sampler) {
#ifdef DEBUG
  ofstream myfile;
  myfile.open("rays.txt", ios::app | ios::out);
#endif
  vec3 final_color(0,0,0);
  vec3 throughput(1,1,1);
  float prev_t = 1;
  ray r1 = r;
  ray r2 = r;
  bool diffuse_bounce = false;
  for(size_t i = 0; i < max_depth; i++) {
    hit_record hrec;
    if(world->hit(r2, 0.001, FLT_MAX, hrec, rng)) { //generated hit record, world space
#ifdef DEBUG
      myfile << i << ", " << r2.A << " ";
      myfile << ", " << hrec.p << "\n ";
#endif
      scatter_record srec;
      final_color += throughput * hrec.mat_ptr->emitted(r2, hrec, hrec.u, hrec.v, hrec.p);
      if(throughput.x() == 0 && throughput.y() == 0 && throughput.z() == 0) {
        return(vec3(0,0,0));
      }
      if(i > roulette_activate) {
        float t = std::max(throughput.x(), std::max(throughput.y(), throughput.z()));
        //From Szecsi, Szirmay-Kalos, and Kelemen
        float prob_continue = std::min(1.0f, std::sqrt(t/prev_t));
        prev_t = t;
        if(rng.unif_rand() > prob_continue) {
          return(final_color);
        }
        throughput *= 1 / prob_continue;
      }
      float pdf_val;
      if(hrec.mat_ptr->scatter(r2, hrec, srec, rng)) { //generates scatter record, world space
        if(srec.is_specular) { //returns specular ray
          r2 = srec.specular_ray;
          throughput *= srec.attenuation;
          continue;
        }
        hitable_pdf p_imp(hlist, hrec.p); //creates pdf of all objects to be sampled
        mixture_pdf p(&p_imp, srec.pdf_ptr); //creates mixture pdf of surface intersected at hrec.p and all sampled objects/lights
        
        //Generates a scatter direction (with origin hrec.p) from the mixture 
        //and saves surface normal from light to use in pdf_value calculation
        //(along with the scatter direction)
        //Translates the world space point into object space point, generates ray assuring intersection, and then translates 
        //ray back into world space
        vec3 offset_p = offset_ray(hrec.p-r2.A, hrec.normal) + r2.A;

        r1 = r2;
        if(!diffuse_bounce) {
          diffuse_bounce = true;
          r2 = ray(offset_p, p.generate(sampler), r2.pri_stack, r2.time()); //scatters a ray from hit point to stratified direction
        } else {
          r2 = ray(offset_p, p.generate(rng), r2.pri_stack, r2.time()); //scatters a ray from hit point to direction
        }
        
        pdf_val = p.value(r2.direction(), rng); //generates a pdf value based the intersection point and the mixture pdf
        throughput *= hrec.mat_ptr->f(r1, hrec, r2) / pdf_val;
      } else {
        return(final_color);
      }
    } else {
      return(final_color);
    }
  }
#ifdef DEBUG
  myfile.close();
#endif
  return(final_color);
}

#ifdef DEBUGBVH
float debug_bvh(const ray& r, hitable *world, random_gen rng) {
  hit_record hrec;
  hrec.bvh_nodes = 0;
  world->hit(r, 0.001, FLT_MAX, hrec, rng);
  return(hrec.bvh_nodes);
}
#endif


float calculate_depth(const ray& r, hitable *world, random_gen rng) {
  hit_record hrec;
  if(world->hit(r, 0.001, FLT_MAX, hrec, rng)) {
    return((r.origin()-hrec.p).length());
  } else {
    return(INFINITY);
  }
}

vec3 calculate_normals(const ray& r, hitable *world, random_gen rng) {
  hit_record hrec;
  if(world->hit(r, 0.001, FLT_MAX, hrec, rng)) {
    hrec.normal.make_unit_vector();
    return((vec3(1,1,1) + hrec.normal)/2);
  } else {
    return(vec3(0,0,0));
  }
}

vec3 calculate_uv(const ray& r, hitable *world, random_gen rng) {
  hit_record hrec;
  if(world->hit(r, 0.001, FLT_MAX, hrec, rng)) {
    return(vec3(hrec.u,hrec.v,1-hrec.u-hrec.v));
  } else {
    return(vec3(0,0,0));
  }
}

// //Does not take into account moving objects
// void calculate_inside(const ray& r_in, hitable *world, random_gen rng) {
//   hit_record hrec;
//   if(world->hit(r_in, 0.001, FLT_MAX, hrec, rng)) {
//     if(hrec.mat_ptr->is_dielectric()) {
//       bool encountered = false;
//       int current_index = -1;
//       for(size_t i = 0; i < r_in.pri_stack->size(); i++) {
//         if(r_in.pri_stack->at(i) == hrec.mat_ptr) {
//           encountered = true;
//           current_index = i;
//         } 
//       }
//       if(dot(r_in.direction(),hrec.normal) > 0) {
//         if(encountered) {
//           r_in.pri_stack->erase(r_in.pri_stack->begin() + current_index);
//         } 
//       } else {
//         r_in.pri_stack->push_back((dielectric*)hrec.mat_ptr);
//       }
//     }
//     ray r = ray(hrec.p, r_in.direction(),  r_in.pri_stack);
//     calculate_inside(r, world, rng);
//   } 
// }

// [[Rcpp::export]]
List render_scene_rcpp(List camera_info, bool ambient_light,
                       IntegerVector type, 
                       NumericVector radius, IntegerVector shape,
                       List position_list,
                       List properties, List velocity, LogicalVector moving,
                       int n,
                       NumericVector& bghigh, NumericVector& bglow,
                       LogicalVector ischeckered, List checkercolors, 
                       List gradient_info,
                       NumericVector noise, LogicalVector isnoise,
                       NumericVector& noisephase, NumericVector& noiseintensity, List noisecolorlist,
                       List& angle,
                       LogicalVector& isimage, CharacterVector& filelocation,
                       List& alphalist,
                       NumericVector& lightintensity,
                       LogicalVector& isflipped,
                       LogicalVector& isvolume, NumericVector& voldensity,
                       LogicalVector& implicit_sample, List& order_rotation_list,
                       float clampval,
                       LogicalVector& isgrouped, List& group_pivot, List& group_translate,
                       List& group_angle, List& group_order_rotation, List& group_scale,
                       LogicalVector& tri_normal_bools, LogicalVector& is_tri_color, List& tri_color_vert,
                       CharacterVector& fileinfo, CharacterVector& filebasedir, 
                       bool progress_bar, int numbercores,
                       bool hasbackground, CharacterVector& background, List& scale_list,
                       NumericVector sigmavec,
                       float rotate_env, float intensity_env, bool verbose, int debug_channel,
                       IntegerVector& shared_id_mat, LogicalVector& is_shared_mat,
                       float min_variance, int min_adaptive_size, List glossyinfo,
                       List image_repeat) {
  auto startfirst = std::chrono::high_resolution_clock::now();
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
  size_t max_depth = as<size_t>(camera_info["max_depth"]);
  size_t roulette_active = as<size_t>(camera_info["roulette_active_depth"]);
  int sample_method = as<int>(camera_info["sample_method"]);
  NumericVector stratified_dim = as<NumericVector>(camera_info["stratified_dim"]);
  
  //Initialize output matrices
  NumericMatrix routput(nx,ny);
  NumericMatrix goutput(nx,ny);
  NumericMatrix boutput(nx,ny);
  
  vec3 lookfrom(lookfromvec[0],lookfromvec[1],lookfromvec[2]);
  vec3 lookat(lookatvec[0],lookatvec[1],lookatvec[2]);
  vec3 backgroundhigh(bghigh[0],bghigh[1],bghigh[2]);
  vec3 backgroundlow(bglow[0],bglow[1],bglow[2]);
  float dist_to_focus = focus_distance;
  CharacterVector alpha_files = as<CharacterVector>(alphalist["alpha_temp_file_names"]);
  LogicalVector has_alpha = as<LogicalVector>(alphalist["alpha_tex_bool"]);
  
  CharacterVector bump_files = as<CharacterVector>(alphalist["bump_temp_file_names"]);
  LogicalVector has_bump = as<LogicalVector>(alphalist["bump_tex_bool"]);
  NumericVector bump_intensity = as<NumericVector>(alphalist["bump_intensity"]);
  
  
  RcppThread::ThreadPool pool(numbercores);
  GetRNGstate();
  random_gen rng(unif_rand() * std::pow(2,32));
  camera cam(lookfrom, lookat, vec3(camera_up(0),camera_up(1),camera_up(2)), fov, float(nx)/float(ny), 
             aperture, dist_to_focus,
             shutteropen, shutterclose);
  ortho_camera ocam(lookfrom, lookat, vec3(camera_up(0),camera_up(1),camera_up(2)),
                    ortho_dimensions(0), ortho_dimensions(1),
                    shutteropen, shutterclose, rng);
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
  std::vector<material* >* shared_materials = new std::vector<material* >;
  
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
  
  
  hitable *worldbvh = build_scene(type, radius, shape, position_list,
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
                                image_repeat, rng);
  auto finish = std::chrono::high_resolution_clock::now();
  if(verbose) {
    std::chrono::duration<double> elapsed = finish - start;
    Rcpp::Rcout << elapsed.count() << " seconds" << "\n";
  }
  
  //Calculate world bounds
  aabb bounding_box_world;
  worldbvh->bounding_box(0,0,bounding_box_world);
  Float world_radius = bounding_box_world.diagonal.length()/2 ;
  vec3 world_center  = bounding_box_world.centroid;
  world_radius = world_radius > (lookfrom - world_center).length() ? world_radius : (lookfrom - world_center).length()*2;
  
  if(fov == 0) {
    Float ortho_diag = sqrt(pow(ortho_dimensions(0),2) + pow(ortho_dimensions(1),2));
    world_radius += ortho_diag;
  }
  //Initialize background
  if(verbose && hasbackground) {
    Rcpp::Rcout << "Loading Environment Image: ";
  }
  start = std::chrono::high_resolution_clock::now();
  texture *background_texture = nullptr;
  material *background_material = nullptr;
  hitable *background_sphere = nullptr;
  Float *background_texture_data = nullptr;
  
  if(hasbackground) {
    background_texture_data = stbi_loadf(background[0], &nx1, &ny1, &nn1, 0);
    background_texture = new image_texture(background_texture_data, nx1, ny1, nn1, 1, 1, intensity_env);
    background_material = new diffuse_light(background_texture);
    background_sphere = new InfiniteAreaLight(nx1, ny1, world_radius*2, world_center,
                                              background_texture, background_material);
    if(rotate_env != 0) {
      background_sphere = new rotate_y(background_sphere, rotate_env);
    }
  } else if(ambient_light) {
    //Check if both high and low are black, and set to FLT_MIN
    if(backgroundhigh.length() == 0 && backgroundlow.length() == 0) {
      backgroundhigh = vec3(FLT_MIN,FLT_MIN,FLT_MIN);
      backgroundlow = vec3(FLT_MIN,FLT_MIN,FLT_MIN);
    }
    background_texture = new gradient_texture(backgroundlow, backgroundhigh, false, false);
    background_material = new diffuse_light(background_texture);
    background_sphere = new InfiniteAreaLight(100, 100, world_radius*2, world_center,
                                              background_texture, background_material);
  } else {
    //Minimum intensity FLT_MIN so the CDF isn't NAN
    background_texture = new constant_texture(vec3(FLT_MIN,FLT_MIN,FLT_MIN));
    background_material = new diffuse_light(background_texture);
    background_sphere = new InfiniteAreaLight(100, 100, world_radius*2, world_center,
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
  hitable *world_full[2];
  world_full[0] = worldbvh;
  world_full[1] = background_sphere;
  int nval = numbertosample == 0 || hasbackground || ambient_light ? 2 : 1;
  
  hitable_list world(world_full, nval);
  bool impl_only_bg = false;
  if(numbertosample == 0 || hasbackground) {
    numbertosample++; 
    impl_only_bg = true;
  }

  std::vector<hitable* > implicit_sample_vector(numbertosample);
  int counter = 0;
  if(verbose) {
    Rcpp::Rcout << "Building Importance Sampling List: ";
  }
  start = std::chrono::high_resolution_clock::now();
  for(int i = 0; i < n; i++)  {
    if(implicit_sample(i)) {
      implicit_sample_vector[counter] = build_imp_sample(type, radius, shape, position_list,
                               properties, velocity,
                               n, shutteropen, shutterclose,
                               angle, i, order_rotation_list,
                               isgrouped, group_pivot, group_translate,
                               group_angle, group_order_rotation, group_scale,
                               fileinfo, filebasedir, scale_list, rng);
      counter++;
    }
  }
  finish = std::chrono::high_resolution_clock::now();
  if(verbose) {
    std::chrono::duration<double> elapsed = finish - start;
    Rcpp::Rcout << elapsed.count() << " seconds" << "\n";
  }
  if(impl_only_bg || hasbackground) {
    implicit_sample_vector[counter] = background_sphere;
  }

  hitable_list hlist;
  hlist = hitable_list(&implicit_sample_vector[0], numbertosample);

  if(verbose && !progress_bar) {
    Rcpp::Rcout << "Starting Raytracing:\n ";
  }
  RProgress::RProgress pb_sampler("Generating Samples [:bar] :percent%");
  pb_sampler.set_width(70);
  RProgress::RProgress pb("Adaptive Raytracing [:bar] :percent%");
  pb.set_width(70);

  if(progress_bar) {
    pb_sampler.set_total(ny);
    if(min_variance != 0.0f) {
      pb.set_total(ns);
    } else {
      pb.set_total(ny);
    }
  }
  if(debug_channel == 1) {
    Float depth_into_scene = 0.0;
    for(int j = ny - 1; j >= 0; j--) {
      for(int i = 0; i < nx; i++) {
        depth_into_scene = 0;
        Float u = Float(i) / Float(nx);
        Float v = Float(j) / Float(ny);
        ray r;
        if(fov != 0) {
          r = cam.get_ray(u,v, vec3(0,0,0), 0);
        } else {
          r = ocam.get_ray(u,v);
        }
        depth_into_scene = calculate_depth(r, &world, rng);
        routput(i,j) = depth_into_scene;
        goutput(i,j) = depth_into_scene;
        boutput(i,j) = depth_into_scene;
      }
    }
  } else if(debug_channel == 2) {
    vec3 normal_map(0,0,0);
    for(int j = ny - 1; j >= 0; j--) {
      for(int i = 0; i < nx; i++) {
        normal_map = vec3(0,0,0);
        Float u = Float(i) / Float(nx);
        Float v = Float(j) / Float(ny);
        ray r;
        if(fov != 0) {
          r = cam.get_ray(u,v, vec3(0,0,0), 0);
        } else {
          r = ocam.get_ray(u,v);
        }
        normal_map = calculate_normals(r, &world, rng);
        routput(i,j) = normal_map.x();
        goutput(i,j) = normal_map.y();
        boutput(i,j) = normal_map.z();
      }
    }
  } else if(debug_channel == 3) {
    vec3 uv_map(0,0,0);
    for(int j = ny - 1; j >= 0; j--) {
      for(int i = 0; i < nx; i++) {
        uv_map = vec3(0,0,0);
        Float u = Float(i) / Float(nx);
        Float v = Float(j) / Float(ny);
        ray r;
        if(fov != 0) {
          r = cam.get_ray(u,v, vec3(0,0,0), 0);
        } else {
          r = ocam.get_ray(u,v);
        }
        uv_map = calculate_uv(r, &world, rng);
        routput(i,j) = uv_map.x();
        goutput(i,j) = uv_map.y();
        boutput(i,j) = uv_map.z();
      }
    }
  } else if(debug_channel == 4) {
#ifdef DEBUGBVH

    Float bvh_intersections = 0.0;
    Float max_intersections = 0.0;
    for(int j = ny - 1; j >= 0; j--) {
      for(int i = 0; i < nx; i++) {
        Float u = Float(i + rng.unif_rand()) / Float(nx);
        Float v = Float(j + rng.unif_rand()) / Float(ny);
        ray r;
        if(fov != 0) {
          r = cam.get_ray(u,v, vec3(0,0,0));
        } else {
          r = ocam.get_ray(u,v);
        }
        bvh_intersections = 0.0;
        bvh_intersections = debug_bvh(r, &world, rng);
        routput(i,j) = bvh_intersections;
        goutput(i,j) = bvh_intersections;
        boutput(i,j) = bvh_intersections;
        max_intersections = bvh_intersections > max_intersections ? bvh_intersections : max_intersections;
      }
    }
    for(int j = ny - 1; j >= 0; j--) {
      for(int i = 0; i < nx; i++) {
      routput(i,j) = routput(i,j) / max_intersections;
      goutput(i,j) = goutput(i,j) / max_intersections;
      boutput(i,j) = boutput(i,j) / max_intersections;
      }
    }
#endif
  } else {
    if(min_variance == 0) {
      std::vector<unsigned int> seeds(ny);
      for(int i = 0; i < ny; i++) {
        seeds[i] = unif_rand() * std::pow(2,32);
      }
      if(verbose) {
        start = std::chrono::high_resolution_clock::now();
        if(sample_method == 0) {
          Rcpp::Rcout << "Allocating random sampler: ";
        } else {
          Rcpp::Rcout << "Allocating stratified (" << 
            stratified_dim(0)<< "x" << stratified_dim(1) << ") sampler: ";
        }
      }
      std::vector<random_gen > rngs;
      std::vector<Sampler* > samplers;
      for(int i = 0; i < nx * ny; i++) {
        random_gen rng_single(unif_rand() * std::pow(2,32));
        rngs.push_back(rng_single);
        Sampler* strat;
        if(sample_method == 0) {
          strat = new RandomSampler(rng_single);
        } else {
          strat = new StratifiedSampler(stratified_dim(0), stratified_dim(1),
                                        true, 5, rng_single);
        }
        strat->StartPixel(vec2(0,0));
        strat->SetSampleNumber(0);
        samplers.push_back(strat);
      }
      RcppThread::ThreadPool pool(numbercores);
      auto worker = [&routput, &goutput, &boutput,
                     nx, ny, ns, &rngs, &samplers, 
                     fov,
                     &cam, &ocam, &world, &hlist,
                     clampval, max_depth, roulette_active, progress_bar, 
                     numbercores] (int j) {
       if(progress_bar && j % numbercores == 0) {
         RcppThread::Rcout << "Progress (" << numbercores << " core): ";
         RcppThread::Rcout << (int)((1-(double)j/double(ny)) * 100) << "%\r";
       }
       std::vector<dielectric*> *mat_stack = new std::vector<dielectric*>;
       
       for(int i = 0; i < nx; i++) {
         vec3 col(0,0,0);
         int index = j + ny * i;
         
         for(int s = 0; s < ns; s++) {
           vec2 u2 = samplers[index]->Get2D();
           Float u = (Float(i) + u2.x()) / Float(nx);
           Float v = (Float(j) + u2.y()) / Float(ny);
           ray r;
           if(fov != 0) {
             r = cam.get_ray(u, v, rand_to_unit(samplers[index]->Get2D()), 
                             samplers[index]->Get1D());
           } else {
             r = ocam.get_ray(u,v);
           }
           r.pri_stack = mat_stack;
           
           col += clamp(de_nan(color(r, &world, &hlist, max_depth, roulette_active, 
                                     rngs[index], samplers[index])),0,clampval);
           mat_stack->clear();
           samplers[index]->StartNextSample();
         }
         col /= Float(ns);
         routput(i,j) = col[0];
         goutput(i,j) = col[1];
         boutput(i,j) = col[2];
       }
       delete mat_stack;
     };
      for(int j = ny - 1; j >= 0; j--) {
        pool.push(worker,j);
      }
      pool.join();
      for(size_t i = 0; i < samplers.size(); i++) {
        delete samplers.at(i);
      }
    } else {
      NumericMatrix routput2(nx,ny);
      NumericMatrix goutput2(nx,ny);
      NumericMatrix boutput2(nx,ny);
      adaptive_sampler adaptive_pixel_sampler(numbercores, nx, ny, ns, debug_channel,
                                              min_variance, min_adaptive_size,
                                              routput, goutput, boutput,
                                              routput2, goutput2, boutput2);
      std::vector<random_gen > rngs;
      std::vector<Sampler* > samplers;
      start = std::chrono::high_resolution_clock::now();
      
      for(int i = 0; i < nx * ny; i++) {
        if(progress_bar && (i % ny == 0)) {
          pb_sampler.tick();
        }
        random_gen rng_single(unif_rand() * std::pow(2,32));
        rngs.push_back(rng_single);
        Sampler* strat;
        if(sample_method == 0) {
          strat = new RandomSampler(rng_single);
        } else {
          strat = new StratifiedSampler(stratified_dim(0), stratified_dim(1),
                                        true, 5, rng_single);
        }
        strat->StartPixel(vec2(0,0));
        strat->SetSampleNumber(0);
        samplers.push_back(strat);
      }

      if(verbose) {
        finish = std::chrono::high_resolution_clock::now();
        if(sample_method == 0) {
          Rcpp::Rcout << "Allocating random sampler: ";
        } else {
          Rcpp::Rcout << "Allocating stratified (" << 
            stratified_dim(0)<< "x" << stratified_dim(1) << ") sampler: ";
        }
        std::chrono::duration<double> elapsed = finish - start;
        Rcpp::Rcout << elapsed.count() << " seconds" << "\n";
      }
      for(size_t s = 0; s < static_cast<size_t>(ns); s++) {
        Rcpp::checkUserInterrupt();
        if(progress_bar) {
          pb.tick();
        }
        RcppThread::ThreadPool pool(numbercores);
        auto worker = [&adaptive_pixel_sampler,
                       nx, ny, s,
                       &rngs, fov, &samplers,
                       &cam, &ocam, &world, &hlist,
                       clampval, max_depth, roulette_active] (int k) {
           int nx_begin = adaptive_pixel_sampler.pixel_chunks[k].startx;
           int ny_begin = adaptive_pixel_sampler.pixel_chunks[k].starty;
           int nx_end = adaptive_pixel_sampler.pixel_chunks[k].endx;
           int ny_end = adaptive_pixel_sampler.pixel_chunks[k].endy;
           
           std::vector<dielectric*> *mat_stack = new std::vector<dielectric*>;
           for(int i = nx_begin; i < nx_end; i++) {
             for(int j = ny_begin; j < ny_end; j++) {
               int index = j + ny * i;
               vec2 u2 = samplers[index]->Get2D();
               Float u = (Float(i) + u2.x()) / Float(nx);
               Float v = (Float(j) + u2.y()) / Float(ny);
               ray r; 
               if(fov != 0) {
                 r = cam.get_ray(u, v, rand_to_unit(samplers[index]->Get2D()), 
                                 samplers[index]->Get1D());
                 
               } else {
                 r = ocam.get_ray(u,v);
               }
               r.pri_stack = mat_stack;

               vec3 col = clamp(de_nan(color(r, &world, &hlist, max_depth, 
                                        roulette_active, rngs[index], samplers[index])),
                           0, clampval);
               mat_stack->clear();
               adaptive_pixel_sampler.add_color_main(i, j, col);
               if(s % 2 == 0) {
                 adaptive_pixel_sampler.add_color_sec(i, j, col);
               }
               samplers[index]->StartNextSample();
             }
           }
           if(s % 2 == 1 && s > 1) {
             adaptive_pixel_sampler.test_for_convergence(k, s, nx_end, nx_begin, ny_end, ny_begin);
           }
           delete mat_stack;
        };
        for(size_t j = 0; j < adaptive_pixel_sampler.size(); j++) {
          pool.push(worker, j);
        }
        pool.join();
        if(s % 2 == 1 && s > 1) {
          adaptive_pixel_sampler.split_remove_chunks(s);
        }
        adaptive_pixel_sampler.max_s++;
      }
      adaptive_pixel_sampler.write_final_pixels();
      for(size_t i = 0; i < samplers.size(); i++) {
        delete samplers.at(i);
      }
    }
  }

  if(verbose) {
    Rcpp::Rcout << "Cleaning up memory..." << "\n";
  }
  delete worldbvh;
  delete background_sphere;
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
  for(size_t i = 0; i < shared_materials->size(); i++) {
    delete shared_materials->at(i);
  }
  delete shared_materials;
  PutRNGstate();
  finish = std::chrono::high_resolution_clock::now();
  if(verbose) {
    std::chrono::duration<double> elapsed = finish - startfirst;
    Rcpp::Rcout << "Total time elapsed: " << elapsed.count() << " seconds" << "\n";
  }
  return(List::create(_["r"] = routput, _["g"] = goutput, _["b"] = boutput));
}


// [[Rcpp::export]]
List tonemap_image(int nx, int ny, 
                   NumericMatrix routput, NumericMatrix goutput, NumericMatrix boutput, 
                   int toneval) {
  for(int j = ny - 1; j >= 0; j--) {
    for(int i = 0; i < nx; i++) {
      if(toneval == 1) {
        routput(i,j) = std::pow(routput(i,j),1.0f/2.2f);
        goutput(i,j) = std::pow(goutput(i,j),1.0f/2.2f);
        boutput(i,j) = std::pow(boutput(i,j),1.0f/2.2f);
      } else if (toneval == 2) {
        float max = (routput(i,j)+goutput(i,j)+boutput(i,j))/3.0f;
        routput(i,j) = reinhard(routput(i,j),max);
        goutput(i,j) = reinhard(goutput(i,j),max);
        boutput(i,j) = reinhard(boutput(i,j),max);
      } else if (toneval == 3) {
        routput(i,j) = hable(routput(i,j));
        goutput(i,j) = hable(goutput(i,j));
        boutput(i,j) = hable(boutput(i,j));
      } else if (toneval == 4) {
        routput(i,j) = hbd(routput(i,j));
        goutput(i,j) = hbd(goutput(i,j));
        boutput(i,j) = hbd(boutput(i,j));
      } else {
        //do nothing
      }
    }
  }
  return(List::create(_["r"] = routput, _["g"] = goutput, _["b"] = boutput));
}
