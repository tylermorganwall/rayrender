#define STB_IMAGE_IMPLEMENTATION 

#ifdef RAY_FLOAT_AS_DOUBLE
typedef double Float;
#else
typedef float Float;
#endif 

#include "vec3.h"
#include "mathinline.h"
#include "camera.h"
#include "float.h"
#include "buildscene.h"
#include "RProgress.h"
#include "rng.h"
#include "tonemap.h"
#include "infinite_area_light.h"
#include "adaptivesampler.h"
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
#include "RcppThread.h"

using namespace std;

inline vec3 de_nan(const vec3& c) {
  vec3 temp = c;
  if(std::isnan(c[0])) temp.e[0] = 0.0f;
  if(std::isnan(c[1])) temp.e[1] = 0.0f;
  if(std::isnan(c[2])) temp.e[2] = 0.0f;
  return(temp);
}

vec3 color(const ray& r, hitable *world, hitable *hlist, int depth, bool tonemap, random_gen& rng) {
  hit_record hrec;
  if(world->hit(r, 0.001, FLT_MAX, hrec, rng)) { //generated hit record, world space
    scatter_record srec;
    vec3 emitted = hrec.mat_ptr->emitted(r, hrec, hrec.u, hrec.v, hrec.p);
    float pdf_val;
    if(depth < 50 && hrec.mat_ptr->scatter(r, hrec, srec, rng)) { //generates scatter record, world space
      if(srec.is_specular) { //returns specular ray
        // srec.specular_ray.origin() = offset_ray(srec.specular_ray.origin(), hrec.normal);
        return(srec.attenuation * 
               color(srec.specular_ray, world, 
                     hlist, depth + 1, tonemap, rng));
      }
      hitable_pdf p_imp(hlist, hrec.p); //creates pdf of all objects to be sampled
      mixture_pdf p(&p_imp, srec.pdf_ptr); //creates mixture pdf of surface intersected at hrec.p and all sampled objects/lights
      
      //Generates a scatter direction (with origin hrec.p) from the mixture 
      //and saves surface normal from light to use in pdf_value calculation
      //(along with the scatter direction)
      //Translates the world space point into object space point, generates ray assuring intersection, and then translates 
      //ray back into world space
      // ray scattered = ray(offset_ray(hrec.p,hrec.normal), p.generate(rng), r.pri_stack, r.time()); //scatters a ray from hit point to direction
      ray scattered = ray(hrec.p, p.generate(rng), r.pri_stack, r.time()); //scatters a ray from hit point to direction
      pdf_val = p.value(scattered.direction(), rng); //generates a pdf value based the intersection point and the mixture pdf
      return(emitted + srec.attenuation *
             hrec.mat_ptr->scattering_pdf(r, hrec, scattered) *
             color(scattered, world,
                  hlist, depth + 1, tonemap, rng) / pdf_val);
    } else {
      return(emitted);
    }
  } else {
    return(vec3(0,0,0));
  }
}

vec3 color_amb(const ray& r, hitable *world, hitable *hlist, int depth,
           const vec3& backgroundhigh, const vec3& backgroundlow, bool tonemap, random_gen& rng) {
  hit_record hrec;
  if(world->hit(r, 0.001, FLT_MAX, hrec, rng)) {
    scatter_record srec;
    vec3 emitted = hrec.mat_ptr->emitted(r, hrec, hrec.u, hrec.v,hrec.p);
    float pdf_val;
    if(depth < 50 && hrec.mat_ptr->scatter(r, hrec, srec, rng)) {
      if(srec.is_specular) {
        return(srec.attenuation * 
               color_amb(srec.specular_ray, world, hlist, depth + 1, 
                         backgroundhigh,backgroundlow,tonemap, rng));
      }
      hitable_pdf p_imp(hlist, hrec.p);
      mixture_pdf p(&p_imp, srec.pdf_ptr);
      ray scattered = ray(hrec.p, p.generate(rng), r.pri_stack, r.time());
      pdf_val = p.value(scattered.direction(), rng);
      return(emitted + srec.attenuation * 
             hrec.mat_ptr->scattering_pdf(r, hrec, scattered) *  
             color_amb(scattered, world, hlist, depth + 1, 
                       backgroundhigh,backgroundlow, tonemap, rng) / pdf_val);
    } else {
      vec3 unit_direction = unit_vector(r.direction());
      float t = 0.5f * (unit_direction.y() + 1.0f);
      return(emitted + (1.0f - t) * backgroundlow + t * backgroundhigh);
    }
  } else {
    vec3 unit_direction = unit_vector(r.direction());
    float t = 0.5f * (unit_direction.y() + 1.0);
    return (1.0 - t) * backgroundlow + t * backgroundhigh;
  }
}

vec3 color_uniform(const ray& r, hitable *world, int depth, bool tonemap, random_gen& rng, texture* background_texture) {
  hit_record hrec;
  if(world->hit(r, 0.001, FLT_MAX, hrec, rng)) {
    scatter_record srec;
    vec3 emitted = hrec.mat_ptr->emitted(r, hrec, hrec.u, hrec.v,hrec.p);
    float pdf_val;
    if(depth < 50 && hrec.mat_ptr->scatter(r, hrec, srec, rng)) {
      if(srec.is_specular) {
        return(srec.attenuation * color_uniform(srec.specular_ray, world, depth + 1, 
                                                tonemap, rng, background_texture));
      }
      cosine_pdf p(hrec.normal);
      ray scattered = ray(hrec.p, p.generate(rng), r.pri_stack, r.time());
      pdf_val = p.value(scattered.direction(), rng);
      return(emitted + srec.attenuation * hrec.mat_ptr->scattering_pdf(r, hrec, scattered) *  
             color_uniform(scattered, world, depth + 1, tonemap, rng, background_texture) / pdf_val);
    } else {
      return(emitted);
    }
  } else {
    vec3 unit_direction = unit_vector(r.direction());
    float phi = atan2(unit_direction.x(),unit_direction.z());
    float u = 0.5 + phi / (2*M_PI);
    float v = 0.5 * (1.0 + unit_direction.y());
    return(background_texture->value(u, v, unit_direction));
  }
}

vec3 color_amb_uniform(const ray& r, hitable *world, int depth,
               const vec3& backgroundhigh, const vec3& backgroundlow, bool tonemap, random_gen& rng) {
  hit_record hrec;
  if(world->hit(r, 0.001, FLT_MAX, hrec, rng)) {
    scatter_record srec;

    vec3 emitted = hrec.mat_ptr->emitted(r, hrec, hrec.u, hrec.v,hrec.p);
    float pdf_val;
    if(depth < 50 && hrec.mat_ptr->scatter(r, hrec, srec, rng)) {
      if(srec.is_specular) {
        return(srec.attenuation * 
               color_amb_uniform(srec.specular_ray, world, depth + 1, 
                                 backgroundhigh,backgroundlow, tonemap, rng));
      }
      cosine_pdf p(hrec.normal);
      ray scattered = ray(hrec.p, p.generate(rng), r.pri_stack, r.time());
      pdf_val = p.value(scattered.direction(), rng);
      return(emitted + srec.attenuation * 
             hrec.mat_ptr->scattering_pdf(r, hrec, scattered) *  
             color_amb_uniform(scattered, world, depth + 1, 
                               backgroundhigh,backgroundlow, tonemap, rng) / pdf_val);
    } else {
      vec3 unit_direction = unit_vector(r.direction());
      float t = 0.5 * (unit_direction.y() + 1.0);
      return(emitted+(1.0 - t) * backgroundlow + t * backgroundhigh);
    }
  } else {
    vec3 unit_direction = unit_vector(r.direction());
    float t = 0.5 * (unit_direction.y() + 1.0);
    return (1.0 - t) * backgroundlow + t * backgroundhigh;
  }
}

#ifdef DEBUGBVH
float debug_bvh(const ray& r, hitable *world, random_gen rng) {
  hit_record hrec;
  hrec.bvh_nodes = 0;
  world->hit(r, 0.000.01, FLT_MAX, hrec, rng);
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
                       float rotate_env, bool verbose, int debug_channel,
                       IntegerVector& shared_id_mat, LogicalVector& is_shared_mat,
                       float min_variance, int min_adaptive_size) {
  auto startfirst = std::chrono::high_resolution_clock::now();
  //Unpack Camera Info
  int nx = as<int>(camera_info["nx"]);
  int ny = as<int>(camera_info["ny"]);
  int ns = as<int>(camera_info["ns"]);
  Float fov = as<int>(camera_info["fov"]);
  NumericVector lookfromvec = as<NumericVector>(camera_info["lookfrom"]);
  NumericVector lookatvec = as<NumericVector>(camera_info["lookat"]);
  Float aperture = as<Float>(camera_info["aperture"]);
  NumericVector camera_up = as<NumericVector>(camera_info["camera_up"]);
  Float shutteropen = as<Float>(camera_info["shutteropen"]);
  Float shutterclose = as<Float>(camera_info["shutterclose"]);
  Float focus_distance = as<Float>(camera_info["focal_distance"]);
  int toneval = as<Float>(camera_info["toneval"]);
  NumericVector ortho_dimensions = as<NumericVector>(camera_info["ortho_dimensions"]);
  
  //Initialize output matrices
  NumericMatrix routput(nx,ny);
  NumericMatrix goutput(nx,ny);
  NumericMatrix boutput(nx,ny);
  
  vec3 lookfrom(lookfromvec[0],lookfromvec[1],lookfromvec[2]);
  vec3 lookat(lookatvec[0],lookatvec[1],lookatvec[2]);
  vec3 backgroundhigh(bghigh[0],bghigh[1],bghigh[2]);
  vec3 backgroundlow(bglow[0],bglow[1],bglow[2]);
  float dist_to_focus = focus_distance;
  bool tonemap = toneval == 1 ? false : true;
  CharacterVector alpha_files = as<CharacterVector>(alphalist["alpha_temp_file_names"]);
  LogicalVector has_alpha = as<LogicalVector>(alphalist["alpha_tex_bool"]);
  RcppThread::ThreadPool pool(numbercores);
  GetRNGstate();
  random_gen rng(unif_rand() * std::pow(2,32));
  camera cam(lookfrom, lookat, vec3(camera_up(0),camera_up(1),camera_up(2)), fov, float(nx)/float(ny), 
             aperture, dist_to_focus,
             shutteropen, shutterclose, rng);
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
  }
  
  
  hitable *worldbvh = build_scene(type, radius, shape, position_list,
                                properties, velocity, moving,
                                n,shutteropen,shutterclose,
                                ischeckered, checkercolors, 
                                gradient_info,
                                noise, isnoise, noisephase, noiseintensity, noisecolorlist,
                                angle, 
                                isimage, has_alpha, alpha_textures, nx_ny_nn_alpha,
                                textures, nx_ny_nn,
                                lightintensity, isflipped,
                                isvolume, voldensity, order_rotation_list, 
                                isgrouped, group_pivot, group_translate,
                                group_angle, group_order_rotation, group_scale,
                                tri_normal_bools, is_tri_color, tri_color_vert, 
                                fileinfo, filebasedir, 
                                scale_list, sigmavec, 
                                shared_id_mat, is_shared_mat, shared_materials, rng);
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
    background_texture = new image_texture(background_texture_data, nx1, ny1, nn1);
    background_material = new diffuse_light(background_texture);
    background_sphere = new InfiniteAreaLight(nx1, ny1, world_radius*2, world_center,
                                              background_texture, background_material);
    if(rotate_env != 0) {
      background_sphere = new rotate_y(background_sphere, rotate_env);
    }
  } else {
    background_texture = new constant_texture(vec3(0,0,0));
    background_material = new lambertian(background_texture);
    background_sphere = new sphere(world_center, world_radius, background_material, nullptr);
  }
  finish = std::chrono::high_resolution_clock::now();
  if(verbose && hasbackground) {
    std::chrono::duration<double> elapsed = finish - start;
    Rcpp::Rcout << elapsed.count() << " seconds" << "\n";
  }
  hitable *world_full[2];
  world_full[0] = worldbvh;
  world_full[1] = background_sphere;
  int nval = hasbackground ? 2 : 1;
  hitable_list world(world_full, nval);
  int numbertosample = 0;
  for(int i = 0; i < implicit_sample.size(); i++) {
    if(implicit_sample(i)) {
      numbertosample++;
    }
  }
  if(hasbackground) {
    numbertosample++;
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
  if(hasbackground) {
    implicit_sample_vector[counter] = background_sphere;
  }
  

  hitable_list hlist;
  if(!implicit_sample_vector.empty()) {
    hlist = hitable_list(&implicit_sample_vector[0], numbertosample);
  } else {
    numbertosample = 0;
  }
  if(verbose && !progress_bar) {
    Rcpp::Rcout << "Starting Raytracing:\n ";
  }
  RProgress::RProgress pb("Adaptive Raytracing [:bar] :percent%");
  pb.set_width(70);

  if(progress_bar) {
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
          r = cam.get_ray(u,v);
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
          r = cam.get_ray(u,v);
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
          r = cam.get_ray(u,v);
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
          r = cam.get_ray(u,v);
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
      RcppThread::ThreadPool pool(numbercores);
      auto worker = [&routput, &goutput, &boutput,
                     ambient_light, nx, ny, ns, seeds, fov,
                     &cam, &ocam, backgroundhigh, backgroundlow, &world, &hlist,
                     numbertosample, clampval, tonemap, progress_bar, 
                     numbercores, background_texture] (int j) {
       if(progress_bar && j % numbercores == 0) {
         RcppThread::Rcout << "Progress (" << numbercores << " core): ";
         RcppThread::Rcout << (int)((1-(double)j/double(ny)) * 100) << "%\r";
       }
       random_gen rng(seeds[j]);
       std::vector<dielectric*> *mat_stack = new std::vector<dielectric*>;
       
       for(int i = 0; i < nx; i++) {
         vec3 col(0,0,0);
         for(int s = 0; s < ns; s++) {
           Float u = Float(i + rng.unif_rand()) / Float(nx);
           Float v = Float(j + rng.unif_rand()) / Float(ny);
           ray r;
           if(fov != 0) {
             r = cam.get_ray(u,v);
           } else {
             r = ocam.get_ray(u,v);
           }
           r.pri_stack = mat_stack;
           
           if(numbertosample) {
             if(ambient_light) {
               col += clamp(de_nan(color_amb(r, &world, &hlist, 0,
                                             backgroundhigh, backgroundlow, tonemap, rng)),0,clampval);
             } else {
               col += clamp(de_nan(color(r, &world, &hlist, 0, tonemap, rng)),0,clampval);
             }
           } else {
             if(ambient_light) {
               col += clamp(de_nan(color_amb_uniform(r, &world, 0, 
                                                     backgroundhigh, backgroundlow, tonemap, rng)),0,clampval);
             } else {
               col += clamp(de_nan(color_uniform(r, &world, 0, tonemap, rng, background_texture)),0,clampval);
             }
           }
           mat_stack->clear();
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
    } else {
      NumericMatrix routput2(nx,ny);
      NumericMatrix goutput2(nx,ny);
      NumericMatrix boutput2(nx,ny);
      adaptive_sampler adaptive_pixel_sampler(numbercores, nx, ny, ns, debug_channel,
                                              min_variance, min_adaptive_size,
                                              routput, goutput, boutput,
                                              routput2, goutput2, boutput2);
      for(int s = 0; s < ns; s++) {
        Rcpp::checkUserInterrupt();
        if(progress_bar) {
          pb.tick();
        }
        std::vector<random_gen > rngs;
        for(int i = 0; i < nx * ny; i++) {
          random_gen rng_single(unif_rand() * std::pow(2,32));
          rngs.push_back(rng_single);
        }
        RcppThread::ThreadPool pool(numbercores);
        auto worker = [&adaptive_pixel_sampler, 
                       ambient_light, nx, ny, s, 
                       &rngs, fov, 
                       &cam, &ocam, backgroundhigh, backgroundlow, &world, &hlist,
                       numbertosample, clampval, tonemap, 
                       background_texture] (int k) {
           int nx_begin = adaptive_pixel_sampler.pixel_chunks[k].startx;
           int ny_begin = adaptive_pixel_sampler.pixel_chunks[k].starty;
           int nx_end = adaptive_pixel_sampler.pixel_chunks[k].endx;
           int ny_end = adaptive_pixel_sampler.pixel_chunks[k].endy;
           random_gen rng = rngs[k];
           std::vector<dielectric*> *mat_stack = new std::vector<dielectric*>;
           for(int i = nx_begin; i < nx_end; i++) {
             for(int j = ny_begin; j < ny_end; j++) {
               vec3 col(0,0,0);
               Float u = Float(i + rng.unif_rand()) / Float(nx);
               Float v = Float(j + rng.unif_rand()) / Float(ny);
               ray r;
               if(fov != 0) {
                 r = cam.get_ray(u,v);
               } else {
                 r = ocam.get_ray(u,v);
               }
               r.pri_stack = mat_stack;
               if(numbertosample) {
                 if(ambient_light) {
                   col = clamp(de_nan(color_amb(r, &world, &hlist, 0,
                                                 backgroundhigh, backgroundlow, tonemap, rng)),0,clampval);
                 } else {
                   col = clamp(de_nan(color(r, &world, &hlist, 0, tonemap, rng)),0,clampval);
                 }
               } else {
                 if(ambient_light) {
                   col = clamp(de_nan(color_amb_uniform(r, &world, 0, 
                                                         backgroundhigh, backgroundlow, tonemap, rng)),0,clampval);
                 } else {
                   col = clamp(de_nan(color_uniform(r, &world, 0, tonemap, rng, background_texture)),0,clampval);
                 }
               }
               mat_stack->clear();
               adaptive_pixel_sampler.add_color_main(i, j, col);
               if(s % 2 == 0) {
                 adaptive_pixel_sampler.add_color_sec(i, j, col);
               }
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
