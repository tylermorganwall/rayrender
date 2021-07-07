#include "debug.h"

void debug_scene(size_t numbercores, size_t nx, size_t ny, size_t ns, int debug_channel,
                Float min_variance, size_t min_adaptive_size, 
                Rcpp::NumericMatrix& routput, Rcpp::NumericMatrix& goutput, Rcpp::NumericMatrix& boutput,
                bool progress_bar, int sample_method, Rcpp::NumericVector& stratified_dim,
                bool verbose, ortho_camera& ocam, camera &cam, environment_camera &ecam, Float fov,
                hitable_list& world, hitable_list& hlist,
                Float clampval, size_t max_depth, size_t roulette_active,
                Rcpp::NumericVector& light_direction, random_gen& rng) {
  if(debug_channel == 1) {
    Float depth_into_scene = 0.0;
    for(unsigned int j = 0; j < ny; j++) {
      for(unsigned int i = 0; i < nx; i++) {
        depth_into_scene = 0;
        Float u = Float(i) / Float(nx);
        Float v = Float(j) / Float(ny);
        ray r;
        if(fov != 0 && fov != 360) {
          r = cam.get_ray(u,v, vec3f(0,0,0), 0);
        } else if (fov == 0){
          r = ocam.get_ray(u,v, rng.unif_rand());
        } else {
          r = ecam.get_ray(u,v, rng.unif_rand());
        }
        depth_into_scene = calculate_depth(r, &world, rng);
        routput(i,j) = depth_into_scene;
        goutput(i,j) = depth_into_scene;
        boutput(i,j) = depth_into_scene;
      }
    }
  } else if(debug_channel == 2) {
    vec3f normal_map(0,0,0);
    for(unsigned int j = 0; j < ny; j++) {
      for(unsigned int i = 0; i < nx; i++) {
        normal_map = vec3f(0,0,0);
        Float u = Float(i) / Float(nx);
        Float v = Float(j) / Float(ny);
        ray r;
        if(fov != 0 && fov != 360) {
          r = cam.get_ray(u,v, vec3f(0,0,0), 0);
        } else if (fov == 0){
          r = ocam.get_ray(u,v, rng.unif_rand());
        } else {
          r = ecam.get_ray(u,v, rng.unif_rand());
        }
        normal_map = calculate_normals(r, &world, rng);
        routput(i,j) = normal_map.x();
        goutput(i,j) = normal_map.y();
        boutput(i,j) = normal_map.z();
      }
    }
  } else if(debug_channel == 3) {
    vec3f uv_map(0,0,0);
    for(unsigned int j = 0; j < ny; j++) {
      for(unsigned int i = 0; i < nx; i++) {
        uv_map = vec3f(0,0,0);
        Float u = Float(i) / Float(nx);
        Float v = Float(j) / Float(ny);
        ray r;
        if(fov != 0 && fov != 360) {
          r = cam.get_ray(u,v, vec3f(0,0,0), 0);
        } else if (fov == 0){
          r = ocam.get_ray(u,v, rng.unif_rand());
        } else {
          r = ecam.get_ray(u,v, rng.unif_rand());
        }
        uv_map = calculate_uv(r, &world, rng);
        routput(i,j) = uv_map.x();
        goutput(i,j) = uv_map.y();
        boutput(i,j) = uv_map.z();
      }
    }
  } else if(debug_channel == 4) {
  #ifdef DEBUGBVH
    for(unsigned int j = 0; j < ny; j++) {
      for(unsigned int i = 0; i < nx; i++) {
        Float u = Float(i) / Float(nx);
        Float v = Float(j) / Float(ny);
        ray r;
        if(fov != 0 && fov != 360) {
          r = cam.get_ray(u,v, vec3f(0,0,0), 0);
        } else if (fov == 0){
          r = ocam.get_ray(u,v, rng.unif_rand());
        } else {
          r = ecam.get_ray(u,v, rng.unif_rand());
        }
        double bvh_intersections = debug_bvh(r, &world, rng);
        routput(i,j) = bvh_intersections;
        goutput(i,j) = bvh_intersections;
        boutput(i,j) = bvh_intersections;
      }
    }
  #endif
  } else if (debug_channel == 6 || debug_channel == 7) {
    for(unsigned int j = 0; j < ny; j++) {
      for(unsigned int i = 0; i < nx; i++) {
        Float u = Float(i) / Float(nx);
        Float v = Float(j) / Float(ny);
        ray r;
        if(fov != 0 && fov != 360) {
          r = cam.get_ray(u,v, vec3f(0,0,0), 0);
        } else if (fov == 0){
          r = ocam.get_ray(u,v, rng.unif_rand());
        } else {
          r = ecam.get_ray(u,v, rng.unif_rand());
        }
        vec3f dpd_val = calculate_dpduv(r, &world, rng, debug_channel == 6);
        routput(i,j) = dpd_val.x();
        goutput(i,j) = dpd_val.y();
        boutput(i,j) = dpd_val.z();
      }
    }
  } else if (debug_channel == 8) {
    std::vector<dielectric*> *mat_stack = new std::vector<dielectric*>;
    
    for(unsigned int j = 0; j < ny; j++) {
      for(unsigned int i = 0; i < nx; i++) {
        Float u = Float(i) / Float(nx);
        Float v = Float(j) / Float(ny);
        ray r;
        if(fov != 0 && fov != 360) {
          r = cam.get_ray(u,v, vec3f(0,0,0), 0);
        } else if (fov == 0){
          r = ocam.get_ray(u,v, rng.unif_rand());
        } else {
          r = ecam.get_ray(u,v, rng.unif_rand());
        }
        r.pri_stack = mat_stack;
        vec3f dpd_val = calculate_color(r, &world, rng);
        mat_stack->clear();
        
        routput(i,j) = dpd_val.x();
        goutput(i,j) = dpd_val.y();
        boutput(i,j) = dpd_val.z();
      }
    }
    delete mat_stack;
  } else if (debug_channel == 9) {
    vec3f light_dir(light_direction(0),light_direction(1),light_direction(2));
    Float n_exp = light_direction(3);
    RcppThread::ThreadPool pool(numbercores);
    auto worker = [&routput, &goutput, &boutput,
                   nx, ny,  fov, light_dir, n_exp,
                   &cam, &ocam, &ecam, &world] (int j) {
                     std::vector<dielectric*> *mat_stack = new std::vector<dielectric*>;
                     random_gen rng(j);
                     for(unsigned int i = 0; i < nx; i++) {
                       Float u = Float(i) / Float(nx);
                       Float v = Float(j) / Float(ny);
                       ray r;
                       if(fov != 0 && fov != 360) {
                         r = cam.get_ray(u,v, vec3f(0,0,0), 0);
                       } else if (fov == 0){
                         r = ocam.get_ray(u,v, rng.unif_rand());
                       } else {
                         r = ecam.get_ray(u,v, rng.unif_rand());
                       }
                       r.pri_stack = mat_stack;
                       vec3f qr = quick_render(r, &world, rng, light_dir, n_exp);
                       mat_stack->clear();
                       
                       routput(i,j) = qr.x();
                       goutput(i,j) = qr.y();
                       boutput(i,j) = qr.z();
                     }
                     delete mat_stack;
                   };
    for(int j = ny - 1; j >= 0; j--) {
      pool.push(worker,j);
    }
    pool.join();
  }
}
