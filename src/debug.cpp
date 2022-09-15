#include "debug.h"
#include "RProgress.h"
#include "adaptivesampler.h"

void debug_scene(size_t numbercores, size_t nx, size_t ny, size_t ns, int debug_channel,
                Float min_variance, size_t min_adaptive_size, 
                RayMatrix& routput, RayMatrix& goutput, RayMatrix& boutput,
                bool progress_bar, int sample_method, Rcpp::NumericVector& stratified_dim,
                bool verbose, 
                RayCamera* cam, 
                Float fov,
                hitable_list& world, hitable_list& hlist,
                Float clampval, size_t max_depth, size_t roulette_active,
                Rcpp::NumericVector& light_direction, random_gen& rng, Float sample_dist,
                bool keep_colors, vec3f backgroundhigh) {
  Environment pkg = Environment::namespace_env("rayrender");
  Function print_time = pkg["print_time"];
  if(debug_channel == 1) {
    Float depth_into_scene = 0.0;
    for(unsigned int j = 0; j < ny; j++) {
      for(unsigned int i = 0; i < nx; i++) {
        depth_into_scene = 0;
        ray r;
        if(fov >= 0) {
          Float u = (Float(i)) / Float(nx);
          Float v = (Float(j)) / Float(ny);
          r = cam->get_ray(u,v, vec3f(0,0,0), rng.unif_rand());
        } else {
          point2f u(rng.unif_rand(),rng.unif_rand());
          point2f u2(rng.unif_rand(),rng.unif_rand());
          CameraSample samp(u, u2, rng.unif_rand());
          cam->GenerateRay(samp, &r);
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
        std::vector<dielectric*> *mat_stack = new std::vector<dielectric*>;
        
        normal_map = vec3f(0,0,0);
        ray r;
        if(fov >= 0) {
          Float u = (Float(i)) / Float(nx);
          Float v = (Float(j)) / Float(ny);
          r = cam->get_ray(u,v, vec3f(0,0,0), rng.unif_rand());
        } else {
          point2f u(rng.unif_rand(),rng.unif_rand());
          point2f u2(rng.unif_rand(),rng.unif_rand());
          CameraSample samp(u, u2, rng.unif_rand());
          cam->GenerateRay(samp, &r);
        }
        r.pri_stack = mat_stack;
        normal_map = calculate_normals(r, &world, max_depth, rng);
        routput(i,j) = normal_map.x();
        goutput(i,j) = normal_map.y();
        boutput(i,j) = normal_map.z();
        delete mat_stack;
      }
    }
  } else if(debug_channel == 3) {
    vec3f uv_map(0,0,0);
    for(unsigned int j = 0; j < ny; j++) {
      for(unsigned int i = 0; i < nx; i++) {
        uv_map = vec3f(0,0,0);
        ray r;
        if(fov >= 0) {
          Float u = (Float(i)) / Float(nx);
          Float v = (Float(j)) / Float(ny);
          r = cam->get_ray(u,v, vec3f(0,0,0), rng.unif_rand());
        } else {
          point2f u(rng.unif_rand(),rng.unif_rand());
          point2f u2(rng.unif_rand(),rng.unif_rand());
          CameraSample samp(u, u2, rng.unif_rand());
          cam->GenerateRay(samp, &r);
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
        if(fov >= 0) {
          Float u = (Float(i)) / Float(nx);
          Float v = (Float(j)) / Float(ny);
          r = cam->get_ray(u,v, vec3f(0,0,0), rng.unif_rand());
        } else {
          point2f u(rng.unif_rand(),rng.unif_rand());
          point2f u2(rng.unif_rand(),rng.unif_rand());
          CameraSample samp(u, u2, rng.unif_rand());
          cam->GenerateRay(samp, &r);
        }
        Float bvh_intersections = debug_bvh(r, &world, rng);
        routput(i,j) = bvh_intersections;
        goutput(i,j) = bvh_intersections;
        boutput(i,j) = bvh_intersections;
      }
    }
  #endif
  } else if (debug_channel == 6 || debug_channel == 7) {
    for(unsigned int j = 0; j < ny; j++) {
      for(unsigned int i = 0; i < nx; i++) {
        ray r;
        if(fov >= 0) {
          Float u = (Float(i)) / Float(nx);
          Float v = (Float(j)) / Float(ny);
          r = cam->get_ray(u,v, vec3f(0,0,0), rng.unif_rand());
        } else {
          point2f u(rng.unif_rand(),rng.unif_rand());
          point2f u2(rng.unif_rand(),rng.unif_rand());
          CameraSample samp(u, u2, rng.unif_rand());
          cam->GenerateRay(samp, &r);
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
        ray r;
        if(fov >= 0) {
          Float u = (Float(i)) / Float(nx);
          Float v = (Float(j)) / Float(ny);
          r = cam->get_ray(u,v, vec3f(0,0,0), rng.unif_rand());
        } else {
          point2f u(rng.unif_rand(),rng.unif_rand());
          point2f u2(rng.unif_rand(),rng.unif_rand());
          CameraSample samp(u, u2, rng.unif_rand());
          cam->GenerateRay(samp, &r);
        }
        r.pri_stack = mat_stack;
        point3f dpd_val = calculate_color(r, &world, rng);
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
                   cam, &world] (int j) {
                     std::vector<dielectric*> *mat_stack = new std::vector<dielectric*>;
                     random_gen rng(j);
                     for(unsigned int i = 0; i < nx; i++) {
                       ray r;
                       if(fov >= 0) {
                         Float u = (Float(i)) / Float(nx);
                         Float v = (Float(j)) / Float(ny);
                         r = cam->get_ray(u,v, vec3f(0,0,0), rng.unif_rand());
                       } else {
                         point2f u(rng.unif_rand(),rng.unif_rand());
                         point2f u2(rng.unif_rand(),rng.unif_rand());
                         CameraSample samp(u, u2, rng.unif_rand());
                         cam->GenerateRay(samp, &r);
                       }
                       r.pri_stack = mat_stack;
                       point3f qr = quick_render(r, &world, rng, light_dir, n_exp);
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
  } else if (debug_channel == 10) {
    RcppThread::ThreadPool pool(numbercores);
    auto worker = [&routput, &goutput, &boutput,
                   nx, ny,  fov, max_depth, &hlist,
                   cam, &world] (int j) {
                     std::vector<dielectric*> *mat_stack = new std::vector<dielectric*>;
                     random_gen rng(j);
                     for(unsigned int i = 0; i < nx; i++) {
                       ray r;
                       if(fov >= 0) {
                         Float u = (Float(i)) / Float(nx);
                         Float v = (Float(j)) / Float(ny);
                         r = cam->get_ray(u,v, vec3f(0,0,0), rng.unif_rand());
                       } else {
                         point2f u(rng.unif_rand(),rng.unif_rand());
                         point2f u2(rng.unif_rand(),rng.unif_rand());
                         CameraSample samp(u, u2, rng.unif_rand());
                         cam->GenerateRay(samp, &r);
                       }
                       r.pri_stack = mat_stack;
                       point3f qr = calculate_position(r, &world, &hlist, max_depth, rng);
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
  } else if (debug_channel == 11) {
    RcppThread::ThreadPool pool(numbercores);
    
    auto worker = [&routput, &goutput, &boutput, 
                   nx, ny,  fov, &hlist, max_depth,ns,
                   cam, &world] (int j) {
                     std::vector<dielectric*> *mat_stack = new std::vector<dielectric*>;
                     random_gen rng(j);
                     for(unsigned int i = 0; i < nx; i++) {
                       for(size_t s = 0; s < static_cast<size_t>(ns); s++) {
                         ray r;
                         if(fov >= 0) {
                           Float u = (Float(i)) / Float(nx);
                           Float v = (Float(j)) / Float(ny);
                           r = cam->get_ray(u,v, vec3f(0,0,0), rng.unif_rand());
                         } else {
                           point2f u(rng.unif_rand(),rng.unif_rand());
                           point2f u2(rng.unif_rand(),rng.unif_rand());
                           CameraSample samp(u, u2, rng.unif_rand());
                           cam->GenerateRay(samp, &r);
                         }
                         r.pri_stack = mat_stack;
                         point3f qr = calculate_bounce_dir(r, &world, &hlist,
                                                           max_depth, rng);
                         mat_stack->clear();
                         
                         routput(i,j) += qr.x()/ns;
                         goutput(i,j) += qr.y()/ns;
                         boutput(i,j) += qr.z()/ns;
                       }
                     }
                     delete mat_stack;
                   };
    for(int j = ny - 1; j >= 0; j--) {
      pool.push(worker,j);
    }
    pool.join();
  } else if (debug_channel == 12) {
    RcppThread::ThreadPool pool(numbercores);
    
    auto worker = [&routput, &goutput, &boutput, 
                   nx, ny,  fov, &hlist, max_depth,
                   cam, &world] (int j) {
                     std::vector<dielectric*> *mat_stack = new std::vector<dielectric*>;
                     random_gen rng(j);
                     for(unsigned int i = 0; i < nx; i++) {
                       ray r;
                       if(fov >= 0) {
                         Float u = (Float(i)) / Float(nx);
                         Float v = (Float(j)) / Float(ny);
                         r = cam->get_ray(u,v, vec3f(0,0,0), rng.unif_rand());
                       } else {
                         point2f u(rng.unif_rand(),rng.unif_rand());
                         point2f u2(rng.unif_rand(),rng.unif_rand());
                         CameraSample samp(u, u2, rng.unif_rand());
                         cam->GenerateRay(samp, &r);
                       }
                       r.pri_stack = mat_stack;
                       Float qr = calculate_time(r, &world, &hlist,
                                                         max_depth, rng);
                       mat_stack->clear();
                       
                       routput(i,j) = qr;
                       goutput(i,j) = qr;
                       boutput(i,j) = qr;
                     }
                     delete mat_stack;
                   };
    for(int j = ny - 1; j >= 0; j--) {
      pool.push(worker,j);
    }
    pool.join();
  } else if (debug_channel == 13) {
    RcppThread::ThreadPool pool(numbercores);
    
    auto worker = [&routput, &goutput, &boutput, 
                   nx, ny,  fov, &hlist, max_depth,
                   cam, &world] (int j) {
                     std::vector<dielectric*> *mat_stack = new std::vector<dielectric*>;
                     random_gen rng(j);
                     for(unsigned int i = 0; i < nx; i++) {
                       ray r;
                       if(fov >= 0) {
                         Float u = (Float(i)) / Float(nx);
                         Float v = (Float(j)) / Float(ny);
                         r = cam->get_ray(u,v, vec3f(0,0,0), rng.unif_rand());
                       } else {
                         point2f u(rng.unif_rand(),rng.unif_rand());
                         point2f u2(rng.unif_rand(),rng.unif_rand());
                         CameraSample samp(u, u2, rng.unif_rand());
                         cam->GenerateRay(samp, &r);
                       }
                       r.pri_stack = mat_stack;
                       point3f qr = calculate_shape(r, &world, &hlist,
                                                    max_depth, rng);
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
  } else if (debug_channel == 14) {
    RcppThread::ThreadPool pool(numbercores);
    
    auto worker = [&routput, &goutput, &boutput,
                   nx, ny,  fov, &hlist, max_depth,ns,
                   cam,&world] (int j) {
                     std::vector<dielectric*> *mat_stack = new std::vector<dielectric*>;
                     random_gen rng(j);
                     for(unsigned int i = 0; i < nx; i++) {
                       for(size_t s = 0; s < static_cast<size_t>(ns); s++) {
                         ray r;
                         if(fov >= 0) {
                           Float u = (Float(i)) / Float(nx);
                           Float v = (Float(j)) / Float(ny);
                           r = cam->get_ray(u,v, vec3f(0,0,0), rng.unif_rand());
                         } else {
                           point2f u(rng.unif_rand(),rng.unif_rand());
                           point2f u2(rng.unif_rand(),rng.unif_rand());
                           CameraSample samp(u, u2, rng.unif_rand());
                           cam->GenerateRay(samp, &r);
                         }
                         r.pri_stack = mat_stack;
                         Float qr = calculate_pdf(r, &world, &hlist,
                                                      max_depth, rng);
                         mat_stack->clear();
                         
                         routput(i,j) += qr/(Float)ns;
                         goutput(i,j) += qr/(Float)ns;
                         boutput(i,j) += qr/(Float)ns;
                       }
                     }
                     delete mat_stack;
                   };
    for(int j = ny - 1; j >= 0; j--) {
      pool.push(worker,j);
    }
    pool.join();
  } else if (debug_channel == 15) {
    RcppThread::ThreadPool pool(numbercores);
    
    auto worker = [&routput, &goutput, &boutput,
                   nx, ny,  fov, &hlist, max_depth, 
                   cam, &world] (int j) {
                     std::vector<dielectric*> *mat_stack = new std::vector<dielectric*>;
                     random_gen rng(j);
                     for(unsigned int i = 0; i < nx; i++) {
                       ray r;
                       if(fov >= 0) {
                         Float u = (Float(i)) / Float(nx);
                         Float v = (Float(j)) / Float(ny);
                         r = cam->get_ray(u,v, vec3f(0,0,0), rng.unif_rand());
                       } else {
                         point2f u(rng.unif_rand(),rng.unif_rand());
                         point2f u2(rng.unif_rand(),rng.unif_rand());
                         CameraSample samp(u, u2, rng.unif_rand());
                         cam->GenerateRay(samp, &r);
                       }
                       r.pri_stack = mat_stack;
                       Float qr = calculate_error(r, &world, &hlist,
                                                max_depth, rng);
                       mat_stack->clear();
                       
                       routput(i,j) = qr;
                       goutput(i,j) = qr;
                       boutput(i,j) = qr;
                     }
                     delete mat_stack;
                   };
    for(int j = ny - 1; j >= 0; j--) {
      pool.push(worker,j);
    }
    pool.join();
  } else if (debug_channel == 16) {
    RcppThread::ThreadPool pool(numbercores);
    
    auto worker = [&routput, &goutput, &boutput, 
                   nx, ny,  fov, &hlist, max_depth, ns,
                   cam,&world] (int j) {
                     std::vector<dielectric*> *mat_stack = new std::vector<dielectric*>;
                     random_gen rng(j);
                     for(unsigned int i = 0; i < nx; i++) {
                       for(size_t s = 0; s < static_cast<size_t>(ns); s++) {
                         ray r;
                         if(fov >= 0) {
                           Float u = (Float(i)) / Float(nx);
                           Float v = (Float(j)) / Float(ny);
                           r = cam->get_ray(u,v, vec3f(0,0,0), rng.unif_rand());
                         } else {
                           point2f u(rng.unif_rand(),rng.unif_rand());
                           point2f u2(rng.unif_rand(),rng.unif_rand());
                           CameraSample samp(u, u2, rng.unif_rand());
                           cam->GenerateRay(samp, &r);
                         }
                         r.pri_stack = mat_stack;
                         Float qr = calculate_bounces(r, &world, &hlist,
                                                    max_depth, rng);
                         mat_stack->clear();
                         
                         routput(i,j) += qr/ns;
                         goutput(i,j) += qr/ns;
                         boutput(i,j) += qr/ns;
                       }
                     }
                     delete mat_stack;
                   };
    for(int j = ny - 1; j >= 0; j--) {
      pool.push(worker,j);
    }
    pool.join();
  } else if (debug_channel == 17) {
    RcppThread::ThreadPool pool(numbercores);
    
    auto worker = [&routput, &goutput, &boutput, 
                   nx, ny,  fov, ns,
                   cam] (int j) {
                     random_gen rng(j);
                     for(unsigned int i = 0; i < nx; i++) {
                       for(size_t s = 0; s < static_cast<size_t>(ns); s++) {
                         ray r;
                         Float weight = 1;
                         if(fov >= 0) {
                           Float u = (Float(i)) / Float(nx);
                           Float v = (Float(j)) / Float(ny);
                           r = cam->get_ray(u,v, vec3f(0,0,0), rng.unif_rand());
                         } else {
                           point2f u(rng.unif_rand(),rng.unif_rand());
                           point2f u2(rng.unif_rand(),rng.unif_rand());
                           CameraSample samp(u, u2, rng.unif_rand());
                           weight = cam->GenerateRay(samp, &r);
                         }
                         vec3f v2 = unit_vector(r.direction());
                         
                         
                         routput(i,j) = weight != 0 ? (v2.x()+1.f)/(2.f) : 0;
                         goutput(i,j) = weight != 0 ? (v2.y()+1.f)/(2.f) : 0;
                         boutput(i,j) = weight != 0 ? (v2.z()+1.f)/(2.f) : 0;
                       }
                     }
                   };
    for(int j = ny - 1; j >= 0; j--) {
      pool.push(worker,j);
    }
    pool.join();
  } else if (debug_channel == 18) {
    RProgress::RProgress pb_sampler("Generating Samples [:bar] :percent%");
    pb_sampler.set_width(70);
    RProgress::RProgress pb("Adaptive AO [:bar] :percent%");
    pb.set_width(70);
    
    if(progress_bar) {
      pb_sampler.set_total(ny);
      pb.set_total(ns);
    }
    RayMatrix routput2(nx,ny);
    RayMatrix goutput2(nx,ny);
    RayMatrix boutput2(nx,ny);
    adaptive_sampler adaptive_pixel_sampler(numbercores, nx, ny, ns, debug_channel,
                                            min_variance, min_adaptive_size,
                                            routput, goutput, boutput,
                                            routput2, goutput2, boutput2);
    std::vector<random_gen > rngs;
    std::vector<std::unique_ptr<Sampler> > samplers;

    for(unsigned int j = 0; j < ny; j++) {
      if(progress_bar) {
        pb_sampler.tick();
      }
      for(unsigned int i = 0; i < nx; i++) {
        random_gen rng_single(unif_rand() * std::pow(2,32));
        rngs.push_back(rng_single);
        if(sample_method == 0) {
          samplers.push_back(std::unique_ptr<Sampler>(new RandomSampler(rng_single)));
          samplers.back()->StartPixel(0,0);
        } else if (sample_method == 1){
          samplers.push_back(std::unique_ptr<Sampler>(new StratifiedSampler(stratified_dim(0), stratified_dim(1),
                                                                            true, 5, rng_single)));
          samplers.back()->StartPixel(0,0);
        } else if (sample_method == 2) {
          samplers.push_back(std::unique_ptr<Sampler>(new SobolSampler(ns, rng_single)));
          samplers.back()->StartPixel(i,j);
        } else {
          samplers.push_back(std::unique_ptr<Sampler>(new SobolBlueNoiseSampler(rng_single)));
          samplers.back()->StartPixel(i,j);
        }
        samplers.back()->SetSampleNumber(0);
      }
    }
    
    print_time(verbose, "Allocating sampler");
    for(size_t s = 0; s < static_cast<size_t>(ns); s++) {
      Rcpp::checkUserInterrupt();
      if(progress_bar) {
        pb.tick();
      }
      RcppThread::ThreadPool pool(numbercores);
      auto worker = [&adaptive_pixel_sampler,
                     nx, ny, s, sample_method,
                     &rngs, fov, &samplers,
                     cam, &world, &hlist, backgroundhigh, 
                     sample_dist, keep_colors] (int k) {
                       int nx_begin = adaptive_pixel_sampler.pixel_chunks[k].startx;
                       int ny_begin = adaptive_pixel_sampler.pixel_chunks[k].starty;
                       int nx_end = adaptive_pixel_sampler.pixel_chunks[k].endx;
                       int ny_end = adaptive_pixel_sampler.pixel_chunks[k].endy;
                       
                       for(int i = nx_begin; i < nx_end; i++) {
                         for(int j = ny_begin; j < ny_end; j++) {
                           int index = j + ny * i;
                           ray r; 
                           vec2f u2 = samplers[index]->Get2D();
                           Float u = (Float(i) + u2.x()) / Float(nx);
                           Float v = (Float(j) + u2.y()) / Float(ny);
                           if(fov >= 0) {
                             Float u = (Float(i)) / Float(nx);
                             Float v = (Float(j)) / Float(ny);
                             r = cam->get_ray(u,v, rand_to_unit(samplers[index]->Get2D()),
                                              samplers[index]->Get1D());
                           } else {
                             CameraSample samp({1-u,1-v},samplers[index]->Get2D(), samplers[index]->Get1D());
                             cam->GenerateRay(samp, &r);
                           }
                           
                           point3f col = clamp_point(calculate_ao(r, &world, &hlist,
                                                      sample_dist, rngs[index], samplers[index].get(),
                                                      keep_colors, backgroundhigh),0.f,1.f);
                           adaptive_pixel_sampler.add_color_main(i, j, col);
                           if(s % 2 == 0) {
                             adaptive_pixel_sampler.add_color_sec(i, j, col);
                           }
                           samplers[index]->StartNextSample();
                         }
                       }
                       if((s % 2 == 1 && s > 3 && sample_method != 2) || (s % 2 == 1 && sample_method == 2 && s > 64)) {
                         adaptive_pixel_sampler.test_for_convergence(k, s, nx_end, nx_begin, ny_end, ny_begin);
                       }
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
  } else if (debug_channel == 19) {
    RcppThread::ThreadPool pool(numbercores);
    
    auto worker = [&routput, &goutput, &boutput, 
                   nx, ny,  fov, &hlist, max_depth,
                   cam, &world] (int j) {
                     std::vector<dielectric*> *mat_stack = new std::vector<dielectric*>;
                     random_gen rng(j);
                     for(unsigned int i = 0; i < nx; i++) {
                       ray r;
                       if(fov >= 0) {
                         Float u = (Float(i)) / Float(nx);
                         Float v = (Float(j)) / Float(ny);
                         r = cam->get_ray(u,v, vec3f(0,0,0), rng.unif_rand());
                       } else {
                         point2f u(rng.unif_rand(),rng.unif_rand());
                         point2f u2(rng.unif_rand(),rng.unif_rand());
                         CameraSample samp(u, u2, rng.unif_rand());
                         cam->GenerateRay(samp, &r);
                       }
                       r.pri_stack = mat_stack;
                       point3f qr = calculate_material(r, &world, &hlist,
                                                    max_depth, rng);
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
