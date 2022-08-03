#include "integrator.h"

#include "float.h"
// #define DEBUG
#include "RcppThread.h"
#include "RProgress.h"
#include "adaptivesampler.h"
#include "color.h"
#include "mathinline.h"
#include "filter.h"
#include "sampler.h"
#include "PreviewDisplay.h"


void pathtracer(std::size_t numbercores, std::size_t nx, std::size_t ny, std::size_t ns, int debug_channel,
                Float min_variance, std::size_t min_adaptive_size, 
                RayMatrix& routput, RayMatrix& goutput, RayMatrix& boutput,
                bool progress_bar, int sample_method, Rcpp::NumericVector& stratified_dim,
                bool verbose, RayCamera* cam, 
                Float fov,
                hitable_list& world, hitable_list& hlist,
                Float clampval, std::size_t max_depth, std::size_t roulette_active,
                PreviewDisplay& display) {
  RProgress::RProgress pb_sampler("Generating Samples [:bar] :percent%");
  pb_sampler.set_width(70);
  RProgress::RProgress pb("Adaptive Raytracing [:bar] :percent%");
  pb.set_width(70);
  
  if(progress_bar) {
    pb_sampler.set_total(ny);
    pb.set_total(ns);
  }
#ifdef DEBUG
  std::remove("rays.txt");
#endif
  RayMatrix routput2(nx,ny);
  RayMatrix goutput2(nx,ny);
  RayMatrix boutput2(nx,ny);
  display.write_fast_output = false;
  adaptive_sampler adaptive_pixel_sampler(numbercores, nx, ny, ns, debug_channel,
                                          min_variance, min_adaptive_size,
                                          routput, goutput, boutput,
                                          routput2, goutput2, boutput2);
  
  int nx_small = nx*0.25;
  int ny_small = ny*0.25;
  RayMatrix routput_small(nx_small,ny_small);
  RayMatrix goutput_small(nx_small,ny_small);
  RayMatrix boutput_small(nx_small,ny_small);
  RayMatrix routput_small2(nx_small,ny_small);
  RayMatrix goutput_small2(nx_small,ny_small);
  RayMatrix boutput_small2(nx_small,ny_small);
  adaptive_sampler adaptive_pixel_sampler_small(numbercores, nx_small, ny_small, 1, debug_channel,
                                                0, 1,
                                                routput_small, goutput_small, boutput_small,
                                                routput_small2, goutput_small2, boutput_small2);
  std::vector<random_gen > rngs;
  std::vector<random_gen > rngs_small;
  
  std::vector<std::unique_ptr<Sampler> > samplers;
  std::vector<std::unique_ptr<Sampler> > samplers_small;
  
  auto start = std::chrono::high_resolution_clock::now();
  
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
  
  for(unsigned int j = 0; j < ny_small; j++) {
    if(progress_bar) {
      pb_sampler.tick();
    }
    for(unsigned int i = 0; i < nx_small; i++) {
      random_gen rng_single(unif_rand() * std::pow(2,32));
      rngs_small.push_back(rng_single);
      if(sample_method == 0) {
        samplers_small.push_back(std::unique_ptr<Sampler>(new RandomSampler(rng_single)));
        samplers_small.back()->StartPixel(0,0);
      } else if (sample_method == 1){
        samplers_small.push_back(std::unique_ptr<Sampler>(new StratifiedSampler(stratified_dim(0), stratified_dim(1),
                                                                          true, 5, rng_single)));
        samplers_small.back()->StartPixel(0,0);
      } else if (sample_method == 2) {
        samplers_small.push_back(std::unique_ptr<Sampler>(new SobolSampler(ns, rng_single)));
        samplers_small.back()->StartPixel(i,j);
      } else {
        samplers_small.push_back(std::unique_ptr<Sampler>(new SobolBlueNoiseSampler(rng_single)));
        samplers_small.back()->StartPixel(i,j);
      }
      samplers_small.back()->SetSampleNumber(0);
    }
  }
  random_gen rng_interactive(unif_rand() * std::pow(2,32));
  
  if(verbose) {
    auto finish = std::chrono::high_resolution_clock::now();
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
    if(progress_bar && !display.preview) {
      pb.tick();
    }
    RcppThread::ThreadPool pool(numbercores);
    if(!display.write_fast_output) {
      auto worker = [&adaptive_pixel_sampler,
                     nx, ny, s, sample_method,
                     &rngs, fov, &samplers,
                     cam, &world, &hlist,
                     clampval, max_depth, roulette_active] (int k) {
                       // MitchellFilter fil(vec2f(1.0),1./3.,1./3.);
                       int nx_begin = adaptive_pixel_sampler.pixel_chunks[k].startx;
                       int ny_begin = adaptive_pixel_sampler.pixel_chunks[k].starty;
                       int nx_end = adaptive_pixel_sampler.pixel_chunks[k].endx;
                       int ny_end = adaptive_pixel_sampler.pixel_chunks[k].endy;
                       
                       std::vector<dielectric*> *mat_stack = new std::vector<dielectric*>;
                       for(int i = nx_begin; i < nx_end; i++) {
                         for(int j = ny_begin; j < ny_end; j++) {
                           int index = j + ny * i;
                           ray r; 
                           vec2f u2 = samplers[index]->Get2D();
                           Float weight(1.0);
                           Float u = (Float(i) + u2.x()) / Float(nx);
                           Float v = (Float(j) + u2.y()) / Float(ny);
                           
                           if(fov >= 0) {
                             r = cam->get_ray(u,v, rand_to_unit(samplers[index]->Get2D()),
                                              samplers[index]->Get1D());
                           } else {
                             CameraSample samp({1-u,1-v},samplers[index]->Get2D(), samplers[index]->Get1D());
                             weight = cam->GenerateRay(samp, &r);
                           }
                           r.pri_stack = mat_stack;
                           point3f col = weight != 0 ? clamp_point(de_nan(color(r, &world, &hlist, max_depth, 
                                                         roulette_active, rngs[index], samplers[index].get())),
                                               0, clampval) * weight * cam->get_iso() : 0;
                           // col = col * fil.Evaluate(u2);
                           mat_stack->clear();
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
    } else {
      auto worker = [&adaptive_pixel_sampler_small,
                     nx_small, ny_small, s, sample_method,
                     &rngs_small, fov, &samplers_small,
                     cam, &world, &hlist,
                     clampval, max_depth, roulette_active] (int k) {
                       // MitchellFilter fil(vec2f(1.0),1./3.,1./3.);
                       int nx_begin = adaptive_pixel_sampler_small.pixel_chunks[k].startx;
                       int ny_begin = adaptive_pixel_sampler_small.pixel_chunks[k].starty;
                       int nx_end = adaptive_pixel_sampler_small.pixel_chunks[k].endx;
                       int ny_end = adaptive_pixel_sampler_small.pixel_chunks[k].endy;
                       
                       std::vector<dielectric*> *mat_stack = new std::vector<dielectric*>;
                       for(int i = nx_begin; i < nx_end; i++) {
                         for(int j = ny_begin; j < ny_end; j++) {
                           int index = j + ny_small * i;
                           ray r; 
                           vec2f u2 = samplers_small[index]->Get2D();
                           Float weight(1.0);
                           Float u = (Float(i) + u2.x()) / Float(nx_small);
                           Float v = (Float(j) + u2.y()) / Float(ny_small);
                           
                           if(fov >= 0) {
                             r = cam->get_ray(u,v, rand_to_unit(samplers_small[index]->Get2D()),
                                              samplers_small[index]->Get1D());
                           } else {
                             CameraSample samp({1-u,1-v},samplers_small[index]->Get2D(), samplers_small[index]->Get1D());
                             weight = cam->GenerateRay(samp, &r);
                           }
                           r.pri_stack = mat_stack;
                           point3f col = weight != 0 ? clamp_point(de_nan(color(r, &world, &hlist, max_depth, 
                                                                                roulette_active, rngs_small[index], samplers_small[index].get())),
                                                                                0, clampval) * weight * cam->get_iso() : 0;
                           // col = col * fil.Evaluate(u2);
                           mat_stack->clear();
                           adaptive_pixel_sampler_small.add_color_main(i, j, col);
                           if(s % 2 == 0) {
                             adaptive_pixel_sampler_small.add_color_sec(i, j, col);
                           }
                           samplers_small[index]->StartNextSample();
                         }
                       }
                       if((s % 2 == 1 && s > 3 && sample_method != 2) || (s % 2 == 1 && sample_method == 2 && s > 64)) {
                         adaptive_pixel_sampler_small.test_for_convergence(k, s, nx_end, nx_begin, ny_end, ny_begin);
                       }
                       delete mat_stack;
                     };
      
      for(size_t j = 0; j < adaptive_pixel_sampler_small.size(); j++) {
        pool.push(worker, j);
      }
      pool.join();
      
      if(s % 2 == 1 && s > 1) {
        adaptive_pixel_sampler_small.split_remove_chunks(s);
      }
      adaptive_pixel_sampler_small.max_s++;
      
      Float ratio_x = (Float)nx_small/(Float)nx;
      Float ratio_y = (Float)ny_small/(Float)ny;
      
      for(int ii = 0; ii < nx; ii++) {
        for(int jj = 0; jj < ny; jj++) {
          adaptive_pixel_sampler.r(ii,jj) = adaptive_pixel_sampler_small.r((Float)ii * ratio_x,(Float)jj * ratio_y);
          adaptive_pixel_sampler.g(ii,jj) = adaptive_pixel_sampler_small.g((Float)ii * ratio_x,(Float)jj * ratio_y);
          adaptive_pixel_sampler.b(ii,jj) = adaptive_pixel_sampler_small.b((Float)ii * ratio_x,(Float)jj * ratio_y);
        }
      }
    }
    display.DrawImage(adaptive_pixel_sampler, adaptive_pixel_sampler_small,
                      s, pb, progress_bar,
                      (Float)s/(Float)ns, &world, rng_interactive);
    if(display.terminate) {
      adaptive_pixel_sampler.ns = s;
      adaptive_pixel_sampler.max_s = s;
      break;
    }
  }
  adaptive_pixel_sampler.write_final_pixels();
}
