#include "../core/integrator.h"

#include "../math/RayMatrix.h"
#include "../math/float.h"
// #define DEBUG
#include "RcppThread.h"
#include "RProgress.h"
#include "../core/adaptivesampler.h"
#include "../core/color.h"
#include "../math/mathinline.h"
#include "../math/filter.h"
#include "../math/sampler.h"
#include "../core/PreviewDisplay.h"
#include "../utils/raylog.h"


void pathtracer(std::size_t numbercores, std::size_t nx, std::size_t ny, std::size_t ns, int debug_channel,
                Float min_variance, std::size_t min_adaptive_size,
                RayMatrix& rgb_output, RayMatrix& normalOutput, RayMatrix& albedoOutput,
                RayMatrix& alpha_output,
                RayMatrix& draw_rgb_output,
                bool progress_bar, int sample_method, int stratified_x, int stratified_y,
                bool verbose, RayCamera* cam,
                Float fov,
                hitable_list& world, hitable_list& hlist,
                Float clampval, std::size_t max_depth, std::size_t roulette_active,
                PreviewDisplay& display, IntegratorType integrator_type, random_gen* rng_override) {
  RProgress::RProgress pb_sampler("Generating Samples [:bar] :percent%");
  pb_sampler.set_width(70);
  RProgress::RProgress pb("Adaptive Raytracing [:bar] :percent%");
  pb.set_width(70);

  Environment pkg = Environment::namespace_env("rayrender");
  Function print_time = pkg["print_time"];

  if(progress_bar) {
    pb_sampler.set_total(ny);
    pb.set_total(ns);
  }
  RayMatrix rgb_output2(nx,ny,3);
  display.write_fast_output = false;
  bool adaptive_on = min_variance > 0;
  adaptive_sampler adaptive_pixel_sampler(numbercores, nx, ny, ns, debug_channel,
                                          min_variance, min_adaptive_size,
                                          rgb_output,
                                          rgb_output2,
                                          normalOutput, albedoOutput,
                                          alpha_output, draw_rgb_output,
                                          adaptive_on);

  size_t nx_small = nx*0.25;
  size_t ny_small = ny*0.25;
  RayMatrix rgb_output_small(nx_small,ny_small,3);
  RayMatrix rgb_output_small2(nx_small,ny_small,3);
  RayMatrix draw_rgb_output_small(nx_small,ny_small,3);

  RayMatrix normal_output_small(nx_small,ny_small,3);
  RayMatrix albedo_output_small(nx_small,ny_small,3);

  RayMatrix alpha_output_small(nx_small,ny_small,1);

  adaptive_sampler adaptive_pixel_sampler_small(numbercores, nx_small, ny_small,
                                                1, debug_channel,
                                                0, 1,
                                                rgb_output_small,
                                                rgb_output_small2,
                                                normal_output_small,
                                                albedo_output_small,
                                                alpha_output_small,
                                                draw_rgb_output_small,
                                                adaptive_on);
  std::vector<random_gen > rngs;
  std::vector<random_gen > rngs_small;

  std::vector<std::unique_ptr<Sampler> > samplers;
  std::vector<std::unique_ptr<Sampler> > samplers_small;
  std::vector<unsigned int> seeds;
  std::vector<unsigned int> seeds_small;

  auto next_seed = [rng_override]() {
    Float rand_unit = rng_override ? rng_override->unif_rand() : unif_rand();
    return static_cast<unsigned int>(rand_unit * std::pow(2, 32));
  };

  seeds.reserve(nx * ny);
  for(unsigned int j = 0; j < ny; j++) {
    if(progress_bar) {
      pb_sampler.tick();
    }
    for(unsigned int i = 0; i < nx; i++) {
      seeds.push_back(next_seed());
    }
  }

  seeds_small.reserve(nx_small * ny_small);
  for(size_t j = 0; j < ny_small; j++) {
    if(progress_bar) {
      pb_sampler.tick();
    }
    for(size_t i = 0; i < nx_small; i++) {
      seeds_small.push_back(next_seed());
    }
  }

  auto reset_sampler_state = [sample_method, ns, stratified_x, stratified_y] (
      size_t width, size_t height, const std::vector<unsigned int>& state_seeds,
      std::vector<random_gen>& state_rngs, std::vector<std::unique_ptr<Sampler> >& state_samplers) {
    state_rngs.clear();
    state_samplers.clear();
    state_rngs.reserve(state_seeds.size());
    state_samplers.reserve(state_seeds.size());
    size_t index = 0;
    for(size_t j = 0; j < height; j++) {
      for(size_t i = 0; i < width; i++) {
        random_gen rng_single(state_seeds[index]);
        state_rngs.push_back(rng_single);
        if(sample_method == 0) {
          state_samplers.push_back(std::unique_ptr<Sampler>(new RandomSampler(rng_single)));
          state_samplers.back()->StartPixel(0,0);
        } else if (sample_method == 1) {
          state_samplers.push_back(std::unique_ptr<Sampler>(new StratifiedSampler(stratified_x, stratified_y,
                                                                                   true, 5, rng_single)));
          state_samplers.back()->StartPixel(0,0);
        } else if (sample_method == 2) {
          state_samplers.push_back(std::unique_ptr<Sampler>(new SobolSampler(ns, rng_single)));
          state_samplers.back()->StartPixel(i,j);
        } else {
          state_samplers.push_back(std::unique_ptr<Sampler>(new SobolBlueNoiseSampler(rng_single)));
          state_samplers.back()->StartPixel(i,j);
        }
        state_samplers.back()->SetSampleNumber(0);
        index++;
      }
    }
  };

  reset_sampler_state(nx, ny, seeds, rngs, samplers);
  reset_sampler_state(nx_small, ny_small, seeds_small, rngs_small, samplers_small);
  random_gen rng_interactive(next_seed());

  auto copy_small_preview = [&adaptive_pixel_sampler, &adaptive_pixel_sampler_small, nx_small, ny_small, nx, ny]() {
    Float ratio_x = (Float)nx_small/(Float)nx;
    Float ratio_y = (Float)ny_small/(Float)ny;
    for(size_t ii = 0; ii < nx; ii++) {
      for(size_t jj = 0; jj < ny; jj++) {
        int iii = (Float)ii * ratio_x;
        int jjj = (Float)jj * ratio_y;
        adaptive_pixel_sampler.rgb(ii,jj,0) = adaptive_pixel_sampler_small.rgb(iii,jjj,0);
        adaptive_pixel_sampler.rgb(ii,jj,1) = adaptive_pixel_sampler_small.rgb(iii,jjj,1);
        adaptive_pixel_sampler.rgb(ii,jj,2) = adaptive_pixel_sampler_small.rgb(iii,jjj,2);

        adaptive_pixel_sampler.normalOutput(ii,jj,0) = adaptive_pixel_sampler_small.normalOutput(iii,jjj,0);
        adaptive_pixel_sampler.normalOutput(ii,jj,1) = adaptive_pixel_sampler_small.normalOutput(iii,jjj,1);
        adaptive_pixel_sampler.normalOutput(ii,jj,2) = adaptive_pixel_sampler_small.normalOutput(iii,jjj,2);

        adaptive_pixel_sampler.albedoOutput(ii,jj,0) = adaptive_pixel_sampler_small.albedoOutput(iii,jjj,0);
        adaptive_pixel_sampler.albedoOutput(ii,jj,1) = adaptive_pixel_sampler_small.albedoOutput(iii,jjj,1);
        adaptive_pixel_sampler.albedoOutput(ii,jj,2) = adaptive_pixel_sampler_small.albedoOutput(iii,jjj,2);

        adaptive_pixel_sampler.a(ii,jj,0) = adaptive_pixel_sampler_small.a(iii,jjj,0);
      }
    }
  };

  auto render_full_sample = [&adaptive_pixel_sampler, numbercores, nx, ny, sample_method,
                             &rngs, fov, &samplers, cam, &world, &hlist,
                             clampval, max_depth, roulette_active, integrator_type] (size_t s) {
    RcppThread::ThreadPool pool(numbercores);
    auto worker = [&adaptive_pixel_sampler,
                   nx, ny, s, sample_method,
                   &rngs, fov, &samplers,
                   cam, &world, &hlist,
                   clampval, max_depth, roulette_active, integrator_type] (int k) {
                     int nx_begin = adaptive_pixel_sampler.pixel_chunks[k].startx;
                     int ny_begin = adaptive_pixel_sampler.pixel_chunks[k].starty;
                     int nx_end = adaptive_pixel_sampler.pixel_chunks[k].endx;
                     int ny_end = adaptive_pixel_sampler.pixel_chunks[k].endy;

                     std::vector<dielectric*> *mat_stack = new std::vector<dielectric*>;
                     for(int i = nx_begin; i < nx_end; i++) {
                       for(int j = ny_begin; j < ny_end; j++) {
                         int index = j + ny * i;
                         Ray r;
                         vec2f u2 = samplers[index]->Get2D();
                         Float weight(1.0);
                         Float u = (Float(i) + u2.xy.x) / Float(nx);
                         Float v = (Float(j) + u2.xy.y) / Float(ny);

                         if(fov >= 0) {
                           r = cam->get_ray(u,v, convert_to_point3(rand_to_unit(samplers[index]->Get2D())),
                                            samplers[index]->Get1D());
                         } else {
                           CameraSample samp({1-u,1-v},samplers[index]->Get2D(), samplers[index]->Get1D());
                           weight = cam->GenerateRay(samp, &r);
                         }
                         r.pri_stack = mat_stack;
                         bool alpha = false;
                         point3f color_sample;
                         normal3f normal_sample;
                         point3f albedo_sample;
                         color(r, &world, &hlist, max_depth, roulette_active, rngs[index],
                               samplers[index].get(), alpha,integrator_type,
                               color_sample, normal_sample, albedo_sample);
                         point3f col = weight != 0 ? clamp_point(de_nan(color_sample),
                                                                 0, clampval) * weight * cam->get_iso() : 0;
                         if(alpha) {
                           adaptive_pixel_sampler.add_alpha_count(i,j);
                         }
                         mat_stack->clear();
                         adaptive_pixel_sampler.add_color_main(i, j, col);
                         if(s % 2 == 0) {
                           adaptive_pixel_sampler.add_color_sec(i, j, col);
                         }
                         adaptive_pixel_sampler.add_albedo(i, j, albedo_sample);
                         adaptive_pixel_sampler.add_normal(i, j, normal_sample);
                         samplers[index]->StartNextSample();
                       }
                     }
                     if (adaptive_pixel_sampler.adaptive_on) {
                      if((s % 2 == 1 && s > 3 && sample_method != 2) || (s % 2 == 1 && sample_method == 2 && s > 64)) {
                        adaptive_pixel_sampler.test_for_convergence(k, s, nx_end, nx_begin, ny_end, ny_begin);
                      }
                     }
                     delete mat_stack;
                   };
    for(size_t j = 0; j < adaptive_pixel_sampler.size(); j++) {
      pool.push(worker, j);
    }
    pool.join();
    if (adaptive_pixel_sampler.adaptive_on) {
      if(s % 2 == 1 && s > 1) {
        adaptive_pixel_sampler.split_remove_chunks(s);
      }
    }
    adaptive_pixel_sampler.max_s++;
  };

  auto render_small_sample = [&adaptive_pixel_sampler_small, numbercores, nx_small, ny_small, sample_method,
                              &rngs_small, fov, &samplers_small, cam, &world, &hlist,
                              clampval, max_depth, roulette_active, integrator_type] (size_t s) {
    RcppThread::ThreadPool pool(numbercores);
    auto worker = [&adaptive_pixel_sampler_small,
                   nx_small, ny_small, s, sample_method,
                   &rngs_small, fov, &samplers_small,
                   cam, &world, &hlist,
                   clampval, max_depth, roulette_active, integrator_type] (int k) {
                     int nx_begin = adaptive_pixel_sampler_small.pixel_chunks[k].startx;
                     int ny_begin = adaptive_pixel_sampler_small.pixel_chunks[k].starty;
                     int nx_end = adaptive_pixel_sampler_small.pixel_chunks[k].endx;
                     int ny_end = adaptive_pixel_sampler_small.pixel_chunks[k].endy;

                     std::vector<dielectric*> *mat_stack = new std::vector<dielectric*>;
                     for(int i = nx_begin; i < nx_end; i++) {
                       for(int j = ny_begin; j < ny_end; j++) {
                         int index = j + ny_small * i;
                         Ray r;
                         vec2f u2 = samplers_small[index]->Get2D();
                         Float weight(1.0);
                         Float u = (Float(i) + u2.xy.x) / Float(nx_small);
                         Float v = (Float(j) + u2.xy.y) / Float(ny_small);

                         if(fov >= 0) {
                           r = cam->get_ray(u,v, convert_to_point3(rand_to_unit(samplers_small[index]->Get2D())),
                                            samplers_small[index]->Get1D());
                         } else {
                           CameraSample samp({1-u,1-v},samplers_small[index]->Get2D(), samplers_small[index]->Get1D());
                           weight = cam->GenerateRay(samp, &r);
                         }
                         r.pri_stack = mat_stack;
                         bool alpha = false;
                         point3f color_sample;
                         normal3f normal_sample;
                         point3f albedo_sample;
                         color(r, &world, &hlist, max_depth, roulette_active,
                               rngs_small[index],
                               samplers_small[index].get(),
                               alpha, integrator_type,
                               color_sample, normal_sample, albedo_sample);
                         point3f col = weight != 0 ? clamp_point(de_nan(color_sample),
                                                                 0, clampval) * weight * cam->get_iso() : 0;
                         if(alpha) {
                           adaptive_pixel_sampler_small.add_alpha_count(i,j);
                         }
                         adaptive_pixel_sampler_small.add_albedo(i, j, albedo_sample);
                         adaptive_pixel_sampler_small.add_normal(i, j, normal_sample);
                         mat_stack->clear();
                         adaptive_pixel_sampler_small.add_color_main(i, j, col);
                         if(s % 2 == 0) {
                           adaptive_pixel_sampler_small.add_color_sec(i, j, col);
                         }
                         samplers_small[index]->StartNextSample();
                       }
                     }
                     if (adaptive_pixel_sampler_small.adaptive_on) {
                        if((s % 2 == 1 && s > 3 && sample_method != 2) || (s % 2 == 1 && sample_method == 2 && s > 64)) {
                          adaptive_pixel_sampler_small.test_for_convergence(k, s, nx_end, nx_begin, ny_end, ny_begin);
                        }
                     }
                     delete mat_stack;
                   };
    for(size_t j = 0; j < adaptive_pixel_sampler_small.size(); j++) {
      pool.push(worker, j);
    }
    pool.join();
    if (adaptive_pixel_sampler_small.adaptive_on) {
      if(s % 2 == 1 && s > 1) {
        adaptive_pixel_sampler_small.split_remove_chunks(s);
      }
    }
    adaptive_pixel_sampler_small.max_s++;
  };

  auto reset_render_state = [&]() {
    adaptive_pixel_sampler.reset();
    adaptive_pixel_sampler_small.reset();
    display.ResetPreviewExposure();
    reset_sampler_state(nx, ny, seeds, rngs, samplers);
    reset_sampler_state(nx_small, ny_small, seeds_small, rngs_small, samplers_small);
  };

  print_time(verbose, "Allocating sampler" );
  if(display.deferred_render) {
    size_t preview_sample = 0;
    bool completed_final_render = false;
    while(!display.terminate && !completed_final_render) {
      while(!display.render_requested && !display.terminate) {
        Rcpp::checkUserInterrupt();
        if(!display.write_fast_output) {
          render_full_sample(preview_sample);
        } else {
          render_small_sample(preview_sample);
          copy_small_preview();
        }
        display.DrawImage(adaptive_pixel_sampler, adaptive_pixel_sampler_small,
                          preview_sample, pb, false,
                          0, &world, rng_interactive);
        if(!display.render_requested && !display.terminate) {
          preview_sample++;
        }
      }
      if(display.terminate) {
        adaptive_pixel_sampler.ns = std::max<size_t>(preview_sample, 1);
        adaptive_pixel_sampler.max_s = adaptive_pixel_sampler.ns;
        break;
      }

      display.write_fast_output = false;
      reset_render_state();

      for(size_t s = 0; s < static_cast<size_t>(ns); s++) {
        Rcpp::checkUserInterrupt();
        if(progress_bar && !display.preview) {
          pb.tick();
        }
        if(!display.write_fast_output) {
          render_full_sample(s);
        } else {
          render_small_sample(s);
          copy_small_preview();
        }
        display.DrawImage(adaptive_pixel_sampler, adaptive_pixel_sampler_small,
                          s, pb, progress_bar,
                          (Float)s/(Float)ns, &world, rng_interactive);
        if(display.terminate && display.preview) {
          adaptive_pixel_sampler.ns = s;
          adaptive_pixel_sampler.max_s = s;
          break;
        }
        if(!display.render_requested) {
          preview_sample = 0;
          reset_render_state();
          break;
        }
      }
      if(display.terminate) {
        break;
      }
      if(display.render_requested) {
        completed_final_render = true;
      }
    }
  } else {
    for(size_t s = 0; s < static_cast<size_t>(ns); s++) {
      Rcpp::checkUserInterrupt();
      if(progress_bar && !display.preview) {
        pb.tick();
      }
      if(!display.write_fast_output) {
        render_full_sample(s);
      } else {
        render_small_sample(s);
        copy_small_preview();
      }
      display.DrawImage(adaptive_pixel_sampler, adaptive_pixel_sampler_small,
                        s, pb, progress_bar,
                        (Float)s/(Float)ns, &world, rng_interactive);
      if(display.terminate && display.preview) {
        adaptive_pixel_sampler.ns = s;
        adaptive_pixel_sampler.max_s = s;
        break;
      }
    }
  }

  #ifdef RAY_COLOR_DEBUG
  Rcpp::Rcout << "ns: " << adaptive_pixel_sampler.ns << " adaptive_pixel_sampler.rgb pre divide:\n";
  adaptive_pixel_sampler.rgb.print();
  #endif
  adaptive_pixel_sampler.write_final_pixels();
}
