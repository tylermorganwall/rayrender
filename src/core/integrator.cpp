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
#include "../core/oidn_aux.h"
#include "../volumes/boundary.h"
#include "../core/oidn_denoiser.h"
#include "../utils/raylog.h"

#include <atomic>
#include <chrono>
#include <exception>
#include <future>
#include <thread>

static const size_t FAST_INTERACTIVE_PREVIEW_SAMPLES = 4;

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
  bool has_media = hlist.volume_scene && hlist.volume_scene->has_media;
  display.volume_scene = hlist.volume_scene;
  display.transparent_volume_background = hlist.volume_scene && hlist.volume_scene->transparent_background;
  if(hlist.volume_scene) {
    hlist.volume_scene->light_sampler = std::make_shared<VolumeLightSampler>(hlist);
    hlist.volume_scene->statistics.Reset();
  }
  adaptive_sampler adaptive_pixel_sampler(numbercores, nx, ny, ns, debug_channel,
                                          min_variance, min_adaptive_size,
                                          rgb_output,
                                          rgb_output2,
                                          normalOutput, albedoOutput,
                                          alpha_output, draw_rgb_output,
                                          adaptive_on, has_media);

  size_t nx_small = nx*0.25;
  size_t ny_small = ny*0.25;
  RayMatrix rgb_output_small(nx_small,ny_small,3);
  RayMatrix rgb_output_small2(nx_small,ny_small,3);
  RayMatrix draw_rgb_output_small(nx_small,ny_small,3);

  RayMatrix normal_output_small(nx_small,ny_small,3);
  RayMatrix albedo_output_small(nx_small,ny_small,3);

  RayMatrix alpha_output_small(nx_small,ny_small,1);

  adaptive_sampler adaptive_pixel_sampler_small(numbercores, nx_small, ny_small,
                                                FAST_INTERACTIVE_PREVIEW_SAMPLES, debug_channel,
                                                0, 1,
                                                rgb_output_small,
                                                rgb_output_small2,
                                                normal_output_small,
                                                albedo_output_small,
                                                alpha_output_small,
                                                draw_rgb_output_small,
                                                adaptive_on, has_media);

#ifdef HAS_OIDN
  RayMatrix oidn_normal_output_small(nx_small, ny_small, 3);
  RayMatrix oidn_albedo_output_small(nx_small, ny_small, 3);
  RayOidnDenoiser fast_preview_denoiser;
  bool denoise_fast_preview = display.denoise;
  if(denoise_fast_preview) {
    fast_preview_denoiser.Setup(rgb_output_small,
                                oidn_albedo_output_small,
                                oidn_normal_output_small,
                                draw_rgb_output_small,
                                nx_small,
                                ny_small,
                                RayOidnQuality::Fast,
                                false,
                                false, !has_media);
  }
#endif

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

  auto reset_sampler_state = [sample_method, ns, stratified_x, stratified_y, integrator_type] (
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
        state_samplers.back()->independent_dimensions = integrator_type == IntegratorType::ShadowRays;
        state_samplers.back()->SetSampleNumber(0);
        index++;
      }
    }
  };

  reset_sampler_state(nx, ny, seeds, rngs, samplers);
  reset_sampler_state(nx_small, ny_small, seeds_small, rngs_small, samplers_small);
  random_gen rng_interactive(next_seed());

#ifdef HAS_OIDN
  auto make_oidn_aux_options = [sample_method,
                                stratified_x,
                                stratified_y,
                                max_depth] (std::size_t samples) {
    OidnAuxRenderOptions options;
    options.samples = samples;
    options.max_depth = max_depth;
    options.max_dielectric_splits = 1;
    options.sample_method = sample_method;
    options.stratified_x = stratified_x;
    options.stratified_y = stratified_y;
    return options;
  };

  auto ensure_full_preview_oidn_aux = [&]() {
    if(has_media) { display.oidn_aux_dirty = false; return; }
    if(!display.preview ||
       !display.denoise ||
       !display.oidn_aux_dirty ||
       display.oidn_albedo_output == nullptr ||
       display.oidn_normal_output == nullptr) {
      return;
    }
    if(display.PollCloseEvent()) {
      return;
    }
    OidnAuxRenderOptions options = make_oidn_aux_options(1);
    render_oidn_aux_features(numbercores,
                             nx,
                             ny,
                             cam,
                             fov,
                             &world,
                             options,
                             *display.oidn_normal_output,
                             *display.oidn_albedo_output,
                             [&display]() { return display.PollCloseEvent(); });
    if(display.terminate) {
      return;
    }
    display.MarkOidnAuxClean(false);
  };

  auto ensure_fast_preview_oidn_aux = [&]() {
    if(has_media) { display.oidn_fast_aux_dirty = false; return; }
    if(!display.preview ||
       !denoise_fast_preview ||
       !display.oidn_fast_aux_dirty) {
      return;
    }
    if(display.PollCloseEvent()) {
      return;
    }
    OidnAuxRenderOptions options = make_oidn_aux_options(1);
    render_oidn_aux_features(numbercores,
                             nx_small,
                             ny_small,
                             cam,
                             fov,
                             &world,
                             options,
                             oidn_normal_output_small,
                             oidn_albedo_output_small,
                             [&display]() { return display.PollCloseEvent(); });
    if(display.terminate) {
      return;
    }
    display.MarkOidnAuxClean(true);
  };
#endif

  auto copy_small_color = [nx_small, ny_small, nx, ny] (RayMatrix& target, RayMatrix& source) {
    Float ratio_x = (Float)nx_small/(Float)nx;
    Float ratio_y = (Float)ny_small/(Float)ny;
    for(size_t ii = 0; ii < nx; ii++) {
      for(size_t jj = 0; jj < ny; jj++) {
        int iii = (Float)ii * ratio_x;
        int jjj = (Float)jj * ratio_y;
        target(ii,jj,0) = source(iii,jjj,0) / (Float)FAST_INTERACTIVE_PREVIEW_SAMPLES;
        target(ii,jj,1) = source(iii,jjj,1) / (Float)FAST_INTERACTIVE_PREVIEW_SAMPLES;
        target(ii,jj,2) = source(iii,jjj,2) / (Float)FAST_INTERACTIVE_PREVIEW_SAMPLES;
      }
    }
  };

  auto copy_small_preview = [&adaptive_pixel_sampler, &adaptive_pixel_sampler_small, nx_small, ny_small, nx, ny]() {
    Float ratio_x = (Float)nx_small/(Float)nx;
    Float ratio_y = (Float)ny_small/(Float)ny;
    for(size_t ii = 0; ii < nx; ii++) {
      for(size_t jj = 0; jj < ny; jj++) {
        int iii = (Float)ii * ratio_x;
        int jjj = (Float)jj * ratio_y;
        adaptive_pixel_sampler.a(ii,jj,0) = adaptive_pixel_sampler_small.a(iii,jjj,0) /
          (Float)FAST_INTERACTIVE_PREVIEW_SAMPLES;
      }
    }
  };

  std::atomic<bool> render_cancelled(false);
  auto wait_for_render_jobs = [&display, &render_cancelled](
      std::vector<std::future<void> >& futures) -> bool {
    std::exception_ptr interrupt_exception = nullptr;
    bool done = false;
    while(!done) {
      done = true;
      for(auto& future : futures) {
        if(future.wait_for(std::chrono::milliseconds(0)) != std::future_status::ready) {
          done = false;
          break;
        }
      }
      if(done) {
        break;
      }
      if(display.PollCloseEvent()) {
        render_cancelled.store(true, std::memory_order_relaxed);
      }
      try {
        RcppThread::checkUserInterrupt();
      } catch(...) {
        render_cancelled.store(true, std::memory_order_relaxed);
        interrupt_exception = std::current_exception();
        break;
      }
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    for(auto& future : futures) {
      future.wait();
    }
    for(auto& future : futures) {
      future.get();
    }
    if(interrupt_exception != nullptr) {
      std::rethrow_exception(interrupt_exception);
    }
    return !render_cancelled.load(std::memory_order_relaxed);
  };

  auto render_full_sample = [&adaptive_pixel_sampler, numbercores, nx, ny, sample_method,
                             &rngs, fov, &samplers, cam, &world, &hlist,
                             clampval, max_depth, roulette_active, integrator_type,
                             &render_cancelled, &wait_for_render_jobs] (size_t s) -> bool {
    render_cancelled.store(false, std::memory_order_relaxed);
    RcppThread::ThreadPool pool(numbercores);
    auto worker = [&adaptive_pixel_sampler,
                   nx, ny, s, sample_method,
                   &rngs, fov, &samplers,
                   cam, &world, &hlist,
                   clampval, max_depth, roulette_active, integrator_type,
                   &render_cancelled] (int k) {
                     int nx_begin = adaptive_pixel_sampler.pixel_chunks[k].startx;
                     int ny_begin = adaptive_pixel_sampler.pixel_chunks[k].starty;
                     int nx_end = adaptive_pixel_sampler.pixel_chunks[k].endx;
                     int ny_end = adaptive_pixel_sampler.pixel_chunks[k].endy;

                     std::vector<dielectric*> *mat_stack = new std::vector<dielectric*>;
                     for(int i = nx_begin; i < nx_end; i++) {
                       if(render_cancelled.load(std::memory_order_relaxed)) {
                         break;
                       }
                       for(int j = ny_begin; j < ny_end; j++) {
                         if(render_cancelled.load(std::memory_order_relaxed)) {
                           break;
                         }
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
                         Float alpha = 0;
                         point3f color_sample;
                         normal3f normal_sample;
                         point3f albedo_sample;
                         color(r, &world, &hlist, max_depth, roulette_active, rngs[index],
                               samplers[index].get(), alpha,integrator_type,
                               color_sample, normal_sample, albedo_sample,
                               &render_cancelled);
                         point3f col = weight != 0 ? clamp_point(de_nan(color_sample),
                                                                 0, clampval) * weight * cam->get_iso() : 0;
                         adaptive_pixel_sampler.add_alpha_count(i,j, alpha);
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
                      if(!render_cancelled.load(std::memory_order_relaxed) &&
                         ((s % 2 == 1 && s > 3 && sample_method != 2) ||
                          (s % 2 == 1 && sample_method == 2 && s > 64))) {
                        adaptive_pixel_sampler.test_for_convergence(k, s, nx_end, nx_begin, ny_end, ny_begin);
                      }
                     }
                     delete mat_stack;
                   };
    std::vector<std::future<void> > futures;
    futures.reserve(adaptive_pixel_sampler.size());
    for(size_t j = 0; j < adaptive_pixel_sampler.size(); j++) {
      futures.push_back(pool.pushReturn(worker, static_cast<int>(j)));
    }
    bool completed_sample = wait_for_render_jobs(futures);
    if (adaptive_pixel_sampler.adaptive_on) {
      if(completed_sample && s % 2 == 1 && s > 1) {
        adaptive_pixel_sampler.split_remove_chunks(s);
      }
    }
    if(completed_sample) {
      adaptive_pixel_sampler.max_s++;
    }
    return completed_sample;
  };

  auto render_small_sample = [&adaptive_pixel_sampler_small, numbercores, nx_small, ny_small, sample_method,
                              &rngs_small, fov, &samplers_small, cam, &world, &hlist,
                              clampval, max_depth, roulette_active, integrator_type,
                              &render_cancelled, &wait_for_render_jobs] (size_t s) -> bool {
    render_cancelled.store(false, std::memory_order_relaxed);
    RcppThread::ThreadPool pool(numbercores);
    auto worker = [&adaptive_pixel_sampler_small,
                   nx_small, ny_small, s, sample_method,
                   &rngs_small, fov, &samplers_small,
                   cam, &world, &hlist,
                   clampval, max_depth, roulette_active, integrator_type,
                   &render_cancelled] (int k) {
                     int nx_begin = adaptive_pixel_sampler_small.pixel_chunks[k].startx;
                     int ny_begin = adaptive_pixel_sampler_small.pixel_chunks[k].starty;
                     int nx_end = adaptive_pixel_sampler_small.pixel_chunks[k].endx;
                     int ny_end = adaptive_pixel_sampler_small.pixel_chunks[k].endy;

                     std::vector<dielectric*> *mat_stack = new std::vector<dielectric*>;
                     for(int i = nx_begin; i < nx_end; i++) {
                       if(render_cancelled.load(std::memory_order_relaxed)) {
                         break;
                       }
                       for(int j = ny_begin; j < ny_end; j++) {
                         if(render_cancelled.load(std::memory_order_relaxed)) {
                           break;
                         }
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
                         Float alpha = 0;
                         point3f color_sample;
                         normal3f normal_sample;
                         point3f albedo_sample;
                         color(r, &world, &hlist, max_depth, roulette_active,
                               rngs_small[index],
                               samplers_small[index].get(),
                               alpha, integrator_type,
                               color_sample, normal_sample, albedo_sample,
                               &render_cancelled);
                         point3f col = weight != 0 ? clamp_point(de_nan(color_sample),
                                                                 0, clampval) * weight * cam->get_iso() : 0;
                         adaptive_pixel_sampler_small.add_alpha_count(i,j, alpha);
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
                        if(!render_cancelled.load(std::memory_order_relaxed) &&
                           ((s % 2 == 1 && s > 3 && sample_method != 2) ||
                            (s % 2 == 1 && sample_method == 2 && s > 64))) {
                          adaptive_pixel_sampler_small.test_for_convergence(k, s, nx_end, nx_begin, ny_end, ny_begin);
                        }
                     }
                     delete mat_stack;
                   };
    std::vector<std::future<void> > futures;
    futures.reserve(adaptive_pixel_sampler_small.size());
    for(size_t j = 0; j < adaptive_pixel_sampler_small.size(); j++) {
      futures.push_back(pool.pushReturn(worker, static_cast<int>(j)));
    }
    bool completed_sample = wait_for_render_jobs(futures);
    if (adaptive_pixel_sampler_small.adaptive_on) {
      if(completed_sample && s % 2 == 1 && s > 1) {
        adaptive_pixel_sampler_small.split_remove_chunks(s);
      }
    }
    if(completed_sample) {
      adaptive_pixel_sampler_small.max_s++;
    }
    return completed_sample;
  };

  auto render_fast_preview_sample = [&render_small_sample,
                                     &copy_small_color,
                                     &copy_small_preview,
                                     &adaptive_pixel_sampler,
                                     &adaptive_pixel_sampler_small,
                                     &display
#ifdef HAS_OIDN
                                     ,
                                     denoise_fast_preview,
                                     &fast_preview_denoiser,
                                     &ensure_fast_preview_oidn_aux
#endif
                                     ] (size_t s) -> bool {
    for(size_t sample_offset = 0;
        sample_offset < FAST_INTERACTIVE_PREVIEW_SAMPLES;
        sample_offset++) {
      if(display.PollCloseEvent()) {
        return false;
      }
      if(!render_small_sample(s * FAST_INTERACTIVE_PREVIEW_SAMPLES + sample_offset)) {
        return false;
      }
    }
    if(display.PollCloseEvent()) {
      return false;
    }
    copy_small_color(adaptive_pixel_sampler.rgb, adaptive_pixel_sampler_small.rgb);
    copy_small_color(adaptive_pixel_sampler.normalOutput, adaptive_pixel_sampler_small.normalOutput);
    copy_small_color(adaptive_pixel_sampler.albedoOutput, adaptive_pixel_sampler_small.albedoOutput);
    copy_small_preview();
    adaptive_pixel_sampler.max_s = std::max(adaptive_pixel_sampler.max_s, s + 1);
#ifdef HAS_OIDN
    if(denoise_fast_preview) {
      if(display.PollCloseEvent()) {
        return true;
      }
      ensure_fast_preview_oidn_aux();
      if(display.PollCloseEvent()) {
        return true;
      }
      fast_preview_denoiser.Execute();
      if(!fast_preview_denoiser.ReportError()) {
        copy_small_color(adaptive_pixel_sampler.draw_rgb_output,
                         adaptive_pixel_sampler_small.draw_rgb_output);
        display.MarkDenoisedPreviewReady(s + 1);
      }
    }
#endif
    return true;
  };

  auto reset_render_state = [&]() {
    adaptive_pixel_sampler.reset();
    adaptive_pixel_sampler_small.reset();
    display.ResetPreviewExposure();
#ifdef HAS_OIDN
    display.InvalidateOidnAux();
#endif
    reset_sampler_state(nx, ny, seeds, rngs, samplers);
    reset_sampler_state(nx_small, ny_small, seeds_small, rngs_small, samplers_small);
  };

  bool termination_sample_count_set = false;
  auto finish_preview_termination = [&](size_t rendered_samples) {
    if(display.preview && !termination_sample_count_set) {
      adaptive_pixel_sampler.ns = std::max<size_t>(rendered_samples, 1);
      adaptive_pixel_sampler.max_s = adaptive_pixel_sampler.ns;
      termination_sample_count_set = true;
    }
  };

  print_time(verbose, "Allocating sampler" );
  if(display.deferred_render) {
    size_t preview_sample = 0;
    bool completed_final_render = false;
    while(!display.terminate && !completed_final_render) {
      while(!display.render_requested && !display.terminate) {
        Rcpp::checkUserInterrupt();
        if(display.PollCloseEvent()) {
          finish_preview_termination(preview_sample);
          break;
        }
        bool rendered_sample = false;
        if(!display.write_fast_output) {
          rendered_sample = render_full_sample(preview_sample);
          if(!rendered_sample) {
            finish_preview_termination(preview_sample);
            break;
          }
          if(display.PollCloseEvent()) {
            finish_preview_termination(preview_sample + 1);
            break;
          }
#ifdef HAS_OIDN
          ensure_full_preview_oidn_aux();
          if(display.PollCloseEvent()) {
            finish_preview_termination(preview_sample + 1);
            break;
          }
#endif
        } else {
          rendered_sample = render_fast_preview_sample(preview_sample);
        }
        if(display.PollCloseEvent()) {
          finish_preview_termination(preview_sample + (rendered_sample ? 1 : 0));
          break;
        }
        display.DrawImage(adaptive_pixel_sampler, adaptive_pixel_sampler_small,
                          preview_sample, pb, false,
                          0, &world, rng_interactive);
        if(!display.render_requested && !display.terminate) {
          preview_sample++;
        }
      }
      if(display.terminate) {
        finish_preview_termination(preview_sample + 1);
        break;
      }

      display.write_fast_output = false;
      reset_render_state();

      for(size_t s = 0; s < static_cast<size_t>(ns); s++) {
        Rcpp::checkUserInterrupt();
        if(display.PollCloseEvent()) {
          finish_preview_termination(s);
          break;
        }
        if(progress_bar && !display.preview) {
          pb.tick();
        }
        bool rendered_sample = false;
        if(!display.write_fast_output) {
          rendered_sample = render_full_sample(s);
          if(!rendered_sample) {
            finish_preview_termination(s);
            break;
          }
          if(display.PollCloseEvent()) {
            finish_preview_termination(s + 1);
            break;
          }
#ifdef HAS_OIDN
          ensure_full_preview_oidn_aux();
          if(display.PollCloseEvent()) {
            finish_preview_termination(s + 1);
            break;
          }
#endif
        } else {
          rendered_sample = render_fast_preview_sample(s);
        }
        if(display.PollCloseEvent()) {
          finish_preview_termination(s + (rendered_sample ? 1 : 0));
          break;
        }
        display.DrawImage(adaptive_pixel_sampler, adaptive_pixel_sampler_small,
                          s, pb, progress_bar,
                          (Float)s/(Float)ns, &world, rng_interactive);
        if(display.terminate && display.preview) {
          finish_preview_termination(s + 1);
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
      if(display.PollCloseEvent()) {
        finish_preview_termination(s);
        break;
      }
      if(progress_bar && !display.preview) {
        pb.tick();
      }
      bool rendered_sample = false;
      if(!display.write_fast_output) {
        rendered_sample = render_full_sample(s);
        if(!rendered_sample) {
          finish_preview_termination(s);
          break;
        }
        if(display.PollCloseEvent()) {
          finish_preview_termination(s + 1);
          break;
        }
#ifdef HAS_OIDN
        ensure_full_preview_oidn_aux();
        if(display.PollCloseEvent()) {
          finish_preview_termination(s + 1);
          break;
        }
#endif
      } else {
        rendered_sample = render_fast_preview_sample(s);
      }
      if(display.PollCloseEvent()) {
        finish_preview_termination(s + (rendered_sample ? 1 : 0));
        break;
      }
      display.DrawImage(adaptive_pixel_sampler, adaptive_pixel_sampler_small,
                        s, pb, progress_bar,
                        (Float)s/(Float)ns, &world, rng_interactive);
      if(display.terminate && display.preview) {
        finish_preview_termination(s + 1);
        break;
      }
    }
  }

  #ifdef RAY_COLOR_DEBUG
  Rcpp::Rcout << "ns: " << adaptive_pixel_sampler.ns << " adaptive_pixel_sampler.rgb pre divide:\n";
  adaptive_pixel_sampler.rgb.print();
  #endif
  // Optional headless capture exercises the same RGBA snapshot path as the UI.
  if(const char* capture=std::getenv("RAYRENDER_VOLUME_SNAPSHOT")) {
    if(*capture && display.transparent_volume_background) {
      display.SetSnapshotFilename(capture);
      display.CaptureVolumeSnapshot(adaptive_pixel_sampler,adaptive_pixel_sampler.rgb,
                                    adaptive_pixel_sampler.ns,&world,rng_interactive);
      display.SavePreviewSnapshot();
    }
  }
  adaptive_pixel_sampler.write_final_pixels();
#ifdef HAS_OIDN
  if(display.terminate && display.HasDenoisedPreview()) {
    adaptive_pixel_sampler.write_final_denoised_pixels(display.DenoisedPreviewSampleCount());
  }
#endif
}
