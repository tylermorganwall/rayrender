#include "../core/oidn_aux.h"

#include <atomic>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <future>
#include <memory>
#include <thread>
#include <vector>

#include "RcppThread.h"

#include "../materials/material.h"
#include "../math/mathinline.h"
#include "../math/sampler.h"

namespace {

struct OidnFeature {
  point3f albedo;
  normal3f normal;
  bool valid;
};

inline point3f ClampOidnAlbedo(const point3f& value) {
  return point3f(clamp(value.xyz.x, static_cast<Float>(0), static_cast<Float>(1)),
                 clamp(value.xyz.y, static_cast<Float>(0), static_cast<Float>(1)),
                 clamp(value.xyz.z, static_cast<Float>(0), static_cast<Float>(1)));
}

inline bool AuxCancelled(const std::atomic<bool>* cancel) {
  return cancel != nullptr && cancel->load(std::memory_order_relaxed);
}

inline unsigned int NextRSeed() {
  return static_cast<unsigned int>(unif_rand() * std::pow(2, 32));
}

inline normal3f FeatureNormal(const hit_record& hrec) {
  return !hrec.has_bump ? hrec.normal : hrec.bump_normal;
}

inline normal3f BlendNormals(const normal3f& reflected,
                             const normal3f& transmitted,
                             Float reflected_weight) {
  vec3f blended = reflected_weight * convert_to_vec3(reflected) +
    (static_cast<Float>(1) - reflected_weight) * convert_to_vec3(transmitted);
  return normal3f(clamp(blended.xyz.x, static_cast<Float>(-1), static_cast<Float>(1)),
                  clamp(blended.xyz.y, static_cast<Float>(-1), static_cast<Float>(1)),
                  clamp(blended.xyz.z, static_cast<Float>(-1), static_cast<Float>(1)));
}

inline bool RefractFeatureRay(const vec3f& wi,
                              const normal3f& n,
                              Float eta,
                              vec3f* wt) {
  if(eta == 1) {
    *wt = -wi;
    return true;
  }
  Float cosThetaI = dot(n, wi);
  Float sin2ThetaI = std::fmax(static_cast<Float>(0),
                               static_cast<Float>(1) - cosThetaI * cosThetaI);
  Float sin2ThetaT = eta * eta * sin2ThetaI;
  if(sin2ThetaT >= 1) {
    return false;
  }
  Float cosThetaT = std::sqrt(static_cast<Float>(1) - sin2ThetaT);
  *wt = eta * -wi + (eta * cosThetaI - cosThetaT) * convert_to_vec3(n);
  return true;
}

std::vector<dielectric*> CopyPriorityStack(const Ray& ray) {
  if(ray.pri_stack == nullptr) {
    return std::vector<dielectric*>();
  }
  return *ray.pri_stack;
}

Ray MakeFeatureRay(const point3f& origin,
                   const vec3f& direction,
                   std::vector<dielectric*>& priority_stack,
                   Float time) {
  return Ray(origin, direction, &priority_stack, time);
}

OidnFeature TraceOidnFeature(const Ray& ray,
                             hitable* world,
                             random_gen& rng,
                             const OidnAuxRenderOptions& options,
                             std::size_t depth,
                             std::size_t dielectric_splits,
                             const std::atomic<bool>* cancel) {
  if(AuxCancelled(cancel) || depth >= options.max_depth) {
    return {point3f(1, 1, 1), normal3f(0, 0, 0), false};
  }

  Ray current_ray = ray;
  hit_record hrec;
  while(!AuxCancelled(cancel) &&
        world->hit(current_ray, static_cast<Float>(0.001), MaxT, hrec, rng)) {
    if(AuxCancelled(cancel)) {
      return {point3f(1, 1, 1), normal3f(0, 0, 0), false};
    }
    if(hrec.alpha_miss) {
      current_ray.o = OffsetRayOrigin(hrec.p,
                                      hrec.pError,
                                      hrec.normal,
                                      current_ray.direction());
      continue;
    }

    normal3f feature_normal = FeatureNormal(hrec);
    material* mat = hrec.mat_ptr;
    if(mat == nullptr) {
      return {point3f(0, 0, 0), feature_normal, true};
    }

    if(mat->is_dielectric()) {
      dielectric* dielectric_mat = static_cast<dielectric*>(mat);
      std::vector<dielectric*> base_stack = CopyPriorityStack(current_ray);

      size_t active_priority_value = dielectric_mat->priority;
      size_t next_down_priority = 100000;
      int current_layer = -1;
      int prev_active = -1;
      bool skip = false;
      for(size_t i = 0; i < base_stack.size(); i++) {
        if(base_stack[i] == dielectric_mat) {
          current_layer = static_cast<int>(i);
          continue;
        }
        if(base_stack[i]->priority < active_priority_value) {
          active_priority_value = base_stack[i]->priority;
          skip = true;
        }
        if(base_stack[i]->priority < next_down_priority &&
           base_stack[i] != dielectric_mat) {
          prev_active = static_cast<int>(i);
          next_down_priority = base_stack[i]->priority;
        }
      }

      bool entering = dot(hrec.normal, current_ray.direction()) < 0;
      if(entering) {
        base_stack.push_back(dielectric_mat);
      }

      point3f offset_p = OffsetRayOrigin(hrec.p,
                                         hrec.pError,
                                         hrec.normal,
                                         current_ray.direction());
      if(skip) {
        if(!entering && current_layer != -1) {
          base_stack.erase(base_stack.begin() + static_cast<size_t>(current_layer));
        }
        Ray pass_ray = MakeFeatureRay(offset_p,
                                      current_ray.direction(),
                                      base_stack,
                                      current_ray.time());
        return TraceOidnFeature(pass_ray,
                                world,
                                rng,
                                options,
                                depth + 1,
                                dielectric_splits,
                                cancel);
      }

      if(dielectric_splits >= options.max_dielectric_splits) {
        return {point3f(1, 1, 1), feature_normal, true};
      }

      Float current_ref_idx = prev_active != -1 ?
        base_stack[static_cast<size_t>(prev_active)]->ref_idx :
        static_cast<Float>(1);
      normal3f outward_normal = entering ? feature_normal : -feature_normal;
      Float ni_over_nt = entering ?
        current_ref_idx / dielectric_mat->ref_idx :
        dielectric_mat->ref_idx / current_ref_idx;
      vec3f wi = -unit_vector(current_ray.direction());
      Float reflected_weight = FrDielectric(dot(wi, outward_normal), ni_over_nt);

      std::vector<dielectric*> reflected_stack = base_stack;
      if(entering && !reflected_stack.empty()) {
        reflected_stack.pop_back();
      }
      vec3f reflected = Reflect(wi, outward_normal);
      Ray reflected_ray = MakeFeatureRay(offset_p,
                                         reflected,
                                         reflected_stack,
                                         current_ray.time());
      OidnFeature reflected_feature = TraceOidnFeature(reflected_ray,
                                                       world,
                                                       rng,
                                                       options,
                                                       depth + 1,
                                                       dielectric_splits + 1,
                                                       cancel);
      if(AuxCancelled(cancel)) {
        return {point3f(1, 1, 1), feature_normal, false};
      }

      std::vector<dielectric*> transmitted_stack = base_stack;
      if(!entering && current_layer != -1) {
        transmitted_stack.erase(transmitted_stack.begin() +
                                static_cast<size_t>(current_layer));
      }
      vec3f refracted(0, 0, 0);
      bool has_refracted = RefractFeatureRay(wi,
                                             outward_normal,
                                             ni_over_nt,
                                             &refracted);
      if(!has_refracted) {
        return reflected_feature.valid ?
          reflected_feature :
          OidnFeature{point3f(1, 1, 1), feature_normal, true};
      }

      Ray transmitted_ray = MakeFeatureRay(offset_p,
                                           refracted,
                                           transmitted_stack,
                                           current_ray.time());
      OidnFeature transmitted_feature = TraceOidnFeature(transmitted_ray,
                                                         world,
                                                         rng,
                                                         options,
                                                         depth + 1,
                                                         dielectric_splits + 1,
                                                         cancel);
      if(AuxCancelled(cancel)) {
        return {point3f(1, 1, 1), feature_normal, false};
      }

      point3f reflected_albedo = reflected_feature.valid ?
        reflected_feature.albedo :
        point3f(1, 1, 1);
      point3f transmitted_albedo = transmitted_feature.valid ?
        transmitted_feature.albedo :
        point3f(1, 1, 1);
      normal3f reflected_normal = reflected_feature.valid ?
        reflected_feature.normal :
        feature_normal;
      normal3f transmitted_normal = transmitted_feature.valid ?
        transmitted_feature.normal :
        feature_normal;

      point3f blended_albedo =
        reflected_weight * reflected_albedo +
        (static_cast<Float>(1) - reflected_weight) * transmitted_albedo;
      normal3f blended_normal = BlendNormals(reflected_normal,
                                             transmitted_normal,
                                             reflected_weight);
      return {ClampOidnAlbedo(blended_albedo), blended_normal, true};
    }

    if(mat->is_delta_specular()) {
      vec3f wi = -unit_vector(current_ray.direction());
      vec3f reflected = Reflect(wi, feature_normal);
      std::vector<dielectric*> reflected_stack = CopyPriorityStack(current_ray);
      Ray reflected_ray = MakeFeatureRay(OffsetRayOrigin(hrec.p,
                                                         hrec.pError,
                                                         hrec.normal,
                                                         reflected),
                                         reflected,
                                         reflected_stack,
                                         current_ray.time());
      OidnFeature reflected_feature = TraceOidnFeature(reflected_ray,
                                                       world,
                                                       rng,
                                                       options,
                                                       depth + 1,
                                                       dielectric_splits,
                                                       cancel);
      if(reflected_feature.valid) {
        return reflected_feature;
      }
    }

    return {ClampOidnAlbedo(mat->get_albedo(hrec)), feature_normal, true};
  }

  return {point3f(0, 0, 0), normal3f(0, 0, 0), false};
}

void BuildAuxSamplers(std::size_t nx,
                      std::size_t ny,
                      const OidnAuxRenderOptions& options,
                      std::vector<random_gen>& rngs,
                      std::vector<std::unique_ptr<Sampler> >& samplers) {
  rngs.clear();
  samplers.clear();
  rngs.reserve(nx * ny);
  samplers.reserve(nx * ny);
  for(std::size_t j = 0; j < ny; j++) {
    for(std::size_t i = 0; i < nx; i++) {
      unsigned int seed = NextRSeed();
      random_gen rng_single(seed);
      rngs.push_back(rng_single);
      if(options.sample_method == 0) {
        samplers.push_back(std::unique_ptr<Sampler>(new RandomSampler(rng_single)));
        samplers.back()->StartPixel(0, 0);
      } else if(options.sample_method == 1) {
        samplers.push_back(std::unique_ptr<Sampler>(
          new StratifiedSampler(options.stratified_x,
                                options.stratified_y,
                                true,
                                5,
                                rng_single)));
        samplers.back()->StartPixel(0, 0);
      } else if(options.sample_method == 2) {
        samplers.push_back(std::unique_ptr<Sampler>(
          new SobolSampler(options.samples, rng_single)));
        samplers.back()->StartPixel(i, j);
      } else {
        samplers.push_back(std::unique_ptr<Sampler>(
          new SobolBlueNoiseSampler(rng_single)));
        samplers.back()->StartPixel(i, j);
      }
      samplers.back()->SetSampleNumber(0);
    }
  }
}

} // namespace

void render_oidn_aux_features(std::size_t numbercores,
                              std::size_t nx,
                              std::size_t ny,
                              RayCamera* cam,
                              Float fov,
                              hitable* world,
                              const OidnAuxRenderOptions& options,
                              RayMatrix& normalOutput,
                              RayMatrix& albedoOutput,
                              std::function<bool()> poll_cancel) {
  if(cam == nullptr || world == nullptr) {
    return;
  }
  if(poll_cancel && poll_cancel()) {
    return;
  }

  normalOutput.reset();
  albedoOutput.reset();

  OidnAuxRenderOptions render_options = options;
  render_options.samples = std::max<std::size_t>(render_options.samples, 1);
  render_options.max_depth = std::max<std::size_t>(render_options.max_depth, 1);
  numbercores = std::max<std::size_t>(numbercores, 1);

  std::vector<random_gen> rngs;
  std::vector<std::unique_ptr<Sampler> > samplers;
  BuildAuxSamplers(nx, ny, render_options, rngs, samplers);

  std::size_t chunk_count = numbercores * numbercores;
  std::size_t nx_chunk = nx / numbercores;
  std::size_t ny_chunk = ny / numbercores;
  std::size_t bonus_x = nx - nx_chunk * numbercores;
  std::size_t bonus_y = ny - ny_chunk * numbercores;
  std::size_t completed_samples = 0;
  std::atomic<bool> render_cancelled(false);

  auto wait_for_aux_jobs = [&render_cancelled, &poll_cancel](
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
      if(poll_cancel && poll_cancel()) {
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

  for(std::size_t s = 0; s < render_options.samples; s++) {
    render_cancelled.store(false, std::memory_order_relaxed);
    RcppThread::ThreadPool pool(numbercores);
    auto worker = [&](int k) {
      std::size_t chunk_x = static_cast<std::size_t>(k) / numbercores;
      std::size_t chunk_y = static_cast<std::size_t>(k) % numbercores;
      std::size_t start_x = chunk_x * nx_chunk;
      std::size_t start_y = chunk_y * ny_chunk;
      std::size_t end_x = (chunk_x + 1) * nx_chunk +
        (chunk_x == numbercores - 1 ? bonus_x : 0);
      std::size_t end_y = (chunk_y + 1) * ny_chunk +
        (chunk_y == numbercores - 1 ? bonus_y : 0);

      for(std::size_t i = start_x; i < end_x; i++) {
        if(render_cancelled.load(std::memory_order_relaxed)) {
          break;
        }
        for(std::size_t j = start_y; j < end_y; j++) {
          if(render_cancelled.load(std::memory_order_relaxed)) {
            break;
          }
          std::size_t index = j + ny * i;
          Ray ray;
          vec2f pixel_sample = samplers[index]->Get2D();
          Float u = (static_cast<Float>(i) + pixel_sample.xy.x) /
            static_cast<Float>(nx);
          Float v = (static_cast<Float>(j) + pixel_sample.xy.y) /
            static_cast<Float>(ny);
          if(fov >= 0) {
            ray = cam->get_ray(u,
                               v,
                               convert_to_point3(rand_to_unit(samplers[index]->Get2D())),
                               samplers[index]->Get1D());
          } else {
            CameraSample camera_sample({1 - u, 1 - v},
                                       samplers[index]->Get2D(),
                                       samplers[index]->Get1D());
            cam->GenerateRay(camera_sample, &ray);
          }

          std::vector<dielectric*> priority_stack;
          ray.pri_stack = &priority_stack;
          OidnFeature feature = TraceOidnFeature(ray,
                                                 world,
                                                  rngs[index],
                                                  render_options,
                                                  0,
                                                  0,
                                                  &render_cancelled);
          if(feature.valid) {
            albedoOutput(i, j, 0) += feature.albedo.xyz.x;
            albedoOutput(i, j, 1) += feature.albedo.xyz.y;
            albedoOutput(i, j, 2) += feature.albedo.xyz.z;
            normalOutput(i, j, 0) += feature.normal.xyz.x;
            normalOutput(i, j, 1) += feature.normal.xyz.y;
            normalOutput(i, j, 2) += feature.normal.xyz.z;
          }
          samplers[index]->StartNextSample();
        }
      }
    };

    std::vector<std::future<void> > futures;
    futures.reserve(chunk_count);
    for(std::size_t k = 0; k < chunk_count; k++) {
      futures.push_back(pool.pushReturn(worker, static_cast<int>(k)));
    }
    if(!wait_for_aux_jobs(futures)) {
      return;
    }
    completed_samples++;
  }

  Float inv_samples = static_cast<Float>(1) /
    static_cast<Float>(std::max<std::size_t>(completed_samples, 1));
  for(std::size_t i = 0; i < nx; i++) {
    for(std::size_t j = 0; j < ny; j++) {
      albedoOutput(i, j, 0) *= inv_samples;
      albedoOutput(i, j, 1) *= inv_samples;
      albedoOutput(i, j, 2) *= inv_samples;
      normalOutput(i, j, 0) *= inv_samples;
      normalOutput(i, j, 1) *= inv_samples;
      normalOutput(i, j, 2) *= inv_samples;
    }
  }
}
