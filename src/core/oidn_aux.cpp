#include "../core/oidn_aux.h"
#include "../core/render_jobs.h"

#include <atomic>
#include <algorithm>
#include <cmath>
#include <future>
#include <memory>
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
    return {point3f(0), normal3f(0), false};
  }

  Ray current_ray = ray;
  hit_record hrec;
  while(!AuxCancelled(cancel) &&
        world->hit(current_ray, static_cast<Float>(0.001), MaxT, hrec, rng)) {
    if(AuxCancelled(cancel)) {
      return {point3f(0), normal3f(0), false};
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
        PropagateRayDifferentials(current_ray,hrec,pass_ray,false,1,true);
        return TraceOidnFeature(pass_ray,
                                world,
                                rng,
                                options,
                                depth + 1,
                                dielectric_splits,
                                cancel);
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
      if(entering && !reflected_stack.empty()) reflected_stack.pop_back();
      std::vector<dielectric*> transmitted_stack = base_stack;
      if(!entering && current_layer != -1)
        transmitted_stack.erase(transmitted_stack.begin() + current_layer);

      vec3f reflected = Reflect(wi, outward_normal), refracted(0);
      bool has_refracted = RefractFeatureRay(wi, outward_normal, ni_over_nt, &refracted);
      if(!has_refracted) reflected_weight = 1;
      // Offset each branch toward its own side of the interface. Using the
      // incoming direction for reflection can put the guide inside the glass.
      auto trace_branch = [&](const vec3f& direction, std::vector<dielectric*>& stack) {
        Ray next = MakeFeatureRay(OffsetRayOrigin(hrec.p, hrec.pError,
                                                   hrec.normal, direction),
                                  direction, stack, current_ray.time());
        PropagateRayDifferentials(current_ray,hrec,next,dot(direction,hrec.normal)*dot(current_ray.d,hrec.normal)>0,1/ni_over_nt);
        return TraceOidnFeature(next, world, rng, options, depth + 1,
                                dielectric_splits + 1, cancel);
      };
      // Split only the first interfaces to bound work. Beyond the split budget,
      // sample a Fresnel branch and keep following it to a non-delta surface.
      // This estimates the same blend without exponential recursion or a white
      // fallback at the back face of every glass object.
      if(reflected_weight >= 1) return trace_branch(reflected, reflected_stack);
      if(reflected_weight <= 0) return trace_branch(refracted, transmitted_stack);
      if(dielectric_splits >= options.max_dielectric_splits) {
        return rng.unif_rand() < reflected_weight ?
          trace_branch(reflected, reflected_stack) :
          trace_branch(refracted, transmitted_stack);
      }
      OidnFeature reflected_feature = trace_branch(reflected, reflected_stack);
      OidnFeature transmitted_feature = trace_branch(refracted, transmitted_stack);
      if(AuxCancelled(cancel)) return {point3f(0), normal3f(0), false};

      point3f reflected_albedo = reflected_feature.valid ?
        reflected_feature.albedo :
        point3f(0);
      point3f transmitted_albedo = transmitted_feature.valid ?
        transmitted_feature.albedo :
        point3f(0);
      normal3f reflected_normal = reflected_feature.valid ?
        reflected_feature.normal :
        normal3f(0);
      normal3f transmitted_normal = transmitted_feature.valid ?
        transmitted_feature.normal :
        normal3f(0);

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
      PropagateRayDifferentials(current_ray,hrec,reflected_ray,false,1);
      OidnFeature reflected_feature = TraceOidnFeature(reflected_ray,
                                                       world,
                                                       rng,
                                                       options,
                                                       depth + 1,
                                                       dielectric_splits,
                                                       cancel);
      return reflected_feature;
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

  // Guide images use the same completion-driven scheduling as the beauty
  // render, with one worker pool for all feature samples.
  RcppThread::ThreadPool pool(numbercores);
  const Float differential_scale = std::max(Float(.125),
    Float(1 / std::sqrt(double(render_options.samples))));

  for(std::size_t s = 0; s < render_options.samples; s++) {
    render_cancelled.store(false, std::memory_order_relaxed);
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
            ray = cam->get_ray_differential(u,
                               v,
                               convert_to_point3(rand_to_unit(samplers[index]->Get2D())),
                               samplers[index]->Get1D(), 1.f/nx, 1.f/ny, differential_scale);
          } else {
            CameraSample camera_sample({1 - u, 1 - v},
                                       samplers[index]->Get2D(),
                                       samplers[index]->Get1D());
            cam->GenerateRayDifferential(camera_sample, &ray, -1.f/nx, -1.f/ny, differential_scale);
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
    if(!wait_for_render_jobs(futures, render_cancelled,
                             [&poll_cancel] { return poll_cancel && poll_cancel(); })) {
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

#ifdef NOT_CRAN
#include <testthat.h>
#include "../hitables/sphere.h"
#include "../core/color.h"
#include "../core/adaptivesampler.h"
#include "../volumes/volpath.h"
#include "../volumes/boundary.h"

context("Denoising features through specular paths") {
  test_that("nested glass reaches a colored non-specular surface beyond the split budget") {
    Transform identity;
    auto matte = std::make_shared<diffuse_material>(
      std::make_shared<constant_texture>(point3f(.2, .4, .7)));
    for (Float ior : {Float(1), Float(1.5)}) {
      auto glass = std::make_shared<dielectric>(point3f(1), ior, point3f(0), 1);
      auto inner = std::make_shared<dielectric>(point3f(1), Float(1.2), point3f(0), 2);
      hitable_list world;
      world.add(std::make_shared<sphere>(10, matte, nullptr, nullptr, &identity, &identity, false));
      world.add(std::make_shared<sphere>(1, glass, nullptr, nullptr, &identity, &identity, false));
      world.add(std::make_shared<sphere>(.5, inner, nullptr, nullptr, &identity, &identity, false));
      random_gen rng(513);
      std::vector<dielectric*> stack;
      Ray ray(point3f(0, 0, -4), vec3f(0, 0, 1), &stack);
      OidnAuxRenderOptions options;
      options.max_depth = 64;
      for (size_t splits : {size_t(0), size_t(1)}) {
        options.max_dielectric_splits = splits;
        for (int i = 0; i < 100; ++i) {
          auto feature = TraceOidnFeature(ray, &world, rng, options, 0, 0, nullptr);
          expect_true(feature.valid);
          expect_true((feature.albedo - point3f(.2, .4, .7)).length() < 1e-5);
          expect_true(std::abs(feature.normal[2]) > .85);
        }
      }
      options.max_depth = 1;
      auto exhausted = TraceOidnFeature(ray, &world, rng, options, 0, 0, nullptr);
      expect_true(exhausted.albedo.length() == 0);
    }
  }

  test_that("mirror guides use the reflected surface and unresolved paths have no guide") {
    Transform identity;
    auto mirror = std::make_shared<metal>(std::make_shared<constant_texture>(point3f(1)),
                                         0, point3f(1), point3f(1));
    hitable_list world;
    world.add(std::make_shared<sphere>(1, mirror, nullptr, nullptr, &identity, &identity, false));
    random_gen rng(93);
    Ray ray(point3f(0, 0, -4), vec3f(0, 0, 1));
    OidnAuxRenderOptions options;
    auto missing = TraceOidnFeature(ray, &world, rng, options, 0, 0, nullptr);
    expect_false(missing.valid);
    auto matte = std::make_shared<diffuse_material>(std::make_shared<constant_texture>(point3f(.3)));
    world.add(std::make_shared<sphere>(10, matte, nullptr, nullptr, &identity, &identity, false));
    auto reflected = TraceOidnFeature(ray, &world, rng, options, 0, 0, nullptr);
    expect_true(reflected.valid);
    expect_true((reflected.albedo - point3f(.3)).length() < 1e-6);
    expect_true(reflected.normal[2] < -.999);
  }

  test_that("all beauty integrators defer guides until after both glass faces") {
    Transform identity;
    auto glass = std::make_shared<dielectric>(point3f(1), 1, point3f(0), 1);
    auto matte = std::make_shared<diffuse_material>(std::make_shared<constant_texture>(point3f(.2, .4, .7)));
    hitable_list world, lights;
    world.add(std::make_shared<sphere>(1, glass, nullptr, nullptr, &identity, &identity, false));
    world.add(std::make_shared<sphere>(10, matte, nullptr, nullptr, &identity, &identity, false));
    lights.add(world.objects.back());
    for (auto integrator : {IntegratorType::Basic, IntegratorType::BasicPathGuiding,
                            IntegratorType::ShadowRays}) {
      random_gen rng(17);
      RandomSampler sampler(rng);
      std::vector<dielectric*> stack;
      Ray ray(point3f(0, 0, -4), vec3f(0, 0, 1), &stack);
      Float alpha;
      point3f radiance, albedo;
      normal3f normal;
      color(ray, &world, &lights, 5, 5, rng, &sampler, alpha, integrator,
            radiance, normal, albedo);
      expect_true((albedo - point3f(.2, .4, .7)).length() < 1e-6);
      expect_true(normal[2] > .999);
    }
  }

  test_that("rough metals remain first-hit guides even though they return a specular ray") {
    Transform identity;
    auto rough = std::make_shared<metal>(std::make_shared<constant_texture>(point3f(.3)),
                                        .2, point3f(1), point3f(1));
    hitable_list world, lights;
    world.add(std::make_shared<sphere>(1, rough, nullptr, nullptr, &identity, &identity, false));
    for (auto integrator : {IntegratorType::Basic, IntegratorType::BasicPathGuiding,
                            IntegratorType::ShadowRays}) {
      random_gen rng(17);
      RandomSampler sampler(rng);
      std::vector<dielectric*> stack;
      Ray ray(point3f(0, 0, -4), vec3f(0, 0, 1), &stack);
      Float alpha;
      point3f radiance, albedo;
      normal3f normal;
      color(ray, &world, &lights, 2, 5, rng, &sampler, alpha, integrator,
            radiance, normal, albedo);
      expect_true((albedo - point3f(.3)).length() < 1e-6);
      expect_true(normal[2] < -.999);
    }
  }

  test_that("volume and diffusion guides survive a preceding priority dielectric") {
    Transform identity, placement = Translate(vec3f(0, 0, 3)), inverse = Inverse(placement);
    Rcpp::NumericMatrix transform(4, 4);
    for (int i = 0; i < 4; ++i) transform(i, i) = 1;
    for (bool diffusion : {false, true}) {
      auto medium = std::make_shared<Medium>(Rcpp::List::create(
        Rcpp::Named("sigma_a") = Rcpp::NumericVector::create(8, 6, 3),
        Rcpp::Named("sigma_s") = Rcpp::NumericVector::create(2, 4, 7),
        Rcpp::Named("density_scale") = 1, Rcpp::Named("g") = 0,
        Rcpp::Named("emission") = Rcpp::NumericVector::create(0, 0, 0),
        Rcpp::Named("temperature") = R_NilValue, Rcpp::Named("emission_scale") = 1,
        Rcpp::Named("temperature_scale") = 1, Rcpp::Named("temperature_offset") = 0,
        Rcpp::Named("medium_transform") = transform));
      medium->subsurface = diffusion;
      medium->subsurface_diffusion = diffusion;
      medium->subsurface_ior = 1;
      medium->diffusion_color = {.2, .4, .7};
      medium->diffusion_radius = {.05, .05, .05};
      auto glass = std::make_shared<dielectric>(point3f(1), 1, point3f(0), 0);
      auto body = std::make_shared<dielectric>(point3f(1), 1, point3f(0), 1);
      auto glass_shape = std::make_shared<sphere>(1, glass, nullptr, nullptr, &identity, &identity, false);
      auto body_shape = std::make_shared<sphere>(.5, body, nullptr, nullptr, &placement, &inverse, false);
      auto scene = std::make_shared<VolumeScene>();
      auto glass_boundary = std::make_shared<MediumBoundary>(glass_shape, nullptr, identity, true, scene->NextBoundaryId());
      auto body_boundary = std::make_shared<MediumBoundary>(body_shape, medium, placement, true, scene->NextBoundaryId());
      scene->boundaries.add(glass_boundary);
      scene->boundaries.add(body_boundary);
      scene->Finish(0, 1);
      hitable_list world, lights;
      world.add(glass_boundary);
      world.add(body_boundary);
      lights.volume_scene = scene;
      random_gen rng(72);
      RandomSampler sampler(rng);
      int scattering_guides = 0;
      for (int i = 0; i < 50; ++i) {
        point3f radiance, albedo;
        normal3f normal;
        Float alpha;
        color_volume(Ray(point3f(0, 0, -4), vec3f(0, 0, 1)), &world, &lights,
                     8, 8, rng, &sampler, alpha, radiance, normal, albedo, nullptr);
        if (albedo.length() > 0) {
          ++scattering_guides;
          expect_true((albedo - point3f(.2, .4, .7)).length() < 1e-6);
          expect_true(normal.length() == 0);
        }
      }
      expect_true(scattering_guides > 0);
    }
  }

  test_that("adaptive guide copies average active samples without rescaling finished pixels") {
    RayMatrix rgb(2, 1, 3), rgb2(2, 1, 3), normals(2, 1, 3), albedos(2, 1, 3),
              alpha(2, 1, 1), draw(2, 1, 3), normal_copy(2, 1, 3), albedo_copy(2, 1, 3);
    adaptive_sampler sampler(1, 2, 1, 8, 0, 0, 1, rgb, rgb2, normals, albedos, alpha, draw, true);
    sampler.max_s = 4;
    sampler.finalized[0] = true;
    normals(0, 0, 2) = 1;
    albedos(0, 0, 0) = .2;
    normals(1, 0, 2) = -4;
    albedos(1, 0, 0) = .8;
    auto* storage = albedo_copy.begin();
    sampler.copy_denoising_features(normal_copy, albedo_copy);
    expect_true(storage == albedo_copy.begin());
    expect_true(std::abs(albedo_copy(0, 0, 0) - .2) < 1e-6);
    expect_true(std::abs(albedo_copy(1, 0, 0) - .2) < 1e-6);
    expect_true(normal_copy(0, 0, 2) == 1);
    expect_true(normal_copy(1, 0, 2) == -1);
    expect_true(normals(1, 0, 2) == -4);
  }
}
#endif
