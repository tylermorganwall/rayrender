#include "../core/color.h"
#include "../volumes/volpath.h"
#include "RcppThread.h"
#include "../math/mathinline.h"
#include "../utils/raylog.h"

// #include "fstream"
// #define DEBUG

static inline bool ColorCancelled(const std::atomic<bool>* cancel) {
  return cancel != nullptr && cancel->load(std::memory_order_relaxed);
}

// Basic path tracing without importance sampling
void color_basic(const Ray &r, hitable *world, size_t max_depth,
                 random_gen &rng, Sampler *sampler, bool &alpha, point3f &color,
                 normal3f &normal, point3f &albedo,
                 const std::atomic<bool>* cancel) {
  point3f final_color(0, 0, 0);
  point3f emit_color(0, 0, 0);
  bool wrote_normal = false;
  bool wrote_albedo = false;
  point3f throughput(1, 1, 1);
  Ray r1 = r;
  Ray r2 = r;
  bool diffuse_bounce = false;
  // To-do: Add logic to detect when rays only go through transmissive surfaces
  // and make those transparent when transparent_background = TRUE
  for (size_t i = 0; i < max_depth; i++) {
    if(ColorCancelled(cancel)) {
      color = final_color;
      return;
    }
    #ifdef RAY_COLOR_DEBUG
    Rcpp::Rcout << i << "th ray: O[" <<  r2.origin() << "] Color: [" << throughput << "]\n";
    #endif
    bool is_invisible = false;
    hit_record hrec;
    START_TIMER("Total Hits");
    if (world->hit(r2, 0.001, MaxT, hrec,
                   rng)) { // generated hit record, world space
      STOP_TIMER("Total Hits");
      if(ColorCancelled(cancel)) {
        color = final_color;
        return;
      }
      scatter_record srec;
      if (hrec.alpha_miss) {
        r2.o =
            OffsetRayOrigin(hrec.p, hrec.pError, hrec.normal, r2.direction());
        continue;
      }
      if (hrec.infinite_area_hit && i == 0) {
        alpha = true;
      }
      emit_color = throughput * hrec.mat_ptr->emitted(r2, hrec, hrec.u, hrec.v,
                                                      hrec.p, is_invisible);
      // Some lights can be invisible until after diffuse bounce
      // If so, generate new ray with intersection point and continue ray
      if (is_invisible && !diffuse_bounce) {
        r2.o =
            OffsetRayOrigin(hrec.p, hrec.pError, hrec.normal, r2.direction());
        continue;
      }
      final_color += emit_color;
      if (throughput.xyz.x == 0 && throughput.xyz.y == 0 && throughput.xyz.z == 0) {
        color = (point3f(0, 0, 0));
        return;
      }
      float pdf_val;
      // generates scatter record and sends out new ray, otherwise exits out
      // with accumulated color
      if (hrec.mat_ptr->scatter(r2, hrec, srec, sampler)) {
        if(!wrote_normal) {
          normal = hrec.normal;
          wrote_normal = true;
        }        
        if(!wrote_albedo) {
          albedo = throughput;
          wrote_albedo = true;
        }
        if (srec.is_specular) { // returns specular ray
          r2 = srec.specular_ray;
          throughput *= srec.attenuation;
          continue;
        }

        // Generates a scatter direction (with origin hrec.p) from the surface
        // and saves surface normal from light to use in pdf_value calculation
        //(along with the scatter direction)
        r1 = r2;
        vec3f dir;
        if (!diffuse_bounce) {
          //`diffuse_bounce` switched by generate()
          dir = srec.pdf_ptr->generate(
              sampler, diffuse_bounce,
              r2.time()); // scatters a ray from hit point to stratified
                          // direction
        } else {
          dir = srec.pdf_ptr->generate(
              rng, diffuse_bounce,
              r2.time()); // scatters a ray from hit point to random direction
        }

        r2 = Ray(OffsetRayOrigin(hrec.p, hrec.pError, hrec.normal, dir), dir,
                 r2.pri_stack, r2.time());
        pdf_val = srec.pdf_ptr->value(
            dir, rng, r2.time()); // generates a pdf value based the
                                  // intersection point and the mixture pdf

        if (pdf_val == 0) {
          break;
        }

        if ((dir.xyz.x == 0 && dir.xyz.y == 0 && dir.xyz.z == 0)) {
          break;
        }

        throughput *= hrec.mat_ptr->f(r1, hrec, r2.direction()) / pdf_val;

      } else {
        color = final_color;
        return;
      }
    } else {
      STOP_TIMER("hit");
      color = final_color;
      return;
    }
  }
  color = final_color;
  return;
}

void color_basic_path_guiding(const Ray &r, hitable *world, hitable_list *hlist,
                              size_t max_depth, size_t roulette_activate,
                              random_gen &rng, Sampler *sampler, bool &alpha,
                              point3f &color, normal3f &normal,
                              point3f &albedo,
                              const std::atomic<bool>* cancel) {
  SCOPED_CONTEXT("Overall");
  SCOPED_TIMER_COUNTER("Color");
  point3f final_color(0, 0, 0);
  point3f emit_color(0, 0, 0);
  bool wrote_normal = false;
  bool wrote_albedo = false;

  point3f throughput(1, 1, 1);
  Ray r1 = r;
  Ray r2 = r;
  bool diffuse_bounce = false;
  // To-do: Add logic to detect when rays only go through transmissive surfaces
  // and make those transparent when transparent_background = TRUE
  for (size_t i = 0; i < max_depth; i++) {
    if(ColorCancelled(cancel)) {
      color = final_color;
      return;
    }
    #ifdef RAY_COLOR_DEBUG
    Rcpp::Rcout << i << "th ray: O[" <<  r2.origin() << "] Color: [" << throughput << "]\n";
    #endif
    bool is_invisible = false;
    hit_record hrec;
    START_TIMER("Total Hits");
    if (world->hit(r2, 0.001, MaxT, hrec,
                   rng)) { // generated hit record, world space
      STOP_TIMER("Total Hits");
      if(ColorCancelled(cancel)) {
        color = final_color;
        return;
      }
      scatter_record srec;
      if (hrec.alpha_miss) {
        r2.o =
            OffsetRayOrigin(hrec.p, hrec.pError, hrec.normal, r2.direction());
        continue;
      }
      if (hrec.infinite_area_hit && i == 0) {
        alpha = true;
      }
      emit_color = throughput * hrec.mat_ptr->emitted(r2, hrec, hrec.u, hrec.v,
                                                      hrec.p, is_invisible);
      // Some lights can be invisible until after diffuse bounce
      // If so, generate new ray with intersection point and continue ray
      if (is_invisible && !diffuse_bounce) {
        r2.o =
            OffsetRayOrigin(hrec.p, hrec.pError, hrec.normal, r2.direction());
        continue;
      }
      final_color += emit_color;
      if (throughput.xyz.x == 0 && throughput.xyz.y == 0 && throughput.xyz.z == 0) {
        if(!wrote_normal) [[unlikely]] {
          normal = hrec.normal;
          wrote_normal = true;
        }
        color = point3f(0, 0, 0);
        return;
      }
      float pdf_val;
      // generates scatter record and sends out new ray, otherwise exits out
      // with accumulated color
      if (hrec.mat_ptr->scatter(r2, hrec, srec, sampler)) {
        if(!wrote_normal) {
          normal = hrec.normal;
          wrote_normal = true;
        }        
        if(!wrote_albedo) {
          albedo = hrec.mat_ptr->get_albedo(hrec);
          wrote_albedo = true;
        }
        if (srec.is_specular) { // returns specular ray
          r2 = srec.specular_ray;
          throughput *= srec.attenuation;
          continue;
        }
        hitable_pdf p_imp(hlist,
                          hrec.p); // creates pdf of all objects to be sampled
        mixture_pdf p(
            &p_imp, srec.pdf_ptr); // creates mixture pdf of surface intersected
                                   // at hrec.p and all sampled objects/lights

        // Generates a scatter direction (with origin hrec.p) from the mixture
        // and saves surface normal from light to use in pdf_value calculation
        //(along with the scatter direction)
        r1 = r2;
        vec3f dir;
        if (!diffuse_bounce) {
          //`diffuse_bounce` switched by generate()
          dir = p.generate(sampler, diffuse_bounce,
                           r2.time()); // scatters a ray from hit point to
                                       // stratified direction
        } else {
          dir = p.generate(
              rng, diffuse_bounce,
              r2.time()); // scatters a ray from hit point to random direction
        }

        r2 = Ray(OffsetRayOrigin(hrec.p, hrec.pError, hrec.normal, dir), dir,
                 r2.pri_stack, r2.time());
        pdf_val = p.value(dir, rng,
                          r2.time()); // generates a pdf value based the
                                      // intersection point and the mixture pdf

        if (pdf_val == 0) [[unlikely]] {
          break;
        }

        if ((dir.xyz.x == 0 && dir.xyz.y == 0 && dir.xyz.z == 0)) [[unlikely]] {
          break;
        }

        throughput *= hrec.mat_ptr->f(r1, hrec, r2.direction()) / pdf_val;
      } else {
        color = final_color;
        return;
      }
    } else {
      STOP_TIMER("hit");
      color = final_color;
      return;
    }
  }
  color = final_color;
  return;
}


void color(const Ray &r, hitable *world, hitable_list *hlist, size_t max_depth,
           size_t roulette_activate, random_gen &rng, Sampler *sampler,
           Float &transparency, IntegratorType type, point3f &color, normal3f &normal,
           point3f &albedo, const std::atomic<bool>* cancel) {
  bool alpha = false;
  switch (type) {
  case IntegratorType::Basic: {
    color_basic(r, world, max_depth, rng, sampler, alpha, color, normal,
                albedo, cancel);
    transparency = alpha ? 1 : 0;
    return;
  }
  case IntegratorType::BasicPathGuiding: {
    color_basic_path_guiding(r, world, hlist, max_depth, roulette_activate, rng,
                             sampler, alpha, color, normal, albedo, cancel);
    transparency = alpha ? 1 : 0;
    return;
  }
  case IntegratorType::ShadowRays: {
    color_volume(r, world, hlist, max_depth, roulette_activate, rng,
                 sampler, transparency, color, normal, albedo, cancel);
    return;
  }
  default: {
    color = point3f(0, 0, 0);
    normal = normal3f(0, 0, 0);
    albedo = point3f(0, 0, 0);
    // Handle error or default case
    return;
  }
  }
}
