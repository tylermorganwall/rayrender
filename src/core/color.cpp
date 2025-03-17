#include "../core/color.h"
#include "../core/RcppThread.h"
#include "../math/mathinline.h"
#include "../utils/raylog.h"

// #include "fstream"
// #define DEBUG

// Basic path tracing without importance sampling
void color_basic(const Ray &r, hitable *world, size_t max_depth,
                 random_gen &rng, Sampler *sampler, bool &alpha, point3f &color,
                 normal3f &normal, point3f &albedo) {
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
    #ifdef RAY_COLOR_DEBUG
    Rcpp::Rcout << i << "th ray: O[" <<  r2.origin() << "] Color: [" << throughput << "]\n";
    #endif
    bool is_invisible = false;
    hit_record hrec;
    START_TIMER("Total Hits");
    if (world->hit(r2, 0.001, MaxT, hrec,
                   rng)) { // generated hit record, world space
      STOP_TIMER("Total Hits");
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
      if (throughput.x == 0 && throughput.y == 0 && throughput.z == 0) {
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

        if ((dir.x == 0 && dir.y == 0 && dir.z == 0)) {
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
                              point3f &albedo) {
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
    #ifdef RAY_COLOR_DEBUG
    Rcpp::Rcout << i << "th ray: O[" <<  r2.origin() << "] Color: [" << throughput << "]\n";
    #endif
    bool is_invisible = false;
    hit_record hrec;
    START_TIMER("Total Hits");
    if (world->hit(r2, 0.001, MaxT, hrec,
                   rng)) { // generated hit record, world space
      STOP_TIMER("Total Hits");
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
      if (throughput.x == 0 && throughput.y == 0 && throughput.z == 0) {
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

        if ((dir.x == 0 && dir.y == 0 && dir.z == 0)) [[unlikely]] {
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

void color_shadow_rays(const Ray &r, hitable *world, hitable_list *hlist,
                       size_t max_depth, size_t roulette_activate,
                       random_gen &rng, Sampler *sampler, bool &alpha,
                       point3f &color, normal3f &normal, point3f &albedo) {
  SCOPED_CONTEXT("Overall");
  SCOPED_TIMER_COUNTER("Color");
  point3f final_color(0, 0, 0);
  point3f emit_color(0, 0, 0);

  bool wrote_normal = false;
  bool wrote_albedo = false;
  point3f throughput(1, 1, 1);
  float prev_t = 1;
  Ray r2 = r;
  bool diffuse_bounce = false;
  // To-do: Add logic to detect when rays only go through transmissive surfaces
  // and make those transparent when transparent_background = TRUE
  for (size_t i = 0; i < max_depth; i++) {
    #ifdef RAY_COLOR_DEBUG
    Rcpp::Rcout << i << "th ray: O[" <<  r2.origin() << "] Color: [" << throughput << "]\n";
    #endif
    bool is_invisible = false;
    hit_record hrec;
    START_TIMER("Total Hits");
    if (world->hit(r2, 0.001, MaxT, hrec,
                   rng)) { // generated hit record, world space
      STOP_TIMER("Total Hits");
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
      if (throughput.x == 0 && throughput.y == 0 && throughput.z == 0) {
        color = point3f(0, 0, 0);
        return;
      }
      if (i > roulette_activate) {
        float t = std::fmax(throughput.x,
                            std::fmax(throughput.y, throughput.z));
        // From Szecsi, Szirmay-Kalos, and Kelemen
        float prob_continue = std::fmin(1.0f, std::sqrt(t / prev_t));
        prev_t = t;
        if (rng.unif_rand() > prob_continue) {
          color = final_color;
          return;
        }
        throughput *= 1 / prob_continue;
      }
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
        if (srec.is_specular) {
          // Handle specular reflection/transmission
          r2 = srec.specular_ray;
          throughput *= srec.attenuation;
          diffuse_bounce = false;
          continue;
        }

        // Compute direct lighting via shadow rays
        point3f direct_light(0, 0, 0);
        int num_lights = hlist->size();
        if (num_lights > 0) {
          int li = rng.UniformUInt32(num_lights);
          // Optionally, sample a subset of lights or implement multiple
          // importance sampling for (int li = 0; li < num_lights; li++) {
          hitable *light = hlist->objects[li].get();

          // Sample a point on the light using the random function
          vec3f wi = light->random(hrec.p, sampler, r2.time());
          float pdf_light = light->pdf_value(hrec.p, wi, sampler, r2.time());
          Float p_light_selection = 1.0f / num_lights;

          // Avoid division by zero
          if (pdf_light <= 0.0f) {
            continue;
          }

          // Create a shadow ray towards the light
          Ray shadow_ray(OffsetRayOrigin(hrec.p, hrec.pError, hrec.normal, wi),
                         wi, r2.pri_stack, r2.time());

          // Check for occlusion using HitP
          if (!world->HitP(shadow_ray, 0.001f, MaxT, rng)) {
            // Unoccluded, compute contribution
            float cos_theta = dot(hrec.normal, wi);
            if (cos_theta > 0) {
              // Evaluate BRDF
              point3f f = hrec.mat_ptr->f(r2, hrec, wi);

              // Get emitted radiance from the light material
              bool is_invisible = false;
              point3f emitted_radiance = light->mat_ptr->emitted(
                  shadow_ray, hrec, hrec.u, hrec.v, hrec.p, is_invisible);

              // Accumulate direct lighting
              direct_light += f * emitted_radiance * cos_theta /
                              (pdf_light * p_light_selection);
            }
          }
          // }
          // Average over all lights
          direct_light /= float(num_lights);
        }
        final_color += throughput * direct_light;

        // Continue with indirect lighting using your mixture PDFs
        // Sample the BRDF to get a new direction
        vec3f scatter_direction;
        float pdf_scatter;

        // Create a mixture PDF that combines your BRDF sampling and light
        // sampling
        hitable_pdf light_pdf(hlist, hrec.p);
        mixture_pdf mixed_pdf(&light_pdf, srec.pdf_ptr);

        // Generate a new direction using the mixture PDF
        scatter_direction =
            mixed_pdf.generate(sampler, diffuse_bounce, r2.time());
        pdf_scatter = mixed_pdf.value(scatter_direction, sampler, r2.time());

        if (pdf_scatter == 0) {
          break;
        }

        // Evaluate the BRDF
        point3f f = hrec.mat_ptr->f(r2, hrec, scatter_direction);

        // Update throughput
        // float cos_theta = dot(scatter_direction, hrec.normal);
        // if (cos_theta <= 0) {
        //     break;
        // }
        throughput *= f / pdf_scatter;

        // Update the ray
        r2 = Ray(OffsetRayOrigin(hrec.p, hrec.pError, hrec.normal,
                                 scatter_direction),
                 scatter_direction, r2.pri_stack, r2.time());

        diffuse_bounce = true;
        continue;
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
           bool &alpha, IntegratorType type, point3f &color, normal3f &normal,
           point3f &albedo) {
  switch (type) {
  case IntegratorType::Basic: {
    color_basic(r, world, max_depth, rng, sampler, alpha, color, normal,
                albedo);
    return;
  }
  case IntegratorType::BasicPathGuiding: {
    color_basic_path_guiding(r, world, hlist, max_depth, roulette_activate, rng,
                             sampler, alpha, color, normal, albedo);
    return;
  }
  case IntegratorType::ShadowRays: {
    color_shadow_rays(r, world, hlist, max_depth, roulette_activate, rng,
                      sampler, alpha, color, normal, albedo);
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