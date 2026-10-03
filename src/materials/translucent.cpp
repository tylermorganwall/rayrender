#include "material.h"
#include <algorithm>

namespace {
vec3f sheet_normal(const Ray& ray, const hit_record& h) {
  vec3f n = unit_vector(convert_to_vec3(h.geometric_normal));
  return dot(ray.direction(), n) > 0 ? -n : n;
}

// Select a hemisphere by its RGB-average energy, then sample cosine-weighted
// solid angle there. Both value() overloads include that hemisphere probability.
// This keeps f/pdf unbiased for colored, reflecting-only and transmitting-only
// sheets, including next-event estimates from lights on the opposite side.
class translucent_pdf final : public pdf {
public:
  translucent_pdf(const vec3f& n, Float reflection_probability)
    : probability(reflection_probability) { frame.build_from_w(n); }
  Float density(const vec3f& direction) const {
    if (!(direction.squared_length() > 0)) return 0;
    Float cosine = dot(unit_vector(direction), frame.w());
    return std::abs(cosine) * M_1_PI * (cosine > 0 ? probability : 1 - probability);
  }
  Float value(const vec3f& d, random_gen&, Float) override { return density(d); }
  Float value(const vec3f& d, Sampler*, Float) override { return density(d); }
  vec3f generate(random_gen& rng, bool& diffuse, Float) override {
    diffuse = true;
    vec3f d = frame.local_to_world(rng.random_cosine_direction());
    return rng.unif_rand() < probability ? d : -d;
  }
  vec3f generate(Sampler* sampler, bool& diffuse, Float) override {
    diffuse = true;
    vec3f d = frame.local_to_world(rand_cosine_direction(sampler->Get2D()));
    return sampler->Get1D() < probability ? d : -d;
  }
private:
  onb frame;
  Float probability;
};
}

std::array<point3f, 2> translucent_material::colors(const hit_record& h) const {
  std::array<point3f, 2> result = {reflection->value(h),
                                  transmission->value(h)};
  for (int c = 0; c < 3; ++c) {
    for (auto& color : result)
      color[c] = std::isfinite(color[c]) ? std::clamp(color[c], Float(0), Float(1)) : 0;
    // Image inputs may exceed unit combined energy even when constants do not.
    Float sum = result[0][c] + result[1][c];
    if (sum > 1) { result[0][c] /= sum; result[1][c] /= sum; }
  }
  return result;
}

point3f translucent_material::get_albedo(const hit_record& h) const {
  auto color = colors(h);
  return color[0] + color[1];
}

point3f translucent_material::f(const Ray& ray, const hit_record& h, const vec3f& wi) const {
  if (!(wi.squared_length() > 0) || !(h.geometric_normal.squared_length() > 0) ||
      dot(ray.direction(), h.geometric_normal) == 0) return point3f(0);
  Float cosine = dot(unit_vector(wi), sheet_normal(ray, h));
  const Float projected_weight = std::abs(cosine) * Float(M_1_PI);
  return colors(h)[cosine > 0 ? 0 : 1] * projected_weight;
}

bool translucent_material::prepare(const Ray& ray, const hit_record& h, scatter_record& s) const {
  if (!(h.geometric_normal.squared_length() > 0) || dot(ray.direction(), h.geometric_normal) == 0)
    return false;
  auto color = colors(h);
  Float reflected = color[0][0] + color[0][1] + color[0][2];
  Float transmitted = color[1][0] + color[1][1] + color[1][2];
  if (!(reflected + transmitted > 0)) return false;
  s.is_specular = false;
  s.attenuation = color[0] + color[1];
  s.pdf_ptr = new translucent_pdf(sheet_normal(ray, h), reflected / (reflected + transmitted));
  return true;
}
bool translucent_material::scatter(const Ray& r, const hit_record& h, scatter_record& s, random_gen&) {
  return prepare(r, h, s);
}
bool translucent_material::scatter(const Ray& r, const hit_record& h, scatter_record& s, Sampler*) {
  return prepare(r, h, s);
}

#ifdef NOT_CRAN
#include <testthat.h>
context("Diffuse sheet transmission") {
  test_that("both hemispheres conserve energy and sample the full BSDF on either side") {
    hit_record h;
    h.p = point3f(0); h.u = h.v = 0;
    h.normal = h.geometric_normal = normal3f(0, 0, 1);
    for (auto r : {point3f(0), point3f(.2, .3, .1), point3f(.8)})
    for (auto t : {point3f(0), point3f(.6, .2, .8), point3f(.8)}) {
      translucent_material material(std::make_shared<constant_texture>(r),
                                    std::make_shared<constant_texture>(t));
      point3f expected = r + t;
      for (int c = 0; c < 3; ++c) expected[c] = std::min(Float(1), expected[c]);
      for (Float side : {-1.f, 1.f}) {
        Ray ray(point3f(0, 0, side), vec3f(0, 0, -side));
        random_gen rng(791);
        RandomSampler sampler(rng);
        scatter_record scatter;
        bool active = material.scatter(ray, h, scatter, &sampler);
        expect_true(active == (expected[0] + expected[1] + expected[2] > 0));
        if (!active) continue;
        point3f integral(0); double pdf_integral = 0;
        for (int i = 0; i < 200; ++i) {
          Float z = -1 + (i + .5f) / 100;
          vec3f wi(std::sqrt(1 - z*z), 0, z);
          integral += material.f(ray, h, wi) * Float(4 * M_PI / 200);
          pdf_integral += scatter.pdf_ptr->value(wi, rng) * (4 * M_PI / 200);
          // Remove the outgoing cosine before checking BSDF reciprocity.
          auto reverse = material.f(Ray(point3f(0), -wi), h, -ray.direction());
          expect_true((material.f(ray, h, wi) / std::abs(z) - reverse).length() < 1e-5);
        }
        expect_true((integral - expected).length() < 1e-5);
        expect_true(std::abs(pdf_integral - 1) < 1e-5);
        for (bool use_sampler : {false, true}) {
          point3f estimate(0); bool valid = true;
          for (int i = 0; i < 20000; ++i) {
            bool diffuse = false;
            vec3f wi = use_sampler ? scatter.pdf_ptr->generate(&sampler, diffuse) :
                                    scatter.pdf_ptr->generate(rng, diffuse);
            Float density = scatter.pdf_ptr->value(wi, rng);
            valid = valid && diffuse && density > 0 &&
                density == scatter.pdf_ptr->value(wi, &sampler);
            if (density > 0) estimate += material.f(ray, h, wi) / density;
          }
          expect_true(valid);
          expect_true((estimate / 20000 - expected).length() < .015);
        }
      }
    }
  }
}
#endif
