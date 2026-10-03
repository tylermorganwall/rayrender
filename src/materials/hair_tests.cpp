#ifdef NOT_CRAN
#include "material.h"
#include <testthat.h>

namespace {
hit_record HairFixture(Float h = .3f) {
  hit_record hit;
  hit.p = point3f(0);
  hit.normal = hit.geometric_normal = normal3f(0, 0, 1);
  hit.dpdu = vec3f(1, 0, 0);
  hit.dpdv = vec3f(0, .005f, 0);
  hit.v = (h + 1) / 2;
  return hit;
}
vec3f HairDirection(Float x, Float phi) {
  Float r = SafeSqrt(1 - x*x);
  return vec3f(x, r*std::cos(phi), r*std::sin(phi));
}
}

context("Hair scattering contracts") {
  test_that("directional PDFs normalize and agree across APIs") {
    for (Float roughness : {.3f, .7f}) for (Float longitudinal : {0.f, .6f, .95f}) {
      hair mat(point3f(.3f, 1.f, 3.f), 1.55f, roughness, roughness, 2);
      auto hit = HairFixture();
      Ray ray(point3f(0), -HairDirection(longitudinal, .4f));
      random_gen rng(724);
      RandomSampler sampler(rng);
      scatter_record scattering;
      mat.scatter(ray, hit, scattering, rng);
      double integral = 0, first_moment = 0;
      Float smallest = Infinity, largest = 0;
      const int nx = 256, np = 512;
      for (int i = 0; i < nx; ++i) for (int j = 0; j < np; ++j) {
        auto direction = HairDirection(-1 + 2 * (i+.5f)/nx, 2*Float(M_PI)*(j+.5f)/np);
        Float density = scattering.pdf_ptr->value(direction, rng);
        integral += density;
        first_moment += density * direction[0];
        smallest = std::min(smallest, density); largest = std::max(largest, density);
        if (i % 32 == 0 && j % 64 == 0)
          expect_true(density == scattering.pdf_ptr->value(direction, &sampler));
      }
      integral *= 4*M_PI/(nx*np); first_moment *= 4*M_PI/(nx*np);
      expect_true(std::abs(integral - 1) < .008);
      expect_true(largest > smallest * 2);
      double sampled_moment = 0;
      for (int i = 0; i < 32768; ++i) {
        bool diffuse;
        auto wi = scattering.pdf_ptr->generate(rng, diffuse);
        sampled_moment += wi[0];
      }
      expect_true(std::abs(sampled_moment/32768 - first_moment) < .02);
    }
  }

  test_that("white furnace remains white including tangent and edge incidence") {
    for (Float roughness : {.1f, .3f, .8f, 1.f})
    for (Float longitudinal : {0.f, .9f, 1.f, -1.f})
    for (Float h : {-.9f, 0.f, 1.f}) {
      hair mat(point3f(0), 1.55f, roughness, roughness, 2);
      auto hit = HairFixture(h);
      Ray ray(point3f(0), -HairDirection(longitudinal, .7f));
      random_gen rng(942);
      RandomSampler sampler(rng);
      scatter_record scattering;
      mat.scatter(ray, hit, scattering, &sampler);
      double energy = 0;
      bool finite = true;
      for (int i = 0; i < 1024; ++i) {
        bool diffuse;
        vec3f wi = i % 2 ? scattering.pdf_ptr->generate(rng, diffuse) :
                          scattering.pdf_ptr->generate(&sampler, diffuse);
        Float density = scattering.pdf_ptr->value(wi, rng);
        point3f value = mat.f(ray, hit, wi);
        finite &= std::isfinite(density) && density > 0 && std::isfinite(value[0]);
        finite &= std::abs(wi.length() - 1) < 2e-5;
        if (density > 0) energy += value[0] / density;
      }
      expect_true(finite);
      expect_true(std::abs(energy/1024 - 1) < .002);
    }
  }

  test_that("pigmented transport has finite weights and a bounded nonblack guide") {
    hair mat(point3f(.4f, 1.f, 4.f), 1.55f, .3f, .3f, 2);
    auto hit = HairFixture();
    random_gen rng(159);
    Ray ray(point3f(0), -HairDirection(.3f, .5f));
    scatter_record scattering;
    mat.scatter(ray, hit, scattering, rng);
    point3f sum(0);
    bool bounded = true;
    for (int i = 0; i < 32768; ++i) {
      bool diffuse;
      auto wi = scattering.pdf_ptr->generate(rng, diffuse);
      Float density = scattering.pdf_ptr->value(wi, rng);
      auto weight = mat.f(ray, hit, wi) / density;
      for (int k = 0; k < 3; ++k)
        bounded &= std::isfinite(weight[k]) && weight[k] >= 0 && weight[k] <= 3.001f;
      sum += weight;
    }
    expect_true(bounded);
    for (int k = 0; k < 3; ++k) {
      expect_true(sum[k]/32768 < 1.01);
      expect_true((mat.get_albedo(hit)[k] > 0 && mat.get_albedo(hit)[k] < 1));
    }
  }

  test_that("derivative scale and normal tilt cannot change the fiber tangent") {
    hair mat(point3f(.4f, 1.f, 4.f), 1.55f, .3f, .3f, 2);
    auto a = HairFixture(), b = a;
    b.dpdu *= 17; b.dpdv *= .01f;
    b.normal = normal3f(.7f, 0, 1);
    random_gen rng(123);
    Ray ray(point3f(0), -HairDirection(.3f, .7f));
    scatter_record sa, sb;
    mat.scatter(ray, a, sa, rng); mat.scatter(ray, b, sb, rng);
    for (int i = 0; i < 100; ++i) {
      auto wi = HairDirection(-.99f + 1.98f*i/99, .6f);
      expect_true(std::abs(sa.pdf_ptr->value(wi,rng) - sb.pdf_ptr->value(wi,rng)) < 2e-5);
      expect_true((mat.f(ray,a,wi) - mat.f(ray,b,wi)).length() < 2e-5);
    }
    b.dpdu = vec3f(0); b.normal = normal3f(0);
    scatter_record fallback;
    expect_true(mat.scatter(ray,b,fallback,rng));
    expect_true(std::isfinite(fallback.pdf_ptr->value(vec3f(0,0,1),rng)));
  }
}
#endif
