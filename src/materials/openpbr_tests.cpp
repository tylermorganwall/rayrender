#ifdef NOT_CRAN
#include "openpbr.h"
#include "../volumes/boundary.h"
#include <testthat.h>

namespace {
Rcpp::List openpbr_test_parameters() {
  Rcpp::Function constructor = Rcpp::Environment::namespace_env("rayrender")["openpbr"];
  Rcpp::List materials = constructor();
  Rcpp::List descriptor = materials[0];
  return Rcpp::clone(Rcpp::as<Rcpp::List>(descriptor["openpbr"]));
}
hit_record openpbr_test_hit() {
  hit_record h;
  h.p = point3f(0);
  h.normal = h.geometric_normal = h.physical_shading_normal = normal3f(0, 0, 1);
  h.dpdu = vec3f(1, 0, 0);
  h.dpdv = vec3f(0, 1, 0);
  h.has_bump = false;
  h.u = h.v = 0;
  return h;
}
}

context("OpenPBR BSDF integration") {
  test_that("joint lobe samples agree with evaluation and directional PDFs") {
    auto parameters = openpbr_test_parameters();
    parameters["specular_roughness"] = .35;
    random_gen rng(30926);
    const auto h = openpbr_test_hit();
    for (int variant = 0; variant < 6; ++variant) {
      auto p = Rcpp::clone(parameters);
      if (variant == 1) p["base_metalness"] = 1.;
      if (variant == 2) { p["coat_weight"] = .7; p["coat_roughness"] = .2; }
      if (variant == 3) p["fuzz_weight"] = .6;
      if (variant == 4) { p["transmission_weight"] = 1.; p["specular_roughness_anisotropy"] = .5; }
      if (variant == 5) { p["thin_film_weight"] = 1.; p["base_metalness"] = 1.; }
      OpenPBRMaterial material(p, std::make_shared<constant_texture>(point3f(.8)), nullptr, false);
      for (Float side : {-1.f, 1.f}) {
        if (side < 0 && variant != 4) continue;
        Ray ray(point3f(0, 0, side), unit_vector(vec3f(.3, 0, -side)));
        auto bsdf = material.Prepare(ray, h);
        int valid = 0;
        for (int i = 0; i < 1000; ++i) {
          Float branch = rng.unif_rand(), u = rng.unif_rand(), v = rng.unif_rand();
          auto sample = bsdf.Sample(branch, u, v);
          if (sample.pdf == 0) continue;
          ++valid;
          expect_true(std::isfinite(sample.pdf));
          expect_true(std::abs(sample.direction.length() - 1) < 1e-5);
          if (sample.specular) continue;
          Float density = bsdf.Pdf(sample.direction);
          auto value = bsdf.Evaluate(sample.direction);
          expect_true(std::abs(density / sample.pdf - 1) < .003);
          for (int c = 0; c < 3; ++c) {
            expect_true((std::isfinite(sample.weight[c]) && sample.weight[c] >= 0));
            expect_true(std::abs(value[c] / sample.pdf - sample.weight[c]) <
                        .003 * std::max(Float(1), sample.weight[c]));
          }
        }
        expect_true(valid > 500);
      }
    }
  }

  test_that("opaque layered white materials conserve furnace energy") {
    auto parameters = openpbr_test_parameters();
    const auto h = openpbr_test_hit();
    random_gen rng(9876);
    for (int variant = 0; variant < 4; ++variant) {
      auto p = Rcpp::clone(parameters);
      if (variant == 1) { p["base_metalness"] = 1.; p["specular_roughness"] = .7; }
      if (variant == 2) { p["coat_weight"] = 1.; p["coat_roughness"] = .3; }
      if (variant == 3) p["fuzz_weight"] = 1.;
      OpenPBRMaterial material(p, std::make_shared<constant_texture>(point3f(1)), nullptr, false);
      auto bsdf = material.Prepare(Ray(point3f(0, 0, 1), vec3f(0, 0, -1)), h);
      point3f energy(0);
      for (int i = 0; i < 10000; ++i) {
        Float branch = rng.unif_rand(), u = rng.unif_rand(), v = rng.unif_rand();
        energy += bsdf.Sample(branch, u, v).weight / 10000;
      }
      for (int c = 0; c < 3; ++c) {
        expect_true(energy[c] > .85);
        expect_true(energy[c] < 1.03);
      }
    }
  }

  test_that("interior coefficients implement Beer attenuation in world units") {
    auto p = openpbr_test_parameters();
    p["transmission_weight"] = 1.;
    p["transmission_depth"] = 2.;
    p["transmission_color"] = Rcpp::NumericVector::create(.8, .5, .2);
    const auto volume = OpenPBRInterior(p);
    const double expected[] = {.8, .5, .2};
    for (int c = 0; c < 3; ++c) {
      expect_true(volume.scattering[c] == 0);
      expect_true(std::abs(std::exp(-2 * volume.absorption[c]) - expected[c]) < 1e-6);
    }
  }

  test_that("raw roughness bypasses the legacy remapping curve") {
    unsigned char pixels[] = {128};
    roughness_texture map(pixels, 1, 1, 1);
    expect_true(std::abs(map.raw_value(.5, .5)[0] - 128.f / 255) < 1e-6);
    expect_true(std::abs(map.value(.5, .5)[0] - map.raw_value(.5, .5)[0]) > .01);
  }
}
#endif
