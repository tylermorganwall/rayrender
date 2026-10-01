#ifdef NOT_CRAN
#include "point.h"
#include <testthat.h>

context("Point and spot light sampling") {
  test_that("point intensity obeys inverse-square falloff") {
    PointLight light;
    light.intensity = point3f(4, 8, 12);
    auto near = light.Sample(point3f(0, 0, 2));
    auto far = light.Sample(point3f(0, 0, 4));
    for (int c = 0; c < 3; ++c) {
      expect_true(std::abs(near.radiance[c] - (c + 1)) < 1e-6);
      expect_true(std::abs(near.radiance[c] - 4 * far.radiance[c]) < 1e-6);
    }
    expect_true(light.Sample(light.position).pmf == 0);
  }
  test_that("spot falloff and integrated power match PBRT v4") {
    PointLight light;
    light.spot = true;
    light.direction = vec3f(0, 0, 1);
    light.intensity = point3f(1);
    light.cos_outer = .5;
    light.cos_inner = .75;
    expect_true(light.Falloff(vec3f(0, 0, 1)) == 1);
    expect_true(light.Falloff(vec3f(1, 0, 0)) == 0);
    const double cosine = .625;
    expect_true(std::abs(light.Falloff(vec3f(std::sqrt(1 - cosine*cosine), 0, cosine)) - .5) < 1e-6);
    double integral = 0;
    const int n = 10000;
    for (int i = 0; i < n; ++i) {
      double z = -1 + 2. * (i + .5) / n;
      integral += light.Falloff(vec3f(std::sqrt(1 - z*z), 0, z)) * 4 * M_PI / n;
    }
    expect_true(std::abs(integral - light.Power()) < 1e-5);
    light.cos_inner = light.cos_outer;
    expect_true(light.Falloff(vec3f(0, 0, 1)) == 1);
    expect_true(light.Falloff(vec3f(1, 0, 0)) == 0);
  }
  test_that("power-weighted selection retains the discrete probability") {
    auto descriptor = [](double intensity) {
      return Rcpp::List::create(Rcpp::_["type"] = "point",
        Rcpp::_["position"] = Rcpp::NumericVector::create(0, 0, 0),
        Rcpp::_["color"] = Rcpp::NumericVector::create(1, 1, 1),
        Rcpp::_["intensity"] = intensity);
    };
    PointLightSet lights(Rcpp::List::create(descriptor(0), descriptor(1), descriptor(3)));
    auto a = lights.Sample(point3f(0, 0, 2), .1);
    auto b = lights.Sample(point3f(0, 0, 2), .5);
    expect_true(std::abs(a.pmf - .25) < 1e-6);
    expect_true(std::abs(b.pmf - .75) < 1e-6);
    expect_true(std::abs(a.radiance[0] / a.pmf - 1) < 1e-6);
    expect_true(std::abs(b.radiance[0] / b.pmf - 1) < 1e-6);
    expect_true(PointLightSet(Rcpp::List::create(descriptor(0))).Empty());
  }
}
#endif
