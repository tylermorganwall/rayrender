#ifdef NOT_CRAN
#include "atmosphere.h"
#include "../core/PreviewDisplay.h"
#include <testthat.h>

namespace {
// Deliberately non-multiplicative fit, like independently fitted finite-distance
// tables. Splitting it must not change the total because null rates are arbitrary.
class FittedAtmosphere final : public Atmosphere {
public:
  AtmosphereSegment Segment(const point3f &, const vec3f &, double d) const override {
    AtmosphereSegment s;
    s.transmission = point3f(std::exp(-.01 * d - .0001 * d * d));
    s.radiance = point3f(2, 3, 4) * (1 - s.transmission[0]);
    return s;
  }
  point3f Transmission(const point3f &p, const vec3f &w, double d) const override {
    return Segment(p, w, d).transmission;
  }
  point3f SkyRadiance(const point3f &, const vec3f &) const override { return point3f(2, 3, 4); }
};
}
context("Atmospheric ray subdivision") {
  test_that("null events and invisible boundaries subdivide one fitted integral") {
    FittedAtmosphere model;
    for (int n : {1, 2, 11, 301}) {
      AtmosphereRay ray(&model);
      point3f tr(1), le(0);
      vec3f direction(0, 0, 1);
      for (int i = 0; i < n; ++i) {
        ray.Start(point3f(0, 0, 100.0 * i / n), direction, true);
        auto s = ray.Advance(point3f(0, 0, 100.0 * (i + 1) / n));
        le += tr * s.radiance;
        tr *= s.transmission;
      }
      auto whole = model.Segment(point3f(0), direction, 100);
      expect_true((tr - whole.transmission).length() < 1e-5);
      expect_true((le - whole.radiance).length() < 1e-5);
      // A distant environment includes the remaining haze, counted just once.
      point3f environment(3, 5, 7);
      expect_true((le + tr * ray.Remaining(environment) - environment).length() < 1e-5);
    }
  }
  test_that("glass pauses the atmosphere and a new direction resets its anchor") {
    FittedAtmosphere model;
    AtmosphereRay ray(&model);
    ray.Start(point3f(0), vec3f(0, 0, 1), true);
    ray.Advance(point3f(0, 0, 20));
    ray.Start(point3f(0, 0, 20), vec3f(0, 0, 1), false);
    auto glass = ray.Advance(point3f(0, 0, 40));
    expect_true(glass.transmission[0] == 1);
    expect_true(glass.radiance[0] == 0);
    ray.Start(point3f(0, 0, 40), vec3f(1, 0, 0), true);
    auto air = ray.Advance(point3f(10, 0, 40));
    auto reference = model.Segment(point3f(0, 0, 40), vec3f(1, 0, 0), 10);
    expect_true((air.radiance - reference.radiance).length() < 1e-6);
  }
}

// Optional end-to-end preview export of an actual rendered RGBA fixture. This
// exercises the preview's native premultiplication, orientation, and snapshot
// code without requiring a window server. Ordinary C++ tests need no fixture.
context("Atmospheric preview export") {
  if (const char* fixture = std::getenv("RAYRENDER_PRAGUE_PREVIEW_FIXTURE")) {
    test_that("a rendered atmospheric image retains its alpha in preview snapshots") {
      Rcpp::Function read = Rcpp::Environment::base_env()["readRDS"];
      Rcpp::List input = read(fixture);
      Rcpp::NumericVector pixels = input["fixed"];
      Rcpp::IntegerVector dim = pixels.attr("dim");
      int height = dim[0], width = dim[1], count = height * width;
      expect_true(dim[2] == 4);
      RayMatrix color(width, height, 3), squared(width, height, 3),
        normal(width, height, 3), albedo(width, height, 3), alpha(width, height, 1),
        draw(width, height, 3);
      adaptive_sampler sampler(1, width, height, 1, 0, 0, 8, color, squared,
                               normal, albedo, alpha, draw, false, true);
      std::fill(sampler.finalized.begin(), sampler.finalized.end(), true);
      for (int y = 0; y < height; ++y) for (int x = 0; x < width; ++x) {
        double a = pixels[y + height * x + 3 * count];
        alpha(width - 1 - x, height - 1 - y, 0) = a;
        for (int c = 0; c < 3; ++c)
          color(width - 1 - x, height - 1 - y, c) = pixels[y + height * x + c * count] * a;
      }
      Transform transform;
#ifdef HAS_OIDN
      PreviewDisplay preview(width, height, false, false, false, 1, nullptr,
                              &transform, &transform, nullptr, nullptr, nullptr, false, false);
#else
      PreviewDisplay preview(width, height, false, false, false, 1, nullptr,
                              &transform, &transform, false);
#endif
      preview.transparent_volume_background = true;
      random_gen rng(0);
      preview.CaptureVolumeSnapshot(sampler, color, 1, nullptr, rng);
      double error = 0;
      for (int y = 0; y < height; ++y) for (int x = 0; x < width; ++x)
        error = std::max(error, std::abs(double(preview.snapshot_alpha[x + width * y]) -
                                            pixels[y + height * x + 3 * count]));
      expect_true(error < 1e-6);
      preview.SetSnapshotFilename(std::string(fixture) + "-preview.png");
      preview.SavePreviewSnapshot();
    }
  }
}

namespace {
class TestAtmosphereLight final : public InfiniteLight, public Atmosphere {
public:
  point3f Radiance(const point3f &, const vec3f &, Float) const override { return point3f(.2); }
  vec3f Sample(const point3f &, vec2f u, Float) const override {
    Float y = 1 - 2 * u[0], r = std::sqrt(std::max(Float(0), 1 - y * y));
    return vec3f(r * std::cos(2 * M_PI * u[1]), y, r * std::sin(2 * M_PI * u[1]));
  }
  Float Pdf(const point3f &, const vec3f &, Float) const override { return 1 / (4 * M_PI); }
  double SamplingWeight() const override { return 1; }
  size_t GetSize() const override { return sizeof(*this); }
  const Atmosphere *GetAtmosphere() const override { return this; }
  AtmosphereSegment Segment(const point3f &, const vec3f &, double) const override { return {}; }
  point3f Transmission(const point3f &, const vec3f &, double) const override { return point3f(.25); }
  point3f CelestialTransmission(const point3f &p, const vec3f &w, InfiniteLightSpectrum) const override {
    return point3f(MaySeeDisk(p, w, 0) ? .5 : 0);
  }
  point3f SkyRadiance(const point3f &, const vec3f &) const override { return point3f(.2); }
  bool MaySeeDisk(const point3f &p, const vec3f &w, double radius) const override {
    return p[1] > 0 || w[1] + radius >= 0;
  }
};
}
context("Atmospheric celestial light mixtures") {
  test_that("sampling and PDFs adapt together when a disk is hidden at one altitude") {
    auto atmosphere = std::make_shared<TestAtmosphereLight>();
    auto texture = std::make_shared<constant_texture>(point3f(100000));
    vec3f center = unit_vector(vec3f(0, -.5, 1));
    auto disk = std::make_shared<DiskInfiniteLight>(texture, 64, 64, center, 20, 0, false,
                                                   InfiniteLightSpectrum::Sun);
    InfiniteLightMixture mixture({atmosphere, disk});
    point3f low(0), high(0, 100, 0);
    double sky_pdf = atmosphere->Pdf(low, center, 0);
    expect_true(std::abs(mixture.Pdf(low, center, 0) - sky_pdf) < 1e-6);
    expect_true(std::abs(mixture.Pdf(high, center, 0) - .5 * (sky_pdf + disk->Pdf(high, center, 0))) < 1e-5);
    expect_true(std::abs(mixture.Radiance(low, center, 0)[0] - .2) < 1e-5);
    expect_true(std::abs(mixture.Radiance(high, center, 0)[0] - 50000.2) < .01);
    int low_hits = 0, high_hits = 0;
    const int n = 128;
    bool valid = true;
    for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x) {
      vec2f u((x + .5f) / n, (y + .5f) / n);
      auto a = mixture.Sample(low, u, 0), b = mixture.Sample(high, u, 0);
      low_hits += disk->Pdf(low, a, 0) > 0;
      high_hits += disk->Pdf(high, b, 0) > 0;
      valid &= mixture.Pdf(low, a, 0) > 0 && mixture.Pdf(high, b, 0) > 0;
    }
    double uniform_hit = (1 - std::cos(10 * M_PI / 180)) / 2;
    expect_true(valid);
    expect_true(std::abs(double(low_hits) / (n * n) - uniform_hit) < .003);
    expect_true(std::abs(double(high_hits) / (n * n) - (.5 + .5 * uniform_hit)) < .003);
  }
}
#endif
