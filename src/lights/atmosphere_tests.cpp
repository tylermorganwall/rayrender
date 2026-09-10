#ifdef NOT_CRAN
#include "atmosphere.h"
#include "../core/PreviewDisplay.h"
#include <testthat.h>
#include <cstdlib>

namespace {
// Deliberately non-multiplicative fit, like independently fitted finite-distance
// tables. Splitting it must not change the total because null rates are arbitrary.
class FittedAtmosphere final : public Atmosphere {
public:
  mutable int calls = 0;
  AtmosphereSegment Segment(const point3f &, const vec3f &, double d) const override {
    ++calls;
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
  test_that("deferred constant weights query only the endpoint and retain the fit anchor") {
    FittedAtmosphere model;
    AtmosphereRay ray(&model);
    ray.Start(point3f(0), vec3f(0, 0, 1), true);
    std::array<double, 3> weight{2, 3, 4};
    for (int i = 1; i <= 100; ++i) {
      ray.Start(point3f(0, 0, i - 1), vec3f(0, 0, 1), true);
      ray.Accumulate(point3f(0, 0, i), weight, .5);
    }
    expect_true(model.calls == 0);
    auto first = ray.Flush();
    expect_true(model.calls == 1);
    auto reference = model.Segment(point3f(0), vec3f(0, 0, 1), 100);
    for (int c = 0; c < 3; ++c) {
      expect_true(std::abs(first.radiance[c] - weight[c] * reference.radiance[c]) < 1e-6);
      expect_true(first.transmission[c] == reference.transmission[c]);
    }
    expect_true(ray.Flush().transmission[0] == 1);
    // A later flush must continue the same cumulative fit, rather than fitting
    // a new local interval and accumulating its zero-distance fitting error.
    for (int c = 0; c < 3; ++c) weight[c] *= first.transmission[c];
    ray.Accumulate(point3f(0, 0, 150), weight, .5);
    auto last = ray.Flush();
    reference = model.Segment(point3f(0), vec3f(0, 0, 1), 150);
    for (int c = 0; c < 3; ++c) {
      expect_true(std::abs(first.radiance[c] + last.radiance[c] - (c + 2) * reference.radiance[c]) < 1e-5);
      expect_true(std::abs(first.transmission[c] * last.transmission[c] - reference.transmission[c]) < 1e-6);
    }
  }

  test_that("a sampled weight correction has the exact weighted integral in expectation") {
    FittedAtmosphere model;
    const double d[] = {10, 20, 45, 60, 100};
    // Include channel-dependent increases and decreases, like colored null
    // weights and compensated roulette. The correction may have either sign.
    const std::array<double, 3> weights[] = {{1, 2, 3}, {.2, 1, 4},
      {.8, .5, 2}, {.3, 2, 1}, {.1, .1, .2}};
    std::array<double, 3> expected{}, mean{}, second{};
    point3f previous(0);
    for (int i = 0; i < 5; ++i) {
      auto value = model.Segment(point3f(0), vec3f(0, 0, 1), d[i]);
      for (int c = 0; c < 3; ++c)
        expected[c] += weights[i][c] * (double(value.radiance[c]) - previous[c]);
      previous = value.radiance;
    }
    random_gen rng(7301);
    const int trials = 100000;
    model.calls = 0;
    for (int trial = 0; trial < trials; ++trial) {
      AtmosphereRay ray(&model);
      ray.Start(point3f(0), vec3f(0, 0, 1), true);
      for (int i = 0; i < 5; ++i)
        ray.Accumulate(point3f(0, 0, d[i]), weights[i], rng.unif_rand());
      auto result = ray.Flush();
      for (int c = 0; c < 3; ++c) {
        mean[c] += result.radiance[c] / trials;
        second[c] += result.radiance[c] * result.radiance[c] / trials;
      }
    }
    expect_true(model.calls == 2 * trials);
    for (int c = 0; c < 3; ++c) {
      double se = std::sqrt(std::max(0.0, second[c] - mean[c] * mean[c]) / trials);
      expect_true(std::abs(mean[c] - expected[c]) < 6 * se + 1e-6);
    }
  }

  test_that("excluded volume intervals restart haze at their exit without bridging the gap") {
    FittedAtmosphere model;
    AtmosphereRay ray(&model);
    vec3f direction(0, 0, 1);
    ray.Start(point3f(0), direction, true);
    auto first = ray.Advance(point3f(0, 0, 10));
    ray.Start(point3f(0, 0, 10), direction, false);
    auto excluded = ray.Advance(point3f(0, 0, 25));
    // Entering/exiting a nested volume keeps haze disabled until the outer exit.
    ray.Start(point3f(0, 0, 25), direction, false);
    auto nested = ray.Advance(point3f(0, 0, 40));
    expect_true(excluded.transmission[0] == 1);
    expect_true(excluded.radiance[0] == 0);
    expect_true(nested.transmission[0] == 1);
    expect_true(nested.radiance[0] == 0);
    ray.Start(point3f(0, 0, 40), direction, true);
    auto last = ray.Advance(point3f(0, 0, 60));
    auto reference = model.Segment(point3f(0, 0, 40), direction, 20);
    expect_true((last.transmission - reference.transmission).length() < 1e-6);
    expect_true((last.radiance - reference.radiance).length() < 1e-6);
    expect_true((ray.Origin(point3f(0)) - point3f(0, 0, 40)).length() == 0);
    point3f environment(3, 5, 7);
    point3f result = first.radiance + first.transmission *
        (last.radiance + last.transmission * ray.Remaining(environment));
    expect_true((result - (first.radiance + first.transmission * environment)).length() < 1e-5);
  }

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
  explicit TestAtmosphereLight(bool finite = true) : finite(finite) {}
  point3f Radiance(const point3f &, const vec3f &, Float) const override { return point3f(.2); }
  vec3f Sample(const point3f &, vec2f u, Float) const override {
    Float y = 1 - 2 * u[0], r = std::sqrt(std::max(Float(0), 1 - y * y));
    return vec3f(r * std::cos(2 * M_PI * u[1]), y, r * std::sin(2 * M_PI * u[1]));
  }
  Float Pdf(const point3f &, const vec3f &, Float) const override { return 1 / (4 * M_PI); }
  double SamplingWeight() const override { return 1; }
  size_t GetSize() const override { return sizeof(*this); }
  const Atmosphere *GetAtmosphere() const override { return this; }
  const Atmosphere *GetTransportAtmosphere() const override { return finite ? this : nullptr; }
  AtmosphereSegment Segment(const point3f &, const vec3f &, double) const override { return {}; }
  point3f Transmission(const point3f &, const vec3f &, double) const override { return point3f(.25); }
  point3f CelestialTransmission(const point3f &p, const vec3f &w, InfiniteLightSpectrum) const override {
    return point3f(MaySeeDisk(p, w, 0) ? .5 : 0);
  }
  point3f SkyRadiance(const point3f &, const vec3f &) const override { return point3f(.2); }
  bool MaySeeDisk(const point3f &p, const vec3f &w, double radius) const override {
    return p[1] > 0 || w[1] + radius >= 0;
  }
private:
  bool finite;
};
}
context("Atmospheric celestial light mixtures") {
  test_that("disabling finite haze preserves celestial filtering and altitude visibility") {
    auto atmosphere = std::make_shared<TestAtmosphereLight>(false);
    auto texture = std::make_shared<constant_texture>(point3f(100000));
    vec3f direction = unit_vector(vec3f(0, -.5, 1));
    auto disk = std::make_shared<DiskInfiniteLight>(texture, 64, 64, direction, 20, 0, false,
                                                   InfiniteLightSpectrum::Sun);
    InfiniteLightMixture mixture({atmosphere, disk});
    expect_true(mixture.GetAtmosphere() == atmosphere.get());
    expect_true(mixture.GetTransportAtmosphere() == nullptr);
    expect_true(std::abs(mixture.Radiance(point3f(0), direction, 0)[0] - .2) < 1e-5);
    expect_true(std::abs(mixture.Radiance(point3f(0, 100, 0), direction, 0)[0] - 50000.2) < .01);
    expect_true(std::abs(mixture.Pdf(point3f(0), direction, 0) - 1 / (4 * M_PI)) < 1e-6);
    InfiniteLightMixture full({std::make_shared<TestAtmosphereLight>()});
    expect_true(full.GetTransportAtmosphere() == full.GetAtmosphere());
  }

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

#ifdef NOT_CRAN
namespace {
Rcpp::List PragueTestDescription(const char *filename) {
  Rcpp::List description;
  description["filename"] = std::string(filename);
  description["meters_per_unit"] = 1.; description["altitude"] = 0.;
  description["visibility"] = 50.; description["albedo"] = .5;
  description["elevation"] = -1.; description["azimuth"] = 250.;
  description["intensity"] = 1.; description["rotation"] = 0.;
  description["origin"] = Rcpp::NumericVector::create(0., 0., 0.);
  description["rgb_gain"] = Rcpp::NumericVector::create(1., 1., 1.);
  description["angular_diameter"] = .533;
  description["render_mode"] = "all"; description["deferred_haze"] = true;
  return description;
}
}

context("Sampled horizon correction") {
  test_that("weighted correction outcomes recover the original clamped source") {
    const char *filename = std::getenv("RAYRENDER_PRAGUE_TEST_FILE");
    if (!filename) return;
    Rcpp::List description = PragueTestDescription(filename);
    for (double probability : {.25, .5, 1.}) {
      description["haze_correction_probability"] = probability;
      PragueInfiniteLight light(description, false);
      expect_true(light.SampleHazeCorrection() == (probability < 1));
      for (double height : {0., 100., 1000., 5000., 15000.})
        for (double slope : {-.1, -.01, 0., .01, .1, .5})
          for (double distance : {.01, 30., 100., 1000., 10000.}) {
            point3f p(0, height, 0);
            vec3f w = unit_vector(vec3f(1, slope, 0));
            AtmosphereSegmentCache exact_cache, low_cache, high_cache;
            auto exact = light.Segment(p, w, distance, exact_cache);
            auto low = light.SampledSegment(p, w, distance, low_cache, probability / 2);
            auto high = light.SampledSegment(p, w, distance, high_cache, (1 + probability) / 2);
            for (int c = 0; c < 3; ++c) {
              double average = probability * double(low.radiance[c]) +
                               (1 - probability) * double(high.radiance[c]);
              expect_true(std::abs(average - exact.radiance[c]) <
                          2e-5 * std::max(1., std::abs(double(exact.radiance[c]))));
              expect_true(low.transmission[c] == exact.transmission[c]);
              expect_true(high.transmission[c] == exact.transmission[c]);
            }
          }
    }
  }

  test_that("changing the table budget reloads the cached model and preserves queries") {
    const char *filename = std::getenv("RAYRENDER_PRAGUE_TEST_FILE");
    if (!filename) return;
    auto description = PragueTestDescription(filename);
    description["transmission_table_max_mb"] = 0.;
    size_t baseline_size;
    AtmosphereSegment expected;
    {
      PragueInfiniteLight light(description, false);
      baseline_size = light.GetSize();
      expected = light.Segment(point3f(0, 5000, 0), unit_vector(vec3f(1, .01, 0)), 10000);
    }
    size_t expanded_size = 0;
    for (double cap : {512., 0., 256., 1., std::numeric_limits<double>::infinity(), 0.}) {
      description["transmission_table_max_mb"] = cap;
      PragueInfiniteLight light(description, false);
      if (cap == 512) {
        expanded_size = light.GetSize();
        expect_true(expanded_size > baseline_size);
      }
      expect_true(light.GetSize() == (cap <= 1 ? baseline_size : expanded_size));
      auto actual = light.Segment(point3f(0, 5000, 0), unit_vector(vec3f(1, .01, 0)), 10000);
      for (int c = 0; c < 3; ++c) {
        expect_true(actual.transmission[c] == expected.transmission[c]);
        expect_true(actual.radiance[c] == expected.radiance[c]);
      }
    }
  }
}
#endif
