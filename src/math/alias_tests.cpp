#ifdef NOT_CRAN
#include "distributions.h"
#include "sampler.h"
#include <random>
#include <set>
#include <testthat.h>

context("CPU alias distributions") {
  test_that("alias samples match independent cell masses, densities and endpoints") {
    const int width = 7, height = 5;
    std::vector<Float> weights(width * height);
    double sum = 0;
    for (int i = 0; i < width * height; ++i) {
      weights[i] = i / width == 2 || i % 4 == 0 ? 0 : 1 + i * i;
      sum += weights[i];
    }
    AliasDistribution2D alias(weights.data(), width, height);
    Distribution2D cdf(weights.data(), width, height);
    double integral = 0, density_error = 0, sample_error = 0;
    for (int i = 0; i < width * height; ++i) {
      vec2f uv((i % width + .5) / width, (i / width + .5) / height);
      double actual = alias.Pdf(uv), expected = weights[i] / sum * width * height;
      integral += actual / (width * height);
      density_error = std::max(density_error, std::abs(actual - expected));
      density_error = std::max(density_error, std::abs(actual - cdf.Pdf(uv)));
    }
    std::mt19937 rng(571);
    const int samples = 1 << 20;
    std::vector<int> counts(weights.size());
    for (int i = 0; i < samples; ++i) {
      vec2f u(Float(rng() >> 8) * 0x1p-24f, Float(rng() >> 8) * 0x1p-24f);
      Float pdf;
      auto uv = alias.SampleContinuous(u, &pdf);
      int cell = int(double(uv[1]) * height) * width + int(double(uv[0]) * width);
      ++counts.at(cell);
      sample_error = std::max(sample_error, double(std::abs(pdf - alias.Pdf(uv))));
    }
    size_t frequency_errors = 0;
    for (size_t i = 0; i < counts.size(); ++i) {
      double expected = samples * weights[i] / sum;
      frequency_errors += std::abs(counts[i] - expected) > 6 * std::sqrt(expected) + 8;
    }
    for (Float u : {Float(0), Float(.25), Float(.5), std::nextafter(Float(1), Float(0)), Float(1)})
      for (Float v : {Float(0), Float(.5), Float(1)}) {
        Float pdf;
        auto uv = alias.SampleContinuous(vec2f(u, v), &pdf);
        sample_error = std::max(sample_error, double(std::abs(pdf - alias.Pdf(uv))));
      }
    expect_true(std::abs(integral - 1) < 1e-6);
    expect_true(density_error < 2e-6);
    expect_true(sample_error < 2e-6);
    expect_true(frequency_errors == 0);
  }
  test_that("black maps and single cells have a uniform finite fallback") {
    for (int count : {1, 17}) {
      std::vector<Float> black(count, 0);
      AliasDistribution2D alias(black.data(), count, 1);
      Float pdf;
      auto uv = alias.SampleContinuous(vec2f(1, 0), &pdf);
      expect_true(pdf == 1);
      expect_true(alias.Pdf(uv) == 1);
      expect_true((uv[0] >= 0 && uv[0] < 1 && uv[1] == 0));
    }
    Float invalid[] = {-1};
    expect_error(AliasDistribution2D(invalid, 1, 1));
    expect_error(AliasDistribution2D(invalid, 0, 1));
  }
  test_that("continuous blue-noise samples retain strata and subsequent dimensions") {
    random_gen rng;
    SobolBlueNoiseSampler regular(rng), continuous(rng);
    regular.independent_dimensions = continuous.independent_dimensions = true;
    std::set<uint32_t> fractions;
    size_t stratum_errors = 0, sequence_errors = 0;
    for (unsigned pixel = 0; pixel < 512; ++pixel) {
      regular.StartPixel(pixel, 21);
      continuous.StartPixel(pixel, 21);
      for (int sample = 0; sample < 8; ++sample) {
        vec2f a = regular.Get2D(), b = continuous.Get2DContinuous();
        for (int axis = 0; axis < 2; ++axis) {
          stratum_errors += int(a[axis] * 256) != int(b[axis] * 256);
          fractions.insert(uint32_t(b[axis] * 0x1p24f) & 65535u);
        }
        sequence_errors += regular.Get1D() != continuous.Get1D();
        regular.StartNextSample();
        continuous.StartNextSample();
      }
    }
    expect_true(stratum_errors == 0);
    expect_true(sequence_errors == 0);
    expect_true(fractions.size() > 4000);
  }
}
#endif
