#include "src/base/base.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <type_traits>
#include <vector>

using namespace rayrender::base;

static_assert(NSpectrumSamples == 4, "pbrt spectral packets use four samples");
static_assert(LambdaMin == 360.f, "pbrt visible wavelength interval starts at 360 nm");
static_assert(LambdaMax == 830.f, "pbrt visible wavelength interval ends at 830 nm");
static_assert(sizeof(SampledSpectrum) == sizeof(Float) * NSpectrumSamples,
              "SampledSpectrum must remain fixed-width inline storage");
static_assert(sizeof(SampledWavelengths) == sizeof(Float) * NSpectrumSamples * 2,
              "SampledWavelengths must remain fixed-width inline storage");

namespace {

bool Approx(Float lhs, Float rhs, Float tolerance = static_cast<Float>(1e-5)) {
  return std::fabs(lhs - rhs) <= tolerance;
}

bool ApproxRel(Float lhs, Float rhs, Float relativeTolerance, Float absoluteTolerance) {
  Float scale = std::max(static_cast<Float>(1), std::max(std::fabs(lhs), std::fabs(rhs)));
  return std::fabs(lhs - rhs) <= std::max(absoluteTolerance, relativeTolerance * scale);
}

void Check(bool condition, const char* message) {
  if (!condition) {
    std::cerr << "PR3 spectrum test failed: " << message << std::endl;
    std::exit(1);
  }
}

Float ReferenceSampleVisibleWavelengths(Float u) {
  return static_cast<Float>(538) -
         static_cast<Float>(138.888889) *
           std::atanh(static_cast<Float>(0.85691062) - static_cast<Float>(1.82750197) * u);
}

Float ReferenceVisibleWavelengthsPDF(Float lambda) {
  if (lambda < static_cast<Float>(360) || lambda > static_cast<Float>(830)) {
    return 0;
  }
  Float c = std::cosh(static_cast<Float>(0.0072) * (lambda - static_cast<Float>(538)));
  return static_cast<Float>(0.0039398042) / (c * c);
}

Float Uniform01(std::uint32_t& state) {
  state = state * 1664525u + 1013904223u;
  return static_cast<Float>((state >> 8) * (1.0 / 16777216.0));
}

void TestSampledSpectrumArithmetic() {
  SampledSpectrum a{1, 2, 3, 4};
  SampledSpectrum b{4, 3, 2, 1};

  SampledSpectrum sum = a + b;
  Check(Approx(sum[0], 5), "SampledSpectrum addition lane 0");
  Check(Approx(sum[3], 5), "SampledSpectrum addition lane 3");

  SampledSpectrum product = a * b;
  Check(Approx(product[0], 4), "SampledSpectrum product lane 0");
  Check(Approx(product[1], 6), "SampledSpectrum product lane 1");
  Check(Approx(product[2], 6), "SampledSpectrum product lane 2");
  Check(Approx(product[3], 4), "SampledSpectrum product lane 3");

  SampledSpectrum scaled = product / 2;
  Check(Approx(scaled[0], 2), "SampledSpectrum scalar division lane 0");
  Check(Approx(scaled[1], 3), "SampledSpectrum scalar division lane 1");

  SampledSpectrum safe = SafeDiv(SampledSpectrum{1, 2, 3, 4}, SampledSpectrum{1, 0, 3, 0});
  Check(Approx(safe[0], 1), "SafeDiv finite lane");
  Check(Approx(safe[1], 0), "SafeDiv zero lane");
  Check(Approx(safe[2], 1), "SafeDiv second finite lane");
  Check(Approx(safe[3], 0), "SafeDiv second zero lane");

  SampledSpectrum clamped = ClampZero(SampledSpectrum{-1, 0, 2, -3});
  Check(Approx(clamped[0], 0), "ClampZero first lane");
  Check(Approx(clamped[2], 2), "ClampZero positive lane");

  SampledSpectrum sqrt = Sqrt(SampledSpectrum{1, 4, 9, 16});
  Check(Approx(sqrt[0], 1), "Sqrt lane 0");
  Check(Approx(sqrt[3], 4), "Sqrt lane 3");

  SampledSpectrum exp = Exp(SampledSpectrum{0, 0, 0, 0});
  Check(Approx(exp[0], 1), "Exp lane 0");
  Check(Approx(exp[3], 1), "Exp lane 3");

  Check(Approx(a.MinComponentValue(), 1), "MinComponentValue");
  Check(Approx(a.MaxComponentValue(), 4), "MaxComponentValue");
  Check(Approx(a.Average(), 2.5f), "Average");
  Check(a.IsPositive(), "IsPositive");
  Check(!a.HasNaNs(), "HasNaNs false");
  Check(!a.IsInf(), "IsInf false");

  SampledSpectrum inf{1, 2, std::numeric_limits<Float>::infinity(), 4};
  Check(inf.IsInf(), "IsInf true");
}

void TestPbrtReferenceValues() {
  struct Reference {
    Float u;
    Float lambda;
    Float pdf;
  };

  const std::array<Reference, 10> references = {{
    {0.000000000f, 360.000002717f, 0.001046822479f},
    {0.001000000f, 360.949692127f, 0.001059148834f},
    {0.010000000f, 369.033757203f, 0.001168901806f},
    {0.100000000f, 424.342922587f, 0.002149193620f},
    {0.250000000f, 479.154062544f, 0.003309324806f},
    {0.500000000f, 545.903013581f, 0.003927075374f},
    {0.750000000f, 616.856224458f, 0.002900074181f},
    {0.900000000f, 686.015891700f, 0.001494392621f},
    {0.990000000f, 795.790474048f, 0.000366770583f},
    {0.999000000f, 825.748958248f, 0.000242284589f},
  }};

  for (const Reference& reference : references) {
    Float lambda = SampleVisibleWavelengths(reference.u);
    Float pdf = VisibleWavelengthsPDF(lambda);
    Check(Approx(lambda, reference.lambda, static_cast<Float>(2e-3)), "SampleVisibleWavelengths reference");
    Check(Approx(pdf, reference.pdf, static_cast<Float>(2e-7)), "VisibleWavelengthsPDF reference");
  }

  for (int i = 0; i <= 1024; ++i) {
    Float u = (static_cast<Float>(i) + static_cast<Float>(0.5)) / static_cast<Float>(1025);
    Float lambda = SampleVisibleWavelengths(u);
    Check(Approx(lambda, ReferenceSampleVisibleWavelengths(u), static_cast<Float>(1e-6)),
          "dense visible inverse comparison");
    Check(Approx(VisibleWavelengthsPDF(lambda), ReferenceVisibleWavelengthsPDF(lambda), static_cast<Float>(1e-8)),
          "dense visible PDF comparison");
  }

  Check(Approx(VisibleWavelengthsPDF(300), 0), "Visible PDF below interval");
  Check(Approx(VisibleWavelengthsPDF(900), 0), "Visible PDF above interval");
}

void TestPacketCorrelations() {
  SampledWavelengths visible = SampledWavelengths::SampleVisible(0.87f);
  const std::array<Float, NSpectrumSamples> expectedU = {0.87f, 0.12f, 0.37f, 0.62f};
  for (int i = 0; i < NSpectrumSamples; ++i) {
    Float expectedLambda = SampleVisibleWavelengths(expectedU[i]);
    Check(Approx(visible[i], expectedLambda, static_cast<Float>(1e-4)), "visible shifted wavelength");
    Check(Approx(visible.PDF(i), VisibleWavelengthsPDF(expectedLambda), static_cast<Float>(1e-8)),
          "visible shifted PDF");
  }
  Check(visible.InvariantsHold(), "visible wavelength invariants");

  SampledWavelengths uniform = SampledWavelengths::SampleUniform(0.9f, 400.f, 700.f);
  const std::array<Float, NSpectrumSamples> expectedUniform = {670.f, 445.f, 520.f, 595.f};
  for (int i = 0; i < NSpectrumSamples; ++i) {
    Check(Approx(uniform[i], expectedUniform[i]), "uniform shifted wavelength");
    Check(Approx(uniform.PDF(i), static_cast<Float>(1.0 / 300.0)), "uniform PDF");
  }
  Check(uniform.InvariantsHold(400.f, 700.f), "uniform wavelength invariants");
}

void TestTermination() {
  SampledWavelengths wavelengths = SampledWavelengths::SampleVisible(0.37f);
  std::array<Float, NSpectrumSamples> originalLambda{};
  std::array<Float, NSpectrumSamples> originalPDF{};
  for (int i = 0; i < NSpectrumSamples; ++i) {
    originalLambda[i] = wavelengths[i];
    originalPDF[i] = wavelengths.PDF(i);
    Check(originalPDF[i] > 0, "pre-termination PDF is positive");
  }

  Check(!wavelengths.SecondaryTerminated(), "secondary wavelengths initially active");
  wavelengths.TerminateSecondary();
  Check(wavelengths.SecondaryTerminated(), "secondary wavelengths terminate");
  Check(Approx(wavelengths.PDF(0), originalPDF[0] / NSpectrumSamples), "primary PDF scaled after termination");

  for (int i = 0; i < NSpectrumSamples; ++i) {
    Check(Approx(wavelengths[i], originalLambda[i]), "termination preserves wavelengths");
    if (i > 0) {
      Check(Approx(wavelengths.PDF(i), 0), "secondary PDF is zero after termination");
    }
  }
  Check(wavelengths.InvariantsHold(), "terminated wavelength invariants");

  wavelengths.TerminateSecondary();
  Check(Approx(wavelengths.PDF(0), originalPDF[0] / NSpectrumSamples), "termination is idempotent");
  for (int i = 1; i < NSpectrumSamples; ++i) {
    Check(Approx(wavelengths.PDF(i), 0), "idempotent secondary PDF remains zero");
  }

  SampledSpectrum pdf = wavelengths.PDF();
  Check(Approx(pdf[0], originalPDF[0] / NSpectrumSamples), "PDF spectrum primary lane");
  Check(Approx(pdf[1], 0), "PDF spectrum secondary lane");
}

void TestVisibleHistogram() {
  constexpr int binCount = 16;
  constexpr int sampleCount = 200000;
  constexpr int integrationStepsPerBin = 128;
  const Float binWidth = (LambdaMax - LambdaMin) / binCount;

  std::array<int, binCount> counts{};
  std::uint32_t state = 0x12345678u;
  for (int i = 0; i < sampleCount; ++i) {
    Float u = Uniform01(state);
    Float lambda = SampleVisibleWavelengths(u);
    int bin = static_cast<int>((lambda - LambdaMin) / binWidth);
    if (bin < 0) {
      bin = 0;
    }
    if (bin >= binCount) {
      bin = binCount - 1;
    }
    ++counts[bin];
  }

  for (int bin = 0; bin < binCount; ++bin) {
    Float start = LambdaMin + bin * binWidth;
    Float expected = 0;
    for (int step = 0; step < integrationStepsPerBin; ++step) {
      Float u = (static_cast<Float>(step) + static_cast<Float>(0.5)) / integrationStepsPerBin;
      Float lambda = start + u * binWidth;
      expected += VisibleWavelengthsPDF(lambda) * binWidth / integrationStepsPerBin;
    }
    Float observed = static_cast<Float>(counts[bin]) / sampleCount;
    Check(ApproxRel(observed, expected, static_cast<Float>(0.09), static_cast<Float>(0.004)),
          "visible wavelength histogram agrees with PDF");
  }
}

} // namespace

int main() {
  TestSampledSpectrumArithmetic();
  TestPbrtReferenceValues();
  TestPacketCorrelations();
  TestTermination();
  TestVisibleHistogram();
  std::cout << "PR3 sampled spectrum tests passed" << std::endl;
  return 0;
}
