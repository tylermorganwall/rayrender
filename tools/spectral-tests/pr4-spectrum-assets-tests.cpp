#include "src/base/base.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

using namespace rayrender::base;

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
    std::cerr << "PR4 spectrum asset test failed: " << message << std::endl;
    std::exit(1);
  }
}

template <typename F>
void CheckThrows(F&& f, const char* message) {
  try {
    f();
  } catch (const std::exception&) {
    return;
  }
  Check(false, message);
}

void TestConstantSpectrum() {
  ConstantSpectrum constant(static_cast<Float>(2.5));
  SampledWavelengths wavelengths = SampledWavelengths::SampleUniform(0.25f);
  SampledSpectrum sampled = constant.Sample(wavelengths);

  Check(Approx(constant(400), static_cast<Float>(2.5)), "constant spectrum evaluation");
  Check(Approx(constant.MaxValue(), static_cast<Float>(2.5)), "constant spectrum max value");
  for (int i = 0; i < NSpectrumSamples; ++i) {
    Check(Approx(sampled[i], static_cast<Float>(2.5)), "constant spectrum sampled lane");
  }
}

void TestPiecewiseLinearSpectrum() {
  PiecewiseLinearSpectrum spectrum({400, 500, 600}, {1, 3, 7});
  Check(Approx(spectrum(400), 1), "piecewise knot value");
  Check(Approx(spectrum(450), 2), "piecewise interpolated value");
  Check(Approx(spectrum(600), 7), "piecewise last knot value");
  Check(Approx(spectrum(399), 0), "piecewise zero extrapolation below range");
  Check(Approx(spectrum(601), 0), "piecewise zero extrapolation above range");
  Check(Approx(spectrum.MaxValue(), 7), "piecewise max value");

  SpectrumDataPolicy constantPolicy;
  constantPolicy.extrapolation = SpectrumExtrapolationPolicy::Constant;
  PiecewiseLinearSpectrum clamped({400, 500, 600}, {1, 3, 7}, constantPolicy);
  Check(Approx(clamped(300), 1), "piecewise constant extrapolation below range");
  Check(Approx(clamped(700), 7), "piecewise constant extrapolation above range");

  SpectrumDataPolicy errorPolicy;
  errorPolicy.extrapolation = SpectrumExtrapolationPolicy::Error;
  PiecewiseLinearSpectrum errorSpectrum({400, 500, 600}, {1, 3, 7}, errorPolicy);
  CheckThrows([&]() { (void)errorSpectrum(300); }, "piecewise error extrapolation below range");

  CheckThrows(
    []() { PiecewiseLinearSpectrum({500, 400}, {1, 2}); },
    "piecewise rejects unsorted wavelengths"
  );
  CheckThrows(
    []() { PiecewiseLinearSpectrum({400, 400}, {1, 2}); },
    "piecewise rejects duplicate wavelengths"
  );
  CheckThrows(
    []() { PiecewiseLinearSpectrum({400, 500}, {1, std::numeric_limits<Float>::infinity()}); },
    "piecewise rejects nonfinite values"
  );
  CheckThrows(
    []() { PiecewiseLinearSpectrum({400, 500}, {1, -1}); },
    "piecewise rejects negative values by default"
  );

  SpectrumDataPolicy sortedPolicy;
  sortedPolicy.order = SpectrumOrderPolicy::Sort;
  PiecewiseLinearSpectrum sorted({600, 400, 500}, {7, 1, 3}, sortedPolicy);
  Check(Approx(sorted(450), 2), "piecewise sorted input policy");

  SpectrumDataPolicy boundedPolicy;
  boundedPolicy.validation = SpectrumValueValidation::Bounded01;
  CheckThrows(
    [&]() { PiecewiseLinearSpectrum({400, 500}, {0, static_cast<Float>(1.2)}, boundedPolicy); },
    "piecewise bounded validation"
  );

  SpectrumDataPolicy finitePolicy;
  finitePolicy.validation = SpectrumValueValidation::Finite;
  PiecewiseLinearSpectrum signedSpectrum({400, 500}, {-1, 1}, finitePolicy);
  Check(Approx(signedSpectrum(450), 0), "piecewise finite-only validation allows signed setup data");
}

void TestDenselySampledSpectrum() {
  DenselySampledSpectrum dense(400, 402, {10, 20, 30});
  Check(Approx(dense(400), 10), "dense spectrum exact sample");
  Check(Approx(dense(static_cast<Float>(401.49)), 20), "dense spectrum lround sample");
  Check(Approx(dense(static_cast<Float>(401.51)), 30), "dense spectrum rounded-up sample");
  Check(Approx(dense(399), 0), "dense spectrum below range");
  Check(Approx(dense(403), 0), "dense spectrum above range");
  Check(Approx(dense.MaxValue(), 30), "dense spectrum max value");

  DenselySampledSpectrum sampled = DenselySampledSpectrum::SampleFunction(
    [](Float lambda) { return lambda - static_cast<Float>(399); },
    400,
    402
  );
  Check(Approx(sampled(400), 1), "dense sampled function first value");
  Check(Approx(sampled(402), 3), "dense sampled function last value");
}

void TestBlackbodySpectrum() {
  struct Reference {
    Float lambda;
    Float temperature;
    Float radiance;
  };
  const Reference references[] = {
    {483, 6000, static_cast<Float>(3.1849e13)},
    {600, 6000, static_cast<Float>(2.86772e13)},
    {500, 3700, static_cast<Float>(1.59845e12)},
    {600, 4500, static_cast<Float>(7.46497e12)}
  };

  for (const Reference& reference : references) {
    Check(
      ApproxRel(Blackbody(reference.lambda, reference.temperature), reference.radiance, static_cast<Float>(0.001), 0),
      "blackbody Planck reference value"
    );
  }

  for (Float temperature : {2700.f, 3000.f, 4500.f, 5600.f, 6000.f}) {
    Float lambdaMax = static_cast<Float>(2.8977721e-3) / temperature * static_cast<Float>(1e9);
    Check(
      Blackbody(lambdaMax * static_cast<Float>(0.99), temperature) < Blackbody(lambdaMax, temperature),
      "blackbody rises before Wien peak"
    );
    Check(
      Blackbody(lambdaMax * static_cast<Float>(1.01), temperature) < Blackbody(lambdaMax, temperature),
      "blackbody falls after Wien peak"
    );
  }

  BlackbodySpectrum normalized(5000);
  Check(Approx(normalized.MaxValue(), 1), "normalized blackbody max value");
  Check(normalized(450) > 0, "normalized blackbody positive visible value");
}

void CheckNamedRegistry(const NamedSpectrumRegistry& registry) {
  Check(registry.Size() == 91, "named spectrum registry size");
  Check(registry.Get("cie-x") != nullptr, "cie-x lookup");
  Check(registry.Get("cie-y") != nullptr, "cie-y lookup");
  Check(registry.Get("cie-z") != nullptr, "cie-z lookup");
  Check(registry.Get("stdillum-D65") != nullptr, "stdillum-D65 lookup");
  Check(registry.Get("metal-Cu-eta") != nullptr, "metal-Cu-eta lookup");
  Check(registry.Get("metal-Cu-k") != nullptr, "metal-Cu-k lookup");
  Check(registry.Get("glass-BK7") != nullptr, "glass-BK7 lookup");
  Check(registry.Get("canon_eos_5d_r") != nullptr, "camera sensor lookup");
  Check(registry.Get("missing-spectrum") == nullptr, "missing named spectrum returns null");

  const Spectrum& x = registry.GetOrThrow("cie-x");
  const Spectrum& y = registry.GetOrThrow("cie-y");
  const Spectrum& z = registry.GetOrThrow("cie-z");

  Check(Approx(x(360), static_cast<Float>(0.0001299), static_cast<Float>(1e-9)), "cie-x source sample");
  Check(Approx(y(555), 1, static_cast<Float>(1e-7)), "cie-y source sample");
  Check(Approx(z(360), static_cast<Float>(0.0006061), static_cast<Float>(1e-9)), "cie-z source sample");
  Check(
    ApproxRel(IntegrateSpectrum(y), CIEYIntegral, static_cast<Float>(1e-5), static_cast<Float>(1e-4)),
    "cie-y dense integral"
  );

  XYZ constantXYZ = SpectrumToXYZ(Spectrum(ConstantSpectrum(1)), x, y, z);
  Check(Approx(constantXYZ.x, 1, static_cast<Float>(0.006)), "constant spectrum x integral");
  Check(Approx(constantXYZ.y, 1, static_cast<Float>(0.006)), "constant spectrum y integral");
  Check(Approx(constantXYZ.z, 1, static_cast<Float>(0.006)), "constant spectrum z integral");

  const Spectrum& d65 = registry.GetOrThrow("stdillum-D65");
  Check(Approx(d65(560), static_cast<Float>(1.01122406507609), static_cast<Float>(1e-6)), "D65 normalized sample");
  Check(
    ApproxRel(InnerProduct(d65, y), CIEYIntegral, static_cast<Float>(1e-5), static_cast<Float>(1e-4)),
    "D65 luminance normalization"
  );

  Check(
    Approx(registry.GetOrThrow("metal-Cu-eta")(551.040771f), static_cast<Float>(0.950375), static_cast<Float>(1e-6)),
    "metal copper eta source sample"
  );
  Check(
    Approx(registry.GetOrThrow("metal-Cu-k")(551.040771f), static_cast<Float>(2.5765), static_cast<Float>(1e-6)),
    "metal copper k source sample"
  );

  const SpectrumMetadata* metadata = registry.Metadata("stdillum-D65");
  Check(metadata != nullptr, "metadata lookup");
  Check(metadata->semantic == SpectrumSemantic::Illuminant, "metadata semantic");
  Check(metadata->wavelengthUnit == "nm", "metadata wavelength unit");
  Check(metadata->sha256.size() == 64, "metadata per-spectrum checksum");
}

void TestNamedSpectrumAssets(const std::string& assetDirectory) {
  NamedSpectrumRegistry registry = NamedSpectrumRegistry::LoadFromDirectory(assetDirectory);
  CheckNamedRegistry(registry);
}

void TestDefaultLookup() {
  std::string assetDirectory = FindSpectralAssetDirectory();
  Check(!assetDirectory.empty(), "default asset directory lookup");
  NamedSpectrumRegistry registry = NamedSpectrumRegistry::LoadDefault();
  CheckNamedRegistry(registry);
}

} // namespace

int main(int argc, char** argv) {
  TestConstantSpectrum();
  TestPiecewiseLinearSpectrum();
  TestDenselySampledSpectrum();
  TestBlackbodySpectrum();

  if (argc > 1) {
    TestNamedSpectrumAssets(argv[1]);
  }
  TestDefaultLookup();

  std::cout << "PR4 spectrum asset tests passed" << std::endl;
  return 0;
}
