#include "src/render/spectral_camera.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

using namespace rayrender::base;
using namespace rayrender::render;

namespace {

constexpr int FilmWidth = 120;
constexpr int FilmHeight = 80;

bool Approx(Float lhs, Float rhs, Float tolerance = static_cast<Float>(1e-5)) {
  return std::fabs(lhs - rhs) <= tolerance;
}

void Check(bool condition, const char* message) {
  if (!condition) {
    std::cerr << "PR21 Camera test failed: " << message << std::endl;
    std::exit(1);
  }
}

void CheckApprox(Float actual, Float expected, Float tolerance, const char* message) {
  if (!Approx(actual, expected, tolerance)) {
    std::cerr << "PR21 Camera test failed: " << message
              << " expected " << expected << " got " << actual << std::endl;
    std::exit(1);
  }
}

void CheckWavelengthsEqual(
  const SampledWavelengths& actual,
  const SampledWavelengths& expected,
  const char* message
) {
  for (int i = 0; i < NSpectrumSamples; ++i) {
    CheckApprox(actual[i], expected[i], static_cast<Float>(1e-7), message);
    CheckApprox(actual.PDF(i), expected.PDF(i), static_cast<Float>(1e-10), message);
  }
  Check(actual.SecondaryTerminated() == expected.SecondaryTerminated(), message);
}

Float DotForward(const CameraRay& ray, const SpectralCamera& camera) {
  return dot(ray.ray.direction(), camera.Forward());
}

CameraSample SampleAt(
  Float sx,
  Float sy,
  point2f lens = point2f(static_cast<Float>(0.5), static_cast<Float>(0.5))
) {
  CameraSample sample;
  sample.pFilm = point2f(sx * FilmWidth, sy * FilmHeight);
  sample.pLens = lens;
  sample.time = static_cast<Float>(0.25);
  return sample;
}

SpectralCameraOptions MakeOptions() {
  SpectralCameraOptions options;
  options.filmWidth = FilmWidth;
  options.filmHeight = FilmHeight;
  options.shutterOpen = static_cast<Float>(0.2);
  options.shutterClose = static_cast<Float>(0.6);
  options.enableDifferentials = true;
  return options;
}

std::vector<RealisticCameraLensElement> SimpleBiconvexLens(EtaSpectrumHandle glass) {
  return {
    RealisticCameraLensElement::Spherical(
      static_cast<Float>(0.06),
      static_cast<Float>(0.006),
      std::move(glass),
      static_cast<Float>(0.018)
    ),
    RealisticCameraLensElement::Spherical(
      static_cast<Float>(-0.06),
      static_cast<Float>(0.04),
      ConstantEtaSpectrum(1),
      static_cast<Float>(0.018)
    )
  };
}

SpectralCamera MakeRealistic(EtaSpectrumHandle glass) {
  RealisticCameraParameters params;
  params.lookfrom = point3f(0, 0, 0);
  params.lookat = point3f(0, 0, -1);
  params.up = vec3f(0, 1, 0);
  params.filmDiagonal = static_cast<Float>(0.035);
  params.lensElements = SimpleBiconvexLens(std::move(glass));
  params.options = MakeOptions();
  return SpectralCamera::Realistic(params);
}

void CheckPositiveFiniteWeight(const CameraRay& ray) {
  for (int i = 0; i < NSpectrumSamples; ++i) {
    Check(ray.weight[i] > 0 && std::isfinite(ray.weight[i]), "realistic camera weight is positive finite");
    CheckApprox(ray.weight[i], ray.weight[0], static_cast<Float>(1e-7), "realistic camera weight is spectral-neutral");
  }
}

void TestNondispersiveRealisticCameraPreservesPacket() {
  SpectralCamera camera = MakeRealistic(ConstantEtaSpectrum(static_cast<Float>(1.5168)));
  SampledWavelengths lambda = SampledWavelengths::SampleUniform(static_cast<Float>(0.31));
  SampledWavelengths before = lambda;

  std::optional<CameraRay> ray = camera.GenerateRay(
    SampleAt(static_cast<Float>(0.5), static_cast<Float>(0.5)),
    lambda
  );
  Check(ray.has_value(), "nondispersive realistic camera generates axial ray");
  Check(DotForward(*ray, camera) > static_cast<Float>(0.99), "axial realistic ray points along camera forward");
  CheckApprox(ray->ray.time(), static_cast<Float>(0.3), static_cast<Float>(1e-6), "realistic camera shutter sampling");
  CheckPositiveFiniteWeight(*ray);
  CheckWavelengthsEqual(lambda, before, "nondispersive realistic camera preserves wavelength packet");
}

void TestDispersiveRealisticCameraTerminatesPacket() {
  EtaSpectrumHandle cauchy = CauchyEtaSpectrum(static_cast<Float>(1.45), static_cast<Float>(0.02));
  Check(EvaluateEtaSpectrum(cauchy, static_cast<Float>(410)) >
          EvaluateEtaSpectrum(cauchy, static_cast<Float>(700)),
        "Cauchy glass eta decreases from blue to red");
  SpectralCamera camera = MakeRealistic(cauchy);
  CameraSample sample = SampleAt(static_cast<Float>(0.72), static_cast<Float>(0.5));

  SampledWavelengths blueLambda = SampledWavelengths::SampleUniform(static_cast<Float>(0.10));
  SampledWavelengths redLambda = SampledWavelengths::SampleUniform(static_cast<Float>(0.75));
  std::optional<CameraRay> blue = camera.GenerateRay(sample, blueLambda);
  std::optional<CameraRay> red = camera.GenerateRay(sample, redLambda);
  Check(blue.has_value(), "dispersive realistic camera generates blue ray");
  Check(red.has_value(), "dispersive realistic camera generates red ray");
  Check(blueLambda.SecondaryTerminated(), "dispersive realistic camera terminates blue packet");
  Check(redLambda.SecondaryTerminated(), "dispersive realistic camera terminates red packet");

  Float blueAngle = std::acos(std::max(static_cast<Float>(-1), std::min(static_cast<Float>(1), DotForward(*blue, camera))));
  Float redAngle = std::acos(std::max(static_cast<Float>(-1), std::min(static_cast<Float>(1), DotForward(*red, camera))));
  Check(std::fabs(blueAngle - redAngle) > static_cast<Float>(1e-5), "dispersive camera changes lens direction by wavelength");
}

void TestLensElementGlassDescriptors(const NamedSpectrumRegistry& registry) {
  RealisticCameraLensElement constant = RealisticCameraLensElement::Spherical(1, 1, static_cast<Float>(1.5), 1);
  RealisticCameraLensElement cauchy = RealisticCameraLensElement::Spherical(
    1,
    1,
    CauchyEtaSpectrum(static_cast<Float>(1.45), static_cast<Float>(0.01)),
    1
  );
  RealisticCameraLensElement sellmeier = RealisticCameraLensElement::Spherical(
    1,
    1,
    SellmeierEtaSpectrum(
      {static_cast<Float>(1.03961212), static_cast<Float>(0.231792344), static_cast<Float>(1.01046945)},
      {static_cast<Float>(0.00600069867), static_cast<Float>(0.0200179144), static_cast<Float>(103.560653)}
    ),
    1
  );
  RealisticCameraLensElement named = RealisticCameraLensElement::Spherical(
    1,
    1,
    NamedEtaSpectrum(registry, "glass-BK7"),
    1
  );
  RealisticCameraLensElement stop = RealisticCameraLensElement::ApertureStop(1, 1);

  Check(!constant.IsDispersive(), "constant glass descriptor is nondispersive");
  Check(cauchy.IsDispersive(), "Cauchy glass descriptor is dispersive");
  Check(sellmeier.IsDispersive(), "Sellmeier glass descriptor is dispersive");
  Check(named.IsDispersive(), "named glass descriptor is dispersive");
  Check(!stop.IsDispersive(), "aperture stop is not dispersive");
  Check(cauchy.Eta(static_cast<Float>(410)) > cauchy.Eta(static_cast<Float>(700)), "Cauchy lens eta is physically ordered");
  Check(named.Eta(static_cast<Float>(410)) > named.Eta(static_cast<Float>(700)), "named BK7 eta is physically ordered");
}

void TestRegionInitializationUsesGeneratedRealisticRay() {
  RealisticCameraParameters params;
  params.lookfrom = point3f(0, 0, 0);
  params.lookat = point3f(0, 0, -1);
  params.filmDiagonal = static_cast<Float>(0.035);
  params.lensElements = SimpleBiconvexLens(ConstantEtaSpectrum(static_cast<Float>(1.5168)));
  params.options = MakeOptions();
  params.options.medium = MediumHandle::FromIndex(2, 3);
  params.options.regionInitialization.mode = CameraRegionInitializationMode::DeferredContainment;
  params.options.regionInitialization.requiresContainmentQuery = true;

  SpectralCamera camera = SpectralCamera::Realistic(params);
  SampledWavelengths lambda = SampledWavelengths::SampleUniform(static_cast<Float>(0.27));
  std::optional<CameraRay> ray = camera.GenerateRay(SampleAt(static_cast<Float>(0.55), static_cast<Float>(0.5)), lambda);
  Check(ray.has_value(), "realistic camera with region initialization generates ray");
  Check(ray->hasInitialMedium, "realistic camera records initial medium on generated ray");
  Check(ray->medium == params.options.medium, "realistic camera carries medium handle");
  Check(ray->regionInitialization.mode == CameraRegionInitializationMode::DeferredContainment, "realistic camera carries region mode");
  Check(ray->regionInitialization.requiresContainmentQuery, "realistic camera carries deferred containment flag");
  Check((ray->ray.origin() - camera.Origin()).length() > static_cast<Float>(0.01),
        "realistic camera initializes state on the final lens-exit ray");
}

void TestMeasuredSensorDeterministicIntegral(
  const NamedSpectrumRegistry& registry,
  const RGBColorSpace& colorSpace
) {
  SensorResponseCurve red = SensorResponseCurve::SampleFunction([](Float) { return static_cast<Float>(1); });
  SensorResponseCurve green = SensorResponseCurve::SampleFunction([](Float) { return static_cast<Float>(2); });
  SensorResponseCurve blue = SensorResponseCurve::SampleFunction([](Float) { return static_cast<Float>(3); });
  PixelSensor sensor = PixelSensor::CreateMeasured(
    red,
    green,
    blue,
    ColorSpaceMatrix3x3::Identity(),
    registry,
    colorSpace
  );
  Spectrum radiance(ConstantSpectrum(static_cast<Float>(2)));
  RGB sensorRGB = sensor.DeterministicSensorRGB(radiance);
  Float samples = LambdaMax - LambdaMin + static_cast<Float>(1);
  CheckApprox(sensorRGB.r, static_cast<Float>(2) * samples, static_cast<Float>(1e-4), "measured sensor red integral");
  CheckApprox(sensorRGB.g, static_cast<Float>(4) * samples, static_cast<Float>(1e-4), "measured sensor green integral");
  CheckApprox(sensorRGB.b, static_cast<Float>(6) * samples, static_cast<Float>(1e-4), "measured sensor blue integral");
  Check(sensor.Mode() == PixelSensorMode::Measured, "measured sensor reports measured mode");
}

} // namespace

int main(int argc, char** argv) {
  try {
    std::string assetDirectory = argc > 1 ? argv[1] : FindSpectralAssetDirectory();
    NamedSpectrumRegistry registry = NamedSpectrumRegistry::LoadFromDirectory(assetDirectory);
    RGBColorSpace colorSpace = LoadSRGBColorSpace(assetDirectory);

    TestNondispersiveRealisticCameraPreservesPacket();
    TestDispersiveRealisticCameraTerminatesPacket();
    TestLensElementGlassDescriptors(registry);
    TestRegionInitializationUsesGeneratedRealisticRay();
    TestMeasuredSensorDeterministicIntegral(registry, colorSpace);
  } catch (const std::exception& error) {
    std::cerr << "PR21 Camera test failed with exception: " << error.what() << std::endl;
    return 1;
  }

  std::cout << "PR21 Camera tests passed" << std::endl;
  return 0;
}
