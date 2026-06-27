#include "src/render/spectral_camera.h"

#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>

using namespace rayrender::base;
using namespace rayrender::render;

namespace {

constexpr int FilmWidth = 100;
constexpr int FilmHeight = 50;

bool Approx(Float lhs, Float rhs, Float tolerance = static_cast<Float>(1e-5)) {
  return std::fabs(lhs - rhs) <= tolerance;
}

void Check(bool condition, const char* message) {
  if (!condition) {
    std::cerr << "PR9 Camera test failed: " << message << std::endl;
    std::exit(1);
  }
}

void CheckApprox(Float actual, Float expected, Float tolerance, const char* message) {
  if (!Approx(actual, expected, tolerance)) {
    std::cerr << "PR9 Camera test failed: " << message
              << " expected " << expected << " got " << actual << std::endl;
    std::exit(1);
  }
}

void CheckPointApprox(const point3f& actual, const point3f& expected, Float tolerance, const char* message) {
  CheckApprox(actual.xyz.x, expected.xyz.x, tolerance, message);
  CheckApprox(actual.xyz.y, expected.xyz.y, tolerance, message);
  CheckApprox(actual.xyz.z, expected.xyz.z, tolerance, message);
}

void CheckVectorApprox(const vec3f& actual, const vec3f& expected, Float tolerance, const char* message) {
  CheckApprox(actual.xyz.x, expected.xyz.x, tolerance, message);
  CheckApprox(actual.xyz.y, expected.xyz.y, tolerance, message);
  CheckApprox(actual.xyz.z, expected.xyz.z, tolerance, message);
}

vec3f Normalize(const vec3f& value) {
  return value / value.length();
}

CameraSample SampleAt(Float sx, Float sy, point2f lens = point2f(static_cast<Float>(0.5), static_cast<Float>(0.5))) {
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
  options.shutterOpen = static_cast<Float>(0.1);
  options.shutterClose = static_cast<Float>(0.5);
  options.enableDifferentials = true;
  return options;
}

SpectralCamera MakePerspective(Float aperture = 0, MediumHandle medium = MediumHandle::Invalid()) {
  PerspectiveCameraParameters params;
  params.lookfrom = point3f(0, 0, 0);
  params.lookat = point3f(0, 0, -1);
  params.up = vec3f(0, 1, 0);
  params.vfov = static_cast<Float>(90);
  params.aspect = static_cast<Float>(2);
  params.aperture = aperture;
  params.focusDistance = static_cast<Float>(1);
  params.options = MakeOptions();
  params.options.medium = medium;
  return SpectralCamera::Perspective(params);
}

SpectralCamera MakeOrthographic() {
  OrthographicCameraParameters params;
  params.lookfrom = point3f(0, 0, 0);
  params.lookat = point3f(0, 0, -1);
  params.up = vec3f(0, 1, 0);
  params.width = static_cast<Float>(4);
  params.height = static_cast<Float>(2);
  params.options = MakeOptions();
  return SpectralCamera::Orthographic(params);
}

void CheckNeutralWeight(const CameraRay& ray, const char* message) {
  for (int i = 0; i < NSpectrumSamples; ++i) {
    CheckApprox(ray.weight[i], 1, static_cast<Float>(1e-7), message);
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

void TestPerspectiveGeometryAndWeight() {
  SpectralCamera camera = MakePerspective();
  SampledWavelengths lambda = SampledWavelengths::SampleUniform(static_cast<Float>(0.25));
  SampledWavelengths before = lambda;

  std::optional<CameraRay> center = camera.GenerateRay(SampleAt(static_cast<Float>(0.5), static_cast<Float>(0.5)), lambda);
  Check(center.has_value(), "perspective center ray generated");
  CheckPointApprox(center->ray.origin(), point3f(0, 0, 0), static_cast<Float>(1e-6), "perspective center origin");
  CheckVectorApprox(center->ray.direction(), vec3f(0, 0, -1), static_cast<Float>(1e-6), "perspective center direction");
  CheckApprox(center->ray.time(), static_cast<Float>(0.2), static_cast<Float>(1e-6), "camera shutter sampling");
  CheckNeutralWeight(*center, "perspective camera weight neutral");

  std::optional<CameraRay> corner = camera.GenerateRay(SampleAt(1, 1), lambda);
  Check(corner.has_value(), "perspective corner ray generated");
  CheckVectorApprox(
    corner->ray.direction(),
    Normalize(vec3f(-2, 1, -1)),
    static_cast<Float>(1e-6),
    "perspective ray matches legacy screen geometry with normalized pbrt direction"
  );
  CheckWavelengthsEqual(lambda, before, "perspective camera leaves wavelength packet unchanged");
}

void TestOrthographicGeometry() {
  SpectralCamera camera = MakeOrthographic();
  SampledWavelengths lambda = SampledWavelengths::SampleVisible(static_cast<Float>(0.37));
  SampledWavelengths before = lambda;

  std::optional<CameraRay> ray = camera.GenerateRay(SampleAt(static_cast<Float>(0.75), static_cast<Float>(0.25)), lambda);
  Check(ray.has_value(), "orthographic ray generated");
  CheckPointApprox(ray->ray.origin(), point3f(-1, static_cast<Float>(-0.5), 0), static_cast<Float>(1e-6), "orthographic origin");
  CheckVectorApprox(ray->ray.direction(), vec3f(0, 0, -1), static_cast<Float>(1e-6), "orthographic direction");
  CheckNeutralWeight(*ray, "orthographic camera weight neutral");
  CheckWavelengthsEqual(lambda, before, "orthographic camera leaves wavelength packet unchanged");
}

void TestThinLensGeometryAndWavelengths() {
  SpectralCamera camera = MakePerspective(static_cast<Float>(0.4));
  CameraSample sample = SampleAt(
    static_cast<Float>(0.25),
    static_cast<Float>(0.75),
    point2f(static_cast<Float>(0.75), static_cast<Float>(0.5))
  );
  SampledWavelengths lambda = SampledWavelengths::SampleUniform(static_cast<Float>(0.6));
  SampledWavelengths before = lambda;

  std::optional<CameraRay> ray = camera.GenerateRay(sample, lambda);
  Check(ray.has_value(), "thin-lens perspective ray generated");

  point3f expectedOrigin(static_cast<Float>(-0.1), 0, 0);
  point3f focusPoint(1, static_cast<Float>(0.5), -1);
  CheckPointApprox(ray->ray.origin(), expectedOrigin, static_cast<Float>(2e-6), "thin-lens origin uses concentric disk sample");
  CheckVectorApprox(
    ray->ray.direction(),
    Normalize(focusPoint - expectedOrigin),
    static_cast<Float>(2e-6),
    "thin-lens direction targets focus plane"
  );
  CheckNeutralWeight(*ray, "thin-lens camera weight neutral");
  CheckWavelengthsEqual(lambda, before, "thin-lens camera leaves wavelength packet unchanged");
}

void TestInitialMediumAndRegionStubs() {
  MediumHandle medium = MediumHandle::FromIndex(3, 7);
  PerspectiveCameraParameters params;
  params.lookfrom = point3f(0, 0, 0);
  params.lookat = point3f(0, 0, -1);
  params.vfov = static_cast<Float>(45);
  params.aspect = 1;
  params.focusDistance = 1;
  params.options = MakeOptions();
  params.options.medium = medium;
  params.options.regionInitialization.mode = CameraRegionInitializationMode::DeferredContainment;
  params.options.regionInitialization.requiresContainmentQuery = true;

  SpectralCamera camera = SpectralCamera::Perspective(params);
  SampledWavelengths lambda = SampledWavelengths::SampleVisible(static_cast<Float>(0.41));
  std::optional<CameraRay> ray = camera.GenerateRay(SampleAt(static_cast<Float>(0.5), static_cast<Float>(0.5)), lambda);
  Check(ray.has_value(), "medium camera ray generated");
  Check(ray->hasInitialMedium, "camera ray records initial medium");
  Check(ray->medium == medium, "camera ray carries medium handle");
  Check(ray->regionInitialization.mode == CameraRegionInitializationMode::DeferredContainment, "camera ray records region mode");
  Check(ray->regionInitialization.requiresContainmentQuery, "camera ray records deferred containment query");
}

void TestDifferentials() {
  SpectralCamera camera = MakePerspective();
  SampledWavelengths lambda = SampledWavelengths::SampleUniform(static_cast<Float>(0.17));
  SampledWavelengths before = lambda;
  std::optional<CameraRay> ray = camera.GenerateRayDifferential(
    SampleAt(static_cast<Float>(0.5), static_cast<Float>(0.5)),
    lambda
  );

  Check(ray.has_value(), "differential ray generated");
  Check(ray->hasDifferentials, "camera ray has differentials");
  CheckPointApprox(ray->rx.origin(), ray->ray.origin(), static_cast<Float>(1e-6), "perspective rx origin");
  CheckPointApprox(ray->ry.origin(), ray->ray.origin(), static_cast<Float>(1e-6), "perspective ry origin");
  Check(ray->rx.direction().xyz.x < ray->ray.direction().xyz.x, "rx direction shifts in image x under legacy frame");
  Check(ray->ry.direction().xyz.y > ray->ray.direction().xyz.y, "ry direction shifts in image y");
  CheckWavelengthsEqual(lambda, before, "differential generation leaves wavelength packet unchanged");
}

void TestCameraTable() {
  SpectralCameraTable table;
  CameraHandle perspective = table.Add(MakePerspective());
  CameraHandle orthographic = table.Add(MakeOrthographic());
  Check(table.Size() == 2, "camera table size");
  Check(perspective != orthographic, "camera handles differ");

  SampledWavelengths lambda = SampledWavelengths::SampleUniform(static_cast<Float>(0.2));
  std::optional<CameraRay> ray = table.GenerateRay(orthographic, SampleAt(static_cast<Float>(0.5), static_cast<Float>(0.5)), lambda);
  Check(ray.has_value(), "camera table dispatch generates ray");
  CheckVectorApprox(ray->ray.direction(), vec3f(0, 0, -1), static_cast<Float>(1e-6), "camera table dispatch direction");

  bool threw = false;
  try {
    (void)table.Get(CameraHandle::Invalid());
  } catch (const std::out_of_range&) {
    threw = true;
  }
  Check(threw, "invalid camera handle is rejected");
}

PixelSensor MakeSensor(const NamedSpectrumRegistry& registry) {
  return PixelSensor::CreateCIE1931(registry, RGBColorSpace::SRGB());
}

void TestFilmCameraOrdering(const std::string& assetDir) {
  NamedSpectrumRegistry registry = NamedSpectrumRegistry::LoadFromDirectory(assetDir);
  PixelSensor sensor = MakeSensor(registry);

  FilmOptions options;
  options.width = FilmWidth;
  options.height = FilmHeight;
  options.sensor = sensor;
  options.outputColorSpace = RGBColorSpace::SRGB();
  options.wavelengthSampling = WavelengthSamplingMode::Uniform;
  Film film(options);

  SpectralCamera camera = MakePerspective();
  CameraSample sample = SampleAt(static_cast<Float>(0.4), static_cast<Float>(0.6));
  Float wavelengthU = static_cast<Float>(0.31);
  CameraFilmSample generated = GenerateCameraRayFromFilm(film, camera, sample, wavelengthU, true);

  SampledWavelengths expected = film.SampleWavelengths(wavelengthU);
  CheckWavelengthsEqual(generated.wavelengths, expected, "Film wavelength sample precedes camera ray generation deterministically");
  Check(generated.cameraRay.has_value(), "Film-camera bridge generated camera ray");
  Check(generated.cameraRay->hasDifferentials, "Film-camera bridge can request differentials");

  SampledWavelengths directLambda = expected;
  std::optional<CameraRay> direct = camera.GenerateRayDifferential(sample, directLambda);
  Check(direct.has_value(), "direct differential ray generated");
  CheckPointApprox(generated.cameraRay->ray.origin(), direct->ray.origin(), static_cast<Float>(1e-6), "Film bridge ray origin deterministic");
  CheckVectorApprox(generated.cameraRay->ray.direction(), direct->ray.direction(), static_cast<Float>(1e-6), "Film bridge ray direction deterministic");
  CheckWavelengthsEqual(generated.wavelengths, directLambda, "camera did not mutate Film wavelength packet");
}

} // namespace

int main(int argc, char** argv) {
  try {
    std::string assetDir;
    if (argc > 1) {
      assetDir = argv[1];
    } else if (const char* env = std::getenv("RAYRENDER_SPECTRAL_ASSET_DIR")) {
      assetDir = env;
    } else {
      std::cerr << "PR9 Camera test failed: missing spectral asset directory" << std::endl;
      return 1;
    }

    TestPerspectiveGeometryAndWeight();
    TestOrthographicGeometry();
    TestThinLensGeometryAndWavelengths();
    TestInitialMediumAndRegionStubs();
    TestDifferentials();
    TestCameraTable();
    TestFilmCameraOrdering(assetDir);
  } catch (const std::exception& error) {
    std::cerr << "PR9 Camera test failed with exception: " << error.what() << std::endl;
    return 1;
  }

  std::cout << "PR9 Camera tests passed" << std::endl;
  return 0;
}
