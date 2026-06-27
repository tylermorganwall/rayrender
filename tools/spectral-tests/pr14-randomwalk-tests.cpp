#include "src/render/spectral_integrator.h"

#include "src/base/rgb_spectrum.h"

#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <string>
#include <utility>

using namespace rayrender::base;
using namespace rayrender::materials;
using namespace rayrender::render;

namespace {

bool Approx(Float lhs, Float rhs, Float tolerance = static_cast<Float>(1e-5)) {
  return std::fabs(lhs - rhs) <= tolerance;
}

void Check(bool condition, const char* message) {
  if (!condition) {
    std::cerr << "PR14 RandomWalk test failed: " << message << std::endl;
    std::exit(1);
  }
}

void CheckApprox(Float actual, Float expected, Float tolerance, const char* message) {
  if (!Approx(actual, expected, tolerance)) {
    std::cerr << "PR14 RandomWalk test failed: " << message
              << " expected " << expected << " got " << actual << std::endl;
    std::exit(1);
  }
}

void CheckSpectrumApprox(
  const SampledSpectrum& actual,
  const SampledSpectrum& expected,
  Float tolerance,
  const char* message
) {
  for (int i = 0; i < NSpectrumSamples; ++i) {
    CheckApprox(actual[i], expected[i], tolerance, message);
  }
}

void CheckRGBApprox(const RGB& actual, const RGB& expected, Float tolerance, const char* message) {
  if (!Approx(actual.r, expected.r, tolerance) ||
      !Approx(actual.g, expected.g, tolerance) ||
      !Approx(actual.b, expected.b, tolerance)) {
    std::cerr << "PR14 RandomWalk test failed: " << message
              << " expected [" << expected.r << ", " << expected.g << ", " << expected.b << "]"
              << " got [" << actual.r << ", " << actual.g << ", " << actual.b << "]"
              << std::endl;
    std::exit(1);
  }
}

Spectrum GaussianSpectrum(Float center, Float width, Float scale = 1) {
  return Spectrum(DenselySampledSpectrum::SampleFunction([=](Float lambda) {
    Float x = (lambda - center) / width;
    return scale * std::exp(static_cast<Float>(-0.5) * x * x);
  }));
}

SensorResponseCurve GaussianResponse(Float center, Float width) {
  return SensorResponseCurve::SampleFunction([=](Float lambda) {
    Float x = (lambda - center) / width;
    return std::exp(static_cast<Float>(-0.5) * x * x);
  });
}

RGBColorSpace IdentityColorSpace() {
  RGBColorSpace colorSpace;
  colorSpace.name = "identity";
  colorSpace.encoding = RGBColorEncoding::Linear();
  colorSpace.rgbToXYZ = ColorSpaceMatrix3x3::Identity();
  colorSpace.xyzToRGB = ColorSpaceMatrix3x3::Identity();
  return colorSpace;
}

PixelSensor MakeBandSensor() {
  NamedSpectrumRegistry registry;
  return PixelSensor::CreateMeasured(
    GaussianResponse(650, 18),
    GaussianResponse(540, 18),
    GaussianResponse(450, 18),
    ColorSpaceMatrix3x3::Identity(),
    registry,
    IdentityColorSpace()
  );
}

PixelSensor MakeFlatSensor() {
  SensorResponseCurve flat = SensorResponseCurve::SampleFunction([](Float) {
    return static_cast<Float>(1);
  });
  NamedSpectrumRegistry registry;
  return PixelSensor::CreateMeasured(
    flat,
    flat,
    flat,
    ColorSpaceMatrix3x3::Identity(),
    registry,
    IdentityColorSpace()
  );
}

Film MakeFilm(const PixelSensor& sensor, int samples = 1) {
  (void)samples;
  FilmOptions options;
  options.width = 1;
  options.height = 1;
  options.sensor = sensor;
  options.outputColorSpace = IdentityColorSpace();
  options.filter = FilmFilter::Triangle(1, 1);
  options.wavelengthSampling = WavelengthSamplingMode::Uniform;
  options.deterministicSingleThread = true;
  return Film(options);
}

SpectralCamera MakeCenterCamera() {
  PerspectiveCameraParameters params;
  params.lookfrom = point3f(0, 0, 0);
  params.lookat = point3f(0, 0, -1);
  params.up = vec3f(0, 1, 0);
  params.vfov = static_cast<Float>(30);
  params.aspect = 1;
  params.aperture = 0;
  params.focusDistance = 1;
  params.options.filmWidth = 1;
  params.options.filmHeight = 1;
  params.options.enableDifferentials = false;
  return SpectralCamera::Perspective(params);
}

struct SpectralTestScene {
  Scene scene;
  SpectralMaterialTable materials;
  SpectralLightTable lights;
};

void AddDiffuseSphere(
  SpectralTestScene& testScene,
  SpectrumTexture reflectance,
  point3f center = point3f(0, 0, -3)
) {
  MaterialHandle material = testScene.materials.Add(Material::Diffuse(std::move(reflectance)));
  ShapeHandle shape = testScene.scene.AddShape(Shape::Sphere(1, center));
  PrimitiveBinding binding;
  binding.material = material;
  testScene.scene.AddPrimitive(shape, binding);
}

void AddUniformInfinite(SpectralTestScene& testScene, LightSpectrum radiance) {
  RegisterLight(testScene.scene, testScene.lights, Light::UniformInfinite(std::move(radiance)));
  Bounds3f bounds;
  testScene.scene.Bounds(&bounds);
  testScene.lights.Preprocess(bounds);
}

RandomWalkRenderOptions OneBounceOptions() {
  RandomWalkRenderOptions options;
  options.maxDepth = 1;
  options.pixelSamples = 1;
  options.seed = 17;
  options.jitterCameraSamples = false;
  return options;
}

void TestDiffuseInfiniteLiClosedForm() {
  SpectralTestScene testScene;
  AddDiffuseSphere(testScene, SpectrumTexture::Constant(static_cast<Float>(0.5)));
  AddUniformInfinite(testScene, LightSpectrum::Constant(static_cast<Float>(2)));

  RandomWalkRenderOptions options = OneBounceOptions();
  RandomWalkIntegrator integrator(testScene.scene, testScene.materials, testScene.lights, options);
  SpectralRandomSampler sampler(options.seed);
  sampler.StartPixelSample({0, 0}, 0);
  ScratchBuffer scratch(4096);
  SampledWavelengths lambda = SampledWavelengths::SampleUniform(static_cast<Float>(0.25));
  SpectralVisibleSurface visibleSurface;
  RandomWalkRenderStats stats;

  SampledSpectrum L = integrator.Li(
    Ray(point3f(0, 0, 0), vec3f(0, 0, -1)),
    lambda,
    sampler,
    scratch,
    &visibleSurface,
    &stats
  );

  CheckSpectrumApprox(L, SampledSpectrum(1), static_cast<Float>(2e-5), "diffuse infinite one-bounce result matches pbrt random-walk algebra");
  Check(visibleSurface.valid, "first visible surface is recorded");
  CheckApprox(visibleSurface.depth, static_cast<Float>(2), static_cast<Float>(1e-5), "visible surface depth is geometric hit distance");
  Check(stats.surfaceHits == 1, "one diffuse surface hit counted");
  Check(stats.materialClosures == 1, "one material closure constructed");
  Check(stats.infiniteLightHits == 1, "infinite light miss counted");
  Check(stats.invalidSamples == 0, "no invalid samples in closed-form scene");
}

void TestAreaLightHit() {
  SpectralTestScene testScene;
  Shape lightShape = Shape::Sphere(1, point3f(0, 0, -3));
  DiffuseAreaLightOptions lightOptions;
  lightOptions.twoSided = false;
  LightHandle light = testScene.lights.Add(
    Light::DiffuseArea(lightShape, LightSpectrum::Constant(static_cast<Float>(3)), lightOptions)
  );
  MaterialHandle material = testScene.materials.Add(Material::Interface());
  ShapeHandle shape = testScene.scene.AddShape(lightShape);
  PrimitiveBinding binding;
  binding.material = material;
  binding.areaLight = light;
  testScene.scene.AddPrimitive(shape, binding);

  RandomWalkRenderOptions options = OneBounceOptions();
  options.maxDepth = 0;
  RandomWalkIntegrator integrator(testScene.scene, testScene.materials, testScene.lights, options);
  SpectralRandomSampler sampler(options.seed);
  sampler.StartPixelSample({0, 0}, 0);
  ScratchBuffer scratch(4096);
  SampledWavelengths lambda = SampledWavelengths::SampleUniform(static_cast<Float>(0.5));
  RandomWalkRenderStats stats;

  SampledSpectrum L = integrator.Li(
    Ray(point3f(0, 0, 0), vec3f(0, 0, -1)),
    lambda,
    sampler,
    scratch,
    nullptr,
    &stats
  );

  CheckSpectrumApprox(L, SampledSpectrum(3), static_cast<Float>(1e-6), "visible area light emission is accumulated before depth termination");
  Check(stats.areaLightHits == 1, "area light hit counted");
  Check(stats.infiniteLightHits == 0, "area light hit does not use infinite light path");
}

void TestRenderRandomWalkFilmOutput() {
  SpectralTestScene testScene;
  AddDiffuseSphere(testScene, SpectrumTexture::Constant(static_cast<Float>(0.5)));
  AddUniformInfinite(testScene, LightSpectrum::Constant(static_cast<Float>(2)));

  Film film = MakeFilm(MakeFlatSensor());
  SpectralCamera camera = MakeCenterCamera();
  RandomWalkRenderOptions options = OneBounceOptions();
  RandomWalkRenderStats stats = RenderRandomWalk(
    testScene.scene,
    testScene.materials,
    testScene.lights,
    camera,
    film,
    options
  );

  SpectralRandomSampler expectedSampler(options.seed);
  expectedSampler.StartPixelSample({0, 0}, 0);
  SampledWavelengths expectedLambda = film.SampleWavelengths(expectedSampler.Get1D());
  Film expectedFilm = MakeFilm(MakeFlatSensor());
  expectedFilm.AddFilteredSample(
    FilmPoint2f{static_cast<Float>(0.5), static_cast<Float>(0.5)},
    SampledSpectrum(1),
    expectedLambda,
    nullptr,
    1
  );

  CheckRGBApprox(
    film.GetPixelLinearRGB({0, 0}),
    expectedFilm.GetPixelLinearRGB({0, 0}),
    static_cast<Float>(2e-4),
    "RenderRandomWalk sends camera, Li, camera weight, and Film conversion through one pipeline"
  );
  Check(stats.pixelSamples == 1, "one pixel sample rendered");
  Check(stats.cameraRays == 1, "one camera ray generated");
}

RGB RenderInfiniteBand(Float center) {
  SpectralTestScene testScene;
  AddUniformInfinite(testScene, LightSpectrum(GaussianSpectrum(center, 8)));
  Film film = MakeFilm(MakeBandSensor());
  SpectralCamera camera = MakeCenterCamera();

  RandomWalkRenderOptions options;
  options.pixelSamples = 4096;
  options.maxDepth = 0;
  options.seed = static_cast<std::uint64_t>(center);
  options.jitterCameraSamples = false;
  RenderRandomWalk(testScene.scene, testScene.materials, testScene.lights, camera, film, options);
  return film.GetPixelLinearRGB({0, 0});
}

void TestNarrowBandEmitterSensorResponses() {
  RGB red = RenderInfiniteBand(650);
  RGB green = RenderInfiniteBand(540);
  RGB blue = RenderInfiniteBand(450);

  Check(red.r > red.g * static_cast<Float>(8) && red.r > red.b * static_cast<Float>(8), "red narrow-band emitter dominates red sensor channel");
  Check(green.g > green.r * static_cast<Float>(8) && green.g > green.b * static_cast<Float>(8), "green narrow-band emitter dominates green sensor channel");
  Check(blue.b > blue.r * static_cast<Float>(8) && blue.b > blue.g * static_cast<Float>(8), "blue narrow-band emitter dominates blue sensor channel");
}

void TestRGBRoleSeparation(const std::string& assetDirectory) {
  RGBColorSpace srgb = LoadSRGBColorSpace(assetDirectory);
  RGB rgb(static_cast<Float>(0.85), static_cast<Float>(0.2), static_cast<Float>(0.08));
  RGBAlbedoSpectrum albedoRole(srgb, rgb);
  RGBIlluminantSpectrum illuminantRole(srgb, rgb);
  Spectrum albedoSpectrum(DenselySampledSpectrum::SampleFunction([&](Float lambda) {
    return albedoRole(lambda);
  }));
  Spectrum illuminantSpectrum(DenselySampledSpectrum::SampleFunction([&](Float lambda) {
    return illuminantRole(lambda);
  }));

  SampledWavelengths lambda = SampledWavelengths::SampleUniform(static_cast<Float>(0.37));
  SampledSpectrum albedo = albedoSpectrum.Sample(lambda);
  SampledSpectrum illuminant = illuminantSpectrum.Sample(lambda);
  Float maxDifference = 0;
  for (int i = 0; i < NSpectrumSamples; ++i) {
    maxDifference = std::max(maxDifference, std::abs(albedo[i] - illuminant[i]));
  }
  Check(maxDifference > static_cast<Float>(0.05), "RGB albedo and RGB illuminant roles produce distinct spectra");

  SpectralTestScene albedoLightScene;
  AddUniformInfinite(albedoLightScene, LightSpectrum(albedoSpectrum));
  Film albedoFilm = MakeFilm(MakeBandSensor());
  RenderRandomWalk(
    albedoLightScene.scene,
    albedoLightScene.materials,
    albedoLightScene.lights,
    MakeCenterCamera(),
    albedoFilm,
    OneBounceOptions()
  );

  SpectralTestScene illuminantLightScene;
  AddUniformInfinite(illuminantLightScene, LightSpectrum::FromRGBIlluminant(srgb, rgb));
  Film illuminantFilm = MakeFilm(MakeBandSensor());
  RenderRandomWalk(
    illuminantLightScene.scene,
    illuminantLightScene.materials,
    illuminantLightScene.lights,
    MakeCenterCamera(),
    illuminantFilm,
    OneBounceOptions()
  );

  RGB albedoRGB = albedoFilm.GetPixelLinearRGB({0, 0});
  RGB illuminantRGB = illuminantFilm.GetPixelLinearRGB({0, 0});
  Float rgbDifference = std::abs(albedoRGB.r - illuminantRGB.r) +
                        std::abs(albedoRGB.g - illuminantRGB.g) +
                        std::abs(albedoRGB.b - illuminantRGB.b);
  Check(rgbDifference > static_cast<Float>(1), "RGB roles remain distinct after random-walk Film conversion");
}

void TestDeltaLightRejectionAndDiagnostics() {
  SpectralTestScene testScene;
  RegisterLight(
    testScene.scene,
    testScene.lights,
    Light::Point(point3f(0, 2, -1), LightSpectrum::Constant(1))
  );

  bool threw = false;
  try {
    RandomWalkIntegrator integrator(testScene.scene, testScene.materials, testScene.lights, OneBounceOptions());
    (void)integrator;
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  Check(threw, "RandomWalkIntegrator rejects delta lights in PR14");

  RandomWalkRenderStats stats;
  Check(!RecordSpectralRadianceDiagnostics(SampledSpectrum{-1, 0, 0, 0}, &stats), "negative radiance diagnostic rejects sample");
  Check(stats.negativeRadiance == 1, "negative radiance diagnostic counted");
  Check(!RecordSpectralRadianceDiagnostics(
          SampledSpectrum{std::numeric_limits<Float>::infinity(), 0, 0, 0},
          &stats
        ),
        "non-finite radiance diagnostic rejects sample");
  Check(stats.nonFiniteRadiance == 1, "non-finite radiance diagnostic counted");
}

} // namespace

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "usage: pr14-randomwalk-tests <spectral-asset-directory>" << std::endl;
    return 1;
  }

  try {
    TestDiffuseInfiniteLiClosedForm();
    TestAreaLightHit();
    TestRenderRandomWalkFilmOutput();
    TestNarrowBandEmitterSensorResponses();
    TestRGBRoleSeparation(argv[1]);
    TestDeltaLightRejectionAndDiagnostics();
  } catch (const std::exception& error) {
    std::cerr << "PR14 RandomWalk test failed with exception: " << error.what() << std::endl;
    return 1;
  }

  std::cout << "PR14 RandomWalk tests passed" << std::endl;
  return 0;
}
