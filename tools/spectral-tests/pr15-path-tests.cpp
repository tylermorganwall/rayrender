#include "src/render/spectral_integrator.h"

#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <optional>
#include <utility>

using namespace rayrender::base;
using namespace rayrender::materials;
using namespace rayrender::render;

namespace {

constexpr Float Pi = static_cast<Float>(3.14159265358979323846264338327950288);

bool Approx(Float lhs, Float rhs, Float tolerance = static_cast<Float>(1e-5)) {
  return std::fabs(lhs - rhs) <= tolerance;
}

void Check(bool condition, const char* message) {
  if (!condition) {
    std::cerr << "PR15 PathIntegrator test failed: " << message << std::endl;
    std::exit(1);
  }
}

void CheckApprox(Float actual, Float expected, Float tolerance, const char* message) {
  if (!Approx(actual, expected, tolerance)) {
    std::cerr << "PR15 PathIntegrator test failed: " << message
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

struct SpectralTestScene {
  Scene scene;
  SpectralMaterialTable materials;
  SpectralLightTable lights;
};

SampledWavelengths TestWavelengths() {
  return SampledWavelengths::SampleUniform(static_cast<Float>(0.42));
}

void PreprocessLights(SpectralTestScene& testScene) {
  Bounds3f bounds;
  testScene.scene.Bounds(&bounds);
  testScene.lights.Preprocess(bounds);
}

void AddDiffuseSphere(
  SpectralTestScene& testScene,
  Float reflectance = static_cast<Float>(0.5),
  point3f center = point3f(0, 0, -3)
) {
  MaterialHandle material =
    testScene.materials.Add(Material::Diffuse(SpectrumTexture::Constant(reflectance)));
  ShapeHandle shape = testScene.scene.AddShape(Shape::Sphere(1, center));
  PrimitiveBinding binding;
  binding.material = material;
  testScene.scene.AddPrimitive(shape, binding);
}

Shape FixedDirectionAreaShape(point3f pLight, normal3f nLight, Float pdf) {
  ShapeCapabilities capabilities;
  capabilities.canIntersect = false;
  capabilities.canBound = false;
  capabilities.canSampleArea = false;
  capabilities.canSampleDirection = true;
  capabilities.supportsAreaLight = true;

  ShapeCallbacks callbacks;
  callbacks.area = []() {
    return static_cast<Float>(1);
  };
  callbacks.sampleDirection = [=](const Interaction&, point2f) -> std::optional<ShapeSample> {
    ShapeSample sample;
    sample.interaction.p = pLight;
    sample.interaction.pError = vec3f(static_cast<Float>(1e-4));
    sample.interaction.n = nLight;
    sample.pdf = pdf;
    return sample;
  };
  callbacks.pdfDirection = [=](const Interaction&, const vec3f&) {
    return pdf;
  };
  return Shape::FromCallbacks("FixedDirectionAreaShape", capabilities, std::move(callbacks));
}

Shape InfinitePlane(Float z, normal3f normal) {
  ShapeCapabilities capabilities;
  capabilities.canIntersect = true;
  capabilities.canBound = false;
  capabilities.hasUV = false;
  capabilities.supportsNormalMap = false;
  capabilities.supportsDisplacement = false;
  capabilities.supportsAreaLight = false;

  ShapeCallbacks callbacks;
  callbacks.intersect = [=](const Ray& ray, Float tMin, Float tMax) -> std::optional<ShapeIntersection> {
    Float dz = ray.direction().xyz.z;
    if (std::fabs(dz) < static_cast<Float>(1e-7)) {
      return std::nullopt;
    }
    Float t = (z - ray.origin().xyz.z) / dz;
    if (t <= tMin || t >= tMax) {
      return std::nullopt;
    }

    LegacySurfaceHit hit;
    hit.p = ray(t);
    hit.t = t;
    hit.normal = normal;
    hit.pError = vec3f(static_cast<Float>(1e-4));
    hit.dpdu = vec3f(1, 0, 0);
    hit.dpdv = vec3f(0, 1, 0);

    ShapeIntersection intersection;
    intersection.tHit = t;
    intersection.interaction = SurfaceInteractionFromLegacyHit(hit, ray);
    return intersection;
  };
  callbacks.hitP = [callbacks](const Ray& ray, Float tMin, Float tMax) {
    return callbacks.intersect(ray, tMin, tMax).has_value();
  };
  return Shape::FromCallbacks("InfinitePlane", capabilities, std::move(callbacks));
}

void AddDiffusePlane(
  SpectralTestScene& testScene,
  Shape plane,
  Float reflectance
) {
  MaterialHandle material =
    testScene.materials.Add(Material::Diffuse(SpectrumTexture::Constant(reflectance)));
  ShapeHandle shape = testScene.scene.AddShape(std::move(plane));
  PrimitiveBinding binding;
  binding.material = material;
  testScene.scene.AddPrimitive(shape, binding);
}

PathRenderOptions OneBouncePathOptions() {
  PathRenderOptions options;
  options.maxDepth = 1;
  options.pixelSamples = 1;
  options.seed = 29;
  options.jitterCameraSamples = false;
  return options;
}

SampledSpectrum EvaluatePathLi(
  const SpectralTestScene& testScene,
  const PathRenderOptions& options,
  PathRenderStats& stats,
  int sampleIndex = 0
) {
  PathIntegrator integrator(testScene.scene, testScene.materials, testScene.lights, options);
  SpectralRandomSampler sampler(options.seed);
  sampler.StartPixelSample({0, 0}, sampleIndex);
  ScratchBuffer scratch(options.scratchBufferBytes);
  SampledWavelengths lambda = TestWavelengths();
  return integrator.Li(
    Ray(point3f(0, 0, 0), vec3f(0, 0, -1)),
    lambda,
    sampler,
    scratch,
    nullptr,
    &stats
  );
}

void TestPowerHeuristic() {
  CheckApprox(
    PowerHeuristic(1, static_cast<Float>(0.25), 1, static_cast<Float>(0.5)),
    static_cast<Float>(0.2),
    static_cast<Float>(1e-6),
    "PowerHeuristic matches pbrt square-weight formula"
  );
}

void TestDeltaDirectLighting() {
  SpectralTestScene testScene;
  AddDiffuseSphere(testScene, static_cast<Float>(0.5));
  RegisterLight(
    testScene.scene,
    testScene.lights,
    Light::Point(point3f(0, 0, 2), LightSpectrum::Constant(static_cast<Float>(32) * Pi))
  );
  PreprocessLights(testScene);

  PathRenderOptions options = OneBouncePathOptions();
  options.sampleBSDF = false;
  PathRenderStats stats;
  SampledSpectrum L = EvaluatePathLi(testScene, options, stats);

  CheckSpectrumApprox(L, SampledSpectrum(1), static_cast<Float>(2e-4), "delta direct-light estimator has pbrt point-light algebra");
  Check(stats.directLightSamples == 1, "delta test sampled one direct light");
  Check(stats.deltaLightSamples == 1, "delta light bypasses competing BSDF PDF");
  Check(stats.directLightContributions == 1, "delta direct light contributes radiance");
}

void TestAreaDirectMIS() {
  SpectralTestScene testScene;
  AddDiffuseSphere(testScene, static_cast<Float>(0.5));

  Float lightPDF = static_cast<Float>(0.25);
  Shape fixedLight = FixedDirectionAreaShape(
    point3f(0, 0, 4),
    normal3f(0, 0, -1),
    lightPDF
  );
  DiffuseAreaLightOptions twoSided;
  twoSided.twoSided = true;
  RegisterLight(
    testScene.scene,
    testScene.lights,
    Light::DiffuseArea(fixedLight, LightSpectrum::Constant(static_cast<Float>(4)), twoSided)
  );
  PreprocessLights(testScene);

  PathRenderOptions options = OneBouncePathOptions();
  options.sampleBSDF = false;
  PathRenderStats stats;
  SampledSpectrum L = EvaluatePathLi(testScene, options, stats);

  Float pLight = lightPDF;
  Float pBSDF = static_cast<Float>(1) / Pi;
  Float expected = static_cast<Float>(4) *
                   (static_cast<Float>(0.5) / Pi) *
                   PowerHeuristic(1, pLight, 1, pBSDF) / pLight;
  CheckSpectrumApprox(L, SampledSpectrum(expected), static_cast<Float>(2e-4), "area direct-light MIS weight matches pbrt formula");
  Check(stats.directLightSamples == 1, "area MIS sampled one direct light");
  Check(stats.directLightContributions == 1, "area MIS contributed direct light");
}

void TestBSDFSampledInfiniteMIS() {
  SpectralTestScene testScene;
  AddDiffuseSphere(testScene, static_cast<Float>(0.5));
  RegisterLight(
    testScene.scene,
    testScene.lights,
    Light::UniformInfinite(LightSpectrum::Constant(static_cast<Float>(2)))
  );
  PreprocessLights(testScene);

  PathRenderOptions options = OneBouncePathOptions();
  options.sampleDirectLighting = false;
  PathRenderStats stats;
  SampledSpectrum L = EvaluatePathLi(testScene, options, stats);

  SpectralRandomSampler expectedSampler(options.seed);
  expectedSampler.StartPixelSample({0, 0}, 0);
  (void)expectedSampler.Get1D(); // alpha
  (void)expectedSampler.Get1D(); // BSDF component sample
  point2f u = expectedSampler.Get2D();
  Float cosTheta = std::sqrt(static_cast<Float>(1) - u.xy.x);
  Float pBSDF = cosTheta / Pi;
  Float pLight = static_cast<Float>(1) / (static_cast<Float>(4) * Pi);
  Float expected = PowerHeuristic(1, pBSDF, 1, pLight);

  CheckSpectrumApprox(L, SampledSpectrum(expected), static_cast<Float>(2e-5), "BSDF-sampled infinite emitter hit uses pbrt MIS weight");
  Check(stats.emitterHitMIS == 1, "infinite emitter hit uses MIS after a non-specular bounce");
  Check(stats.infiniteLightHits == 1, "infinite emitter hit counted");
}

void TestCombinedUniformInfiniteMISConverges() {
  SpectralTestScene testScene;
  AddDiffuseSphere(testScene, static_cast<Float>(0.5));
  RegisterLight(
    testScene.scene,
    testScene.lights,
    Light::UniformInfinite(LightSpectrum::Constant(static_cast<Float>(2)))
  );
  PreprocessLights(testScene);

  PathRenderOptions options = OneBouncePathOptions();
  options.seed = 991;
  options.russianRoulette = false;
  PathIntegrator integrator(testScene.scene, testScene.materials, testScene.lights, options);
  ScratchBuffer scratch(options.scratchBufferBytes);
  SpectralRandomSampler sampler(options.seed);

  SampledSpectrum sum(0);
  constexpr int sampleCount = 8192;
  PathRenderStats stats;
  for (int i = 0; i < sampleCount; ++i) {
    sampler.StartPixelSample({0, 0}, i);
    scratch.Reset();
    SampledWavelengths lambda = TestWavelengths();
    sum += integrator.Li(
      Ray(point3f(0, 0, 0), vec3f(0, 0, -1)),
      lambda,
      sampler,
      scratch,
      nullptr,
      &stats
    );
  }
  SampledSpectrum mean = sum / static_cast<Float>(sampleCount);
  CheckSpectrumApprox(mean, SampledSpectrum(1), static_cast<Float>(0.04), "combined direct and BSDF MIS converges to diffuse infinite-light reference");
  Check(stats.directLightContributions > 0, "combined MIS includes direct-light contributions");
  Check(stats.infiniteLightHits > 0, "combined MIS includes BSDF-sampled environment hits");
}

void TestMaxDepthEmissionOrdering() {
  SpectralTestScene testScene;
  Shape lightShape = Shape::Sphere(1, point3f(0, 0, -3));
  LightHandle light = testScene.lights.Add(
    Light::DiffuseArea(lightShape, LightSpectrum::Constant(static_cast<Float>(3)))
  );
  MaterialHandle material =
    testScene.materials.Add(Material::Diffuse(SpectrumTexture::Constant(static_cast<Float>(0.5))));
  ShapeHandle shape = testScene.scene.AddShape(lightShape);
  PrimitiveBinding binding;
  binding.material = material;
  binding.areaLight = light;
  testScene.scene.AddPrimitive(shape, binding);
  AttachLightToScene(testScene.scene, light, testScene.lights.Get(light));
  PreprocessLights(testScene);

  PathRenderOptions options = OneBouncePathOptions();
  options.maxDepth = 0;
  PathRenderStats stats;
  SampledSpectrum L = EvaluatePathLi(testScene, options, stats);

  CheckSpectrumApprox(L, SampledSpectrum(3), static_cast<Float>(1e-6), "area emission is added before max-depth termination");
  Check(stats.areaLightHits == 1, "area emission counted at max depth");
  Check(stats.maxDepthTerminations == 1, "max-depth termination follows pbrt surface ordering");
  Check(stats.directLightSamples == 0, "max-depth termination prevents direct-light sampling");
}

void TestRussianRouletteAccounting() {
  SpectralTestScene testScene;
  AddDiffusePlane(testScene, InfinitePlane(static_cast<Float>(-1), normal3f(0, 0, 1)), static_cast<Float>(0.05));
  AddDiffusePlane(testScene, InfinitePlane(static_cast<Float>(1), normal3f(0, 0, -1)), static_cast<Float>(0.05));
  PreprocessLights(testScene);

  PathRenderOptions rouletteOptions;
  rouletteOptions.maxDepth = 10;
  rouletteOptions.seed = 410;
  rouletteOptions.sampleDirectLighting = false;
  rouletteOptions.russianRoulette = true;

  bool sawTermination = false;
  for (int i = 0; i < 64 && !sawTermination; ++i) {
    PathRenderStats stats;
    (void)EvaluatePathLi(testScene, rouletteOptions, stats, i);
    sawTermination = stats.russianRouletteTerminations > 0;
  }
  Check(sawTermination, "Russian roulette terminates low-throughput multi-bounce paths");

  PathRenderOptions disabled = rouletteOptions;
  disabled.russianRoulette = false;
  PathRenderStats disabledStats;
  (void)EvaluatePathLi(testScene, disabled, disabledStats, 0);
  Check(disabledStats.russianRouletteChecks == 0, "roulette disabled skips roulette checks");
  Check(disabledStats.russianRouletteTerminations == 0, "roulette disabled has no roulette terminations");
}

} // namespace

int main() {
  try {
    TestPowerHeuristic();
    TestDeltaDirectLighting();
    TestAreaDirectMIS();
    TestBSDFSampledInfiniteMIS();
    TestCombinedUniformInfiniteMISConverges();
    TestMaxDepthEmissionOrdering();
    TestRussianRouletteAccounting();
  } catch (const std::exception& error) {
    std::cerr << "PR15 PathIntegrator test failed with exception: " << error.what() << std::endl;
    return 1;
  }

  std::cout << "PR15 PathIntegrator tests passed" << std::endl;
  return 0;
}
