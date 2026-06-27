#include "src/render/spectral_integrator.h"

#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <optional>
#include <string>

using namespace rayrender::base;
using namespace rayrender::materials;
using namespace rayrender::render;

namespace {

bool Approx(Float lhs, Float rhs, Float tolerance = static_cast<Float>(1e-5)) {
  return std::fabs(lhs - rhs) <= tolerance;
}

void Check(bool condition, const char* message) {
  if (!condition) {
    std::cerr << "PR16 Conductor test failed: " << message << std::endl;
    std::exit(1);
  }
}

void CheckApprox(Float actual, Float expected, Float tolerance, const char* message) {
  if (!Approx(actual, expected, tolerance)) {
    std::cerr << "PR16 Conductor test failed: " << message
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
    if (!Approx(actual[i], expected[i], tolerance)) {
      std::cerr << "PR16 Conductor test failed: " << message
                << " component " << i << " expected " << expected[i]
                << " got " << actual[i] << std::endl;
      std::exit(1);
    }
  }
}

vec3f Unit(vec3f value) {
  return value / value.length();
}

Float ReferenceFrComplexNormal(Float eta, Float k) {
  Float numerator = (eta - 1) * (eta - 1) + k * k;
  Float denominator = (eta + 1) * (eta + 1) + k * k;
  return numerator / denominator;
}

SampledWavelengths TestWavelengths() {
  return SampledWavelengths::SampleUniform(static_cast<Float>(0.42));
}

MaterialEvalContext BaseContext() {
  MaterialEvalContext ctx;
  ctx.p = point3f(0, 0, 0);
  ctx.uv = point2f(static_cast<Float>(0.5), static_cast<Float>(0.5));
  ctx.dpdu = vec3f(1, 0, 0);
  ctx.dpdv = vec3f(0, 1, 0);
  ctx.dpdus = ctx.dpdu;
  ctx.dpdvs = ctx.dpdv;
  ctx.n = normal3f(0, 0, 1);
  ctx.ns = normal3f(0, 0, 1);
  ctx.dndus = normal3f(0, 0, 0);
  ctx.dndvs = normal3f(0, 0, 0);
  ctx.wo = vec3f(0, 0, 1);
  return ctx;
}

struct CountingTextureEvaluator {
  mutable int floatCount = 0;
  mutable int spectrumCount = 0;
  UniversalTextureEvaluator evaluator;

  Float operator()(const FloatTexture& texture, const TextureEvalContext& ctx) const {
    ++floatCount;
    return evaluator(texture, ctx);
  }

  SampledSpectrum operator()(
    const SpectrumTexture& texture,
    const TextureEvalContext& ctx,
    const SampledWavelengths& lambda
  ) const {
    ++spectrumCount;
    return evaluator(texture, ctx, lambda);
  }
};

struct SpectralTestScene {
  Scene scene;
  SpectralMaterialTable materials;
  SpectralLightTable lights;
};

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

void AddConductorPlane(
  SpectralTestScene& testScene,
  const Spectrum& eta,
  const Spectrum& k
) {
  MaterialHandle material = testScene.materials.Add(Material::Conductor(
    ConductorMaterial::FromEtaK(
      SpectrumTexture::Constant(eta),
      SpectrumTexture::Constant(k),
      FloatTexture::Constant(0),
      false
    )
  ));
  ShapeHandle shape = testScene.scene.AddShape(InfinitePlane(static_cast<Float>(-1), normal3f(0, 0, 1)));
  PrimitiveBinding binding;
  binding.material = material;
  testScene.scene.AddPrimitive(shape, binding);
}

SampledSpectrum EvaluateSmoothConductorPlane(
  const Spectrum& eta,
  const Spectrum& k,
  SampledWavelengths& lambda,
  PathRenderStats& stats
) {
  SpectralTestScene testScene;
  AddConductorPlane(testScene, eta, k);
  RegisterLight(
    testScene.scene,
    testScene.lights,
    Light::UniformInfinite(LightSpectrum::Constant(1))
  );
  Bounds3f bounds;
  testScene.scene.Bounds(&bounds);
  testScene.lights.Preprocess(bounds);

  PathRenderOptions options;
  options.maxDepth = 1;
  options.seed = 16;
  options.jitterCameraSamples = false;
  options.sampleDirectLighting = true;
  options.russianRoulette = false;

  PathIntegrator integrator(testScene.scene, testScene.materials, testScene.lights, options);
  SpectralRandomSampler sampler(options.seed);
  sampler.StartPixelSample({0, 0}, 0);
  ScratchBuffer scratch(options.scratchBufferBytes);
  return integrator.Li(
    Ray(point3f(0, 0, 0), vec3f(0, 0, -1)),
    lambda,
    sampler,
    scratch,
    nullptr,
    &stats
  );
}

void TestSmoothConductorBxDF() {
  SampledSpectrum eta({static_cast<Float>(0.2), static_cast<Float>(0.5), 1, 2});
  SampledSpectrum k({static_cast<Float>(3), static_cast<Float>(2), static_cast<Float>(1.5), static_cast<Float>(4)});
  ConductorBxDF conductor(TrowbridgeReitzDistribution(0, 0), eta, k);

  Check(conductor.Flags() == BxDFFlags::SpecularReflection, "smooth conductor is specular reflection");
  vec3f wo = Unit(vec3f(static_cast<Float>(0.2), static_cast<Float>(0.3), static_cast<Float>(0.9327379)));
  std::optional<BSDFSample> sample = conductor.Sample_f(
    wo,
    0,
    point2f(static_cast<Float>(0.37), static_cast<Float>(0.61)),
    TransportMode::Radiance
  );
  Check(sample.has_value(), "smooth conductor samples a reflection");
  Check(sample->IsSpecular(), "smooth conductor sample is specular");
  Check(sample->IsReflection(), "smooth conductor sample is reflective");
  CheckApprox(sample->wi.xyz.x, -wo.xyz.x, static_cast<Float>(1e-6), "smooth conductor reflected x");
  CheckApprox(sample->wi.xyz.y, -wo.xyz.y, static_cast<Float>(1e-6), "smooth conductor reflected y");
  CheckApprox(sample->wi.xyz.z, wo.xyz.z, static_cast<Float>(1e-6), "smooth conductor reflected z");
  CheckApprox(sample->pdf, 1, static_cast<Float>(1e-6), "smooth conductor PDF is delta mass");

  SampledSpectrum expected =
    FrComplex(rayrender::render::AbsCosTheta(sample->wi), eta, k) /
    rayrender::render::AbsCosTheta(sample->wi);
  CheckSpectrumApprox(sample->f, expected, static_cast<Float>(1e-6), "smooth conductor sample f matches pbrt formula");
  CheckSpectrumApprox(conductor.f(wo, sample->wi, TransportMode::Radiance), SampledSpectrum(0), static_cast<Float>(1e-6), "smooth conductor f is zero off the delta path");
  CheckApprox(conductor.PDF(wo, sample->wi, TransportMode::Radiance), 0, static_cast<Float>(1e-6), "smooth conductor PDF() is zero off the delta path");

  SampledSpectrum rho = conductor.rho();
  for (int i = 0; i < NSpectrumSamples; ++i) {
    CheckApprox(
      rho[i],
      ReferenceFrComplexNormal(eta[i], k[i]),
      static_cast<Float>(1e-6),
      "smooth conductor rho normal-incidence estimate"
    );
  }
}

void TestRoughConductorBxDF() {
  SampledSpectrum eta(static_cast<Float>(0.6));
  SampledSpectrum k(static_cast<Float>(2.4));
  ConductorBxDF conductor(
    TrowbridgeReitzDistribution(static_cast<Float>(0.35), static_cast<Float>(0.7)),
    eta,
    k
  );
  vec3f wo = Unit(vec3f(static_cast<Float>(0.25), static_cast<Float>(0.15), 1));
  std::optional<BSDFSample> sample = conductor.Sample_f(
    wo,
    0,
    point2f(static_cast<Float>(0.31), static_cast<Float>(0.73)),
    TransportMode::Radiance
  );
  Check(sample.has_value(), "rough conductor samples a reflection");
  Check(sample->IsGlossy(), "rough conductor sample is glossy");
  Check(sample->IsReflection(), "rough conductor sample is reflective");
  Check(rayrender::render::SameHemisphere(wo, sample->wi), "rough conductor sample remains in reflection hemisphere");
  Check(sample->pdf > 0 && std::isfinite(sample->pdf), "rough conductor sample PDF is finite positive");

  Float pdf = conductor.PDF(wo, sample->wi, TransportMode::Radiance);
  CheckApprox(sample->pdf, pdf, static_cast<Float>(1e-6), "rough conductor sampled PDF agrees with PDF()");
  CheckSpectrumApprox(sample->f, conductor.f(wo, sample->wi, TransportMode::Radiance), static_cast<Float>(1e-5), "rough conductor sampled f agrees with f()");
  Check(!conductor.Sample_f(
    wo,
    0,
    point2f(static_cast<Float>(0.2), static_cast<Float>(0.4)),
    TransportMode::Radiance,
    BxDFReflTransFlags::Transmission
  ).has_value(), "rough conductor respects transmission-only sample flags");

  ConductorBxDF regularized(
    TrowbridgeReitzDistribution(static_cast<Float>(0.05), static_cast<Float>(0.2)),
    eta,
    k
  );
  regularized.Regularize();
  std::optional<BSDFSample> regularizedSample = regularized.Sample_f(
    wo,
    0,
    point2f(static_cast<Float>(0.41), static_cast<Float>(0.19)),
    TransportMode::Radiance
  );
  Check(regularizedSample.has_value(), "regularized conductor remains sampleable");
}

void TestConductorMaterialEtaKAndPacket() {
  Material material = Material::Conductor(ConductorMaterial::FromEtaK(
    SpectrumTexture::Constant(static_cast<Float>(0.7)),
    SpectrumTexture::Constant(static_cast<Float>(2.1)),
    FloatTexture::Constant(static_cast<Float>(0.25)),
    FloatTexture::Constant(static_cast<Float>(0.64)),
    true
  ));
  MaterialEvalContext ctx = BaseContext();
  SampledWavelengths lambda = TestWavelengths();
  CountingTextureEvaluator evaluator;
  ScratchBuffer scratch(1024);

  BSDF bsdf = material.GetBSDF(evaluator, ctx, lambda, scratch);
  Check(bsdf.Flags() == BxDFFlags::GlossyReflection, "rough conductor material creates glossy reflection BSDF");
  Check(evaluator.spectrumCount == 2, "eta/k conductor evaluates two spectrum textures");
  Check(evaluator.floatCount == 2, "anisotropic conductor evaluates two roughness textures");
  Check(!lambda.SecondaryTerminated(), "conductor material leaves all wavelength packet components active");
  CheckSpectrumApprox(
    bsdf.rho(),
    FrComplex(1, SampledSpectrum(static_cast<Float>(0.7)), SampledSpectrum(static_cast<Float>(2.1))),
    static_cast<Float>(1e-6),
    "conductor BSDF rho exposes spectral Fresnel estimate"
  );
}

void TestReflectanceCompatibilityMaterial() {
  Material material = Material::Conductor(ConductorMaterial::FromReflectance(
    SpectrumTexture::Constant(static_cast<Float>(0.25)),
    FloatTexture::Constant(0),
    true
  ));
  MaterialEvalContext ctx = BaseContext();
  SampledWavelengths lambda = TestWavelengths();
  CountingTextureEvaluator evaluator;
  ScratchBuffer scratch(1024);

  BSDF bsdf = material.GetBSDF(evaluator, ctx, lambda, scratch);
  std::optional<BSDFSample> sample = bsdf.Sample_f(
    vec3f(0, 0, 1),
    0,
    point2f(static_cast<Float>(0.5), static_cast<Float>(0.5))
  );
  Check(sample.has_value(), "reflectance compatibility conductor samples smooth reflection");
  CheckSpectrumApprox(sample->f, SampledSpectrum(static_cast<Float>(0.25)), static_cast<Float>(2e-5), "reflectance compatibility produces requested normal reflectance");
  CheckSpectrumApprox(bsdf.rho(), SampledSpectrum(static_cast<Float>(0.25)), static_cast<Float>(2e-5), "reflectance compatibility rho matches requested reflectance");
  Check(evaluator.spectrumCount == 1, "reflectance compatibility evaluates one spectrum texture");
  Check(evaluator.floatCount == 2, "reflectance compatibility evaluates roughness as u/v textures");
}

void TestNamedMetalReferenceScene(
  const NamedSpectrumRegistry& registry,
  const char* etaName,
  const char* kName
) {
  const Spectrum& etaSpectrum = registry.GetOrThrow(etaName);
  const Spectrum& kSpectrum = registry.GetOrThrow(kName);
  SampledWavelengths lambda = TestWavelengths();
  SampledWavelengths expectedLambda = lambda;
  PathRenderStats stats;
  SampledSpectrum L = EvaluateSmoothConductorPlane(etaSpectrum, kSpectrum, lambda, stats);
  SampledSpectrum expected = FrComplex(
    1,
    etaSpectrum.Sample(expectedLambda),
    kSpectrum.Sample(expectedLambda)
  );

  CheckSpectrumApprox(L, expected, static_cast<Float>(2e-5), "named smooth conductor plane matches pbrt Fresnel reference");
  Check(stats.bsdfSamples == 1, "named conductor scene samples one specular BSDF event");
  Check(stats.infiniteLightHits == 1, "named conductor scene reflects to the infinite light");
  Check(stats.directLightSamples == 0, "specular conductor is not sampled by direct-light MIS");
  Check(!lambda.SecondaryTerminated(), "named conductor path leaves wavelength packet unterminated");
}

void TestNamedMetalReferenceScenes() {
  NamedSpectrumRegistry registry = NamedSpectrumRegistry::LoadDefault();
  TestNamedMetalReferenceScene(registry, "metal-Cu-eta", "metal-Cu-k");
  TestNamedMetalReferenceScene(registry, "metal-Ag-eta", "metal-Ag-k");
  TestNamedMetalReferenceScene(registry, "metal-Au-eta", "metal-Au-k");
}

} // namespace

int main() {
  try {
    TestSmoothConductorBxDF();
    TestRoughConductorBxDF();
    TestConductorMaterialEtaKAndPacket();
    TestReflectanceCompatibilityMaterial();
    TestNamedMetalReferenceScenes();
  } catch (const std::exception& error) {
    std::cerr << "PR16 Conductor test failed with exception: " << error.what() << std::endl;
    return 1;
  }

  std::cout << "PR16 Conductor tests passed" << std::endl;
  return 0;
}
