#include "src/render/spectral_integrator.h"

#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <optional>
#include <string>
#include <utility>
#include <vector>

using namespace rayrender::base;
using namespace rayrender::materials;
using namespace rayrender::render;

namespace {

bool Approx(Float lhs, Float rhs, Float tolerance = static_cast<Float>(1e-5)) {
  return std::fabs(lhs - rhs) <= tolerance;
}

void Check(bool condition, const char* message) {
  if (!condition) {
    std::cerr << "PR18 Dielectric test failed: " << message << std::endl;
    std::exit(1);
  }
}

void CheckApprox(Float actual, Float expected, Float tolerance, const char* message) {
  if (!Approx(actual, expected, tolerance)) {
    std::cerr << "PR18 Dielectric test failed: " << message
              << " expected " << expected << " got " << actual << std::endl;
    std::exit(1);
  }
}

vec3f Unit(vec3f value) {
  return value / value.length();
}

SampledWavelengths TestWavelengths() {
  return SampledWavelengths::SampleUniform(static_cast<Float>(0.37));
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

void PreprocessLights(SpectralTestScene& testScene) {
  Bounds3f bounds;
  testScene.scene.Bounds(&bounds);
  testScene.lights.Preprocess(bounds);
}

void AddSphereBoundary(
  SpectralTestScene& testScene,
  RegionId region,
  Material material,
  point3f center = point3f(0, 0, -3)
) {
  MaterialHandle materialHandle = testScene.materials.Add(std::move(material));
  ShapeHandle shape = testScene.scene.AddShape(Shape::Sphere(static_cast<Float>(1), center));
  PrimitiveBinding binding;
  binding.material = materialHandle;
  binding.dielectricBoundaries.push_back({region, RegionSide::NegativeNormal});
  testScene.scene.AddPrimitive(shape, binding);
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
  callbacks.intersect = [=](const Ray& ray, Float tMin, Float tMax)
    -> std::optional<ShapeIntersection> {
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
  callbacks.bounds = [](Bounds3f*) { return false; };

  return Shape::FromCallbacks("infinite-plane", capabilities, std::move(callbacks));
}

void TestEtaSpectrumResolution() {
  EtaSpectrumHandle cauchy = CauchyEtaSpectrum(
    static_cast<Float>(1.5),
    static_cast<Float>(0.004)
  );
  Check(!EtaSpectrumIsConstant(cauchy), "Cauchy dispersion is nonconstant");
  Check(
    EvaluateEtaSpectrum(cauchy, static_cast<Float>(450)) >
      EvaluateEtaSpectrum(cauchy, static_cast<Float>(650)),
    "Cauchy eta decreases toward longer wavelengths"
  );

  DielectricRegionTable table;
  RegionId glass = table.AddRegion(cauchy, 0, "cauchy glass");
  DielectricPathState state = DielectricPathState::FromInitialRegions(&table, {});
  ResolvedDielectricTransition enterGlass =
    state.Analyze({{glass, RegionSide::NegativeNormal}}, normal3f(0, 0, 1), vec3f(0, 0, -1));
  Check(
    enterGlass.kind == DielectricBoundaryKind::ScatteringInterface,
    "dispersive exterior-to-glass boundary scatters"
  );
  Check(!enterGlass.interface.ratioIsConstant, "exterior-to-Cauchy ratio is dispersive");
  Check(!enterGlass.interface.IsUnity(), "exterior-to-Cauchy ratio is not unity");

  RegionId matched = table.AddRegion(cauchy, -1, "same eta region");
  state.Commit(enterGlass.token);
  ResolvedDielectricTransition enterMatched =
    state.Analyze({{matched, RegionSide::NegativeNormal}}, normal3f(0, 0, 1), vec3f(0, 0, -1));
  Check(
    enterMatched.kind == DielectricBoundaryKind::IndexMatchedNull,
    "shared dispersive eta handle is exact index-matched null"
  );
  Check(enterMatched.interface.ratioIsConstant, "shared eta handle proves constant ratio");
  Check(enterMatched.interface.IsUnity(), "shared eta handle proves unity");

  DielectricRegionTable scaledTable;
  RegionId outer = scaledTable.AddRegion(
    CauchyEtaSpectrum(static_cast<Float>(1.5), static_cast<Float>(0.006)),
    0,
    "outer"
  );
  RegionId inner = scaledTable.AddRegion(
    CauchyEtaSpectrum(static_cast<Float>(3.0), static_cast<Float>(0.012)),
    -1,
    "inner"
  );
  DielectricPathState scaledState =
    DielectricPathState::FromInitialRegions(&scaledTable, {outer});
  ResolvedDielectricTransition enterInner =
    scaledState.Analyze({{inner, RegionSide::NegativeNormal}}, normal3f(0, 0, 1), vec3f(0, 0, -1));
  Check(enterInner.interface.ratioIsConstant, "scaled Cauchy coefficients prove constant ratio");
  CheckApprox(
    enterInner.interface.Eta(static_cast<Float>(470)),
    static_cast<Float>(2),
    static_cast<Float>(1e-6),
    "scaled Cauchy ratio evaluates to two"
  );

  NamedSpectrumRegistry registry = NamedSpectrumRegistry::LoadDefault();
  EtaSpectrumHandle bk7 = NamedEtaSpectrum(registry, "glass-BK7");
  Check(EvaluateEtaSpectrum(bk7, static_cast<Float>(550)) > 1, "named glass BK7 eta loads");
  Check(!EtaSpectrumIsConstant(bk7), "named glass BK7 is dispersive");
}

void TestRoughDielectricBxDF() {
  DielectricBxDF dielectric(
    static_cast<Float>(1.5),
    TrowbridgeReitzDistribution(static_cast<Float>(0.35), static_cast<Float>(0.22))
  );
  vec3f wo = Unit(vec3f(static_cast<Float>(0.22), static_cast<Float>(-0.16), static_cast<Float>(0.96)));

  std::optional<BSDFSample> reflection = dielectric.Sample_f(
    wo,
    static_cast<Float>(0),
    point2f(static_cast<Float>(0.31), static_cast<Float>(0.77)),
    TransportMode::Radiance,
    BxDFReflTransFlags::Reflection
  );
  Check(reflection.has_value(), "rough dielectric reflection samples");
  Check(reflection->IsGlossy(), "rough dielectric reflection is glossy");
  Check(reflection->IsReflection(), "rough dielectric reflection flag");
  Check(
    rayrender::render::SameHemisphere(wo, reflection->wi),
    "rough reflection stays in same hemisphere"
  );
  Check(reflection->pdf > 0, "rough reflection PDF is positive");
  CheckApprox(
    dielectric.PDF(wo, reflection->wi, TransportMode::Radiance, BxDFReflTransFlags::Reflection),
    reflection->pdf,
    static_cast<Float>(1e-4),
    "rough reflection PDF matches Sample_f"
  );
  CheckApprox(
    dielectric.f(wo, reflection->wi, TransportMode::Radiance)[0],
    reflection->f[0],
    static_cast<Float>(1e-4),
    "rough reflection f matches Sample_f"
  );

  std::optional<BSDFSample> transmission;
  const point2f candidates[] = {
    point2f(static_cast<Float>(0.11), static_cast<Float>(0.19)),
    point2f(static_cast<Float>(0.47), static_cast<Float>(0.53)),
    point2f(static_cast<Float>(0.82), static_cast<Float>(0.29)),
    point2f(static_cast<Float>(0.26), static_cast<Float>(0.91))
  };
  for (point2f u : candidates) {
    transmission = dielectric.Sample_f(
      wo,
      static_cast<Float>(1),
      u,
      TransportMode::Radiance,
      BxDFReflTransFlags::Transmission
    );
    if (transmission) {
      break;
    }
  }
  Check(transmission.has_value(), "rough dielectric transmission samples");
  Check(transmission->IsGlossy(), "rough dielectric transmission is glossy");
  Check(transmission->IsTransmission(), "rough dielectric transmission flag");
  Check(
    !rayrender::render::SameHemisphere(wo, transmission->wi),
    "rough transmission crosses hemisphere"
  );
  Check(transmission->eta > 0, "rough transmission carries eta");
  CheckApprox(
    dielectric.PDF(wo, transmission->wi, TransportMode::Radiance, BxDFReflTransFlags::Transmission),
    transmission->pdf,
    static_cast<Float>(1e-4),
    "rough transmission PDF matches Sample_f"
  );
  CheckApprox(
    dielectric.f(wo, transmission->wi, TransportMode::Radiance)[0],
    transmission->f[0],
    static_cast<Float>(1e-4),
    "rough transmission f matches Sample_f"
  );
}

void TestDispersiveMaterialTermination() {
  DielectricRegionTable table;
  RegionId glass = table.AddRegion(
    CauchyEtaSpectrum(static_cast<Float>(1.5), static_cast<Float>(0.004)),
    0,
    "cauchy"
  );
  DielectricPathState state = DielectricPathState::FromInitialRegions(&table, {});
  ResolvedDielectricTransition transition =
    state.Analyze({{glass, RegionSide::NegativeNormal}}, normal3f(0, 0, 1), vec3f(0, 0, -1));

  Material material = Material::Dielectric();
  MaterialEvalContext ctx = BaseContext();
  ctx.dielectric = &transition.interface;
  SampledWavelengths lambda = TestWavelengths();
  Check(!lambda.SecondaryTerminated(), "test wavelengths begin active");
  ScratchBuffer scratch(1024);
  CountingTextureEvaluator evaluator;
  BSDF bsdf = material.GetBSDF(evaluator, ctx, lambda, scratch);
  Check(lambda.SecondaryTerminated(), "dispersive dielectric terminates before sampling");
  std::optional<BSDFSample> reflected = bsdf.Sample_f(
    ctx.wo,
    static_cast<Float>(0),
    point2f(static_cast<Float>(0.1), static_cast<Float>(0.2)),
    TransportMode::Radiance,
    BxDFReflTransFlags::Reflection
  );
  Check(reflected.has_value(), "reflection can be sampled after termination");
  Check(lambda.SecondaryTerminated(), "reflection branch keeps pbrt early termination");

  std::optional<BSDFSample> roughTransmission = DielectricBxDF(
    transition.interface.Eta(lambda[0]),
    TrowbridgeReitzDistribution(static_cast<Float>(0.3), static_cast<Float>(0.3))
  ).Sample_f(
    ctx.wo,
    static_cast<Float>(1),
    point2f(static_cast<Float>(0.2), static_cast<Float>(0.7)),
    TransportMode::Radiance,
    BxDFReflTransFlags::Transmission
  );
  Check(roughTransmission.has_value(), "rough dielectric can produce transmission");
  Check(!state.Contains(glass), "region state is unchanged before rough transmission commit");
  state.Commit(transition.token);
  Check(state.Contains(glass), "region state commits only after transmission decision");
}

void TestThinDielectricMaterial() {
  MaterialEvalContext ctx = BaseContext();
  CountingTextureEvaluator evaluator;

  Material constantThin = Material::ThinDielectric(
    ThinDielectricMaterial(ConstantEtaSpectrum(static_cast<Float>(1.5)))
  );
  SampledWavelengths constantLambda = TestWavelengths();
  ScratchBuffer constantScratch(1024);
  BSDF constantBSDF = constantThin.GetBSDF(evaluator, ctx, constantLambda, constantScratch);
  Check(!constantLambda.SecondaryTerminated(), "constant thin dielectric keeps wavelengths");
  Check(IsSpecular(constantBSDF.Flags()), "thin dielectric is specular");
  std::optional<BSDFSample> thinTransmission = constantBSDF.Sample_f(
    ctx.wo,
    static_cast<Float>(1),
    point2f(static_cast<Float>(0.3), static_cast<Float>(0.4)),
    TransportMode::Radiance,
    BxDFReflTransFlags::Transmission
  );
  Check(thinTransmission.has_value(), "thin dielectric transmission samples");
  Check(thinTransmission->IsTransmission(), "thin dielectric transmission flag");
  CheckApprox(
    thinTransmission->eta,
    static_cast<Float>(1),
    static_cast<Float>(1e-6),
    "thin dielectric does not update etaScale"
  );

  Material dispersiveThin = Material::ThinDielectric(
    ThinDielectricMaterial(CauchyEtaSpectrum(static_cast<Float>(1.5), static_cast<Float>(0.004)))
  );
  SampledWavelengths dispersiveLambda = TestWavelengths();
  ScratchBuffer dispersiveScratch(1024);
  (void)dispersiveThin.GetBSDF(evaluator, ctx, dispersiveLambda, dispersiveScratch);
  Check(dispersiveLambda.SecondaryTerminated(), "dispersive thin dielectric terminates wavelengths");
}

void TestPathTerminationDiagnostics() {
  DielectricRegionTable regions;
  RegionId glass = regions.AddRegion(
    CauchyEtaSpectrum(static_cast<Float>(1.5), static_cast<Float>(0.004)),
    0,
    "cauchy"
  );
  SpectralTestScene dielectricScene;
  AddSphereBoundary(dielectricScene, glass, Material::Dielectric());
  PreprocessLights(dielectricScene);

  PathRenderOptions options;
  options.maxDepth = 1;
  options.sampleDirectLighting = false;
  options.sampleBSDF = false;
  PathIntegrator integrator(
    dielectricScene.scene,
    dielectricScene.materials,
    dielectricScene.lights,
    options,
    &regions
  );
  SampledWavelengths lambda = TestWavelengths();
  SpectralRandomSampler sampler(9);
  sampler.StartPixelSample({0, 0}, 0);
  ScratchBuffer scratch(4096);
  PathRenderStats stats;
  (void)integrator.Li(
    Ray(point3f(0, 0, 0), Unit(vec3f(0, 0, -1))),
    lambda,
    sampler,
    scratch,
    nullptr,
    &stats
  );
  Check(lambda.SecondaryTerminated(), "path dielectric terminated wavelengths");
  Check(stats.wavelengthTerminations == 1, "path records one wavelength termination");
  Check(
    stats.dielectricWavelengthTerminations == 1,
    "path records dielectric termination location"
  );
  Check(
    stats.thinDielectricWavelengthTerminations == 0,
    "path does not misclassify dielectric as thin"
  );

  DielectricRegionTable nullRegions;
  EtaSpectrumHandle sharedEta =
    CauchyEtaSpectrum(static_cast<Float>(1.5), static_cast<Float>(0.004));
  RegionId shell = nullRegions.AddRegion(sharedEta, 0, "shared shell");
  RegionId matched = nullRegions.AddRegion(sharedEta, -1, "shared insert");
  DielectricPathState initialNullState =
    DielectricPathState::FromInitialRegions(&nullRegions, {shell});
  SpectralTestScene nullScene;
  AddSphereBoundary(nullScene, matched, Material::Dielectric());
  PreprocessLights(nullScene);

  PathIntegrator nullIntegrator(
    nullScene.scene,
    nullScene.materials,
    nullScene.lights,
    options,
    &nullRegions
  );
  SampledWavelengths nullLambda = TestWavelengths();
  SpectralRandomSampler nullSampler(11);
  nullSampler.StartPixelSample({0, 0}, 0);
  ScratchBuffer nullScratch(4096);
  PathRenderStats nullStats;
  (void)nullIntegrator.Li(
    Ray(point3f(0, 0, 0), Unit(vec3f(0, 0, -1))),
    nullLambda,
    nullSampler,
    nullScratch,
    nullptr,
    &nullStats,
    &initialNullState
  );
  Check(!nullLambda.SecondaryTerminated(), "index-matched dispersive null keeps wavelengths");
  Check(
    nullStats.dielectricIndexMatchedSkips > 0,
    "path records dispersive index-matched null traversal"
  );
  Check(nullStats.wavelengthTerminations == 0, "null traversal records no wavelength termination");
  Check(nullStats.materialClosures == 0, "null traversal skips dielectric material evaluation");

  SpectralTestScene thinScene;
  MaterialHandle material = thinScene.materials.Add(Material::ThinDielectric(
    ThinDielectricMaterial(CauchyEtaSpectrum(static_cast<Float>(1.5), static_cast<Float>(0.004)))
  ));
  ShapeHandle plane = thinScene.scene.AddShape(InfinitePlane(static_cast<Float>(-1), normal3f(0, 0, 1)));
  PrimitiveBinding binding;
  binding.material = material;
  thinScene.scene.AddPrimitive(plane, binding);
  PreprocessLights(thinScene);

  PathIntegrator thinIntegrator(thinScene.scene, thinScene.materials, thinScene.lights, options);
  SampledWavelengths thinLambda = TestWavelengths();
  SpectralRandomSampler thinSampler(10);
  thinSampler.StartPixelSample({0, 0}, 0);
  ScratchBuffer thinScratch(4096);
  PathRenderStats thinStats;
  (void)thinIntegrator.Li(
    Ray(point3f(0, 0, 0), Unit(vec3f(0, 0, -1))),
    thinLambda,
    thinSampler,
    thinScratch,
    nullptr,
    &thinStats
  );
  Check(thinLambda.SecondaryTerminated(), "path thin dielectric terminated wavelengths");
  Check(
    thinStats.thinDielectricWavelengthTerminations == 1,
    "path records thin dielectric termination location"
  );
}

} // namespace

int main() {
  try {
    TestEtaSpectrumResolution();
    TestRoughDielectricBxDF();
    TestDispersiveMaterialTermination();
    TestThinDielectricMaterial();
    TestPathTerminationDiagnostics();
  } catch (const std::exception& ex) {
    std::cerr << "PR18 Dielectric test threw exception: " << ex.what() << std::endl;
    return 1;
  }
  std::cout << "PR18 Dielectric tests passed" << std::endl;
  return 0;
}
