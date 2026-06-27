#include "src/render/spectral_integrator.h"

#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <optional>
#include <string>
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
    std::cerr << "PR17 Dielectric test failed: " << message << std::endl;
    std::exit(1);
  }
}

void CheckApprox(Float actual, Float expected, Float tolerance, const char* message) {
  if (!Approx(actual, expected, tolerance)) {
    std::cerr << "PR17 Dielectric test failed: " << message
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

vec3f Unit(vec3f value) {
  return value / value.length();
}

SampledWavelengths TestWavelengths() {
  return SampledWavelengths::SampleUniform(static_cast<Float>(0.41));
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

template <typename Callback>
void ExpectThrows(Callback callback, const char* message) {
  bool threw = false;
  try {
    callback();
  } catch (const std::exception&) {
    threw = true;
  }
  Check(threw, message);
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

void AddDielectricSphereBoundary(
  SpectralTestScene& testScene,
  RegionId region
) {
  MaterialHandle material = testScene.materials.Add(Material::Dielectric());
  ShapeHandle shape =
    testScene.scene.AddShape(Shape::Sphere(static_cast<Float>(1), point3f(0, 0, -3)));
  PrimitiveBinding binding;
  binding.material = material;
  binding.dielectricBoundaries.push_back({region, RegionSide::NegativeNormal});
  testScene.scene.AddPrimitive(shape, binding);
}

void PreprocessLights(SpectralTestScene& testScene) {
  Bounds3f bounds;
  testScene.scene.Bounds(&bounds);
  testScene.lights.Preprocess(bounds);
}

void TestDielectricStateMachine() {
  DielectricRegionTable table;
  RegionId glass = table.AddRegion(static_cast<Float>(1.5), 0, "glass");
  RegionId matched = table.AddRegion(static_cast<Float>(1.5), -1, "matched");
  RegionId liquid = table.AddRegion(static_cast<Float>(1.33), -2, "liquid");
  RegionId ignoredShell = table.AddRegion(static_cast<Float>(1.1), 10, "ignored shell");

  DielectricPathState state = DielectricPathState::FromInitialRegions(&table, {});
  ResolvedDielectricTransition enterGlass = state.Analyze(
    {{glass, RegionSide::NegativeNormal}},
    normal3f(0, 0, 1),
    vec3f(0, 0, -1)
  );
  Check(enterGlass.kind == DielectricBoundaryKind::ScatteringInterface, "air to glass scatters");
  Check(enterGlass.activeBefore == ExteriorRegionId, "glass entry starts in exterior");
  Check(enterGlass.activeAfter == glass, "glass entry activates glass");
  CheckApprox(enterGlass.interface.Eta(), static_cast<Float>(1.5), static_cast<Float>(1e-6), "glass interface eta");
  Check(!state.Contains(glass), "reflection-side state remains exterior before commit");

  DielectricPathState copy = state;
  copy.Commit(enterGlass.token);
  Check(copy.Contains(glass), "visibility-style copy can commit independently");
  Check(!state.Contains(glass), "copy commit leaves original state unchanged");

  state.Commit(enterGlass.token);
  Check(state.ActiveRegionId() == glass, "transmission commit enters glass");
  ExpectThrows([&state, &enterGlass]() { state.Commit(enterGlass.token); }, "stale transition token is rejected");

  ResolvedDielectricTransition enterMatched = state.Analyze(
    {{matched, RegionSide::NegativeNormal}},
    normal3f(0, 0, 1),
    vec3f(0, 0, -1)
  );
  Check(enterMatched.kind == DielectricBoundaryKind::IndexMatchedNull, "eta-matched active change is null traversal");
  state.Commit(enterMatched.token);
  Check(state.ActiveRegionId() == matched, "index-matched traversal still updates membership");

  ResolvedDielectricTransition enterIgnoredShell = state.Analyze(
    {{ignoredShell, RegionSide::NegativeNormal}},
    normal3f(0, 0, 1),
    vec3f(0, 0, -1)
  );
  Check(enterIgnoredShell.kind == DielectricBoundaryKind::PrioritySkipped, "lower-priority overlap is skipped");
  state.Commit(enterIgnoredShell.token);
  Check(state.Contains(ignoredShell), "priority-skipped boundary still commits membership");
  Check(state.ActiveRegionId() == matched, "priority-skipped boundary does not change active region");

  ResolvedDielectricTransition enterLiquid = state.Analyze(
    {{liquid, RegionSide::NegativeNormal}},
    normal3f(0, 0, 1),
    vec3f(0, 0, -1)
  );
  Check(enterLiquid.kind == DielectricBoundaryKind::ScatteringInterface, "higher-priority nested liquid scatters");
  state.Commit(enterLiquid.token);
  Check(state.ActiveRegionId() == liquid, "nested lower numeric priority wins");

  ResolvedDielectricTransition exitLiquid = state.Analyze(
    {{liquid, RegionSide::NegativeNormal}},
    normal3f(0, 0, -1),
    vec3f(0, 0, -1)
  );
  state.Commit(exitLiquid.token);
  Check(state.ActiveRegionId() == matched, "leaving nested liquid reveals previous active region");

  DielectricRegionTable equalPriority;
  RegionId a = equalPriority.AddRegion(static_cast<Float>(1.2), 0, "a");
  RegionId b = equalPriority.AddRegion(static_cast<Float>(1.4), 0, "b");
  DielectricPathState equalState = DielectricPathState::FromInitialRegions(&equalPriority, {});
  bool equalPriorityFailed = false;
  try {
    (void)equalState.Analyze(
      {{a, RegionSide::NegativeNormal}, {b, RegionSide::NegativeNormal}},
      normal3f(0, 0, 1),
      vec3f(0, 0, -1)
    );
  } catch (const std::exception&) {
    equalPriorityFailed = true;
  }
  Check(equalPriorityFailed, "equal-priority active overlaps are invalid");

  DielectricPathState repeated = DielectricPathState::FromInitialRegions(&table, {glass, glass});
  Check(repeated.Members().size() == 1, "repeated RegionId memberships are stored once");
}

void TestPointContainmentInitialization() {
  DielectricRegionTable table;
  RegionId glass = table.AddRegion(static_cast<Float>(1.5), 0, "glass");
  Scene scene;
  ShapeHandle shape = scene.AddShape(Shape::Sphere(static_cast<Float>(1), point3f(0, 0, 0)));
  PrimitiveBinding binding;
  binding.dielectricBoundaries.push_back({glass, RegionSide::NegativeNormal});
  scene.AddPrimitive(shape, binding);

  DielectricPathState inside =
    DielectricPathState::FromPointContainment(&table, scene, point3f(0, 0, 0));
  Check(inside.Contains(glass), "camera inside sphere initializes contained region");
  Check(inside.ActiveRegionId() == glass, "camera inside sphere initializes active glass");

  DielectricPathState outside =
    DielectricPathState::FromPointContainment(&table, scene, point3f(2, 0, 0));
  Check(!outside.Contains(glass), "camera outside sphere initializes exterior");
  Check(outside.ActiveRegionId() == ExteriorRegionId, "outside camera active region is exterior");

  Scene csgScene;
  Shape csg = Shape::CSGSignedDistance(
    [](const point3f& p) {
      return (p - point3f(0, 0, 0)).length() - static_cast<Float>(1);
    },
    CSGShapeOptions{Bounds3f(point3f(-2, -2, -2), point3f(2, 2, 2))},
    "CSG sphere"
  );
  ShapeHandle csgShape = csgScene.AddShape(std::move(csg));
  PrimitiveBinding csgBinding;
  csgBinding.dielectricBoundaries.push_back({glass, RegionSide::Inside});
  csgScene.AddPrimitive(csgShape, csgBinding);
  DielectricPathState csgInside =
    DielectricPathState::FromPointContainment(&table, csgScene, point3f(0, 0, 0));
  Check(csgInside.Contains(glass), "CSG inside side initializes contained region");
}

void TestSmoothDielectricBxDF() {
  DielectricBxDF dielectric(static_cast<Float>(1.5), TrowbridgeReitzDistribution(0, 0));
  BxDFFlags flags = dielectric.Flags();
  Check(IsReflective(flags), "smooth dielectric is reflective");
  Check(IsTransmissive(flags), "smooth dielectric is transmissive");
  Check(IsSpecular(flags), "smooth dielectric is specular");

  vec3f wo(0, 0, 1);
  Float R = rayrender::render::FrDielectric(1, static_cast<Float>(1.5));
  Float T = static_cast<Float>(1) - R;
  std::optional<BSDFSample> reflected = dielectric.Sample_f(
    wo,
    0,
    point2f(static_cast<Float>(0.25), static_cast<Float>(0.75)),
    TransportMode::Radiance
  );
  Check(reflected.has_value(), "smooth dielectric samples reflection");
  Check(reflected->IsReflection() && reflected->IsSpecular(), "reflection sample flags");
  CheckApprox(reflected->wi.xyz.z, static_cast<Float>(1), static_cast<Float>(1e-6), "reflection direction");
  CheckApprox(reflected->pdf, R, static_cast<Float>(1e-6), "reflection probability");
  CheckSpectrumApprox(reflected->f, SampledSpectrum(R), static_cast<Float>(1e-6), "reflection f at normal incidence");

  std::optional<BSDFSample> transmitted = dielectric.Sample_f(
    wo,
    static_cast<Float>(0.99),
    point2f(static_cast<Float>(0.25), static_cast<Float>(0.75)),
    TransportMode::Radiance
  );
  Check(transmitted.has_value(), "smooth dielectric samples transmission");
  Check(transmitted->IsTransmission() && transmitted->IsSpecular(), "transmission sample flags");
  CheckApprox(transmitted->wi.xyz.z, static_cast<Float>(-1), static_cast<Float>(1e-6), "transmission direction");
  CheckApprox(transmitted->pdf, T, static_cast<Float>(1e-6), "transmission probability");
  CheckApprox(transmitted->eta, static_cast<Float>(1.5), static_cast<Float>(1e-6), "transmission eta");
  CheckSpectrumApprox(
    transmitted->f,
    SampledSpectrum(T / (static_cast<Float>(1.5) * static_cast<Float>(1.5))),
    static_cast<Float>(1e-6),
    "radiance transmission divides by eta squared"
  );

  std::optional<BSDFSample> transmissionOnly = dielectric.Sample_f(
    wo,
    0,
    point2f(static_cast<Float>(0.25), static_cast<Float>(0.75)),
    TransportMode::Radiance,
    BxDFReflTransFlags::Transmission
  );
  Check(transmissionOnly.has_value() && transmissionOnly->IsTransmission(), "transmission-only sampling respects flags");
  CheckApprox(transmissionOnly->pdf, static_cast<Float>(1), static_cast<Float>(1e-6), "transmission-only delta PDF");

  vec3f insideWo = Unit(vec3f(static_cast<Float>(0.9), 0, static_cast<Float>(-0.4358899)));
  std::optional<BSDFSample> tir = dielectric.Sample_f(
    insideWo,
    static_cast<Float>(0.5),
    point2f(static_cast<Float>(0.1), static_cast<Float>(0.2)),
    TransportMode::Radiance
  );
  Check(tir.has_value(), "TIR samples reflection");
  Check(tir->IsReflection(), "TIR is reflective");
  CheckApprox(tir->pdf, static_cast<Float>(1), static_cast<Float>(1e-6), "TIR reflection probability is one");
}

void TestDielectricMaterial() {
  Material material = Material::Dielectric();
  MaterialEvalContext ctx = BaseContext();
  ResolvedDielectricInterface resolved;
  resolved.outsideRegionId = 0;
  resolved.insideRegionId = 1;
  resolved.etaOutside = 1;
  resolved.etaInside = static_cast<Float>(1.5);
  resolved.ratioIsConstant = true;
  ctx.dielectric = &resolved;

  SampledWavelengths lambda = TestWavelengths();
  CountingTextureEvaluator evaluator;
  ScratchBuffer scratch(1024);
  BSDF bsdf = material.GetBSDF(evaluator, ctx, lambda, scratch);
  Check(IsReflective(bsdf.Flags()), "dielectric material creates reflective BSDF");
  Check(IsTransmissive(bsdf.Flags()), "dielectric material creates transmissive BSDF");
  Check(IsSpecular(bsdf.Flags()), "dielectric material creates specular BSDF");
  Check(evaluator.floatCount == 2, "dielectric material evaluates u/v roughness textures");
  Check(!lambda.SecondaryTerminated(), "constant-ratio dielectric keeps wavelength packet active");

  bool missingInterfaceFailed = false;
  try {
    MaterialEvalContext missing = BaseContext();
    SampledWavelengths missingLambda = TestWavelengths();
    ScratchBuffer missingScratch(1024);
    (void)material.GetBSDF(evaluator, missing, missingLambda, missingScratch);
  } catch (const std::exception&) {
    missingInterfaceFailed = true;
  }
  Check(missingInterfaceFailed, "dielectric material requires resolved interface");

  bool roughFailed = false;
  try {
    Material rough = Material::Dielectric(DielectricMaterial::FromRoughness(
      FloatTexture::Constant(static_cast<Float>(0.2)),
      false
    ));
    SampledWavelengths roughLambda = TestWavelengths();
    ScratchBuffer roughScratch(1024);
    (void)rough.GetBSDF(evaluator, ctx, roughLambda, roughScratch);
  } catch (const std::exception&) {
    roughFailed = true;
  }
  Check(roughFailed, "rough dielectric is deferred");
}

void TestIndexMatchedPathTraversal() {
  DielectricRegionTable table;
  RegionId matchedAir = table.AddRegion(static_cast<Float>(1), -1, "matched air");
  SpectralTestScene testScene;
  AddDielectricSphereBoundary(testScene, matchedAir);
  RegisterLight(
    testScene.scene,
    testScene.lights,
    Light::UniformInfinite(LightSpectrum::Constant(static_cast<Float>(1)))
  );
  PreprocessLights(testScene);

  PathRenderOptions options;
  options.maxDepth = 1;
  options.seed = 17;
  options.sampleDirectLighting = false;
  options.russianRoulette = false;
  PathIntegrator integrator(testScene.scene, testScene.materials, testScene.lights, options, &table);
  SpectralRandomSampler sampler(options.seed);
  sampler.StartPixelSample({0, 0}, 0);
  ScratchBuffer scratch(options.scratchBufferBytes);
  SampledWavelengths lambda = TestWavelengths();
  PathRenderStats stats;
  SampledSpectrum L = integrator.Li(
    Ray(point3f(0, 0, 0), vec3f(0, 0, -1)),
    lambda,
    sampler,
    scratch,
    nullptr,
    &stats
  );

  CheckSpectrumApprox(L, SampledSpectrum(1), static_cast<Float>(1e-6), "index-matched region traversal reaches infinite light");
  Check(stats.dielectricIndexMatchedSkips == 2, "entry and exit are index-matched null traversals");
  Check(stats.dielectricScatteringInterfaces == 0, "index-matched traversal does not scatter");
  Check(stats.bsdfSamples == 0, "index-matched traversal does not sample a BSDF");
  Check(stats.materialClosures == 0, "index-matched traversal does not construct material closures");
  Check(stats.infiniteLightHits == 1, "path exits to one infinite light");
}

} // namespace

int main() {
  try {
    TestDielectricStateMachine();
    TestPointContainmentInitialization();
    TestSmoothDielectricBxDF();
    TestDielectricMaterial();
    TestIndexMatchedPathTraversal();
  } catch (const std::exception& error) {
    std::cerr << "PR17 Dielectric test failed with exception: " << error.what() << std::endl;
    return 1;
  }

  std::cout << "PR17 Dielectric tests passed" << std::endl;
  return 0;
}
