#include "src/render/spectral_medium.h"

#include "src/render/spectral_integrator.h"

#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <optional>
#include <string>

using namespace rayrender::base;
using namespace rayrender::materials;
using namespace rayrender::render;

namespace {

constexpr Float Pi = static_cast<Float>(3.14159265358979323846264338327950288);
constexpr Float Inv4Pi = static_cast<Float>(1) / (static_cast<Float>(4) * Pi);

bool Approx(Float lhs, Float rhs, Float tolerance = static_cast<Float>(1e-5)) {
  return std::fabs(lhs - rhs) <= tolerance;
}

void Check(bool condition, const char* message) {
  if (!condition) {
    std::cerr << "PR20 Medium test failed: " << message << std::endl;
    std::exit(1);
  }
}

void CheckApprox(Float actual, Float expected, Float tolerance, const char* message) {
  if (!Approx(actual, expected, tolerance)) {
    std::cerr << "PR20 Medium test failed: " << message
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

void CheckSpectrumFiniteNonNegative(const SampledSpectrum& spectrum, const char* message) {
  for (int i = 0; i < NSpectrumSamples; ++i) {
    if (!std::isfinite(spectrum[i]) || spectrum[i] < 0) {
      std::cerr << "PR20 Medium test failed: " << message
                << " component " << i << " got " << spectrum[i] << std::endl;
      std::exit(1);
    }
  }
}

template <typename Fn>
void CheckThrows(Fn&& fn, const char* message) {
  try {
    fn();
  } catch (const std::exception&) {
    return;
  }
  Check(false, message);
}

vec3f Unit(vec3f value) {
  return value / value.length();
}

SampledWavelengths TestWavelengths() {
  return SampledWavelengths::SampleUniform(static_cast<Float>(0.25));
}

Spectrum Constant(Float value) {
  return Spectrum(ConstantSpectrum(value));
}

Spectrum LinearSpectrum(Float startValue, Float endValue) {
  SpectrumDataPolicy policy;
  policy.order = SpectrumOrderPolicy::RequireSorted;
  policy.extrapolation = SpectrumExtrapolationPolicy::Constant;
  policy.validation = SpectrumValueValidation::NonNegative;
  return Spectrum(PiecewiseLinearSpectrum(
    {LambdaMin, LambdaMax},
    {startValue, endValue},
    policy
  ));
}

void TestHenyeyGreensteinPhaseFunction() {
  HenyeyGreensteinPhaseFunction isotropic(0);
  vec3f wo = vec3f(0, 0, 1);
  vec3f wi = Unit(vec3f(static_cast<Float>(0.2), static_cast<Float>(0.3), static_cast<Float>(0.93)));

  CheckApprox(
    isotropic.p(wo, wi),
    Inv4Pi,
    static_cast<Float>(1e-6),
    "isotropic HG evaluates to 1 / 4pi"
  );
  CheckApprox(
    isotropic.PDF(wo, wi),
    isotropic.p(wo, wi),
    static_cast<Float>(1e-7),
    "HG PDF matches phase value"
  );

  std::optional<PhaseFunctionSample> sample =
    isotropic.Sample_p(wo, point2f(static_cast<Float>(0.2), static_cast<Float>(0.37)));
  Check(sample.has_value(), "HG sampling returns a sample");
  CheckApprox(sample->wi.length(), static_cast<Float>(1), static_cast<Float>(1e-5), "HG sampled direction is normalized");
  CheckApprox(sample->pdf, sample->p, static_cast<Float>(1e-7), "HG sample pdf equals p");
  CheckApprox(sample->pdf, isotropic.PDF(wo, sample->wi), static_cast<Float>(1e-6), "HG sampled PDF is self-consistent");

  HenyeyGreensteinPhaseFunction forward(static_cast<Float>(0.45));
  Check(forward.p(wo, -wo) > forward.p(wo, wo), "positive g follows pbrt wo/wi convention");
  CheckThrows([]() { HenyeyGreensteinPhaseFunction invalid(1); }, "HG rejects g at one");
  CheckThrows([]() { HenyeyGreensteinPhaseFunction invalid(-1); }, "HG rejects g at negative one");
}

void TestPhaseFunctionTable() {
  PhaseFunctionTable table;
  PhaseFunctionHandle handle = table.Add(HenyeyGreensteinPhaseFunction(static_cast<Float>(0.2)));
  Check(table.Size() == 1, "phase table records size");
  Check(table.IsValid(handle), "phase table accepts live handle");
  PhaseFunction phase = table.Get(handle);
  Check(static_cast<bool>(phase), "phase table returns dispatch wrapper");
  CheckApprox(
    phase.PDF(vec3f(0, 0, 1), vec3f(0, 0, 1)),
    table.GetHenyeyGreenstein(handle).PDF(vec3f(0, 0, 1), vec3f(0, 0, 1)),
    static_cast<Float>(1e-7),
    "phase dispatch forwards PDF"
  );

  PhaseFunctionHandle stale =
    PhaseFunctionHandle::FromIndex(handle.Index(), handle.Generation() + 1);
  Check(!table.IsValid(stale), "phase table rejects stale generation");
  CheckThrows([&]() { (void)table.Get(stale); }, "phase table throws on stale handle");
}

void TestHomogeneousBeerLambert() {
  HomogeneousMedium medium(
    Constant(static_cast<Float>(0.25)),
    Constant(static_cast<Float>(0.15)),
    static_cast<Float>(2)
  );
  SampledWavelengths lambda = TestWavelengths();
  Ray ray(point3f(0, 0, 0), vec3f(0, 0, 2));
  SampledSpectrum tr = medium.Tr(ray, static_cast<Float>(3), lambda);
  Float expected = std::exp(-static_cast<Float>(0.8) * static_cast<Float>(6));
  CheckSpectrumApprox(
    tr,
    SampledSpectrum(expected),
    static_cast<Float>(2e-6),
    "homogeneous Beer-Lambert transmittance uses sigma_t, scale, and ray length"
  );
}

void TestHomogeneousSpectralPropertiesAndMajorant() {
  HomogeneousMedium medium(
    LinearSpectrum(static_cast<Float>(0.1), static_cast<Float>(0.2)),
    Constant(static_cast<Float>(0.05)),
    static_cast<Float>(3),
    Constant(static_cast<Float>(4)),
    static_cast<Float>(0.5),
    static_cast<Float>(0.1)
  );
  SampledWavelengths lambda = TestWavelengths();
  MediumProperties properties = medium.SamplePoint(point3f(1, 2, 3), lambda);
  Spectrum sigmaA = LinearSpectrum(static_cast<Float>(0.1), static_cast<Float>(0.2));

  for (int i = 0; i < NSpectrumSamples; ++i) {
    CheckApprox(
      properties.sigmaA[i],
      sigmaA(lambda[i]) * static_cast<Float>(3),
      static_cast<Float>(1e-6),
      "homogeneous medium samples spectral sigma_a"
    );
    CheckApprox(
      properties.sigmaS[i],
      static_cast<Float>(0.15),
      static_cast<Float>(1e-6),
      "homogeneous medium scales sigma_s"
    );
    CheckApprox(
      properties.emission[i],
      static_cast<Float>(2),
      static_cast<Float>(1e-6),
      "homogeneous medium scales emission"
    );
  }
  Check(properties.HasScattering(), "homogeneous medium reports scattering phase");
  Check(medium.IsEmissive(), "homogeneous medium reports emission");

  RayMajorantIterator iterator =
    medium.SampleRay(Ray(point3f(0, 0, 0), vec3f(0, 0, -1)), static_cast<Float>(2), lambda);
  std::optional<RayMajorantSegment> segment = iterator.Next();
  Check(segment.has_value(), "homogeneous majorant iterator returns one segment");
  CheckApprox(segment->tMin, 0, static_cast<Float>(1e-7), "majorant segment starts at zero");
  CheckApprox(segment->tMax, static_cast<Float>(2), static_cast<Float>(1e-7), "majorant segment preserves tMax");
  CheckSpectrumApprox(segment->sigmaMaj, properties.SigmaT(), static_cast<Float>(1e-6), "majorant is sigma_t");
  Check(!iterator.Next().has_value(), "homogeneous majorant iterator is exhausted after one segment");

  RayMajorantIterator scaledIterator =
    medium.SampleRay(Ray(point3f(0, 0, 0), vec3f(0, 0, -2)), static_cast<Float>(2), lambda);
  std::optional<RayMajorantSegment> scaledSegment = scaledIterator.Next();
  Check(scaledSegment.has_value(), "non-unit ray direction still returns a majorant segment");
  CheckSpectrumApprox(
    scaledSegment->sigmaMaj,
    properties.SigmaT() * static_cast<Float>(2),
    static_cast<Float>(1e-6),
    "majorant scales with non-unit ray direction"
  );
}

void TestHomogeneousDistanceSampling() {
  HomogeneousMedium medium(
    Constant(static_cast<Float>(0.1)),
    Constant(static_cast<Float>(0.3))
  );
  SampledWavelengths lambda = TestWavelengths();
  Ray ray(
    point3f(0, 0, 0),
    vec3f(0, 0, -1),
    static_cast<Float>(0),
    static_cast<Float>(10)
  );
  std::optional<MediumSample> sample = medium.SampleDistance(
    ray,
    static_cast<Float>(10),
    lambda,
    static_cast<Float>(0.1),
    static_cast<Float>(0.2)
  );
  Check(sample.has_value(), "homogeneous medium distance sampling returns a sample");
  Float expectedDistance = -std::log(static_cast<Float>(0.8)) / static_cast<Float>(0.4);
  Check(sample->sampledMedium, "distance sample lands inside scattering medium");
  CheckApprox(sample->t, expectedDistance, static_cast<Float>(1e-5), "sampled medium t matches exponential distance");
  Check(sample->pdf > 0 && std::isfinite(sample->pdf), "medium sample PDF is valid");
  Check(static_cast<bool>(sample->phase), "medium sample carries phase function");
  CheckSpectrumFiniteNonNegative(sample->beta, "medium sample beta is finite non-negative");

  HomogeneousMedium vacuum(Constant(0), Constant(0));
  std::optional<MediumSample> escaped = vacuum.SampleDistance(
    ray,
    static_cast<Float>(3),
    lambda,
    static_cast<Float>(0.4),
    static_cast<Float>(0.7)
  );
  Check(escaped.has_value(), "vacuum distance sampling returns escaped segment");
  Check(!escaped->sampledMedium, "vacuum distance sampling does not create a medium event");
  CheckSpectrumApprox(escaped->beta, SampledSpectrum(1), static_cast<Float>(1e-6), "vacuum escaped beta is one");

  HomogeneousMedium absorbing(Constant(static_cast<Float>(0.5)), Constant(0));
  MediumProperties absorbingProperties = absorbing.SamplePoint(point3f(0, 0, 0), lambda);
  Check(!absorbingProperties.HasScattering(), "pure absorption medium has no scattering phase event");
}

void TestMediumTable() {
  MediumTable table;
  MediumHandle handle = table.Add(HomogeneousMedium(
    Constant(static_cast<Float>(0.2)),
    Constant(static_cast<Float>(0.1))
  ));
  Check(table.Size() == 1, "medium table records size");
  Check(table.IsValid(handle), "medium table accepts live handle");
  Medium medium = table.Get(handle);
  Check(static_cast<bool>(medium), "medium table returns dispatch wrapper");

  SampledWavelengths lambda = TestWavelengths();
  Ray ray(point3f(0, 0, 0), vec3f(0, 0, -1));
  CheckSpectrumApprox(
    medium.Tr(ray, static_cast<Float>(2), lambda),
    table.GetHomogeneous(handle).Tr(ray, static_cast<Float>(2), lambda),
    static_cast<Float>(1e-7),
    "medium dispatch forwards transmittance"
  );

  MediumHandle stale = MediumHandle::FromIndex(handle.Index(), handle.Generation() + 1);
  Check(!table.IsValid(stale), "medium table rejects stale generation");
  CheckThrows([&]() { (void)table.Get(stale); }, "medium table throws on stale handle");
}

void TestDielectricRegionMediumBinding() {
  MediumTable media;
  MediumHandle absorbing = media.Add(HomogeneousMedium(
    Constant(static_cast<Float>(0.2)),
    Constant(0)
  ));
  MediumHandle scattering = media.Add(HomogeneousMedium(
    Constant(static_cast<Float>(0.05)),
    Constant(static_cast<Float>(0.15))
  ));

  DielectricRegionTable regions;
  RegionId glass = regions.AddRegion(static_cast<Float>(1), 5, "absorbing glass", absorbing);
  RegionId cloud = regions.AddRegion(static_cast<Float>(1), 1, "scattering inclusion", scattering);
  Check(!regions.VacuumOnly(), "region table with media is not vacuum-only");

  DielectricPathState state(&regions);
  std::vector<RegionBoundaryAttachment> boundary = {
    {glass, RegionSide::NegativeNormal}
  };
  ResolvedDielectricTransition enter = state.Analyze(
    boundary,
    normal3f(0, 0, 1),
    vec3f(0, 0, -1)
  );
  Check(!state.ActiveMedium().IsValid(), "active-before medium remains exterior before commit");
  Check(enter.activeBefore == ExteriorRegionId, "enter transition active-before is exterior");
  Check(enter.activeAfter == glass, "enter transition active-after is glass");
  state.Commit(enter.token);
  Check(state.ActiveMedium() == absorbing, "committed glass transition activates region medium");

  DielectricPathState nested = DielectricPathState::FromInitialRegions(&regions, {glass, cloud});
  Check(nested.ActiveMedium() == scattering, "higher-priority nested region supplies active medium");

  ResolvedDielectricTransition exit = state.Analyze(
    boundary,
    normal3f(0, 0, 1),
    vec3f(0, 0, 1)
  );
  Check(state.ActiveMedium() == absorbing, "exit transition keeps active-before medium before commit");
  Check(exit.activeBefore == glass, "exit transition active-before is glass");
  Check(exit.activeAfter == ExteriorRegionId, "exit transition active-after is exterior");
  state.Commit(exit.token);
  Check(!state.ActiveMedium().IsValid(), "committed exit clears active medium");
}

void TestSurfaceIntegratorsRejectActiveMedium() {
  MediumTable media;
  MediumHandle fog = media.Add(HomogeneousMedium(
    Constant(static_cast<Float>(0.1)),
    Constant(static_cast<Float>(0.2))
  ));
  DielectricRegionTable regions;
  RegionId fogRegion = regions.AddRegion(static_cast<Float>(1), 1, "fog", fog);
  DielectricPathState initial =
    DielectricPathState::FromInitialRegions(&regions, {fogRegion});

  Scene scene;
  SpectralMaterialTable materials;
  SpectralLightTable lights;
  Ray ray(point3f(0, 0, 0), vec3f(0, 0, -1));
  SampledWavelengths lambda = TestWavelengths();

  RandomWalkRenderOptions randomOptions;
  RandomWalkIntegrator randomWalk(scene, materials, lights, randomOptions, &regions);
  SpectralRandomSampler randomSampler(7);
  randomSampler.StartPixelSample({0, 0}, 0);
  ScratchBuffer randomScratch(randomOptions.scratchBufferBytes);
  RandomWalkRenderStats randomStats;
  (void)randomWalk.Li(
    ray,
    lambda,
    randomSampler,
    randomScratch,
    nullptr,
    &randomStats,
    &initial
  );
  Check(randomStats.unsupportedMediumInteractions == 1, "RandomWalk rejects active region medium");
  Check(randomStats.raysTraced == 0, "RandomWalk rejects active medium before tracing a ray");

  PathRenderOptions pathOptions;
  PathIntegrator path(scene, materials, lights, pathOptions, &regions);
  SpectralRandomSampler pathSampler(11);
  pathSampler.StartPixelSample({0, 0}, 0);
  ScratchBuffer pathScratch(pathOptions.scratchBufferBytes);
  PathRenderStats pathStats;
  (void)path.Li(
    ray,
    lambda,
    pathSampler,
    pathScratch,
    nullptr,
    &pathStats,
    &initial
  );
  Check(pathStats.unsupportedMediumInteractions == 1, "Path rejects active region medium");
  Check(pathStats.raysTraced == 0, "Path rejects active medium before tracing a ray");
}

} // namespace

int main() {
  TestHenyeyGreensteinPhaseFunction();
  TestPhaseFunctionTable();
  TestHomogeneousBeerLambert();
  TestHomogeneousSpectralPropertiesAndMajorant();
  TestHomogeneousDistanceSampling();
  TestMediumTable();
  TestDielectricRegionMediumBinding();
  TestSurfaceIntegratorsRejectActiveMedium();
  return 0;
}
