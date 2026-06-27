#include "spectral_integrator.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace rayrender {
namespace render {
namespace {

constexpr Float RayEpsilon = static_cast<Float>(1e-4);
constexpr std::uint64_t SplitMixIncrement = 0x9e3779b97f4a7c15ULL;

Float AbsDot(const vec3f& v, const vec3f& w) {
  return std::abs(dot(v, w));
}

bool IsValidPositivePDF(Float pdf) {
  return pdf > 0 && std::isfinite(pdf);
}

bool IsNullBSDF(const BSDF& bsdf) {
  return !base::HasAny(bsdf.Flags(), base::BxDFFlags::All);
}

bool HasNegativeComponent(const base::SampledSpectrum& spectrum) {
  for (int i = 0; i < base::NSpectrumSamples; ++i) {
    if (spectrum[i] < 0) {
      return true;
    }
  }
  return false;
}

bool IsFiniteNonNegative(const base::SampledSpectrum& spectrum) {
  return spectrum.IsFinite() && !HasNegativeComponent(spectrum);
}

void AddContribution(
  base::SampledSpectrum& target,
  const base::SampledSpectrum& contribution,
  RandomWalkRenderStats* stats
) {
  if (!RecordSpectralRadianceDiagnostics(contribution, stats)) {
    if (stats) {
      ++stats->invalidSamples;
    }
    return;
  }
  target += contribution;
}

void SetFirstVisibleSurface(
  const SurfaceInteraction& interaction,
  SpectralVisibleSurface* visibleSurface
) {
  if (visibleSurface == nullptr || visibleSurface->valid) {
    return;
  }
  visibleSurface->valid = true;
  visibleSurface->depth = interaction.tHit;
  visibleSurface->primitiveId = interaction.primitive.IsValid()
                                  ? static_cast<int>(interaction.primitive.Index())
                                  : -1;
}

void ValidateOptions(const RandomWalkRenderOptions& options) {
  if (options.pixelSamples <= 0) {
    throw std::invalid_argument("RandomWalkRenderOptions pixelSamples must be positive");
  }
  if (options.maxDepth < 0) {
    throw std::invalid_argument("RandomWalkRenderOptions maxDepth must be non-negative");
  }
  if (options.scratchBufferBytes == 0) {
    throw std::invalid_argument("RandomWalkRenderOptions scratchBufferBytes must be positive");
  }
  if (options.maxAlphaSkips < 0) {
    throw std::invalid_argument("RandomWalkRenderOptions maxAlphaSkips must be non-negative");
  }
}

void ValidateOptions(const PathRenderOptions& options) {
  if (options.pixelSamples <= 0) {
    throw std::invalid_argument("PathRenderOptions pixelSamples must be positive");
  }
  if (options.maxDepth < 0) {
    throw std::invalid_argument("PathRenderOptions maxDepth must be non-negative");
  }
  if (options.scratchBufferBytes == 0) {
    throw std::invalid_argument("PathRenderOptions scratchBufferBytes must be positive");
  }
  if (options.maxNullSkips < 0) {
    throw std::invalid_argument("PathRenderOptions maxNullSkips must be non-negative");
  }
}

std::vector<base::LightHandle> SceneLightHandles(const Scene& scene) {
  std::vector<base::LightHandle> lights;
  lights.reserve(scene.Lights().size() + scene.InfiniteLights().size());
  lights.insert(lights.end(), scene.Lights().begin(), scene.Lights().end());
  lights.insert(lights.end(), scene.InfiniteLights().begin(), scene.InfiniteLights().end());
  return lights;
}

std::vector<RegionBoundaryAttachment> DielectricBoundaries(
  const SurfaceInteraction& interaction
) {
  std::vector<RegionBoundaryAttachment> boundaries = interaction.dielectricBoundaries;
  if (interaction.hasDielectricRegion) {
    boundaries.push_back({
      static_cast<RegionId>(interaction.dielectricRegionId),
      RegionSide::NegativeNormal
    });
  }
  return boundaries;
}

bool HasDielectricBoundary(const SurfaceInteraction& interaction) {
  return interaction.hasDielectricRegion || !interaction.dielectricBoundaries.empty();
}

bool HasUnsupportedMedium(const SurfaceInteraction& interaction) {
  return interaction.hasMedium || interaction.hasMediumInterface;
}

void RecordDielectricKind(
  DielectricBoundaryKind kind,
  RandomWalkRenderStats* stats
) {
  if (stats == nullptr) {
    return;
  }
  switch (kind) {
  case DielectricBoundaryKind::PrioritySkipped:
    ++stats->dielectricPrioritySkips;
    return;
  case DielectricBoundaryKind::IndexMatchedNull:
    ++stats->dielectricIndexMatchedSkips;
    return;
  case DielectricBoundaryKind::ScatteringInterface:
    ++stats->dielectricScatteringInterfaces;
    return;
  }
}

void RecordDielectricKind(
  DielectricBoundaryKind kind,
  PathRenderStats* stats
) {
  if (stats == nullptr) {
    return;
  }
  switch (kind) {
  case DielectricBoundaryKind::PrioritySkipped:
    ++stats->dielectricPrioritySkips;
    return;
  case DielectricBoundaryKind::IndexMatchedNull:
    ++stats->dielectricIndexMatchedSkips;
    return;
  case DielectricBoundaryKind::ScatteringInterface:
    ++stats->dielectricScatteringInterfaces;
    return;
  }
}

template <typename Stats>
std::optional<ResolvedDielectricTransition> AnalyzeDielectricTransition(
  const DielectricRegionTable* regionTable,
  const DielectricPathState& regions,
  const SurfaceInteraction& interaction,
  const Ray& ray,
  Stats* stats
) {
  if (!HasDielectricBoundary(interaction)) {
    return std::nullopt;
  }
  if (regionTable == nullptr) {
    if (stats) {
      ++stats->unsupportedMediumInteractions;
    }
    return std::nullopt;
  }
  try {
    std::vector<RegionBoundaryAttachment> boundaries = DielectricBoundaries(interaction);
    ResolvedDielectricTransition transition =
      regions.Analyze(boundaries, interaction.n, ray.direction());
    RecordDielectricKind(transition.kind, stats);
    return transition;
  } catch (const std::exception&) {
    if (stats) {
      ++stats->dielectricStateErrors;
    }
    return std::nullopt;
  }
}

template <typename Stats>
bool CommitDielectricTransition(
  DielectricPathState& regions,
  const DielectricTransitionToken& token,
  Stats* stats
) {
  try {
    regions.Commit(token);
    return true;
  } catch (const std::exception&) {
    if (stats) {
      ++stats->dielectricStateErrors;
    }
    return false;
  }
}

DielectricPathState InitialDielectricState(
  const DielectricRegionTable* regionTable,
  const Scene& scene,
  const Ray& ray,
  const DielectricPathState* explicitState
) {
  if (explicitState != nullptr) {
    return *explicitState;
  }
  if (regionTable == nullptr) {
    return DielectricPathState();
  }
  return DielectricPathState::FromPointContainment(regionTable, scene, ray.origin());
}

bool RegionAwareUnoccluded(
  const VisibilityTester& visibility,
  const Scene& scene,
  const DielectricRegionTable* regionTable,
  DielectricPathState regions,
  PathRenderStats* stats
) {
  if (regionTable == nullptr) {
    return visibility.Unoccluded(scene.GetAggregate());
  }

  SpawnedRay spawned = visibility.SpawnRay();
  Ray ray = spawned.ray;
  for (int i = 0; i < 64; ++i) {
    std::optional<PrimitiveIntersection> hit =
      scene.Intersect(ray, RayEpsilon, ray.tMax);
    if (!hit) {
      return true;
    }
    const SurfaceInteraction& interaction = hit->interaction;
    if (HasUnsupportedMedium(interaction)) {
      if (stats) {
        ++stats->unsupportedMediumInteractions;
      }
      return false;
    }
    std::optional<ResolvedDielectricTransition> transition =
      AnalyzeDielectricTransition(regionTable, regions, interaction, ray, stats);
    if (!transition || transition->kind == DielectricBoundaryKind::ScatteringInterface) {
      return false;
    }
    if (!CommitDielectricTransition(regions, transition->token, stats)) {
      return false;
    }
    if (stats) {
      ++stats->copiedVisibilityRegionTraversals;
    }
    SpawnedRay next = interaction.SpawnRayTo(visibility.p1.p);
    ray = next.ray;
    ray.tMax = static_cast<Float>(1) - RayEpsilon;
  }
  if (stats) {
    ++stats->dielectricStateErrors;
  }
  return false;
}

LightSampleContext MakeDirectLightSampleContext(
  const SurfaceInteraction& interaction,
  const BSDF& bsdf
) {
  LightSampleContext ctx(interaction);
  BxDFFlags flags = bsdf.Flags();
  if (base::IsReflective(flags) && !base::IsTransmissive(flags)) {
    ctx.p = interaction.OffsetRayOrigin(interaction.wo);
  } else if (base::IsTransmissive(flags) && !base::IsReflective(flags)) {
    ctx.p = interaction.OffsetRayOrigin(-interaction.wo);
  }
  return ctx;
}

Float AbsShadingCosine(const BSDF& bsdf, const vec3f& w) {
  return AbsDot(w, bsdf.Frame().Z());
}

} // namespace

SpectralRandomSampler::SpectralRandomSampler(std::uint64_t seed)
  : seed_(seed), state_(Mix(seed)) {}

void SpectralRandomSampler::StartPixelSample(
  FilmPoint2i pixel,
  int sampleIndex,
  std::uint64_t stream
) {
  std::uint64_t value = seed_;
  value ^= SplitMixIncrement + static_cast<std::uint64_t>(pixel.x) +
           (static_cast<std::uint64_t>(pixel.y) << 32);
  value = Mix(value);
  value ^= static_cast<std::uint64_t>(sampleIndex) + SplitMixIncrement;
  value = Mix(value);
  value ^= stream + SplitMixIncrement;
  state_ = Mix(value);
  dimension_ = 0;
}

Float SpectralRandomSampler::Get1D() {
  ++dimension_;
  return FloatFromBits(NextUInt());
}

point2f SpectralRandomSampler::Get2D() {
  return point2f(Get1D(), Get1D());
}

std::uint64_t SpectralRandomSampler::SamplesIssued() const {
  return dimension_;
}

std::uint64_t SpectralRandomSampler::Mix(std::uint64_t value) {
  value += SplitMixIncrement;
  value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31);
}

Float SpectralRandomSampler::FloatFromBits(std::uint64_t value) {
  constexpr double scale = 1.0 / 9007199254740992.0;
  double sample = static_cast<double>(value >> 11) * scale;
  Float result = static_cast<Float>(sample);
  return std::min(result, std::nextafter(static_cast<Float>(1), static_cast<Float>(0)));
}

std::uint64_t SpectralRandomSampler::NextUInt() {
  state_ += SplitMixIncrement;
  return Mix(state_);
}

RandomWalkWorkerState::RandomWalkWorkerState(
  std::uint64_t seed,
  std::size_t scratchBufferBytes
)
  : sampler_(seed), scratch_(scratchBufferBytes) {}

SpectralRandomSampler& RandomWalkWorkerState::Sampler() {
  return sampler_;
}

base::ScratchBuffer& RandomWalkWorkerState::Scratch() {
  return scratch_;
}

bool RecordSpectralRadianceDiagnostics(
  const base::SampledSpectrum& spectrum,
  RandomWalkRenderStats* stats
) {
  if (!spectrum.IsFinite()) {
    if (stats) {
      ++stats->nonFiniteRadiance;
    }
    return false;
  }
  if (HasNegativeComponent(spectrum)) {
    if (stats) {
      ++stats->negativeRadiance;
    }
    return false;
  }
  return true;
}

bool RecordSpectralRadianceDiagnostics(
  const base::SampledSpectrum& spectrum,
  PathRenderStats* stats
) {
  if (!spectrum.IsFinite()) {
    if (stats) {
      ++stats->nonFiniteRadiance;
    }
    return false;
  }
  if (HasNegativeComponent(spectrum)) {
    if (stats) {
      ++stats->negativeRadiance;
    }
    return false;
  }
  return true;
}

Float PowerHeuristic(int nf, Float fPdf, int ng, Float gPdf) {
  Float f = static_cast<Float>(nf) * fPdf;
  Float g = static_cast<Float>(ng) * gPdf;
  Float f2 = f * f;
  if (std::isinf(f2)) {
    return 1;
  }
  Float g2 = g * g;
  Float denom = f2 + g2;
  return denom > 0 ? f2 / denom : 0;
}

void ValidateRandomWalkScene(
  const Scene& scene,
  const SpectralLightTable& lightTable,
  const RandomWalkRenderOptions& options
) {
  ValidateOptions(options);
  if (!options.rejectDeltaLights) {
    return;
  }
  for (base::LightHandle handle : scene.Lights()) {
    const Light& light = lightTable.Get(handle);
    if (IsDeltaLight(light.Flags())) {
      throw std::invalid_argument(
        "RandomWalkIntegrator supports area and infinite lights only; delta lights are deferred"
      );
    }
  }
}

RandomWalkIntegrator::RandomWalkIntegrator(
  const Scene& scene,
  const materials::SpectralMaterialTable& materialTable,
  const SpectralLightTable& lightTable,
  RandomWalkRenderOptions options,
  const DielectricRegionTable* dielectricRegions
)
  : scene_(&scene),
    materialTable_(&materialTable),
    lightTable_(&lightTable),
    dielectricRegions_(dielectricRegions),
    options_(options) {
  ValidateRandomWalkScene(scene, lightTable, options_);
}

base::SampledSpectrum RandomWalkIntegrator::Li(
  const Ray& ray,
  base::SampledWavelengths& lambda,
  SpectralRandomSampler& sampler,
  base::ScratchBuffer& scratch,
  SpectralVisibleSurface* visibleSurface,
  RandomWalkRenderStats* stats,
  const DielectricPathState* initialRegions
) const {
  materials::UniversalTextureEvaluator textureEvaluator;
  base::SampledSpectrum L(0);
  base::SampledSpectrum beta(1);
  Ray currentRay = ray;
  DielectricPathState regions =
    InitialDielectricState(dielectricRegions_, *scene_, ray, initialRegions);
  int depth = 0;
  int alphaSkips = 0;

  while (true) {
    if (stats) {
      ++stats->raysTraced;
    }
    std::optional<PrimitiveIntersection> hit =
      scene_->Intersect(currentRay, RayEpsilon, currentRay.tMax);

    if (!hit) {
      for (base::LightHandle handle : scene_->InfiniteLights()) {
        const Light& light = lightTable_->Get(handle);
        base::SampledSpectrum emitted = light.Le(currentRay, lambda);
        if (emitted) {
          AddContribution(L, beta * emitted, stats);
          if (stats) {
            ++stats->infiniteLightHits;
          }
        }
      }
      return L;
    }

    SurfaceInteraction interaction = hit->interaction;
    if (stats) {
      ++stats->surfaceHits;
    }

    if (options_.rejectNonVacuum && HasUnsupportedMedium(interaction)) {
      if (stats) {
        ++stats->unsupportedMediumInteractions;
      }
      return L;
    }

    const materials::Material* material = nullptr;
    materials::MaterialEvalContext materialCtx(interaction);
    if (interaction.material.IsValid()) {
      material = &materialTable_->Get(interaction.material);
      materials::MaterialAlphaResult alpha =
        material->EvaluateAlpha(textureEvaluator, materialCtx, sampler.Get1D());
      if (!alpha.accepted) {
        if (++alphaSkips > options_.maxAlphaSkips) {
          if (stats) {
            ++stats->invalidSamples;
          }
          return L;
        }
        if (stats) {
          ++stats->alphaSkips;
        }
        currentRay = interaction.SpawnRay(currentRay.direction()).ray;
        continue;
      }
    }

    std::optional<ResolvedDielectricTransition> dielectricTransition =
      AnalyzeDielectricTransition(dielectricRegions_, regions, interaction, currentRay, stats);
    if (HasDielectricBoundary(interaction) && !dielectricTransition) {
      return L;
    }
    if (dielectricTransition && dielectricTransition->IsNullTraversal()) {
      if (!CommitDielectricTransition(regions, dielectricTransition->token, stats)) {
        return L;
      }
      currentRay = interaction.SpawnRay(currentRay.direction()).ray;
      currentRay.tMax = Infinity;
      continue;
    }
    if (dielectricTransition) {
      materialCtx.dielectric = &dielectricTransition->interface;
    }

    SetFirstVisibleSurface(interaction, visibleSurface);

    if (interaction.hasAreaLight) {
      const Light& light = lightTable_->Get(interaction.areaLight);
      base::SampledSpectrum emitted = light.L(interaction, interaction.wo, lambda);
      if (emitted) {
        AddContribution(L, beta * emitted, stats);
        if (stats) {
          ++stats->areaLightHits;
        }
      }
    }

    if (depth >= options_.maxDepth) {
      if (stats) {
        ++stats->maxDepthTerminations;
      }
      return L;
    }
    if (material == nullptr || !material->IsValid()) {
      return L;
    }

    scratch.Reset();
    render::BSDF bsdf = material->GetBSDF(textureEvaluator, materialCtx, lambda, scratch);
    if (!bsdf) {
      if (stats) {
        ++stats->nullBSDFs;
      }
      return L;
    }
    if (stats) {
      ++stats->materialClosures;
    }

    std::optional<BSDFSample> bs = bsdf.Sample_f(
      interaction.wo,
      sampler.Get1D(),
      sampler.Get2D(),
      base::TransportMode::Radiance
    );
    if (!bs) {
      if (stats) {
        ++stats->nullBSDFs;
      }
      return L;
    }
    if (!IsValidPositivePDF(bs->pdf)) {
      if (stats) {
        ++stats->invalidPDFs;
      }
      return L;
    }

    Float cosTerm = AbsDot(bs->wi, bsdf.Frame().Z());
    if (!(cosTerm > 0) || !std::isfinite(cosTerm)) {
      if (stats) {
        ++stats->invalidSamples;
      }
      return L;
    }
    base::SampledSpectrum factor = bs->f * (cosTerm / bs->pdf);
    if (!IsFiniteNonNegative(factor)) {
      RecordSpectralRadianceDiagnostics(factor, stats);
      if (stats) {
        ++stats->invalidSamples;
      }
      return L;
    }
    beta *= factor;
    if (!IsFiniteNonNegative(beta)) {
      RecordSpectralRadianceDiagnostics(beta, stats);
      if (stats) {
        ++stats->invalidSamples;
      }
      return L;
    }
    if (!beta) {
      return L;
    }

    if (dielectricTransition && bs->IsTransmission()) {
      if (!CommitDielectricTransition(regions, dielectricTransition->token, stats)) {
        return L;
      }
      if (stats) {
        ++stats->dielectricTransmissionCommits;
      }
    }

    currentRay = interaction.SpawnRay(bs->wi).ray;
    currentRay.tMax = Infinity;
    ++depth;
  }
}

RandomWalkRenderStats RenderRandomWalk(
  const Scene& scene,
  const materials::SpectralMaterialTable& materialTable,
  SpectralLightTable& lightTable,
  const SpectralCamera& camera,
  Film& film,
  const RandomWalkRenderOptions& options,
  const DielectricRegionTable* dielectricRegions
) {
  ValidateOptions(options);
  Bounds3f sceneBounds;
  scene.Bounds(&sceneBounds);
  lightTable.Preprocess(sceneBounds);
  ValidateRandomWalkScene(scene, lightTable, options);

  RandomWalkRenderStats stats;
  RandomWalkIntegrator integrator(scene, materialTable, lightTable, options, dielectricRegions);
  RandomWalkWorkerState worker(options.seed, options.scratchBufferBytes);

  for (int y = 0; y < film.Height(); ++y) {
    for (int x = 0; x < film.Width(); ++x) {
      FilmPoint2i pixel{x, y};
      for (int sampleIndex = 0; sampleIndex < options.pixelSamples; ++sampleIndex) {
        worker.Sampler().StartPixelSample(pixel, sampleIndex);
        ++stats.pixelSamples;

        CameraSample sample;
        if (options.jitterCameraSamples) {
          sample.pFilm = point2f(
            static_cast<Float>(x) + worker.Sampler().Get1D(),
            static_cast<Float>(y) + worker.Sampler().Get1D()
          );
          sample.pLens = worker.Sampler().Get2D();
          sample.time = worker.Sampler().Get1D();
        } else {
          sample.pFilm = point2f(
            static_cast<Float>(x) + static_cast<Float>(0.5),
            static_cast<Float>(y) + static_cast<Float>(0.5)
          );
          sample.pLens = point2f(static_cast<Float>(0.5), static_cast<Float>(0.5));
          sample.time = static_cast<Float>(0.5);
        }
        sample.filterWeight = 1;

        Float wavelengthSample = worker.Sampler().Get1D();
        CameraFilmSample cameraSample = GenerateCameraRayFromFilm(
          film,
          camera,
          sample,
          wavelengthSample,
          options.generateDifferentials
        );
        if (!cameraSample.cameraRay) {
          ++stats.invalidSamples;
          continue;
        }
        ++stats.cameraRays;

        if (options.rejectNonVacuum && cameraSample.cameraRay->hasInitialMedium) {
          ++stats.unsupportedMediumInteractions;
          continue;
        }

        worker.Scratch().Reset();
        SpectralVisibleSurface visibleSurface;
        base::SampledSpectrum L = integrator.Li(
          cameraSample.cameraRay->ray,
          cameraSample.wavelengths,
          worker.Sampler(),
          worker.Scratch(),
          &visibleSurface,
          &stats
        );
        L *= cameraSample.cameraRay->weight;
        if (RecordSpectralRadianceDiagnostics(L, &stats)) {
          film.AddFilteredSample(
            FilmPoint2f{sample.pFilm.xy.x, sample.pFilm.xy.y},
            L,
            cameraSample.wavelengths,
            &visibleSurface,
            sample.filterWeight
          );
        } else {
          ++stats.invalidSamples;
        }
        worker.Scratch().Reset();
      }
    }
  }

  return stats;
}

PathIntegrator::PathIntegrator(
  const Scene& scene,
  const materials::SpectralMaterialTable& materialTable,
  const SpectralLightTable& lightTable,
  PathRenderOptions options,
  const DielectricRegionTable* dielectricRegions
)
  : scene_(&scene),
    materialTable_(&materialTable),
    lightTable_(&lightTable),
    dielectricRegions_(dielectricRegions),
    lightSampler_(lightTable, SceneLightHandles(scene)),
    options_(options) {
  ValidateOptions(options_);
}

base::SampledSpectrum PathIntegrator::SampleLd(
  const SurfaceInteraction& interaction,
  const BSDF& bsdf,
  base::SampledWavelengths& lambda,
  SpectralRandomSampler& sampler,
  PathRenderStats* stats,
  const DielectricPathState& regions
) const {
  if (lightSampler_.Size() == 0) {
    return base::SampledSpectrum(0);
  }
  LightSampleContext ctx = MakeDirectLightSampleContext(interaction, bsdf);
  std::optional<SampledLight> sampledLight = lightSampler_.Sample(sampler.Get1D());
  point2f uLight = sampler.Get2D();
  if (!sampledLight || sampledLight->light == nullptr ||
      !IsValidPositivePDF(sampledLight->pmf)) {
    if (stats) {
      ++stats->invalidPDFs;
    }
    return base::SampledSpectrum(0);
  }
  if (stats) {
    ++stats->directLightSamples;
  }

  std::optional<LightLiSample> ls = sampledLight->light->SampleLi(
    ctx,
    uLight,
    lambda,
    LightSampleMode::CompletePDF
  );
  if (!ls || !ls->L || !IsValidPositivePDF(ls->pdf)) {
    if (ls && !IsValidPositivePDF(ls->pdf) && stats) {
      ++stats->invalidPDFs;
    }
    return base::SampledSpectrum(0);
  }
  if (ls->delta && stats) {
    ++stats->deltaLightSamples;
  }

  base::SampledSpectrum f =
    bsdf.f(interaction.wo, ls->wi, base::TransportMode::Radiance) *
    AbsShadingCosine(bsdf, ls->wi);
  if (!f || !RecordSpectralRadianceDiagnostics(f, stats)) {
    return base::SampledSpectrum(0);
  }
  if (!RegionAwareUnoccluded(ls->visibility, *scene_, dielectricRegions_, regions, stats)) {
    if (stats) {
      ++stats->occludedShadowRays;
    }
    return base::SampledSpectrum(0);
  }

  Float pLight = sampledLight->pmf * ls->pdf;
  if (!IsValidPositivePDF(pLight)) {
    if (stats) {
      ++stats->invalidPDFs;
    }
    return base::SampledSpectrum(0);
  }
  Float weight = 1;
  if (!ls->delta) {
    Float pBSDF = bsdf.PDF(interaction.wo, ls->wi, base::TransportMode::Radiance);
    if (pBSDF < 0 || !std::isfinite(pBSDF)) {
      if (stats) {
        ++stats->invalidPDFs;
      }
      return base::SampledSpectrum(0);
    }
    weight = PowerHeuristic(1, pLight, 1, pBSDF);
  }

  base::SampledSpectrum contribution = ls->L * f * (weight / pLight);
  if (!RecordSpectralRadianceDiagnostics(contribution, stats)) {
    if (stats) {
      ++stats->invalidSamples;
    }
    return base::SampledSpectrum(0);
  }
  if (contribution && stats) {
    ++stats->directLightContributions;
  }
  return contribution;
}

base::SampledSpectrum PathIntegrator::Li(
  const Ray& ray,
  base::SampledWavelengths& lambda,
  SpectralRandomSampler& sampler,
  base::ScratchBuffer& scratch,
  SpectralVisibleSurface* visibleSurface,
  PathRenderStats* stats,
  const DielectricPathState* initialRegions
) const {
  materials::UniversalTextureEvaluator textureEvaluator;
  base::SampledSpectrum L(0);
  base::SampledSpectrum beta(1);
  Ray currentRay = ray;
  DielectricPathState regions =
    InitialDielectricState(dielectricRegions_, *scene_, ray, initialRegions);
  int depth = 0;
  int nullSkips = 0;
  Float pBSDF = 0;
  Float etaScale = 1;
  bool specularBounce = false;
  bool anyNonSpecularBounces = false;
  LightSampleContext previousLightContext;

  while (true) {
    if (stats) {
      ++stats->raysTraced;
    }
    std::optional<PrimitiveIntersection> hit =
      scene_->Intersect(currentRay, RayEpsilon, currentRay.tMax);

    if (!hit) {
      for (base::LightHandle handle : scene_->InfiniteLights()) {
        const Light& light = lightTable_->Get(handle);
        base::SampledSpectrum emitted = light.Le(currentRay, lambda);
        if (!emitted) {
          continue;
        }
        Float weight = 1;
        if (depth != 0 && !specularBounce) {
          Float pLight = lightSampler_.PMF(handle) *
                         light.PDF_Li(
                           previousLightContext,
                           currentRay.direction(),
                           LightSampleMode::CompletePDF
                         );
          if (pLight < 0 || !std::isfinite(pLight)) {
            if (stats) {
              ++stats->invalidPDFs;
            }
            continue;
          }
          weight = PowerHeuristic(1, pBSDF, 1, pLight);
          if (stats) {
            ++stats->emitterHitMIS;
          }
        }
        base::SampledSpectrum contribution = beta * emitted * weight;
        if (RecordSpectralRadianceDiagnostics(contribution, stats)) {
          L += contribution;
          if (stats) {
            ++stats->infiniteLightHits;
          }
        } else if (stats) {
          ++stats->invalidSamples;
        }
      }
      return L;
    }

    SurfaceInteraction interaction = hit->interaction;
    if (stats) {
      ++stats->surfaceHits;
    }

    if (options_.rejectNonVacuum && HasUnsupportedMedium(interaction)) {
      if (stats) {
        ++stats->unsupportedMediumInteractions;
      }
      return L;
    }

    const materials::Material* material = nullptr;
    materials::MaterialEvalContext materialCtx(interaction);
    if (interaction.material.IsValid()) {
      material = &materialTable_->Get(interaction.material);
      materials::MaterialAlphaResult alpha =
        material->EvaluateAlpha(textureEvaluator, materialCtx, sampler.Get1D());
      if (!alpha.accepted) {
        if (++nullSkips > options_.maxNullSkips) {
          if (stats) {
            ++stats->invalidSamples;
          }
          return L;
        }
        if (stats) {
          ++stats->nullSkips;
        }
        currentRay = interaction.SpawnRay(currentRay.direction()).ray;
        currentRay.tMax = Infinity;
        continue;
      }
    }

    std::optional<ResolvedDielectricTransition> dielectricTransition =
      AnalyzeDielectricTransition(dielectricRegions_, regions, interaction, currentRay, stats);
    if (HasDielectricBoundary(interaction) && !dielectricTransition) {
      return L;
    }
    if (dielectricTransition && dielectricTransition->IsNullTraversal()) {
      if (++nullSkips > options_.maxNullSkips) {
        if (stats) {
          ++stats->invalidSamples;
        }
        return L;
      }
      if (!CommitDielectricTransition(regions, dielectricTransition->token, stats)) {
        return L;
      }
      if (stats) {
        ++stats->nullSkips;
      }
      currentRay = interaction.SpawnRay(currentRay.direction()).ray;
      currentRay.tMax = Infinity;
      continue;
    }
    if (dielectricTransition) {
      materialCtx.dielectric = &dielectricTransition->interface;
    }

    if (interaction.hasAreaLight) {
      const Light& light = lightTable_->Get(interaction.areaLight);
      base::SampledSpectrum emitted = light.L(interaction, interaction.wo, lambda);
      if (emitted) {
        Float weight = 1;
        if (depth != 0 && !specularBounce) {
          Float pLight = lightSampler_.PMF(interaction.areaLight) *
                         light.PDF_Li(
                           previousLightContext,
                           currentRay.direction(),
                           LightSampleMode::CompletePDF
                         );
          if (pLight < 0 || !std::isfinite(pLight)) {
            if (stats) {
              ++stats->invalidPDFs;
            }
            return L;
          }
          weight = PowerHeuristic(1, pBSDF, 1, pLight);
          if (stats) {
            ++stats->emitterHitMIS;
          }
        }
        base::SampledSpectrum contribution = beta * emitted * weight;
        if (RecordSpectralRadianceDiagnostics(contribution, stats)) {
          L += contribution;
          if (stats) {
            ++stats->areaLightHits;
          }
        } else if (stats) {
          ++stats->invalidSamples;
        }
      }
    }

    if (material == nullptr || !material->IsValid()) {
      return L;
    }

    scratch.Reset();
    BSDF bsdf = material->GetBSDF(textureEvaluator, materialCtx, lambda, scratch);
    if (stats) {
      ++stats->materialClosures;
    }

    if (!bsdf || IsNullBSDF(bsdf)) {
      if (++nullSkips > options_.maxNullSkips) {
        if (stats) {
          ++stats->invalidSamples;
        }
        return L;
      }
      if (stats) {
        ++stats->nullSkips;
      }
      currentRay = interaction.SpawnRay(currentRay.direction()).ray;
      currentRay.tMax = Infinity;
      continue;
    }

    SetFirstVisibleSurface(interaction, visibleSurface);

    if (options_.regularize && anyNonSpecularBounces) {
      bsdf.Regularize();
    }

    if (depth++ == options_.maxDepth) {
      if (stats) {
        ++stats->maxDepthTerminations;
      }
      return L;
    }

    if (options_.sampleDirectLighting && base::IsNonSpecular(bsdf.Flags())) {
      base::SampledSpectrum Ld = SampleLd(interaction, bsdf, lambda, sampler, stats, regions);
      if (Ld) {
        base::SampledSpectrum contribution = beta * Ld;
        if (RecordSpectralRadianceDiagnostics(contribution, stats)) {
          L += contribution;
        } else if (stats) {
          ++stats->invalidSamples;
        }
      }
    }

    if (!options_.sampleBSDF) {
      return L;
    }

    std::optional<BSDFSample> bs = bsdf.Sample_f(
      interaction.wo,
      sampler.Get1D(),
      sampler.Get2D(),
      base::TransportMode::Radiance
    );
    if (!bs) {
      return L;
    }
    if (!IsValidPositivePDF(bs->pdf)) {
      if (stats) {
        ++stats->invalidPDFs;
      }
      return L;
    }

    Float cosTerm = AbsShadingCosine(bsdf, bs->wi);
    if (!(cosTerm > 0) || !std::isfinite(cosTerm)) {
      if (stats) {
        ++stats->invalidSamples;
      }
      return L;
    }
    base::SampledSpectrum factor = bs->f * (cosTerm / bs->pdf);
    if (!IsFiniteNonNegative(factor)) {
      RecordSpectralRadianceDiagnostics(factor, stats);
      if (stats) {
        ++stats->invalidSamples;
      }
      return L;
    }

    beta *= factor;
    if (!IsFiniteNonNegative(beta)) {
      RecordSpectralRadianceDiagnostics(beta, stats);
      if (stats) {
        ++stats->invalidSamples;
      }
      return L;
    }
    if (!beta) {
      return L;
    }

    pBSDF = bs->pdfIsProportional
              ? bsdf.PDF(interaction.wo, bs->wi, base::TransportMode::Radiance)
              : bs->pdf;
    if (pBSDF < 0 || !std::isfinite(pBSDF)) {
      if (stats) {
        ++stats->invalidPDFs;
      }
      return L;
    }
    specularBounce = bs->IsSpecular();
    anyNonSpecularBounces = anyNonSpecularBounces || !specularBounce;
    if (bs->IsTransmission()) {
      etaScale *= bs->eta * bs->eta;
      if (dielectricTransition) {
        if (!CommitDielectricTransition(regions, dielectricTransition->token, stats)) {
          return L;
        }
        if (stats) {
          ++stats->dielectricTransmissionCommits;
        }
      }
    }
    previousLightContext = LightSampleContext(interaction);
    if (stats) {
      ++stats->bsdfSamples;
    }

    currentRay = interaction.SpawnRay(bs->wi).ray;
    currentRay.tMax = Infinity;

    if (options_.russianRoulette) {
      base::SampledSpectrum rrBeta = beta * etaScale;
      if (rrBeta.MaxComponentValue() < 1 && depth > 1) {
        if (stats) {
          ++stats->russianRouletteChecks;
        }
        Float q = std::max(static_cast<Float>(0), static_cast<Float>(1) - rrBeta.MaxComponentValue());
        if (sampler.Get1D() < q) {
          if (stats) {
            ++stats->russianRouletteTerminations;
          }
          return L;
        }
        beta /= static_cast<Float>(1) - q;
        if (!IsFiniteNonNegative(beta)) {
          RecordSpectralRadianceDiagnostics(beta, stats);
          if (stats) {
            ++stats->invalidSamples;
          }
          return L;
        }
      }
    }
  }
}

PathRenderStats RenderPath(
  const Scene& scene,
  const materials::SpectralMaterialTable& materialTable,
  SpectralLightTable& lightTable,
  const SpectralCamera& camera,
  Film& film,
  const PathRenderOptions& options,
  const DielectricRegionTable* dielectricRegions
) {
  ValidateOptions(options);
  Bounds3f sceneBounds;
  scene.Bounds(&sceneBounds);
  lightTable.Preprocess(sceneBounds);

  PathRenderStats stats;
  PathIntegrator integrator(scene, materialTable, lightTable, options, dielectricRegions);
  RandomWalkWorkerState worker(options.seed, options.scratchBufferBytes);

  for (int y = 0; y < film.Height(); ++y) {
    for (int x = 0; x < film.Width(); ++x) {
      FilmPoint2i pixel{x, y};
      for (int sampleIndex = 0; sampleIndex < options.pixelSamples; ++sampleIndex) {
        worker.Sampler().StartPixelSample(pixel, sampleIndex);
        ++stats.pixelSamples;

        CameraSample sample;
        if (options.jitterCameraSamples) {
          sample.pFilm = point2f(
            static_cast<Float>(x) + worker.Sampler().Get1D(),
            static_cast<Float>(y) + worker.Sampler().Get1D()
          );
          sample.pLens = worker.Sampler().Get2D();
          sample.time = worker.Sampler().Get1D();
        } else {
          sample.pFilm = point2f(
            static_cast<Float>(x) + static_cast<Float>(0.5),
            static_cast<Float>(y) + static_cast<Float>(0.5)
          );
          sample.pLens = point2f(static_cast<Float>(0.5), static_cast<Float>(0.5));
          sample.time = static_cast<Float>(0.5);
        }
        sample.filterWeight = 1;

        Float wavelengthSample = worker.Sampler().Get1D();
        CameraFilmSample cameraSample = GenerateCameraRayFromFilm(
          film,
          camera,
          sample,
          wavelengthSample,
          options.generateDifferentials
        );
        if (!cameraSample.cameraRay) {
          ++stats.invalidSamples;
          continue;
        }
        ++stats.cameraRays;

        if (options.rejectNonVacuum && cameraSample.cameraRay->hasInitialMedium) {
          ++stats.unsupportedMediumInteractions;
          continue;
        }

        worker.Scratch().Reset();
        SpectralVisibleSurface visibleSurface;
        base::SampledSpectrum L = integrator.Li(
          cameraSample.cameraRay->ray,
          cameraSample.wavelengths,
          worker.Sampler(),
          worker.Scratch(),
          &visibleSurface,
          &stats
        );
        L *= cameraSample.cameraRay->weight;
        if (RecordSpectralRadianceDiagnostics(L, &stats)) {
          film.AddFilteredSample(
            FilmPoint2f{sample.pFilm.xy.x, sample.pFilm.xy.y},
            L,
            cameraSample.wavelengths,
            &visibleSurface,
            sample.filterWeight
          );
        } else {
          ++stats.invalidSamples;
        }
        worker.Scratch().Reset();
      }
    }
  }

  return stats;
}

} // namespace render
} // namespace rayrender
