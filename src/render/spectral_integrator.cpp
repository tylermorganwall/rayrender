#include "spectral_integrator.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

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
  RandomWalkRenderOptions options
)
  : scene_(&scene),
    materialTable_(&materialTable),
    lightTable_(&lightTable),
    options_(options) {
  ValidateRandomWalkScene(scene, lightTable, options_);
}

base::SampledSpectrum RandomWalkIntegrator::Li(
  const Ray& ray,
  base::SampledWavelengths& lambda,
  SpectralRandomSampler& sampler,
  base::ScratchBuffer& scratch,
  SpectralVisibleSurface* visibleSurface,
  RandomWalkRenderStats* stats
) const {
  materials::UniversalTextureEvaluator textureEvaluator;
  base::SampledSpectrum L(0);
  base::SampledSpectrum beta(1);
  Ray currentRay = ray;
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

    if (
      options_.rejectNonVacuum &&
      (interaction.hasMedium || interaction.hasMediumInterface || interaction.hasDielectricRegion)
    ) {
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
  const RandomWalkRenderOptions& options
) {
  ValidateOptions(options);
  Bounds3f sceneBounds;
  scene.Bounds(&sceneBounds);
  lightTable.Preprocess(sceneBounds);
  ValidateRandomWalkScene(scene, lightTable, options);

  RandomWalkRenderStats stats;
  RandomWalkIntegrator integrator(scene, materialTable, lightTable, options);
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
