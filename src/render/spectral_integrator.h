#ifndef RAYRENDER_RENDER_SPECTRAL_INTEGRATOR_H
#define RAYRENDER_RENDER_SPECTRAL_INTEGRATOR_H

#include "spectral_camera.h"
#include "spectral_film.h"
#include "spectral_light.h"
#include "spectral_scene.h"

#include "../base/base.h"
#include "../materials/spectral_material.h"

#include <cstddef>
#include <cstdint>

namespace rayrender {
namespace render {

class SpectralRandomSampler {
public:
  explicit SpectralRandomSampler(std::uint64_t seed = 0);

  void StartPixelSample(FilmPoint2i pixel, int sampleIndex, std::uint64_t stream = 0);
  Float Get1D();
  point2f Get2D();
  std::uint64_t SamplesIssued() const;

private:
  static std::uint64_t Mix(std::uint64_t value);
  static Float FloatFromBits(std::uint64_t value);
  std::uint64_t NextUInt();

  std::uint64_t seed_ = 0;
  std::uint64_t state_ = 0;
  std::uint64_t dimension_ = 0;
};

struct RandomWalkRenderOptions {
  int pixelSamples = 1;
  int maxDepth = 50;
  std::uint64_t seed = 0;
  std::size_t scratchBufferBytes = 4096;
  bool jitterCameraSamples = true;
  bool generateDifferentials = false;
  bool rejectNonVacuum = true;
  bool rejectDeltaLights = true;
  int maxAlphaSkips = 64;
};

struct RandomWalkRenderStats {
  std::uint64_t pixelSamples = 0;
  std::uint64_t cameraRays = 0;
  std::uint64_t raysTraced = 0;
  std::uint64_t surfaceHits = 0;
  std::uint64_t materialClosures = 0;
  std::uint64_t areaLightHits = 0;
  std::uint64_t infiniteLightHits = 0;
  std::uint64_t alphaSkips = 0;
  std::uint64_t nullBSDFs = 0;
  std::uint64_t maxDepthTerminations = 0;
  std::uint64_t invalidSamples = 0;
  std::uint64_t invalidPDFs = 0;
  std::uint64_t nonFiniteRadiance = 0;
  std::uint64_t negativeRadiance = 0;
  std::uint64_t unsupportedMediumInteractions = 0;
  std::uint64_t dielectricPrioritySkips = 0;
  std::uint64_t dielectricIndexMatchedSkips = 0;
  std::uint64_t dielectricScatteringInterfaces = 0;
  std::uint64_t dielectricTransmissionCommits = 0;
  std::uint64_t dielectricStateErrors = 0;
  std::uint64_t wavelengthTerminations = 0;
  std::uint64_t dielectricWavelengthTerminations = 0;
  std::uint64_t thinDielectricWavelengthTerminations = 0;
};

class RandomWalkWorkerState {
public:
  RandomWalkWorkerState(std::uint64_t seed, std::size_t scratchBufferBytes);

  SpectralRandomSampler& Sampler();
  base::ScratchBuffer& Scratch();

private:
  SpectralRandomSampler sampler_;
  base::ScratchBuffer scratch_;
};

bool RecordSpectralRadianceDiagnostics(
  const base::SampledSpectrum& spectrum,
  RandomWalkRenderStats* stats
);

Float PowerHeuristic(int nf, Float fPdf, int ng, Float gPdf);

void ValidateRandomWalkScene(
  const Scene& scene,
  const SpectralLightTable& lightTable,
  const RandomWalkRenderOptions& options
);

class RandomWalkIntegrator {
public:
  RandomWalkIntegrator(
    const Scene& scene,
    const materials::SpectralMaterialTable& materialTable,
    const SpectralLightTable& lightTable,
    RandomWalkRenderOptions options = {},
    const DielectricRegionTable* dielectricRegions = nullptr
  );

  base::SampledSpectrum Li(
    const Ray& ray,
    base::SampledWavelengths& lambda,
    SpectralRandomSampler& sampler,
    base::ScratchBuffer& scratch,
    SpectralVisibleSurface* visibleSurface = nullptr,
    RandomWalkRenderStats* stats = nullptr,
    const DielectricPathState* initialRegions = nullptr
  ) const;

private:
  const Scene* scene_ = nullptr;
  const materials::SpectralMaterialTable* materialTable_ = nullptr;
  const SpectralLightTable* lightTable_ = nullptr;
  const DielectricRegionTable* dielectricRegions_ = nullptr;
  RandomWalkRenderOptions options_;
};

RandomWalkRenderStats RenderRandomWalk(
  const Scene& scene,
  const materials::SpectralMaterialTable& materialTable,
  SpectralLightTable& lightTable,
  const SpectralCamera& camera,
  Film& film,
  const RandomWalkRenderOptions& options = {},
  const DielectricRegionTable* dielectricRegions = nullptr
);

struct PathRenderOptions {
  int pixelSamples = 1;
  int maxDepth = 5;
  std::uint64_t seed = 0;
  std::size_t scratchBufferBytes = 4096;
  bool jitterCameraSamples = true;
  bool generateDifferentials = false;
  bool rejectNonVacuum = true;
  bool regularize = false;
  bool russianRoulette = true;
  bool sampleDirectLighting = true;
  bool sampleBSDF = true;
  int maxNullSkips = 64;
};

struct PathRenderStats {
  std::uint64_t pixelSamples = 0;
  std::uint64_t cameraRays = 0;
  std::uint64_t raysTraced = 0;
  std::uint64_t surfaceHits = 0;
  std::uint64_t materialClosures = 0;
  std::uint64_t directLightSamples = 0;
  std::uint64_t directLightContributions = 0;
  std::uint64_t deltaLightSamples = 0;
  std::uint64_t bsdfSamples = 0;
  std::uint64_t areaLightHits = 0;
  std::uint64_t infiniteLightHits = 0;
  std::uint64_t emitterHitMIS = 0;
  std::uint64_t nullSkips = 0;
  std::uint64_t maxDepthTerminations = 0;
  std::uint64_t russianRouletteChecks = 0;
  std::uint64_t russianRouletteTerminations = 0;
  std::uint64_t invalidSamples = 0;
  std::uint64_t invalidPDFs = 0;
  std::uint64_t occludedShadowRays = 0;
  std::uint64_t nonFiniteRadiance = 0;
  std::uint64_t negativeRadiance = 0;
  std::uint64_t unsupportedMediumInteractions = 0;
  std::uint64_t dielectricPrioritySkips = 0;
  std::uint64_t dielectricIndexMatchedSkips = 0;
  std::uint64_t dielectricScatteringInterfaces = 0;
  std::uint64_t dielectricTransmissionCommits = 0;
  std::uint64_t dielectricStateErrors = 0;
  std::uint64_t copiedVisibilityRegionTraversals = 0;
  std::uint64_t wavelengthTerminations = 0;
  std::uint64_t dielectricWavelengthTerminations = 0;
  std::uint64_t thinDielectricWavelengthTerminations = 0;
};

bool RecordSpectralRadianceDiagnostics(
  const base::SampledSpectrum& spectrum,
  PathRenderStats* stats
);

class PathIntegrator {
public:
  PathIntegrator(
    const Scene& scene,
    const materials::SpectralMaterialTable& materialTable,
    const SpectralLightTable& lightTable,
    PathRenderOptions options = {},
    const DielectricRegionTable* dielectricRegions = nullptr
  );

  base::SampledSpectrum Li(
    const Ray& ray,
    base::SampledWavelengths& lambda,
    SpectralRandomSampler& sampler,
    base::ScratchBuffer& scratch,
    SpectralVisibleSurface* visibleSurface = nullptr,
    PathRenderStats* stats = nullptr,
    const DielectricPathState* initialRegions = nullptr
  ) const;

private:
  base::SampledSpectrum SampleLd(
    const SurfaceInteraction& interaction,
    const BSDF& bsdf,
    base::SampledWavelengths& lambda,
    SpectralRandomSampler& sampler,
    PathRenderStats* stats,
    const DielectricPathState& regions
  ) const;

  const Scene* scene_ = nullptr;
  const materials::SpectralMaterialTable* materialTable_ = nullptr;
  const SpectralLightTable* lightTable_ = nullptr;
  const DielectricRegionTable* dielectricRegions_ = nullptr;
  UniformLightSampler lightSampler_;
  PathRenderOptions options_;
};

PathRenderStats RenderPath(
  const Scene& scene,
  const materials::SpectralMaterialTable& materialTable,
  SpectralLightTable& lightTable,
  const SpectralCamera& camera,
  Film& film,
  const PathRenderOptions& options = {},
  const DielectricRegionTable* dielectricRegions = nullptr
);

} // namespace render
} // namespace rayrender

#endif
