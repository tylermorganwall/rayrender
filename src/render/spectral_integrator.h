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
    RandomWalkRenderOptions options = {}
  );

  base::SampledSpectrum Li(
    const Ray& ray,
    base::SampledWavelengths& lambda,
    SpectralRandomSampler& sampler,
    base::ScratchBuffer& scratch,
    SpectralVisibleSurface* visibleSurface = nullptr,
    RandomWalkRenderStats* stats = nullptr
  ) const;

private:
  const Scene* scene_ = nullptr;
  const materials::SpectralMaterialTable* materialTable_ = nullptr;
  const SpectralLightTable* lightTable_ = nullptr;
  RandomWalkRenderOptions options_;
};

RandomWalkRenderStats RenderRandomWalk(
  const Scene& scene,
  const materials::SpectralMaterialTable& materialTable,
  SpectralLightTable& lightTable,
  const SpectralCamera& camera,
  Film& film,
  const RandomWalkRenderOptions& options = {}
);

} // namespace render
} // namespace rayrender

#endif
