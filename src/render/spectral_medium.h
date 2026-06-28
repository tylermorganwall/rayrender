#ifndef RAYRENDER_RENDER_SPECTRAL_MEDIUM_H
#define RAYRENDER_RENDER_SPECTRAL_MEDIUM_H

#include "../base/base.h"
#include "../core/ray.h"
#include "../math/vectypes.h"

#include <cstddef>
#include <optional>
#include <vector>

namespace rayrender {
namespace render {

struct PhaseFunctionSample {
  Float p = 0;
  vec3f wi;
  Float pdf = 0;
};

class HenyeyGreensteinPhaseFunction {
public:
  HenyeyGreensteinPhaseFunction() = default;
  explicit HenyeyGreensteinPhaseFunction(Float g);

  Float p(const vec3f& wo, const vec3f& wi) const;
  std::optional<PhaseFunctionSample> Sample_p(const vec3f& wo, point2f u) const;
  Float PDF(const vec3f& wo, const vec3f& wi) const;
  Float G() const;

private:
  Float g_ = 0;
};

class PhaseFunction {
public:
  PhaseFunction() = default;
  explicit PhaseFunction(const HenyeyGreensteinPhaseFunction* phase);

  explicit operator bool() const;
  Float p(const vec3f& wo, const vec3f& wi) const;
  std::optional<PhaseFunctionSample> Sample_p(const vec3f& wo, point2f u) const;
  Float PDF(const vec3f& wo, const vec3f& wi) const;

private:
  enum class Kind {
    None,
    HenyeyGreenstein
  };

  Kind kind_ = Kind::None;
  const void* ptr_ = nullptr;
};

struct RayMajorantSegment {
  Float tMin = 0;
  Float tMax = 0;
  base::SampledSpectrum sigmaMaj = base::SampledSpectrum(0);
};

class HomogeneousMajorantIterator {
public:
  HomogeneousMajorantIterator() = default;
  HomogeneousMajorantIterator(Float tMin, Float tMax, base::SampledSpectrum sigmaMaj);

  std::optional<RayMajorantSegment> Next();

private:
  RayMajorantSegment segment_;
  bool called_ = true;
};

class RayMajorantIterator {
public:
  RayMajorantIterator() = default;
  explicit RayMajorantIterator(HomogeneousMajorantIterator iterator);

  std::optional<RayMajorantSegment> Next();

private:
  std::optional<HomogeneousMajorantIterator> homogeneous_;
};

struct MediumProperties {
  base::SampledSpectrum sigmaA = base::SampledSpectrum(0);
  base::SampledSpectrum sigmaS = base::SampledSpectrum(0);
  PhaseFunction phase;
  base::SampledSpectrum emission = base::SampledSpectrum(0);

  base::SampledSpectrum SigmaT() const;
  bool HasScattering() const;
};

struct MediumSample {
  base::SampledSpectrum beta = base::SampledSpectrum(0);
  base::SampledSpectrum transmittance = base::SampledSpectrum(0);
  point3f p;
  Float t = 0;
  Float pdf = 0;
  PhaseFunction phase;
  bool sampledMedium = false;
};

class HomogeneousMedium {
public:
  HomogeneousMedium() = default;
  HomogeneousMedium(
    base::Spectrum sigmaA,
    base::Spectrum sigmaS,
    Float scale = 1,
    base::Spectrum emission = base::Spectrum(base::ConstantSpectrum(0)),
    Float emissionScale = 1,
    Float g = 0
  );

  bool IsEmissive() const;
  MediumProperties SamplePoint(
    const point3f& p,
    const base::SampledWavelengths& lambda
  ) const;
  RayMajorantIterator SampleRay(
    const Ray& ray,
    Float tMax,
    const base::SampledWavelengths& lambda
  ) const;
  base::SampledSpectrum Tr(
    const Ray& ray,
    Float tMax,
    const base::SampledWavelengths& lambda
  ) const;
  std::optional<MediumSample> SampleDistance(
    const Ray& ray,
    Float tMax,
    const base::SampledWavelengths& lambda,
    Float uChannel,
    Float uDistance
  ) const;

  const base::Spectrum& SigmaASpectrum() const;
  const base::Spectrum& SigmaSSpectrum() const;
  const base::Spectrum& EmissionSpectrum() const;
  Float Scale() const;
  Float EmissionScale() const;
  const HenyeyGreensteinPhaseFunction& Phase() const;

private:
  base::Spectrum sigmaA_ = base::Spectrum(base::ConstantSpectrum(0));
  base::Spectrum sigmaS_ = base::Spectrum(base::ConstantSpectrum(0));
  base::Spectrum emission_ = base::Spectrum(base::ConstantSpectrum(0));
  Float scale_ = 1;
  Float emissionScale_ = 1;
  HenyeyGreensteinPhaseFunction phase_;
};

class Medium {
public:
  Medium() = default;
  explicit Medium(const HomogeneousMedium* medium);

  explicit operator bool() const;
  bool IsEmissive() const;
  MediumProperties SamplePoint(
    const point3f& p,
    const base::SampledWavelengths& lambda
  ) const;
  RayMajorantIterator SampleRay(
    const Ray& ray,
    Float tMax,
    const base::SampledWavelengths& lambda
  ) const;
  base::SampledSpectrum Tr(
    const Ray& ray,
    Float tMax,
    const base::SampledWavelengths& lambda
  ) const;
  std::optional<MediumSample> SampleDistance(
    const Ray& ray,
    Float tMax,
    const base::SampledWavelengths& lambda,
    Float uChannel,
    Float uDistance
  ) const;

private:
  enum class Kind {
    None,
    Homogeneous
  };

  Kind kind_ = Kind::None;
  const void* ptr_ = nullptr;
};

class PhaseFunctionTable {
public:
  base::PhaseFunctionHandle Add(HenyeyGreensteinPhaseFunction phase);
  PhaseFunction Get(base::PhaseFunctionHandle handle) const;
  const HenyeyGreensteinPhaseFunction& GetHenyeyGreenstein(
    base::PhaseFunctionHandle handle
  ) const;
  bool IsValid(base::PhaseFunctionHandle handle) const;
  std::size_t Size() const;

private:
  std::vector<HenyeyGreensteinPhaseFunction> phases_;
  std::vector<base::PhaseFunctionHandle::GenerationType> generations_;
};

class MediumTable {
public:
  base::MediumHandle Add(HomogeneousMedium medium);
  Medium Get(base::MediumHandle handle) const;
  const HomogeneousMedium& GetHomogeneous(base::MediumHandle handle) const;
  bool IsValid(base::MediumHandle handle) const;
  std::size_t Size() const;

private:
  std::vector<HomogeneousMedium> media_;
  std::vector<base::MediumHandle::GenerationType> generations_;
};

} // namespace render
} // namespace rayrender

#endif
