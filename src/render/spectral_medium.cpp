#include "spectral_medium.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace rayrender {
namespace render {
namespace {

constexpr Float Pi = static_cast<Float>(3.14159265358979323846264338327950288);
constexpr Float Inv4Pi = static_cast<Float>(1) / (static_cast<Float>(4) * Pi);

Float Clamp(Float value, Float low, Float high) {
  return std::min(std::max(value, low), high);
}

Float Sqr(Float value) {
  return value * value;
}

Float SafeSqrt(Float value) {
  return std::sqrt(std::max(static_cast<Float>(0), value));
}

Float OneMinusEpsilon() {
  return std::nextafter(static_cast<Float>(1), static_cast<Float>(0));
}

Float HenyeyGreenstein(Float cosTheta, Float g) {
  g = Clamp(g, static_cast<Float>(-0.99), static_cast<Float>(0.99));
  Float denom = static_cast<Float>(1) + Sqr(g) + static_cast<Float>(2) * g * cosTheta;
  return Inv4Pi * (static_cast<Float>(1) - Sqr(g)) / (denom * SafeSqrt(denom));
}

vec3f UnitOrFallback(const vec3f& value, const vec3f& fallback) {
  Float length = value.length();
  if (!(length > 0) || !std::isfinite(length)) {
    return fallback;
  }
  return value / length;
}

vec3f OrthogonalFallback(const vec3f& z) {
  if (std::abs(z.xyz.x) > std::abs(z.xyz.y)) {
    return UnitOrFallback(vec3f(-z.xyz.z, 0, z.xyz.x), vec3f(1, 0, 0));
  }
  return UnitOrFallback(vec3f(0, z.xyz.z, -z.xyz.y), vec3f(1, 0, 0));
}

Float RayDirectionLength(const Ray& ray) {
  Float length = ray.direction().length();
  return length > 0 && std::isfinite(length) ? length : 0;
}

Float EffectiveTMax(const Ray& ray, Float tMax) {
  Float effective = std::min(tMax, ray.tMax);
  if (effective < 0) {
    return 0;
  }
  return effective;
}

base::SampledSpectrum SafeExpNeg(const base::SampledSpectrum& sigmaT, Float distance) {
  base::SampledSpectrum result;
  for (int i = 0; i < base::NSpectrumSamples; ++i) {
    if (distance == std::numeric_limits<Float>::infinity()) {
      result[i] = sigmaT[i] > 0 ? 0 : 1;
    } else {
      result[i] = std::exp(-sigmaT[i] * distance);
    }
  }
  return result;
}

void ValidateFiniteNonNegative(Float value, const char* name) {
  if (!std::isfinite(value) || value < 0) {
    throw std::invalid_argument(std::string(name) + " must be finite and non-negative");
  }
}

void ValidateNonNegativeSpectrum(const base::Spectrum& spectrum, const char* name) {
  if (!spectrum.IsValid()) {
    return;
  }
  const Float samples[] = {
    base::LambdaMin,
    static_cast<Float>(550),
    base::LambdaMax
  };
  for (Float lambda : samples) {
    Float value = spectrum(lambda);
    if (!std::isfinite(value) || value < 0) {
      throw std::invalid_argument(std::string(name) + " must be non-negative over visible wavelengths");
    }
  }
  Float maxValue = spectrum.MaxValue();
  if (!std::isfinite(maxValue)) {
    throw std::invalid_argument(std::string(name) + " maximum must be finite");
  }
}

template <typename Handle>
bool IsValidHandle(
  Handle handle,
  std::size_t size,
  const std::vector<typename Handle::GenerationType>& generations
) {
  return handle.IsValid() && handle.Index() < size &&
         generations[handle.Index()] == handle.Generation();
}

} // namespace

HenyeyGreensteinPhaseFunction::HenyeyGreensteinPhaseFunction(Float g) : g_(g) {
  if (!(g_ > -1 && g_ < 1) || !std::isfinite(g_)) {
    throw std::invalid_argument("Henyey-Greenstein g must be finite and in (-1, 1)");
  }
}

Float HenyeyGreensteinPhaseFunction::p(const vec3f& wo, const vec3f& wi) const {
  return HenyeyGreenstein(dot(wo, wi), g_);
}

std::optional<PhaseFunctionSample> HenyeyGreensteinPhaseFunction::Sample_p(
  const vec3f& wo,
  point2f u
) const {
  Float g = Clamp(g_, static_cast<Float>(-0.99), static_cast<Float>(0.99));
  u.e[0] = Clamp(u[0], 0, OneMinusEpsilon());
  u.e[1] = Clamp(u[1], 0, OneMinusEpsilon());

  Float cosTheta;
  if (std::abs(g) < static_cast<Float>(1e-3)) {
    cosTheta = static_cast<Float>(1) - static_cast<Float>(2) * u[0];
  } else {
    Float sqrTerm =
      (static_cast<Float>(1) - Sqr(g)) /
      (static_cast<Float>(1) + g - static_cast<Float>(2) * g * u[0]);
    cosTheta =
      -static_cast<Float>(1) / (static_cast<Float>(2) * g) *
      (static_cast<Float>(1) + Sqr(g) - Sqr(sqrTerm));
  }

  Float sinTheta = SafeSqrt(static_cast<Float>(1) - Sqr(cosTheta));
  Float phi = static_cast<Float>(2) * Pi * u[1];
  vec3f z = UnitOrFallback(wo, vec3f(0, 0, 1));
  vec3f x = OrthogonalFallback(z);
  vec3f y = cross(z, x);
  vec3f wi =
    sinTheta * std::cos(phi) * x +
    sinTheta * std::sin(phi) * y +
    cosTheta * z;
  Float pdf = HenyeyGreenstein(cosTheta, g);
  if (!(pdf > 0) || !std::isfinite(pdf) || wi.squared_length() == 0) {
    return {};
  }
  return PhaseFunctionSample{pdf, wi, pdf};
}

Float HenyeyGreensteinPhaseFunction::PDF(const vec3f& wo, const vec3f& wi) const {
  return p(wo, wi);
}

Float HenyeyGreensteinPhaseFunction::G() const {
  return g_;
}

PhaseFunction::PhaseFunction(const HenyeyGreensteinPhaseFunction* phase) {
  if (phase != nullptr) {
    kind_ = Kind::HenyeyGreenstein;
    ptr_ = phase;
  }
}

PhaseFunction::operator bool() const {
  return ptr_ != nullptr;
}

Float PhaseFunction::p(const vec3f& wo, const vec3f& wi) const {
  if (kind_ == Kind::HenyeyGreenstein) {
    return static_cast<const HenyeyGreensteinPhaseFunction*>(ptr_)->p(wo, wi);
  }
  return 0;
}

std::optional<PhaseFunctionSample> PhaseFunction::Sample_p(
  const vec3f& wo,
  point2f u
) const {
  if (kind_ == Kind::HenyeyGreenstein) {
    return static_cast<const HenyeyGreensteinPhaseFunction*>(ptr_)->Sample_p(wo, u);
  }
  return {};
}

Float PhaseFunction::PDF(const vec3f& wo, const vec3f& wi) const {
  if (kind_ == Kind::HenyeyGreenstein) {
    return static_cast<const HenyeyGreensteinPhaseFunction*>(ptr_)->PDF(wo, wi);
  }
  return 0;
}

HomogeneousMajorantIterator::HomogeneousMajorantIterator(
  Float tMin,
  Float tMax,
  base::SampledSpectrum sigmaMaj
) : segment_{tMin, tMax, sigmaMaj}, called_(false) {}

std::optional<RayMajorantSegment> HomogeneousMajorantIterator::Next() {
  if (called_) {
    return {};
  }
  called_ = true;
  return segment_;
}

RayMajorantIterator::RayMajorantIterator(HomogeneousMajorantIterator iterator)
  : homogeneous_(std::move(iterator)) {}

std::optional<RayMajorantSegment> RayMajorantIterator::Next() {
  if (!homogeneous_) {
    return {};
  }
  return homogeneous_->Next();
}

base::SampledSpectrum MediumProperties::SigmaT() const {
  return sigmaA + sigmaS;
}

bool MediumProperties::HasScattering() const {
  return static_cast<bool>(sigmaS) && static_cast<bool>(phase);
}

HomogeneousMedium::HomogeneousMedium(
  base::Spectrum sigmaA,
  base::Spectrum sigmaS,
  Float scale,
  base::Spectrum emission,
  Float emissionScale,
  Float g
) : sigmaA_(std::move(sigmaA)),
    sigmaS_(std::move(sigmaS)),
    emission_(std::move(emission)),
    scale_(scale),
    emissionScale_(emissionScale),
    phase_(g) {
  ValidateFiniteNonNegative(scale_, "homogeneous medium scale");
  ValidateFiniteNonNegative(emissionScale_, "homogeneous medium emission scale");
  ValidateNonNegativeSpectrum(sigmaA_, "homogeneous medium sigma_a");
  ValidateNonNegativeSpectrum(sigmaS_, "homogeneous medium sigma_s");
  ValidateNonNegativeSpectrum(emission_, "homogeneous medium emission");
}

bool HomogeneousMedium::IsEmissive() const {
  return emissionScale_ > 0 && emission_.MaxValue() > 0;
}

MediumProperties HomogeneousMedium::SamplePoint(
  const point3f&,
  const base::SampledWavelengths& lambda
) const {
  MediumProperties properties;
  properties.sigmaA = sigmaA_.Sample(lambda) * scale_;
  properties.sigmaS = sigmaS_.Sample(lambda) * scale_;
  properties.phase = PhaseFunction(&phase_);
  properties.emission = emission_.Sample(lambda) * emissionScale_;
  return properties;
}

RayMajorantIterator HomogeneousMedium::SampleRay(
  const Ray& ray,
  Float tMax,
  const base::SampledWavelengths& lambda
) const {
  Float effectiveTMax = EffectiveTMax(ray, tMax);
  MediumProperties properties = SamplePoint(ray.origin(), lambda);
  base::SampledSpectrum sigmaMaj = properties.SigmaT() * RayDirectionLength(ray);
  return RayMajorantIterator(HomogeneousMajorantIterator(0, effectiveTMax, sigmaMaj));
}

base::SampledSpectrum HomogeneousMedium::Tr(
  const Ray& ray,
  Float tMax,
  const base::SampledWavelengths& lambda
) const {
  Float effectiveTMax = EffectiveTMax(ray, tMax);
  Float rayLength = RayDirectionLength(ray);
  if (effectiveTMax == 0 || rayLength == 0) {
    return base::SampledSpectrum(1);
  }
  Float distance = effectiveTMax == std::numeric_limits<Float>::infinity()
                     ? std::numeric_limits<Float>::infinity()
                     : effectiveTMax * rayLength;
  MediumProperties properties = SamplePoint(ray.origin(), lambda);
  return SafeExpNeg(properties.SigmaT(), distance);
}

std::optional<MediumSample> HomogeneousMedium::SampleDistance(
  const Ray& ray,
  Float tMax,
  const base::SampledWavelengths& lambda,
  Float uChannel,
  Float uDistance
) const {
  Float effectiveTMax = EffectiveTMax(ray, tMax);
  Float rayLength = RayDirectionLength(ray);
  if (rayLength == 0) {
    return {};
  }
  MediumProperties properties = SamplePoint(ray.origin(), lambda);
  base::SampledSpectrum sigmaT = properties.SigmaT();
  int channel = std::min(
    static_cast<int>(Clamp(uChannel, 0, OneMinusEpsilon()) * base::NSpectrumSamples),
    base::NSpectrumSamples - 1
  );
  Float sigmaTChannel = sigmaT[channel];
  Float distance = std::numeric_limits<Float>::infinity();
  if (sigmaTChannel > 0) {
    Float u = Clamp(uDistance, 0, OneMinusEpsilon());
    distance = -std::log(static_cast<Float>(1) - u) / sigmaTChannel;
  }

  Float sampledT = distance / rayLength;
  bool sampledMedium = sampledT < effectiveTMax;
  Float t = sampledMedium ? sampledT : effectiveTMax;
  base::SampledSpectrum tr = Tr(ray, t, lambda);
  base::SampledSpectrum density = sampledMedium ? sigmaT * tr : tr;
  Float pdf = density.Average();
  if (!(pdf > 0) || !std::isfinite(pdf)) {
    return {};
  }

  MediumSample sample;
  sample.transmittance = tr;
  sample.p = ray(t);
  sample.t = t;
  sample.pdf = pdf;
  sample.phase = properties.phase;
  sample.sampledMedium = sampledMedium;
  sample.beta = sampledMedium ? tr * properties.sigmaS / pdf : tr / pdf;
  return sample;
}

const base::Spectrum& HomogeneousMedium::SigmaASpectrum() const {
  return sigmaA_;
}

const base::Spectrum& HomogeneousMedium::SigmaSSpectrum() const {
  return sigmaS_;
}

const base::Spectrum& HomogeneousMedium::EmissionSpectrum() const {
  return emission_;
}

Float HomogeneousMedium::Scale() const {
  return scale_;
}

Float HomogeneousMedium::EmissionScale() const {
  return emissionScale_;
}

const HenyeyGreensteinPhaseFunction& HomogeneousMedium::Phase() const {
  return phase_;
}

Medium::Medium(const HomogeneousMedium* medium) {
  if (medium != nullptr) {
    kind_ = Kind::Homogeneous;
    ptr_ = medium;
  }
}

Medium::operator bool() const {
  return ptr_ != nullptr;
}

bool Medium::IsEmissive() const {
  if (kind_ == Kind::Homogeneous) {
    return static_cast<const HomogeneousMedium*>(ptr_)->IsEmissive();
  }
  return false;
}

MediumProperties Medium::SamplePoint(
  const point3f& p,
  const base::SampledWavelengths& lambda
) const {
  if (kind_ == Kind::Homogeneous) {
    return static_cast<const HomogeneousMedium*>(ptr_)->SamplePoint(p, lambda);
  }
  return {};
}

RayMajorantIterator Medium::SampleRay(
  const Ray& ray,
  Float tMax,
  const base::SampledWavelengths& lambda
) const {
  if (kind_ == Kind::Homogeneous) {
    return static_cast<const HomogeneousMedium*>(ptr_)->SampleRay(ray, tMax, lambda);
  }
  return {};
}

base::SampledSpectrum Medium::Tr(
  const Ray& ray,
  Float tMax,
  const base::SampledWavelengths& lambda
) const {
  if (kind_ == Kind::Homogeneous) {
    return static_cast<const HomogeneousMedium*>(ptr_)->Tr(ray, tMax, lambda);
  }
  return base::SampledSpectrum(1);
}

std::optional<MediumSample> Medium::SampleDistance(
  const Ray& ray,
  Float tMax,
  const base::SampledWavelengths& lambda,
  Float uChannel,
  Float uDistance
) const {
  if (kind_ == Kind::Homogeneous) {
    return static_cast<const HomogeneousMedium*>(ptr_)->SampleDistance(
      ray,
      tMax,
      lambda,
      uChannel,
      uDistance
    );
  }
  return {};
}

base::PhaseFunctionHandle PhaseFunctionTable::Add(HenyeyGreensteinPhaseFunction phase) {
  base::PhaseFunctionHandle handle = base::PhaseFunctionHandle::FromIndex(
    static_cast<base::PhaseFunctionHandle::IndexType>(phases_.size()),
    1
  );
  phases_.push_back(std::move(phase));
  generations_.push_back(handle.Generation());
  return handle;
}

PhaseFunction PhaseFunctionTable::Get(base::PhaseFunctionHandle handle) const {
  return PhaseFunction(&GetHenyeyGreenstein(handle));
}

const HenyeyGreensteinPhaseFunction& PhaseFunctionTable::GetHenyeyGreenstein(
  base::PhaseFunctionHandle handle
) const {
  if (!IsValid(handle)) {
    throw std::out_of_range("invalid phase function handle");
  }
  return phases_[handle.Index()];
}

bool PhaseFunctionTable::IsValid(base::PhaseFunctionHandle handle) const {
  return IsValidHandle(handle, phases_.size(), generations_);
}

std::size_t PhaseFunctionTable::Size() const {
  return phases_.size();
}

base::MediumHandle MediumTable::Add(HomogeneousMedium medium) {
  base::MediumHandle handle = base::MediumHandle::FromIndex(
    static_cast<base::MediumHandle::IndexType>(media_.size()),
    1
  );
  media_.push_back(std::move(medium));
  generations_.push_back(handle.Generation());
  return handle;
}

Medium MediumTable::Get(base::MediumHandle handle) const {
  return Medium(&GetHomogeneous(handle));
}

const HomogeneousMedium& MediumTable::GetHomogeneous(base::MediumHandle handle) const {
  if (!IsValid(handle)) {
    throw std::out_of_range("invalid medium handle");
  }
  return media_[handle.Index()];
}

bool MediumTable::IsValid(base::MediumHandle handle) const {
  return IsValidHandle(handle, media_.size(), generations_);
}

std::size_t MediumTable::Size() const {
  return media_.size();
}

} // namespace render
} // namespace rayrender
