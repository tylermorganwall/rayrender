#include "spectral_bsdf.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>

namespace rayrender {
namespace render {

namespace {

constexpr Float Pi = static_cast<Float>(3.14159265358979323846264338327950288);
constexpr Float InvPi = static_cast<Float>(1) / Pi;

Float Clamp(Float value, Float low, Float high) {
  return std::min(std::max(value, low), high);
}

Float Sqr(Float value) {
  return value * value;
}

Float SafeSqrt(Float value) {
  return std::sqrt(std::max(static_cast<Float>(0), value));
}

Float Lerp(Float t, Float a, Float b) {
  return (static_cast<Float>(1) - t) * a + t * b;
}

bool IsInf(Float value) {
  return std::isinf(value);
}

Float AbsDot(const vec3f& v, const vec3f& w) {
  return std::abs(dot(v, w));
}

vec3f UnitOrFallback(const vec3f& v, const vec3f& fallback) {
  Float length = v.length();
  if (length == 0 || !std::isfinite(length)) {
    return fallback;
  }
  return v / length;
}

vec3f OrthogonalFallback(const vec3f& z) {
  if (std::abs(z.xyz.x) > std::abs(z.xyz.y)) {
    return UnitOrFallback(vec3f(-z.xyz.z, 0, z.xyz.x), vec3f(1, 0, 0));
  }
  return UnitOrFallback(vec3f(0, z.xyz.z, -z.xyz.y), vec3f(1, 0, 0));
}

vec3f FaceForward(vec3f v, const vec3f& reference) {
  return dot(v, reference) < 0 ? -v : v;
}

} // namespace

bool BSDFSample::IsReflection() const {
  return base::HasFlag(flags, BxDFFlags::Reflection);
}

bool BSDFSample::IsTransmission() const {
  return base::HasFlag(flags, BxDFFlags::Transmission);
}

bool BSDFSample::IsDiffuse() const {
  return base::HasFlag(flags, BxDFFlags::Diffuse);
}

bool BSDFSample::IsGlossy() const {
  return base::HasFlag(flags, BxDFFlags::Glossy);
}

bool BSDFSample::IsSpecular() const {
  return base::HasFlag(flags, BxDFFlags::Specular);
}

Float CosTheta(const vec3f& w) {
  return w.xyz.z;
}

Float Cos2Theta(const vec3f& w) {
  return w.xyz.z * w.xyz.z;
}

Float AbsCosTheta(const vec3f& w) {
  return std::abs(w.xyz.z);
}

Float Sin2Theta(const vec3f& w) {
  return std::max(static_cast<Float>(0), static_cast<Float>(1) - Cos2Theta(w));
}

Float Tan2Theta(const vec3f& w) {
  Float cos2Theta = Cos2Theta(w);
  if (cos2Theta == 0) {
    return std::numeric_limits<Float>::infinity();
  }
  return Sin2Theta(w) / cos2Theta;
}

Float CosPhi(const vec3f& w) {
  Float sin2Theta = Sin2Theta(w);
  if (sin2Theta == 0) {
    return 1;
  }
  return Clamp(w.xyz.x / std::sqrt(sin2Theta), -1, 1);
}

Float SinPhi(const vec3f& w) {
  Float sin2Theta = Sin2Theta(w);
  if (sin2Theta == 0) {
    return 0;
  }
  return Clamp(w.xyz.y / std::sqrt(sin2Theta), -1, 1);
}

bool SameHemisphere(const vec3f& w, const vec3f& wp) {
  return w.xyz.z * wp.xyz.z > 0;
}

vec3f Reflect(const vec3f& wo, const vec3f& n) {
  return -wo + static_cast<Float>(2) * dot(wo, n) * n;
}

bool Refract(const vec3f& wi, normal3f n, Float eta, Float* etap, vec3f* wt) {
  Float cosThetaI = dot(convert_to_vec3(n), wi);
  if (cosThetaI < 0) {
    eta = static_cast<Float>(1) / eta;
    cosThetaI = -cosThetaI;
    n = -n;
  }

  Float sin2ThetaI = std::max(static_cast<Float>(0), static_cast<Float>(1) - Sqr(cosThetaI));
  Float sin2ThetaT = sin2ThetaI / Sqr(eta);
  if (sin2ThetaT >= 1) {
    return false;
  }
  Float cosThetaT = SafeSqrt(static_cast<Float>(1) - sin2ThetaT);
  *wt = -wi / eta + (cosThetaI / eta - cosThetaT) * convert_to_vec3(n);
  if (etap != nullptr) {
    *etap = eta;
  }
  return true;
}

vec3f SampleUniformDiskPolar(point2f u) {
  Float r = std::sqrt(u[0]);
  Float theta = static_cast<Float>(2) * Pi * u[1];
  return vec3f(r * std::cos(theta), r * std::sin(theta), 0);
}

vec3f SampleCosineHemisphere(point2f u) {
  vec3f d = SampleUniformDiskPolar(u);
  Float z = SafeSqrt(static_cast<Float>(1) - d.xyz.x * d.xyz.x - d.xyz.y * d.xyz.y);
  return vec3f(d.xyz.x, d.xyz.y, z);
}

Float CosineHemispherePDF(Float cosTheta) {
  return cosTheta * InvPi;
}

Float FrDielectric(Float cosThetaI, Float eta) {
  cosThetaI = Clamp(cosThetaI, -1, 1);
  if (cosThetaI < 0) {
    eta = static_cast<Float>(1) / eta;
    cosThetaI = -cosThetaI;
  }

  Float sin2ThetaI = static_cast<Float>(1) - Sqr(cosThetaI);
  Float sin2ThetaT = sin2ThetaI / Sqr(eta);
  if (sin2ThetaT >= 1) {
    return 1;
  }
  Float cosThetaT = SafeSqrt(static_cast<Float>(1) - sin2ThetaT);

  Float rParallel = (eta * cosThetaI - cosThetaT) / (eta * cosThetaI + cosThetaT);
  Float rPerpendicular = (cosThetaI - eta * cosThetaT) / (cosThetaI + eta * cosThetaT);
  return (Sqr(rParallel) + Sqr(rPerpendicular)) / static_cast<Float>(2);
}

Float FrComplex(Float cosThetaI, Float eta, Float k) {
  using Complex = std::complex<Float>;

  cosThetaI = Clamp(cosThetaI, 0, 1);
  Float sin2ThetaI = static_cast<Float>(1) - Sqr(cosThetaI);
  Complex etaComplex(eta, k);
  Complex sin2ThetaT = sin2ThetaI / (etaComplex * etaComplex);
  Complex cosThetaT = std::sqrt(Complex(1, 0) - sin2ThetaT);

  Complex rParallel =
    (etaComplex * cosThetaI - cosThetaT) / (etaComplex * cosThetaI + cosThetaT);
  Complex rPerpendicular =
    (cosThetaI - etaComplex * cosThetaT) / (cosThetaI + etaComplex * cosThetaT);
  return static_cast<Float>((std::norm(rParallel) + std::norm(rPerpendicular)) / 2);
}

base::SampledSpectrum FrComplex(
  Float cosThetaI,
  const base::SampledSpectrum& eta,
  const base::SampledSpectrum& k
) {
  base::SampledSpectrum result;
  for (int i = 0; i < base::NSpectrumSamples; ++i) {
    result[i] = FrComplex(cosThetaI, eta[i], k[i]);
  }
  return result;
}

TrowbridgeReitzDistribution::TrowbridgeReitzDistribution(Float alphaX, Float alphaY)
  : alphaX_(alphaX), alphaY_(alphaY) {
  if (!EffectivelySmooth()) {
    alphaX_ = std::max(alphaX_, static_cast<Float>(1e-4));
    alphaY_ = std::max(alphaY_, static_cast<Float>(1e-4));
  }
}

Float TrowbridgeReitzDistribution::D(const vec3f& wm) const {
  if (alphaX_ == 0 || alphaY_ == 0) {
    return 0;
  }
  Float tan2Theta = Tan2Theta(wm);
  if (IsInf(tan2Theta)) {
    return 0;
  }
  Float cos4Theta = Sqr(Cos2Theta(wm));
  if (cos4Theta < static_cast<Float>(1e-16)) {
    return 0;
  }
  Float e = tan2Theta * (Sqr(CosPhi(wm) / alphaX_) + Sqr(SinPhi(wm) / alphaY_));
  return static_cast<Float>(1) /
         (Pi * alphaX_ * alphaY_ * cos4Theta * Sqr(static_cast<Float>(1) + e));
}

Float TrowbridgeReitzDistribution::G1(const vec3f& w) const {
  return static_cast<Float>(1) / (static_cast<Float>(1) + Lambda(w));
}

Float TrowbridgeReitzDistribution::Lambda(const vec3f& w) const {
  Float tan2Theta = Tan2Theta(w);
  if (IsInf(tan2Theta)) {
    return 0;
  }
  Float alpha2 = Sqr(CosPhi(w) * alphaX_) + Sqr(SinPhi(w) * alphaY_);
  return (std::sqrt(static_cast<Float>(1) + alpha2 * tan2Theta) - static_cast<Float>(1)) /
         static_cast<Float>(2);
}

Float TrowbridgeReitzDistribution::G(const vec3f& wo, const vec3f& wi) const {
  return static_cast<Float>(1) /
         (static_cast<Float>(1) + Lambda(wo) + Lambda(wi));
}

Float TrowbridgeReitzDistribution::D(const vec3f& w, const vec3f& wm) const {
  Float absCosTheta = AbsCosTheta(w);
  if (absCosTheta == 0) {
    return 0;
  }
  return G1(w) / absCosTheta * D(wm) * std::abs(dot(w, wm));
}

Float TrowbridgeReitzDistribution::PDF(const vec3f& w, const vec3f& wm) const {
  return D(w, wm);
}

vec3f TrowbridgeReitzDistribution::Sample_wm(const vec3f& w, point2f u) const {
  if (EffectivelySmooth()) {
    return vec3f(0, 0, 1);
  }

  vec3f wh = UnitOrFallback(vec3f(alphaX_ * w.xyz.x, alphaY_ * w.xyz.y, w.xyz.z), vec3f(0, 0, 1));
  if (wh.xyz.z < 0) {
    wh = -wh;
  }

  vec3f t1 = wh.xyz.z < static_cast<Float>(0.99999)
                ? UnitOrFallback(cross(vec3f(0, 0, 1), wh), vec3f(1, 0, 0))
                : vec3f(1, 0, 0);
  vec3f t2 = cross(wh, t1);

  vec3f p = SampleUniformDiskPolar(u);
  Float h = SafeSqrt(static_cast<Float>(1) - Sqr(p.xyz.x));
  p.xyz.y = Lerp((static_cast<Float>(1) + wh.xyz.z) / static_cast<Float>(2), h, p.xyz.y);

  Float pz = SafeSqrt(static_cast<Float>(1) - Sqr(p.xyz.x) - Sqr(p.xyz.y));
  vec3f nh = p.xyz.x * t1 + p.xyz.y * t2 + pz * wh;
  return UnitOrFallback(
    vec3f(alphaX_ * nh.xyz.x, alphaY_ * nh.xyz.y, std::max(static_cast<Float>(1e-6), nh.xyz.z)),
    vec3f(0, 0, 1)
  );
}

bool TrowbridgeReitzDistribution::EffectivelySmooth() const {
  return std::max(alphaX_, alphaY_) < static_cast<Float>(1e-3);
}

void TrowbridgeReitzDistribution::Regularize() {
  if (alphaX_ < static_cast<Float>(0.3)) {
    alphaX_ = Clamp(static_cast<Float>(2) * alphaX_, static_cast<Float>(0.1), static_cast<Float>(0.3));
  }
  if (alphaY_ < static_cast<Float>(0.3)) {
    alphaY_ = Clamp(static_cast<Float>(2) * alphaY_, static_cast<Float>(0.1), static_cast<Float>(0.3));
  }
}

Float TrowbridgeReitzDistribution::AlphaX() const {
  return alphaX_;
}

Float TrowbridgeReitzDistribution::AlphaY() const {
  return alphaY_;
}

Float TrowbridgeReitzDistribution::RoughnessToAlpha(Float roughness) {
  return std::sqrt(roughness);
}

DiffuseBxDF::DiffuseBxDF(base::SampledSpectrum reflectance) : reflectance_(reflectance) {}

BxDFFlags DiffuseBxDF::Flags() const {
  return reflectance_ ? BxDFFlags::DiffuseReflection : BxDFFlags::Unset;
}

base::SampledSpectrum DiffuseBxDF::f(
  const vec3f& wo,
  const vec3f& wi,
  TransportMode mode
) const {
  (void)mode;
  if (!SameHemisphere(wo, wi)) {
    return base::SampledSpectrum(0);
  }
  return reflectance_ * InvPi;
}

std::optional<BSDFSample> DiffuseBxDF::Sample_f(
  const vec3f& wo,
  Float uc,
  point2f u,
  TransportMode mode,
  BxDFReflTransFlags sampleFlags
) const {
  (void)uc;
  (void)mode;
  if (!base::HasFlag(sampleFlags, BxDFReflTransFlags::Reflection)) {
    return {};
  }

  vec3f wi = SampleCosineHemisphere(u);
  if (wo.xyz.z < 0) {
    wi.xyz.z *= -1;
  }
  Float pdf = CosineHemispherePDF(AbsCosTheta(wi));
  return BSDFSample{reflectance_ * InvPi, wi, pdf, BxDFFlags::DiffuseReflection};
}

Float DiffuseBxDF::PDF(
  const vec3f& wo,
  const vec3f& wi,
  TransportMode mode,
  BxDFReflTransFlags sampleFlags
) const {
  (void)mode;
  if (!base::HasFlag(sampleFlags, BxDFReflTransFlags::Reflection) || !SameHemisphere(wo, wi)) {
    return 0;
  }
  return CosineHemispherePDF(AbsCosTheta(wi));
}

base::SampledSpectrum DiffuseBxDF::rho() const {
  return reflectance_;
}

void DiffuseBxDF::Regularize() {}

ConductorBxDF::ConductorBxDF(
  TrowbridgeReitzDistribution distribution,
  base::SampledSpectrum eta,
  base::SampledSpectrum k
)
  : distribution_(distribution),
    eta_(eta),
    k_(k) {}

BxDFFlags ConductorBxDF::Flags() const {
  return distribution_.EffectivelySmooth()
           ? BxDFFlags::SpecularReflection
           : BxDFFlags::GlossyReflection;
}

base::SampledSpectrum ConductorBxDF::f(
  const vec3f& wo,
  const vec3f& wi,
  TransportMode mode
) const {
  (void)mode;
  if (!SameHemisphere(wo, wi) || distribution_.EffectivelySmooth()) {
    return base::SampledSpectrum(0);
  }

  Float cosThetaO = AbsCosTheta(wo);
  Float cosThetaI = AbsCosTheta(wi);
  if (cosThetaI == 0 || cosThetaO == 0) {
    return base::SampledSpectrum(0);
  }

  vec3f wm = wi + wo;
  if (wm.squared_length() == 0) {
    return base::SampledSpectrum(0);
  }
  wm = UnitOrFallback(wm, vec3f(0, 0, 1));

  base::SampledSpectrum F = FrComplex(AbsDot(wo, wm), eta_, k_);
  return distribution_.D(wm) * F * distribution_.G(wo, wi) /
         (static_cast<Float>(4) * cosThetaI * cosThetaO);
}

std::optional<BSDFSample> ConductorBxDF::Sample_f(
  const vec3f& wo,
  Float uc,
  point2f u,
  TransportMode mode,
  BxDFReflTransFlags sampleFlags
) const {
  (void)uc;
  (void)mode;
  if (!base::HasFlag(sampleFlags, BxDFReflTransFlags::Reflection)) {
    return {};
  }

  if (distribution_.EffectivelySmooth()) {
    vec3f wi(-wo.xyz.x, -wo.xyz.y, wo.xyz.z);
    Float absCosTheta = AbsCosTheta(wi);
    if (absCosTheta == 0) {
      return {};
    }
    base::SampledSpectrum f = FrComplex(absCosTheta, eta_, k_) / absCosTheta;
    return BSDFSample{f, wi, 1, BxDFFlags::SpecularReflection};
  }

  if (wo.xyz.z == 0) {
    return {};
  }
  vec3f wm = distribution_.Sample_wm(wo, u);
  vec3f wi = Reflect(wo, wm);
  if (!SameHemisphere(wo, wi)) {
    return {};
  }

  Float absDotWoWm = AbsDot(wo, wm);
  if (absDotWoWm == 0) {
    return {};
  }
  Float pdf = distribution_.PDF(wo, wm) /
              (static_cast<Float>(4) * absDotWoWm);
  Float cosThetaO = AbsCosTheta(wo);
  Float cosThetaI = AbsCosTheta(wi);
  if (cosThetaI == 0 || cosThetaO == 0 || pdf == 0) {
    return {};
  }

  base::SampledSpectrum F = FrComplex(absDotWoWm, eta_, k_);
  base::SampledSpectrum f =
    distribution_.D(wm) * F * distribution_.G(wo, wi) /
    (static_cast<Float>(4) * cosThetaI * cosThetaO);
  return BSDFSample{f, wi, pdf, BxDFFlags::GlossyReflection};
}

Float ConductorBxDF::PDF(
  const vec3f& wo,
  const vec3f& wi,
  TransportMode mode,
  BxDFReflTransFlags sampleFlags
) const {
  (void)mode;
  if (!base::HasFlag(sampleFlags, BxDFReflTransFlags::Reflection) ||
      !SameHemisphere(wo, wi) ||
      distribution_.EffectivelySmooth()) {
    return 0;
  }

  vec3f wm = wo + wi;
  if (wm.squared_length() == 0) {
    return 0;
  }
  wm = FaceForward(UnitOrFallback(wm, vec3f(0, 0, 1)), vec3f(0, 0, 1));
  Float absDotWoWm = AbsDot(wo, wm);
  if (absDotWoWm == 0) {
    return 0;
  }
  return distribution_.PDF(wo, wm) / (static_cast<Float>(4) * absDotWoWm);
}

base::SampledSpectrum ConductorBxDF::rho() const {
  return FrComplex(1, eta_, k_);
}

void ConductorBxDF::Regularize() {
  distribution_.Regularize();
}

DielectricBxDF::DielectricBxDF(Float eta, TrowbridgeReitzDistribution distribution)
  : eta_(eta), distribution_(distribution) {}

BxDFFlags DielectricBxDF::Flags() const {
  BxDFFlags flags = eta_ == 1
                       ? BxDFFlags::Transmission
                       : (BxDFFlags::Reflection | BxDFFlags::Transmission);
  return flags |
         (distribution_.EffectivelySmooth() ? BxDFFlags::Specular : BxDFFlags::Glossy);
}

base::SampledSpectrum DielectricBxDF::f(
  const vec3f& wo,
  const vec3f& wi,
  TransportMode mode
) const {
  (void)wo;
  (void)wi;
  (void)mode;
  return base::SampledSpectrum(0);
}

std::optional<BSDFSample> DielectricBxDF::Sample_f(
  const vec3f& wo,
  Float uc,
  point2f u,
  TransportMode mode,
  BxDFReflTransFlags sampleFlags
) const {
  (void)u;
  if (!distribution_.EffectivelySmooth()) {
    return {};
  }

  Float R = FrDielectric(CosTheta(wo), eta_);
  Float T = static_cast<Float>(1) - R;
  Float pr = R;
  Float pt = T;
  if (!base::HasFlag(sampleFlags, BxDFReflTransFlags::Reflection)) {
    pr = 0;
  }
  if (!base::HasFlag(sampleFlags, BxDFReflTransFlags::Transmission)) {
    pt = 0;
  }
  if (pr == 0 && pt == 0) {
    return {};
  }

  if (uc < pr / (pr + pt)) {
    vec3f wi(-wo.xyz.x, -wo.xyz.y, wo.xyz.z);
    Float absCosTheta = AbsCosTheta(wi);
    if (absCosTheta == 0) {
      return {};
    }
    return BSDFSample{
      base::SampledSpectrum(R / absCosTheta),
      wi,
      pr / (pr + pt),
      BxDFFlags::SpecularReflection
    };
  }

  vec3f wi;
  Float etap = 1;
  if (!Refract(wo, normal3f(0, 0, 1), eta_, &etap, &wi)) {
    return {};
  }
  Float absCosTheta = AbsCosTheta(wi);
  if (absCosTheta == 0) {
    return {};
  }
  base::SampledSpectrum ft(T / absCosTheta);
  if (mode == TransportMode::Radiance) {
    ft /= (etap * etap);
  }
  return BSDFSample{
    ft,
    wi,
    pt / (pr + pt),
    BxDFFlags::SpecularTransmission,
    etap
  };
}

Float DielectricBxDF::PDF(
  const vec3f& wo,
  const vec3f& wi,
  TransportMode mode,
  BxDFReflTransFlags sampleFlags
) const {
  (void)wo;
  (void)wi;
  (void)mode;
  (void)sampleFlags;
  return 0;
}

base::SampledSpectrum DielectricBxDF::rho() const {
  return base::SampledSpectrum(1);
}

void DielectricBxDF::Regularize() {
  distribution_.Regularize();
}

BxDFFlags NullBxDF::Flags() const {
  return BxDFFlags::Unset;
}

base::SampledSpectrum NullBxDF::f(
  const vec3f& wo,
  const vec3f& wi,
  TransportMode mode
) const {
  (void)wo;
  (void)wi;
  (void)mode;
  return base::SampledSpectrum(0);
}

std::optional<BSDFSample> NullBxDF::Sample_f(
  const vec3f& wo,
  Float uc,
  point2f u,
  TransportMode mode,
  BxDFReflTransFlags sampleFlags
) const {
  (void)wo;
  (void)uc;
  (void)u;
  (void)mode;
  (void)sampleFlags;
  return {};
}

Float NullBxDF::PDF(
  const vec3f& wo,
  const vec3f& wi,
  TransportMode mode,
  BxDFReflTransFlags sampleFlags
) const {
  (void)wo;
  (void)wi;
  (void)mode;
  (void)sampleFlags;
  return 0;
}

void NullBxDF::Regularize() {}

BxDF::BxDF(DiffuseBxDF* bxdf) : kind_(bxdf ? Kind::Diffuse : Kind::None), ptr_(bxdf) {}

BxDF::BxDF(ConductorBxDF* bxdf) : kind_(bxdf ? Kind::Conductor : Kind::None), ptr_(bxdf) {}

BxDF::BxDF(DielectricBxDF* bxdf) : kind_(bxdf ? Kind::Dielectric : Kind::None), ptr_(bxdf) {}

BxDF::BxDF(NullBxDF* bxdf) : kind_(bxdf ? Kind::Null : Kind::None), ptr_(bxdf) {}

BxDF::operator bool() const {
  return ptr_ != nullptr && kind_ != Kind::None;
}

BxDFFlags BxDF::Flags() const {
  switch (kind_) {
  case Kind::Diffuse:
    return static_cast<const DiffuseBxDF*>(ptr_)->Flags();
  case Kind::Conductor:
    return static_cast<const ConductorBxDF*>(ptr_)->Flags();
  case Kind::Dielectric:
    return static_cast<const DielectricBxDF*>(ptr_)->Flags();
  case Kind::Null:
    return static_cast<const NullBxDF*>(ptr_)->Flags();
  case Kind::None:
    return BxDFFlags::Unset;
  }
  return BxDFFlags::Unset;
}

base::SampledSpectrum BxDF::f(
  const vec3f& wo,
  const vec3f& wi,
  TransportMode mode
) const {
  switch (kind_) {
  case Kind::Diffuse:
    return static_cast<const DiffuseBxDF*>(ptr_)->f(wo, wi, mode);
  case Kind::Conductor:
    return static_cast<const ConductorBxDF*>(ptr_)->f(wo, wi, mode);
  case Kind::Dielectric:
    return static_cast<const DielectricBxDF*>(ptr_)->f(wo, wi, mode);
  case Kind::Null:
    return static_cast<const NullBxDF*>(ptr_)->f(wo, wi, mode);
  case Kind::None:
    return base::SampledSpectrum(0);
  }
  return base::SampledSpectrum(0);
}

std::optional<BSDFSample> BxDF::Sample_f(
  const vec3f& wo,
  Float uc,
  point2f u,
  TransportMode mode,
  BxDFReflTransFlags sampleFlags
) const {
  switch (kind_) {
  case Kind::Diffuse:
    return static_cast<const DiffuseBxDF*>(ptr_)->Sample_f(wo, uc, u, mode, sampleFlags);
  case Kind::Conductor:
    return static_cast<const ConductorBxDF*>(ptr_)->Sample_f(wo, uc, u, mode, sampleFlags);
  case Kind::Dielectric:
    return static_cast<const DielectricBxDF*>(ptr_)->Sample_f(wo, uc, u, mode, sampleFlags);
  case Kind::Null:
    return static_cast<const NullBxDF*>(ptr_)->Sample_f(wo, uc, u, mode, sampleFlags);
  case Kind::None:
    return {};
  }
  return {};
}

Float BxDF::PDF(
  const vec3f& wo,
  const vec3f& wi,
  TransportMode mode,
  BxDFReflTransFlags sampleFlags
) const {
  switch (kind_) {
  case Kind::Diffuse:
    return static_cast<const DiffuseBxDF*>(ptr_)->PDF(wo, wi, mode, sampleFlags);
  case Kind::Conductor:
    return static_cast<const ConductorBxDF*>(ptr_)->PDF(wo, wi, mode, sampleFlags);
  case Kind::Dielectric:
    return static_cast<const DielectricBxDF*>(ptr_)->PDF(wo, wi, mode, sampleFlags);
  case Kind::Null:
    return static_cast<const NullBxDF*>(ptr_)->PDF(wo, wi, mode, sampleFlags);
  case Kind::None:
    return 0;
  }
  return 0;
}

base::SampledSpectrum BxDF::rho() const {
  switch (kind_) {
  case Kind::Diffuse:
    return static_cast<const DiffuseBxDF*>(ptr_)->rho();
  case Kind::Conductor:
    return static_cast<const ConductorBxDF*>(ptr_)->rho();
  case Kind::Dielectric:
    return static_cast<const DielectricBxDF*>(ptr_)->rho();
  case Kind::Null:
  case Kind::None:
    return base::SampledSpectrum(0);
  }
  return base::SampledSpectrum(0);
}

void BxDF::Regularize() {
  switch (kind_) {
  case Kind::Diffuse:
    static_cast<DiffuseBxDF*>(ptr_)->Regularize();
    return;
  case Kind::Conductor:
    static_cast<ConductorBxDF*>(ptr_)->Regularize();
    return;
  case Kind::Dielectric:
    static_cast<DielectricBxDF*>(ptr_)->Regularize();
    return;
  case Kind::Null:
    static_cast<NullBxDF*>(ptr_)->Regularize();
    return;
  case Kind::None:
    return;
  }
}

BSDFFrame::BSDFFrame(normal3f ns, vec3f dpdus) {
  z_ = UnitOrFallback(convert_to_vec3(ns), vec3f(0, 0, 1));
  x_ = dpdus - dot(dpdus, z_) * z_;
  x_ = UnitOrFallback(x_, OrthogonalFallback(z_));
  y_ = cross(z_, x_);
}

vec3f BSDFFrame::ToLocal(const vec3f& v) const {
  return vec3f(dot(v, x_), dot(v, y_), dot(v, z_));
}

vec3f BSDFFrame::FromLocal(const vec3f& v) const {
  return v.xyz.x * x_ + v.xyz.y * y_ + v.xyz.z * z_;
}

const vec3f& BSDFFrame::X() const {
  return x_;
}

const vec3f& BSDFFrame::Y() const {
  return y_;
}

const vec3f& BSDFFrame::Z() const {
  return z_;
}

BSDF::BSDF(normal3f geometricNormal, normal3f shadingNormal, vec3f dpdus, BxDF bxdf)
  : bxdf_(bxdf), shadingFrame_(shadingNormal, dpdus) {
  geometricNormal_ =
    convert_to_normal3(UnitOrFallback(convert_to_vec3(geometricNormal), vec3f(0, 0, 1)));
}

BSDF::operator bool() const {
  return static_cast<bool>(bxdf_);
}

BxDFFlags BSDF::Flags() const {
  return bxdf_.Flags();
}

vec3f BSDF::RenderToLocal(const vec3f& v) const {
  return shadingFrame_.ToLocal(v);
}

vec3f BSDF::LocalToRender(const vec3f& v) const {
  return shadingFrame_.FromLocal(v);
}

base::SampledSpectrum BSDF::f(
  const vec3f& woRender,
  const vec3f& wiRender,
  TransportMode mode
) const {
  vec3f wi = RenderToLocal(wiRender);
  vec3f wo = RenderToLocal(woRender);
  if (wo.xyz.z == 0) {
    return base::SampledSpectrum(0);
  }
  return bxdf_.f(wo, wi, mode);
}

std::optional<BSDFSample> BSDF::Sample_f(
  const vec3f& woRender,
  Float uc,
  point2f u,
  TransportMode mode,
  BxDFReflTransFlags sampleFlags
) const {
  vec3f wo = RenderToLocal(woRender);
  if (wo.xyz.z == 0 || !base::HasAny(bxdf_.Flags(), sampleFlags)) {
    return {};
  }

  std::optional<BSDFSample> bs = bxdf_.Sample_f(wo, uc, u, mode, sampleFlags);
  if (!bs || !bs->f || bs->pdf == 0 || bs->wi.xyz.z == 0) {
    return {};
  }
  bs->wi = LocalToRender(bs->wi);
  return bs;
}

Float BSDF::PDF(
  const vec3f& woRender,
  const vec3f& wiRender,
  TransportMode mode,
  BxDFReflTransFlags sampleFlags
) const {
  vec3f wo = RenderToLocal(woRender);
  vec3f wi = RenderToLocal(wiRender);
  if (wo.xyz.z == 0) {
    return 0;
  }
  return bxdf_.PDF(wo, wi, mode, sampleFlags);
}

base::SampledSpectrum BSDF::rho() const {
  return bxdf_.rho();
}

void BSDF::Regularize() {
  bxdf_.Regularize();
}

const BSDFFrame& BSDF::Frame() const {
  return shadingFrame_;
}

normal3f BSDF::GeometricNormal() const {
  return geometricNormal_;
}

} // namespace render
} // namespace rayrender
