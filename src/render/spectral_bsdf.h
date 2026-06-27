#ifndef RAYRENDER_RENDER_SPECTRAL_BSDF_H
#define RAYRENDER_RENDER_SPECTRAL_BSDF_H

#include "../base/base.h"
#include "../math/vectypes.h"

#include <optional>

namespace rayrender {
namespace render {

using base::BxDFFlags;
using base::BxDFReflTransFlags;
using base::TransportMode;

struct BSDFSample {
  base::SampledSpectrum f;
  vec3f wi;
  Float pdf = 0;
  BxDFFlags flags = BxDFFlags::Unset;
  Float eta = 1;
  bool pdfIsProportional = false;

  bool IsReflection() const;
  bool IsTransmission() const;
  bool IsDiffuse() const;
  bool IsGlossy() const;
  bool IsSpecular() const;
};

Float CosTheta(const vec3f& w);
Float Cos2Theta(const vec3f& w);
Float AbsCosTheta(const vec3f& w);
Float Sin2Theta(const vec3f& w);
Float Tan2Theta(const vec3f& w);
Float CosPhi(const vec3f& w);
Float SinPhi(const vec3f& w);
bool SameHemisphere(const vec3f& w, const vec3f& wp);
vec3f Reflect(const vec3f& wo, const vec3f& n);

vec3f SampleUniformDiskPolar(point2f u);
vec3f SampleCosineHemisphere(point2f u);
Float CosineHemispherePDF(Float cosTheta);

Float FrDielectric(Float cosThetaI, Float eta);
Float FrComplex(Float cosThetaI, Float eta, Float k);
base::SampledSpectrum FrComplex(
  Float cosThetaI,
  const base::SampledSpectrum& eta,
  const base::SampledSpectrum& k
);

class TrowbridgeReitzDistribution {
public:
  TrowbridgeReitzDistribution() = default;
  TrowbridgeReitzDistribution(Float alphaX, Float alphaY);

  Float D(const vec3f& wm) const;
  Float G1(const vec3f& w) const;
  Float Lambda(const vec3f& w) const;
  Float G(const vec3f& wo, const vec3f& wi) const;
  Float D(const vec3f& w, const vec3f& wm) const;
  Float PDF(const vec3f& w, const vec3f& wm) const;
  vec3f Sample_wm(const vec3f& w, point2f u) const;
  bool EffectivelySmooth() const;
  void Regularize();

  Float AlphaX() const;
  Float AlphaY() const;
  static Float RoughnessToAlpha(Float roughness);

private:
  Float alphaX_ = 0;
  Float alphaY_ = 0;
};

class DiffuseBxDF {
public:
  DiffuseBxDF() = default;
  explicit DiffuseBxDF(base::SampledSpectrum reflectance);

  BxDFFlags Flags() const;
  base::SampledSpectrum f(const vec3f& wo, const vec3f& wi, TransportMode mode) const;
  std::optional<BSDFSample> Sample_f(
    const vec3f& wo,
    Float uc,
    point2f u,
    TransportMode mode,
    BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All
  ) const;
  Float PDF(
    const vec3f& wo,
    const vec3f& wi,
    TransportMode mode,
    BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All
  ) const;
  base::SampledSpectrum rho() const;
  void Regularize();

private:
  base::SampledSpectrum reflectance_;
};

class NullBxDF {
public:
  BxDFFlags Flags() const;
  base::SampledSpectrum f(const vec3f& wo, const vec3f& wi, TransportMode mode) const;
  std::optional<BSDFSample> Sample_f(
    const vec3f& wo,
    Float uc,
    point2f u,
    TransportMode mode,
    BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All
  ) const;
  Float PDF(
    const vec3f& wo,
    const vec3f& wi,
    TransportMode mode,
    BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All
  ) const;
  void Regularize();
};

class BxDF {
public:
  BxDF() = default;
  explicit BxDF(DiffuseBxDF* bxdf);
  explicit BxDF(NullBxDF* bxdf);

  explicit operator bool() const;
  BxDFFlags Flags() const;
  base::SampledSpectrum f(const vec3f& wo, const vec3f& wi, TransportMode mode) const;
  std::optional<BSDFSample> Sample_f(
    const vec3f& wo,
    Float uc,
    point2f u,
    TransportMode mode = TransportMode::Radiance,
    BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All
  ) const;
  Float PDF(
    const vec3f& wo,
    const vec3f& wi,
    TransportMode mode,
    BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All
  ) const;
  void Regularize();

private:
  enum class Kind {
    None,
    Diffuse,
    Null
  };

  Kind kind_ = Kind::None;
  void* ptr_ = nullptr;
};

class BSDFFrame {
public:
  BSDFFrame() = default;
  BSDFFrame(normal3f ns, vec3f dpdus);

  vec3f ToLocal(const vec3f& v) const;
  vec3f FromLocal(const vec3f& v) const;
  const vec3f& X() const;
  const vec3f& Y() const;
  const vec3f& Z() const;

private:
  vec3f x_ = vec3f(1, 0, 0);
  vec3f y_ = vec3f(0, 1, 0);
  vec3f z_ = vec3f(0, 0, 1);
};

class BSDF {
public:
  BSDF() = default;
  BSDF(normal3f geometricNormal, normal3f shadingNormal, vec3f dpdus, BxDF bxdf);

  explicit operator bool() const;
  BxDFFlags Flags() const;
  vec3f RenderToLocal(const vec3f& v) const;
  vec3f LocalToRender(const vec3f& v) const;
  base::SampledSpectrum f(
    const vec3f& woRender,
    const vec3f& wiRender,
    TransportMode mode = TransportMode::Radiance
  ) const;
  std::optional<BSDFSample> Sample_f(
    const vec3f& woRender,
    Float uc,
    point2f u,
    TransportMode mode = TransportMode::Radiance,
    BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All
  ) const;
  Float PDF(
    const vec3f& woRender,
    const vec3f& wiRender,
    TransportMode mode = TransportMode::Radiance,
    BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All
  ) const;
  void Regularize();
  const BSDFFrame& Frame() const;
  normal3f GeometricNormal() const;

private:
  BxDF bxdf_;
  BSDFFrame shadingFrame_;
  normal3f geometricNormal_ = normal3f(0, 0, 1);
};

} // namespace render
} // namespace rayrender

#endif
