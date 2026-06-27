#include "src/render/spectral_bsdf.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

using namespace rayrender::base;
using namespace rayrender::render;
namespace base = rayrender::base;

namespace {

constexpr Float Pi = static_cast<Float>(3.14159265358979323846264338327950288);
constexpr Float InvPi = static_cast<Float>(1) / Pi;

bool Approx(Float lhs, Float rhs, Float tolerance = static_cast<Float>(1e-5)) {
  return std::fabs(lhs - rhs) <= tolerance;
}

void Check(bool condition, const char* message) {
  if (!condition) {
    std::cerr << "PR12 BSDF test failed: " << message << std::endl;
    std::exit(1);
  }
}

void CheckApprox(Float actual, Float expected, Float tolerance, const char* message) {
  if (!Approx(actual, expected, tolerance)) {
    std::cerr << "PR12 BSDF test failed: " << message
              << " expected " << expected << " got " << actual << std::endl;
    std::exit(1);
  }
}

void CheckSpectrumConstant(const SampledSpectrum& spectrum, Float expected, const char* message) {
  for (int i = 0; i < NSpectrumSamples; ++i) {
    if (!Approx(spectrum[i], expected, static_cast<Float>(1e-4))) {
      std::cerr << "PR12 BSDF test failed: " << message
                << " component " << i << " expected " << expected
                << " got " << spectrum[i] << std::endl;
      std::exit(1);
    }
  }
}

vec3f Unit(vec3f v) {
  return v / v.length();
}

Float ReferenceFrDielectric(Float cosThetaIValue, Float etaValue) {
  long double cosThetaI = std::min(
    std::max(static_cast<long double>(cosThetaIValue), static_cast<long double>(-1)),
    static_cast<long double>(1)
  );
  long double eta = etaValue;
  if (cosThetaI < 0) {
    eta = static_cast<long double>(1) / eta;
    cosThetaI = -cosThetaI;
  }

  long double sin2ThetaI = static_cast<long double>(1) - cosThetaI * cosThetaI;
  long double sin2ThetaT = sin2ThetaI / (eta * eta);
  if (sin2ThetaT >= 1) {
    return 1;
  }
  long double cosThetaT = std::sqrt(static_cast<double>(1 - sin2ThetaT));
  long double rParallel = (eta * cosThetaI - cosThetaT) / (eta * cosThetaI + cosThetaT);
  long double rPerpendicular = (cosThetaI - eta * cosThetaT) / (cosThetaI + eta * cosThetaT);
  return static_cast<Float>((rParallel * rParallel + rPerpendicular * rPerpendicular) / 2);
}

Float ReferenceFrComplexNormal(Float eta, Float k) {
  Float numerator = (eta - 1) * (eta - 1) + k * k;
  Float denominator = (eta + 1) * (eta + 1) + k * k;
  return numerator / denominator;
}

void TestFlagsAndTransportMode() {
  Check(base::HasFlag(BxDFFlags::DiffuseReflection, BxDFFlags::Reflection), "diffuse reflection flag reflects");
  Check(base::HasFlag(BxDFFlags::DiffuseReflection, BxDFFlags::Diffuse), "diffuse reflection flag is diffuse");
  Check(base::HasAny(BxDFFlags::DiffuseReflection, BxDFReflTransFlags::Reflection), "BxDF/sample flag overlap works");
  Check(!base::HasAny(BxDFFlags::DiffuseReflection, BxDFReflTransFlags::Transmission), "reflection does not match transmission sample flag");
  Check(base::IsNonSpecular(BxDFFlags::GlossyReflection), "glossy flag is non-specular");
  Check(!base::IsNonSpecular(BxDFFlags::SpecularReflection), "specular flag is not non-specular");
  Check(!TransportMode::Radiance == TransportMode::Importance, "transport mode negation flips radiance");
  Check(!TransportMode::Importance == TransportMode::Radiance, "transport mode negation flips importance");
}

void TestFresnel() {
  CheckApprox(
    FrDielectric(1, static_cast<Float>(1.5)),
    static_cast<Float>(0.04),
    static_cast<Float>(1e-6),
    "dielectric Fresnel normal incidence"
  );
  CheckApprox(
    FrDielectric(static_cast<Float>(0.5), static_cast<Float>(1.5)),
    ReferenceFrDielectric(static_cast<Float>(0.5), static_cast<Float>(1.5)),
    static_cast<Float>(1e-6),
    "dielectric Fresnel oblique incidence"
  );
  CheckApprox(
    FrDielectric(static_cast<Float>(-0.35), static_cast<Float>(1.5)),
    ReferenceFrDielectric(static_cast<Float>(-0.35), static_cast<Float>(1.5)),
    static_cast<Float>(1e-6),
    "dielectric Fresnel flips entering side"
  );
  CheckApprox(
    FrDielectric(static_cast<Float>(0.5), static_cast<Float>(0.5)),
    1,
    static_cast<Float>(1e-6),
    "dielectric Fresnel detects total internal reflection"
  );

  Float conductorReference = ReferenceFrComplexNormal(static_cast<Float>(2), static_cast<Float>(3));
  CheckApprox(
    FrComplex(1, static_cast<Float>(2), static_cast<Float>(3)),
    conductorReference,
    static_cast<Float>(1e-6),
    "conductor Fresnel normal incidence"
  );
  SampledSpectrum conductor = FrComplex(
    1,
    SampledSpectrum(2),
    SampledSpectrum(3)
  );
  CheckSpectrumConstant(conductor, conductorReference, "sampled conductor Fresnel is componentwise");
}

void TestDiffuseBxDF() {
  DiffuseBxDF diffuse(SampledSpectrum(static_cast<Float>(0.6)));
  vec3f wo(0, 0, 1);
  vec3f wi(Unit(vec3f(static_cast<Float>(0.5), 0, static_cast<Float>(0.5))));

  CheckSpectrumConstant(diffuse.f(wo, wi, TransportMode::Radiance), static_cast<Float>(0.6) * InvPi, "diffuse f excludes cosine");
  CheckApprox(
    diffuse.PDF(wo, wi, TransportMode::Radiance),
    AbsCosTheta(wi) * InvPi,
    static_cast<Float>(1e-6),
    "diffuse PDF is cosine-weighted"
  );
  CheckSpectrumConstant(diffuse.rho(), static_cast<Float>(0.6), "diffuse rho is reflectance");

  Float averageReflectance = 0;
  int samples = 0;
  for (int y = 0; y < 16; ++y) {
    for (int x = 0; x < 16; ++x) {
      point2f u((x + static_cast<Float>(0.5)) / 16, (y + static_cast<Float>(0.5)) / 16);
      std::optional<BSDFSample> sample = diffuse.Sample_f(wo, 0, u, TransportMode::Radiance);
      Check(sample.has_value(), "diffuse Sample_f returns reflection sample");
      Check(sample->IsReflection(), "diffuse sample reports reflection");
      Check(sample->IsDiffuse(), "diffuse sample reports diffuse");
      Check(sample->wi.xyz.z > 0, "diffuse sample follows positive hemisphere");
      Float pdf = diffuse.PDF(wo, sample->wi, TransportMode::Radiance);
      CheckApprox(sample->pdf, pdf, static_cast<Float>(1e-6), "diffuse sampled PDF agrees with PDF()");
      averageReflectance += sample->f[0] * AbsCosTheta(sample->wi) / sample->pdf;
      ++samples;
    }
  }
  averageReflectance /= static_cast<Float>(samples);
  CheckApprox(averageReflectance, static_cast<Float>(0.6), static_cast<Float>(1e-5), "diffuse cosine integral returns reflectance");

  std::optional<BSDFSample> negativeWo = diffuse.Sample_f(
    vec3f(0, 0, -1),
    0,
    point2f(static_cast<Float>(0.25), static_cast<Float>(0.5)),
    TransportMode::Radiance
  );
  Check(negativeWo.has_value(), "diffuse samples below the surface");
  Check(negativeWo->wi.xyz.z < 0, "diffuse sampling flips to outgoing hemisphere");

  std::optional<BSDFSample> noTransmission = diffuse.Sample_f(
    wo,
    0,
    point2f(static_cast<Float>(0.25), static_cast<Float>(0.25)),
    TransportMode::Radiance,
    BxDFReflTransFlags::Transmission
  );
  Check(!noTransmission.has_value(), "diffuse reflection respects transmission-only sample flags");

  vec3f grazing(Unit(vec3f(1, 0, static_cast<Float>(0.01))));
  CheckSpectrumConstant(
    diffuse.f(wo, grazing, TransportMode::Radiance),
    static_cast<Float>(0.6) * InvPi,
    "diffuse f has no hidden cosine at grazing angles"
  );
}

void TestTrowbridgeReitzDistribution() {
  CheckApprox(
    TrowbridgeReitzDistribution::RoughnessToAlpha(static_cast<Float>(0.25)),
    static_cast<Float>(0.5),
    static_cast<Float>(1e-6),
    "roughness remaps to pbrt sqrt alpha"
  );

  TrowbridgeReitzDistribution distribution(static_cast<Float>(0.5), static_cast<Float>(0.8));
  vec3f w = Unit(vec3f(static_cast<Float>(0.3), static_cast<Float>(0.4), static_cast<Float>(0.8660254)));
  vec3f wm = Unit(vec3f(static_cast<Float>(0.2), static_cast<Float>(0.1), static_cast<Float>(1)));
  Check(distribution.D(wm) > 0, "microfacet D is positive for upper-hemisphere normal");
  Check(distribution.G1(w) > 0 && distribution.G1(w) <= 1, "microfacet G1 is bounded");
  CheckApprox(distribution.PDF(w, wm), distribution.D(w, wm), static_cast<Float>(1e-7), "microfacet PDF delegates visible-normal density");

  Float integral = 0;
  int sampleCount = 0;
  constexpr int zSamples = 96;
  constexpr int phiSamples = 192;
  Float uniformHemispherePdf = static_cast<Float>(1) / (static_cast<Float>(2) * Pi);
  for (int zi = 0; zi < zSamples; ++zi) {
    Float z = (zi + static_cast<Float>(0.5)) / zSamples;
    Float r = std::sqrt(std::max(static_cast<Float>(0), static_cast<Float>(1) - z * z));
    for (int pi = 0; pi < phiSamples; ++pi) {
      Float phi = static_cast<Float>(2) * Pi * (pi + static_cast<Float>(0.5)) / phiSamples;
      vec3f sampleWm(r * std::cos(phi), r * std::sin(phi), z);
      if (dot(w, sampleWm) > 0) {
        integral += distribution.PDF(w, sampleWm) / uniformHemispherePdf;
      }
      ++sampleCount;
    }
  }
  integral /= static_cast<Float>(sampleCount);
  CheckApprox(integral, 1, static_cast<Float>(0.015), "microfacet visible-normal PDF integrates to one");

  for (int i = 0; i < 256; ++i) {
    Float u0 = (i + static_cast<Float>(0.5)) / 256;
    Float u1 = std::fmod(static_cast<Float>(0.37) * i + static_cast<Float>(0.19), static_cast<Float>(1));
    vec3f sampled = distribution.Sample_wm(w, point2f(u0, u1));
    Check(sampled.xyz.z > 0, "microfacet sampling returns upper-hemisphere normal");
    CheckApprox(sampled.length(), 1, static_cast<Float>(2e-5), "microfacet sampled normal is normalized");
    Check(distribution.PDF(w, sampled) > 0, "microfacet sampled normal has positive PDF");
    Check(std::isfinite(distribution.PDF(w, sampled)), "microfacet sampled PDF is finite");
  }

  TrowbridgeReitzDistribution smooth(0, 0);
  Check(smooth.EffectivelySmooth(), "zero roughness is effectively smooth");
  CheckApprox(
    smooth.Sample_wm(vec3f(0, 0, -1), point2f(static_cast<Float>(0.2), static_cast<Float>(0.8))).xyz.z,
    1,
    static_cast<Float>(1e-6),
    "smooth microfacet sampling returns pbrt canonical normal"
  );

  TrowbridgeReitzDistribution regularized(static_cast<Float>(0.05), static_cast<Float>(0.2));
  regularized.Regularize();
  CheckApprox(regularized.AlphaX(), static_cast<Float>(0.1), static_cast<Float>(1e-6), "regularize clamps small alpha x");
  CheckApprox(regularized.AlphaY(), static_cast<Float>(0.3), static_cast<Float>(1e-6), "regularize clamps alpha y");
}

void TestBSDFFrameWrapper() {
  DiffuseBxDF diffuse(SampledSpectrum(static_cast<Float>(0.7)));
  BSDF bsdf(
    normal3f(0, 1, 0),
    normal3f(0, 1, 0),
    vec3f(1, 0, 0),
    BxDF(&diffuse)
  );

  Check(static_cast<bool>(bsdf), "BSDF wrapper holds BxDF");
  CheckApprox(bsdf.RenderToLocal(vec3f(0, 1, 0)).xyz.z, 1, static_cast<Float>(1e-6), "shading normal maps to local z");
  CheckApprox(bsdf.LocalToRender(vec3f(0, 0, 1)).xyz.y, 1, static_cast<Float>(1e-6), "local z maps to shading normal");
  CheckApprox(dot(convert_to_vec3(bsdf.GeometricNormal()), vec3f(0, 1, 0)), 1, static_cast<Float>(1e-6), "geometric normal is stored normalized");

  vec3f woRender(0, 1, 0);
  vec3f wiRender = bsdf.LocalToRender(Unit(vec3f(static_cast<Float>(0.3), 0, static_cast<Float>(0.7))));
  CheckSpectrumConstant(bsdf.f(woRender, wiRender), static_cast<Float>(0.7) * InvPi, "BSDF delegates f in local frame");
  CheckApprox(
    bsdf.PDF(woRender, wiRender),
    AbsCosTheta(bsdf.RenderToLocal(wiRender)) * InvPi,
    static_cast<Float>(1e-6),
    "BSDF PDF delegates in local frame"
  );
  CheckSpectrumConstant(bsdf.f(wiRender, woRender), static_cast<Float>(0.7) * InvPi, "diffuse BSDF is reciprocal");

  std::optional<BSDFSample> sample = bsdf.Sample_f(
    woRender,
    0,
    point2f(static_cast<Float>(0.125), static_cast<Float>(0.625))
  );
  Check(sample.has_value(), "BSDF wrapper samples child BxDF");
  Check(dot(sample->wi, vec3f(0, 1, 0)) > 0, "BSDF sample is returned in render frame");
  CheckApprox(
    sample->pdf,
    bsdf.PDF(woRender, sample->wi),
    static_cast<Float>(1e-6),
    "BSDF sample PDF agrees after world/local transform"
  );

  NullBxDF nullBxDF;
  BSDF nullBsdf(
    normal3f(0, 0, 1),
    normal3f(0, 0, 1),
    vec3f(1, 0, 0),
    BxDF(&nullBxDF)
  );
  Check(nullBsdf.Flags() == BxDFFlags::Unset, "null BxDF has no flags");
  Check(!nullBsdf.Sample_f(vec3f(0, 0, 1), 0, point2f(static_cast<Float>(0.5), static_cast<Float>(0.5))).has_value(), "null BxDF produces no sample");

  BxDF empty(static_cast<DiffuseBxDF*>(nullptr));
  Check(!static_cast<bool>(empty), "null BxDF pointer constructs empty dispatch handle");
  Check(empty.Flags() == BxDFFlags::Unset, "empty dispatch handle has no flags");
}

} // namespace

int main() {
  TestFlagsAndTransportMode();
  TestFresnel();
  TestDiffuseBxDF();
  TestTrowbridgeReitzDistribution();
  TestBSDFFrameWrapper();
  std::cout << "PR12 BSDF tests passed" << std::endl;
  return 0;
}
