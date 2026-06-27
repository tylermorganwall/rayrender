#include "src/materials/spectral_material.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <string>

using namespace rayrender::base;
using namespace rayrender::materials;
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
    std::cerr << "PR19 Material test failed: " << message << std::endl;
    std::exit(1);
  }
}

void CheckApprox(Float actual, Float expected, Float tolerance, const char* message) {
  if (!Approx(actual, expected, tolerance)) {
    std::cerr << "PR19 Material test failed: " << message
              << " expected " << expected << " got " << actual << std::endl;
    std::exit(1);
  }
}

void CheckSpectrumConstant(const SampledSpectrum& spectrum, Float expected, const char* message) {
  for (int i = 0; i < NSpectrumSamples; ++i) {
    if (!Approx(spectrum[i], expected, static_cast<Float>(1e-4))) {
      std::cerr << "PR19 Material test failed: " << message
                << " component " << i << " expected " << expected
                << " got " << spectrum[i] << std::endl;
      std::exit(1);
    }
  }
}

void CheckSpectrumFiniteNonNegative(const SampledSpectrum& spectrum, const char* message) {
  for (int i = 0; i < NSpectrumSamples; ++i) {
    if (!std::isfinite(spectrum[i]) || spectrum[i] < 0) {
      std::cerr << "PR19 Material test failed: " << message
                << " component " << i << " got " << spectrum[i] << std::endl;
      std::exit(1);
    }
  }
}

vec3f Unit(vec3f value) {
  return value / value.length();
}

SampledWavelengths TestWavelengths() {
  return SampledWavelengths::SampleUniform(static_cast<Float>(0.41));
}

MaterialEvalContext BaseContext() {
  MaterialEvalContext ctx;
  ctx.p = point3f(0, 0, 0);
  ctx.uv = point2f(static_cast<Float>(0.5), static_cast<Float>(0.5));
  ctx.dpdu = vec3f(1, 0, 0);
  ctx.dpdv = vec3f(0, 1, 0);
  ctx.dpdus = ctx.dpdu;
  ctx.dpdvs = ctx.dpdv;
  ctx.n = normal3f(0, 0, 1);
  ctx.ns = normal3f(0, 0, 1);
  ctx.dndus = normal3f(0, 0, 0);
  ctx.dndvs = normal3f(0, 0, 0);
  ctx.wo = vec3f(0, 0, 1);
  return ctx;
}

struct CountingTextureEvaluator {
  mutable int floatCount = 0;
  mutable int spectrumCount = 0;
  UniversalTextureEvaluator evaluator;

  Float operator()(const FloatTexture& texture, const TextureEvalContext& ctx) const {
    ++floatCount;
    return evaluator(texture, ctx);
  }

  SampledSpectrum operator()(
    const SpectrumTexture& texture,
    const TextureEvalContext& ctx,
    const SampledWavelengths& lambda
  ) const {
    ++spectrumCount;
    return evaluator(texture, ctx, lambda);
  }
};

void TestDiffuseTransmissionBxDF() {
  DiffuseTransmissionBxDF bxdf(
    SampledSpectrum(static_cast<Float>(0.25)),
    SampledSpectrum(static_cast<Float>(0.5))
  );
  vec3f wo = Unit(vec3f(static_cast<Float>(0.2), static_cast<Float>(0.1), static_cast<Float>(0.97)));
  vec3f wiReflect =
    Unit(vec3f(static_cast<Float>(-0.1), static_cast<Float>(0.3), static_cast<Float>(0.94)));
  vec3f wiTransmit =
    Unit(vec3f(static_cast<Float>(0.1), static_cast<Float>(0.2), static_cast<Float>(-0.97)));

  Check(base::HasFlag(bxdf.Flags(), BxDFFlags::Reflection), "diffuse transmission reflects");
  Check(base::HasFlag(bxdf.Flags(), BxDFFlags::Transmission), "diffuse transmission transmits");
  Check(base::HasFlag(bxdf.Flags(), BxDFFlags::Diffuse), "diffuse transmission is diffuse");
  CheckSpectrumConstant(
    bxdf.f(wo, wiReflect, TransportMode::Radiance),
    static_cast<Float>(0.25) * InvPi,
    "diffuse transmission same-side f"
  );
  CheckSpectrumConstant(
    bxdf.f(wo, wiTransmit, TransportMode::Radiance),
    static_cast<Float>(0.5) * InvPi,
    "diffuse transmission opposite-side f"
  );

  std::optional<BSDFSample> reflection = bxdf.Sample_f(
    wo,
    0,
    point2f(static_cast<Float>(0.3), static_cast<Float>(0.7)),
    TransportMode::Radiance,
    BxDFReflTransFlags::Reflection
  );
  Check(reflection.has_value(), "reflection-only diffuse transmission samples");
  Check(reflection->IsReflection(), "reflection-only sample is reflection");
  CheckApprox(
    reflection->pdf,
    bxdf.PDF(wo, reflection->wi, TransportMode::Radiance, BxDFReflTransFlags::Reflection),
    static_cast<Float>(1e-6),
    "reflection-only PDF matches"
  );

  std::optional<BSDFSample> transmission = bxdf.Sample_f(
    wo,
    static_cast<Float>(0.99),
    point2f(static_cast<Float>(0.2), static_cast<Float>(0.4)),
    TransportMode::Radiance,
    BxDFReflTransFlags::Transmission
  );
  Check(transmission.has_value(), "transmission-only diffuse transmission samples");
  Check(transmission->IsTransmission(), "transmission-only sample is transmission");
  CheckApprox(
    transmission->pdf,
    bxdf.PDF(wo, transmission->wi, TransportMode::Radiance, BxDFReflTransFlags::Transmission),
    static_cast<Float>(1e-6),
    "transmission-only PDF matches"
  );
}

void TestCoatedDiffuseBxDF() {
  CoatedDiffuseBxDF coated(
    DielectricBxDF(
      static_cast<Float>(1.5),
      TrowbridgeReitzDistribution(static_cast<Float>(0.35), static_cast<Float>(0.24))
    ),
    DiffuseBxDF(SampledSpectrum(static_cast<Float>(0.55))),
    static_cast<Float>(0.04),
    SampledSpectrum(0),
    0,
    8,
    4
  );
  vec3f wo =
    Unit(vec3f(static_cast<Float>(0.22), static_cast<Float>(-0.12), static_cast<Float>(0.97)));

  Check(base::HasFlag(coated.Flags(), BxDFFlags::Reflection), "coated diffuse reflects");
  Check(base::HasFlag(coated.Flags(), BxDFFlags::Diffuse), "coated diffuse has diffuse flag");

  std::optional<BSDFSample> sample = coated.Sample_f(
    wo,
    0,
    point2f(static_cast<Float>(0.43), static_cast<Float>(0.71)),
    TransportMode::Radiance,
    BxDFReflTransFlags::Reflection
  );
  Check(sample.has_value(), "coated diffuse samples reflected top layer");
  Check(sample->pdf > 0 && std::isfinite(sample->pdf), "coated diffuse sample PDF is finite");
  Check(sample->IsReflection(), "coated diffuse sample is reflection");
  Check(sample->pdfIsProportional, "coated diffuse sample reports proportional PDF");
  CheckSpectrumFiniteNonNegative(sample->f, "coated diffuse sample f is valid");
  Check(
    coated.PDF(wo, sample->wi, TransportMode::Radiance, BxDFReflTransFlags::Reflection) > 0,
    "coated diffuse PDF is positive for sampled direction"
  );
  CheckSpectrumFiniteNonNegative(
    coated.f(wo, sample->wi, TransportMode::Radiance),
    "coated diffuse f is valid"
  );

  coated.Regularize();
  std::optional<BSDFSample> regularized = coated.Sample_f(
    wo,
    0,
    point2f(static_cast<Float>(0.11), static_cast<Float>(0.83)),
    TransportMode::Radiance,
    BxDFReflTransFlags::Reflection
  );
  Check(regularized.has_value(), "regularized coated diffuse remains sampleable");
}

void TestCoatedConductorBxDF() {
  CoatedConductorBxDF coated(
    DielectricBxDF(
      static_cast<Float>(1.45),
      TrowbridgeReitzDistribution(static_cast<Float>(0.28), static_cast<Float>(0.31))
    ),
    ConductorBxDF(
      TrowbridgeReitzDistribution(static_cast<Float>(0.42), static_cast<Float>(0.37)),
      SampledSpectrum(static_cast<Float>(0.18)),
      SampledSpectrum(static_cast<Float>(3.1))
    ),
    static_cast<Float>(0.03),
    SampledSpectrum(0),
    0,
    8,
    4
  );
  vec3f wo =
    Unit(vec3f(static_cast<Float>(-0.18), static_cast<Float>(0.19), static_cast<Float>(0.97)));
  vec3f wi =
    Unit(vec3f(static_cast<Float>(0.2), static_cast<Float>(0.14), static_cast<Float>(0.96)));

  Check(base::HasFlag(coated.Flags(), BxDFFlags::Reflection), "coated conductor reflects");
  Check(base::HasFlag(coated.Flags(), BxDFFlags::Glossy), "coated conductor is glossy");
  CheckSpectrumFiniteNonNegative(
    coated.f(wo, wi, TransportMode::Radiance),
    "coated conductor f is valid"
  );

  std::optional<BSDFSample> sample = coated.Sample_f(
    wo,
    0,
    point2f(static_cast<Float>(0.37), static_cast<Float>(0.61)),
    TransportMode::Radiance,
    BxDFReflTransFlags::Reflection
  );
  Check(sample.has_value(), "coated conductor samples reflected top layer");
  Check(sample->IsReflection(), "coated conductor sample is reflection");
  Check(sample->pdf > 0 && std::isfinite(sample->pdf), "coated conductor sample PDF is finite");
  CheckSpectrumFiniteNonNegative(sample->f, "coated conductor sample f is valid");
  Check(
    coated.PDF(wo, sample->wi, TransportMode::Radiance, BxDFReflTransFlags::Reflection) > 0,
    "coated conductor PDF is positive for sampled direction"
  );
}

void TestDiffuseTransmissionMaterial() {
  Material material = Material::DiffuseTransmission(DiffuseTransmissionMaterial(
    SpectrumTexture::Constant(static_cast<Float>(0.8)),
    SpectrumTexture::Constant(static_cast<Float>(0.6)),
    static_cast<Float>(0.5),
    std::nullopt,
    std::nullopt
  ));
  Check(material.Type() == MaterialType::DiffuseTransmission, "diffuse transmission material type");
  Check(
    std::string(MaterialTypeName(material.Type())) == "DiffuseTransmission",
    "diffuse transmission material type name"
  );

  MaterialEvalContext ctx = BaseContext();
  SampledWavelengths lambda = TestWavelengths();
  ScratchBuffer scratch(4096);
  CountingTextureEvaluator evaluator;
  BSDF bsdf = material.GetBSDF(evaluator, ctx, lambda, scratch);
  Check(static_cast<bool>(bsdf), "diffuse transmission material creates BSDF");
  Check(base::HasFlag(bsdf.Flags(), BxDFFlags::Transmission), "material BSDF transmits");
  CheckSpectrumConstant(
    bsdf.f(vec3f(0, 0, 1), vec3f(0, 0, 1)),
    static_cast<Float>(0.4) * InvPi,
    "material reflectance scale"
  );
  CheckSpectrumConstant(
    bsdf.f(vec3f(0, 0, 1), vec3f(0, 0, -1)),
    static_cast<Float>(0.3) * InvPi,
    "material transmittance scale"
  );
  Check(evaluator.spectrumCount == 2, "diffuse transmission evaluates two spectrum textures");
  Check(evaluator.floatCount == 0, "diffuse transmission evaluates no float textures");
}

void TestCoatedDiffuseMaterial() {
  CoatedDiffuseMaterial coated(
    SpectrumTexture::Constant(static_cast<Float>(0.5)),
    FloatTexture::Constant(static_cast<Float>(0.16)),
    FloatTexture::Constant(static_cast<Float>(0.02)),
    SpectrumTexture::Constant(0),
    FloatTexture::Constant(0),
    ConstantEtaSpectrum(static_cast<Float>(1.45)),
    true,
    6,
    3
  );
  Material material = Material::CoatedDiffuse(coated);
  Check(material.Type() == MaterialType::CoatedDiffuse, "coated diffuse material type");

  MaterialEvalContext ctx = BaseContext();
  SampledWavelengths lambda = TestWavelengths();
  ScratchBuffer scratch(8192);
  CountingTextureEvaluator evaluator;
  BSDF bsdf = material.GetBSDF(evaluator, ctx, lambda, scratch);
  Check(static_cast<bool>(bsdf), "coated diffuse material creates BSDF");
  Check(!lambda.SecondaryTerminated(), "constant coated diffuse eta preserves wavelengths");
  Check(base::HasFlag(bsdf.Flags(), BxDFFlags::Reflection), "coated diffuse material reflects");
  Check(base::HasFlag(bsdf.Flags(), BxDFFlags::Diffuse), "coated diffuse material is diffuse");
  Check(evaluator.spectrumCount == 2, "coated diffuse evaluates reflectance and albedo");
  Check(evaluator.floatCount == 4, "coated diffuse evaluates roughness, thickness, and g");

  Material dispersive = Material::CoatedDiffuse(CoatedDiffuseMaterial(
    SpectrumTexture::Constant(static_cast<Float>(0.5)),
    FloatTexture::Constant(static_cast<Float>(0.1)),
    FloatTexture::Constant(static_cast<Float>(0.01)),
    SpectrumTexture::Constant(0),
    FloatTexture::Constant(0),
    CauchyEtaSpectrum(static_cast<Float>(1.5), static_cast<Float>(0.004)),
    true,
    4,
    2
  ));
  SampledWavelengths dispersiveLambda = TestWavelengths();
  ScratchBuffer dispersiveScratch(8192);
  CountingTextureEvaluator dispersiveEvaluator;
  BSDF dispersiveBsdf =
    dispersive.GetBSDF(dispersiveEvaluator, ctx, dispersiveLambda, dispersiveScratch);
  Check(static_cast<bool>(dispersiveBsdf), "dispersive coated diffuse creates BSDF");
  Check(
    dispersiveLambda.SecondaryTerminated(),
    "dispersive coated diffuse eta terminates secondary wavelengths"
  );
}

void TestCoatedConductorMaterial() {
  Material material = Material::CoatedConductor(CoatedConductorMaterial::FromReflectance(
    SpectrumTexture::Constant(static_cast<Float>(0.82)),
    FloatTexture::Constant(static_cast<Float>(0.12)),
    FloatTexture::Constant(static_cast<Float>(0.14)),
    FloatTexture::Constant(static_cast<Float>(0.19)),
    FloatTexture::Constant(static_cast<Float>(0.21)),
    FloatTexture::Constant(static_cast<Float>(0.02)),
    SpectrumTexture::Constant(0),
    FloatTexture::Constant(0),
    ConstantEtaSpectrum(static_cast<Float>(1.5)),
    true,
    6,
    3
  ));
  Check(material.Type() == MaterialType::CoatedConductor, "coated conductor material type");

  MaterialEvalContext ctx = BaseContext();
  SampledWavelengths lambda = TestWavelengths();
  ScratchBuffer scratch(8192);
  CountingTextureEvaluator evaluator;
  BSDF bsdf = material.GetBSDF(evaluator, ctx, lambda, scratch);
  Check(static_cast<bool>(bsdf), "coated conductor material creates BSDF");
  Check(!lambda.SecondaryTerminated(), "constant coated conductor eta preserves wavelengths");
  Check(base::HasFlag(bsdf.Flags(), BxDFFlags::Reflection), "coated conductor material reflects");
  Check(base::HasFlag(bsdf.Flags(), BxDFFlags::Glossy), "coated conductor material is glossy");
  Check(evaluator.spectrumCount == 2, "coated conductor evaluates reflectance and albedo");
  Check(evaluator.floatCount == 6, "coated conductor evaluates roughnesses, thickness, and g");

  Material dispersive = Material::CoatedConductor(CoatedConductorMaterial::FromEtaK(
    SpectrumTexture::Constant(static_cast<Float>(0.25)),
    SpectrumTexture::Constant(static_cast<Float>(3)),
    FloatTexture::Constant(static_cast<Float>(0.1)),
    FloatTexture::Constant(static_cast<Float>(0.1)),
    FloatTexture::Constant(static_cast<Float>(0.2)),
    FloatTexture::Constant(static_cast<Float>(0.2)),
    FloatTexture::Constant(static_cast<Float>(0.01)),
    SpectrumTexture::Constant(0),
    FloatTexture::Constant(0),
    CauchyEtaSpectrum(static_cast<Float>(1.5), static_cast<Float>(0.004)),
    true,
    4,
    2
  ));
  SampledWavelengths dispersiveLambda = TestWavelengths();
  ScratchBuffer dispersiveScratch(8192);
  CountingTextureEvaluator dispersiveEvaluator;
  BSDF dispersiveBsdf =
    dispersive.GetBSDF(dispersiveEvaluator, ctx, dispersiveLambda, dispersiveScratch);
  Check(static_cast<bool>(dispersiveBsdf), "dispersive coated conductor creates BSDF");
  Check(
    dispersiveLambda.SecondaryTerminated(),
    "dispersive coated conductor interface eta terminates secondary wavelengths"
  );
}

} // namespace

int main() {
  TestDiffuseTransmissionBxDF();
  TestCoatedDiffuseBxDF();
  TestCoatedConductorBxDF();
  TestDiffuseTransmissionMaterial();
  TestCoatedDiffuseMaterial();
  TestCoatedConductorMaterial();
  std::cout << "PR19 Material tests passed" << std::endl;
  return 0;
}
