#include "src/materials/spectral_material.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <new>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

using namespace rayrender::base;
using namespace rayrender::materials;
using namespace rayrender::render;

namespace {

constexpr Float Pi = static_cast<Float>(3.14159265358979323846264338327950288);
constexpr Float InvPi = static_cast<Float>(1) / Pi;

bool g_countAllocations = false;
std::size_t g_allocationCount = 0;

bool Approx(Float lhs, Float rhs, Float tolerance = static_cast<Float>(1e-5)) {
  return std::fabs(lhs - rhs) <= tolerance;
}

void Check(bool condition, const char* message) {
  if (!condition) {
    std::cerr << "PR13 Material test failed: " << message << std::endl;
    std::exit(1);
  }
}

void CheckApprox(Float actual, Float expected, Float tolerance, const char* message) {
  if (!Approx(actual, expected, tolerance)) {
    std::cerr << "PR13 Material test failed: " << message
              << " expected " << expected << " got " << actual << std::endl;
    std::exit(1);
  }
}

void CheckSpectrumConstant(const SampledSpectrum& spectrum, Float expected, const char* message) {
  for (int i = 0; i < NSpectrumSamples; ++i) {
    if (!Approx(spectrum[i], expected, static_cast<Float>(1e-4))) {
      std::cerr << "PR13 Material test failed: " << message
                << " component " << i << " expected " << expected
                << " got " << spectrum[i] << std::endl;
      std::exit(1);
    }
  }
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
  ctx.dudx = static_cast<Float>(0.01);
  ctx.dvdy = static_cast<Float>(0.01);
  return ctx;
}

SampledWavelengths TestWavelengths() {
  return SampledWavelengths::SampleUniform(static_cast<Float>(0.35));
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

void TestDiffuseMaterialMatchesDirectBxDF() {
  Material material = Material::Diffuse(SpectrumTexture::Constant(static_cast<Float>(0.6)));
  MaterialEvalContext ctx = BaseContext();
  SampledWavelengths lambda = TestWavelengths();
  ScratchBuffer scratch(1024);
  CountingTextureEvaluator evaluator;

  BSDF materialBsdf = material.GetBSDF(evaluator, ctx, lambda, scratch);
  DiffuseBxDF directBxDF(SampledSpectrum(static_cast<Float>(0.6)));
  BSDF directBsdf(ctx.n, ctx.ns, ctx.dpdus, BxDF(&directBxDF));

  vec3f wo(0, 0, 1);
  vec3f wi(vec3f(static_cast<Float>(0.3), static_cast<Float>(0.4), static_cast<Float>(0.8660254)));
  wi = wi / wi.length();
  CheckSpectrumConstant(materialBsdf.f(wo, wi), static_cast<Float>(0.6) * InvPi, "material diffuse f");
  CheckSpectrumConstant(directBsdf.f(wo, wi), static_cast<Float>(0.6) * InvPi, "direct diffuse f");
  CheckApprox(materialBsdf.PDF(wo, wi), directBsdf.PDF(wo, wi), static_cast<Float>(1e-6), "material diffuse PDF equals direct PDF");
  Check(evaluator.spectrumCount == 1, "diffuse reflectance texture is evaluated once");
  Check(evaluator.floatCount == 0, "diffuse material without alpha or bump evaluates no float textures");
  Check(scratch.Used() >= sizeof(DiffuseBxDF), "diffuse BxDF is scratch allocated");
}

void TestReflectanceClamping() {
  Material high = Material::Diffuse(SpectrumTexture::Constant(static_cast<Float>(2)));
  Material low = Material::Diffuse(SpectrumTexture::Constant(static_cast<Float>(-1)));
  MaterialEvalContext ctx = BaseContext();
  SampledWavelengths lambda = TestWavelengths();
  UniversalTextureEvaluator evaluator;

  ScratchBuffer highScratch(1024);
  BSDF highBsdf = high.GetBSDF(evaluator, ctx, lambda, highScratch);
  CheckSpectrumConstant(highBsdf.f(vec3f(0, 0, 1), vec3f(0, 0, 1)), InvPi, "diffuse reflectance clamps high values");

  ScratchBuffer lowScratch(1024);
  BSDF lowBsdf = low.GetBSDF(evaluator, ctx, lambda, lowScratch);
  CheckSpectrumConstant(lowBsdf.f(vec3f(0, 0, 1), vec3f(0, 0, 1)), 0, "diffuse reflectance clamps negative values");
  Check(lowBsdf.Flags() == BxDFFlags::Unset, "zero diffuse reflectance produces unset flags");
}

void TestNoHeapAllocationPerDiffuseHit() {
  Material material = Material::Diffuse(SpectrumTexture::Constant(static_cast<Float>(0.5)));
  MaterialEvalContext ctx = BaseContext();
  SampledWavelengths lambda = TestWavelengths();
  ScratchBuffer scratch(1024);
  CountingTextureEvaluator evaluator;

  g_allocationCount = 0;
  g_countAllocations = true;
  BSDF bsdf = material.GetBSDF(evaluator, ctx, lambda, scratch);
  g_countAllocations = false;

  Check(static_cast<bool>(bsdf), "diffuse material produces a BSDF");
  Check(g_allocationCount == 0, "diffuse GetBSDF performs no heap allocation per hit");
  Check(evaluator.spectrumCount == 1, "no-heap diffuse path still evaluates reflectance once");
}

void TestAlphaAndBumpInvariants() {
  MaterialEvalContext ctx = BaseContext();
  SampledWavelengths lambda = TestWavelengths();
  CountingTextureEvaluator alphaEvaluator;
  ScratchBuffer scratch(1024);

  Material transparent = Material::Diffuse(DiffuseMaterial(
    SpectrumTexture::Constant(static_cast<Float>(0.8)),
    FloatTexture::Constant(0),
    std::nullopt
  ));
  MaterialAlphaResult transparentAlpha = transparent.EvaluateAlpha(alphaEvaluator, ctx, static_cast<Float>(0.5));
  Check(!transparentAlpha.accepted, "zero alpha rejects hit before closure construction");
  CheckApprox(transparentAlpha.alpha, 0, static_cast<Float>(1e-6), "zero alpha value");
  Check(alphaEvaluator.floatCount == 1, "alpha texture evaluated once");
  Check(scratch.Used() == 0, "alpha rejection does not allocate a closure");
  CheckApprox(dot(convert_to_vec3(ctx.n), vec3f(0, 0, 1)), 1, static_cast<Float>(1e-6), "alpha test preserves geometric normal");

  Material partiallyOpaque = Material::Diffuse(DiffuseMaterial(
    SpectrumTexture::Constant(static_cast<Float>(0.8)),
    FloatTexture::Constant(static_cast<Float>(0.75)),
    std::nullopt
  ));
  CountingTextureEvaluator opaqueEvaluator;
  MaterialAlphaResult opaqueAlpha = partiallyOpaque.EvaluateAlpha(
    opaqueEvaluator,
    ctx,
    static_cast<Float>(0.5)
  );
  Check(opaqueAlpha.accepted, "alpha sample below alpha accepts hit");
  CheckApprox(opaqueAlpha.alpha, static_cast<Float>(0.75), static_cast<Float>(1e-6), "alpha value is clamped and reported");

  ctx.uv = point2f(static_cast<Float>(0.499), static_cast<Float>(0.5));
  FloatTexture bumpTexture = FloatTexture::Checker(
    FloatTexture::Constant(0),
    FloatTexture::Constant(10),
    2,
    static_cast<Float>(1000),
    1
  );
  Material bumped = Material::Diffuse(DiffuseMaterial(
    SpectrumTexture::Constant(static_cast<Float>(0.4)),
    std::nullopt,
    std::move(bumpTexture)
  ));
  CountingTextureEvaluator bumpEvaluator;
  BumpMapResult bump = bumped.EvaluateBump(bumpEvaluator, ctx);
  Check(bump.hasBump, "bump material reports bump result");
  Check(bumpEvaluator.floatCount == 3, "bump texture is evaluated at center, du, and dv");
  CheckApprox(dot(convert_to_vec3(bump.geometricNormal), vec3f(0, 0, 1)), 1, static_cast<Float>(1e-6), "bump preserves geometric normal");
  CheckApprox(bump.shadingNormal.length(), 1, static_cast<Float>(1e-5), "bump shading normal is normalized");
  Check(dot(convert_to_vec3(bump.shadingNormal), vec3f(0, 0, 1)) < static_cast<Float>(0.9999), "bump changes shading normal only");

  CountingTextureEvaluator bumpedBsdfEvaluator;
  BSDF bumpedBsdf = bumped.GetBSDF(bumpedBsdfEvaluator, ctx, lambda, scratch);
  CheckApprox(dot(convert_to_vec3(bumpedBsdf.GeometricNormal()), vec3f(0, 0, 1)), 1, static_cast<Float>(1e-6), "bumped BSDF stores geometric normal");
  Check(dot(bumpedBsdf.Frame().Z(), vec3f(0, 0, 1)) < static_cast<Float>(0.9999), "bumped BSDF frame uses bumped shading normal");
  Check(bumpedBsdfEvaluator.spectrumCount == 1, "bumped material evaluates reflectance once");
  Check(bumpedBsdfEvaluator.floatCount == 3, "bumped material evaluates bump texture three times");
}

void TestInterfaceMaterialAndTable() {
  Material material = Material::Interface();
  MaterialEvalContext ctx = BaseContext();
  SampledWavelengths lambda = TestWavelengths();
  ScratchBuffer scratch(1024);
  CountingTextureEvaluator evaluator;

  BSDF bsdf = material.GetBSDF(evaluator, ctx, lambda, scratch);
  Check(bsdf.Flags() == BxDFFlags::Unset, "interface material has no scattering flags");
  Check(!bsdf.Sample_f(vec3f(0, 0, 1), 0, point2f(static_cast<Float>(0.5), static_cast<Float>(0.5))).has_value(), "interface material samples no scattering");
  Check(evaluator.floatCount == 0 && evaluator.spectrumCount == 0, "interface material evaluates no textures");

  SpectralMaterialTable table;
  Material diffuse = Material::Diffuse(SpectrumTexture::Constant(static_cast<Float>(0.25)));
  MaterialHandle handle = table.Add(diffuse);
  Check(table.Size() == 1, "material table stores one material");
  Check(table.Get(handle).Type() == MaterialType::Diffuse, "material table returns stored diffuse material");
  MaterialHandle stale = MaterialHandle::FromIndex(handle.Index(), 2);
  bool threw = false;
  try {
    (void)table.Get(stale);
  } catch (const std::out_of_range&) {
    threw = true;
  }
  Check(threw, "material table rejects stale handles");

  MaterialTextureRequirements requirements = diffuse.TextureRequirements();
  Check(requirements.spectrumTextures, "diffuse material reports spectrum texture requirement");
  Check(!requirements.alphaTexture && !requirements.bumpTexture, "plain diffuse material reports no alpha or bump requirements");
  Check(material.CanEvaluateTextures(UniversalTextureEvaluator()), "interface material can evaluate with universal evaluator");
}

void TestSurfaceInteractionContextAndScratchReset() {
  SurfaceInteraction interaction;
  interaction.p = point3f(1, 2, 3);
  interaction.uv = point2f(static_cast<Float>(0.25), static_cast<Float>(0.75));
  interaction.n = normal3f(0, 0, 1);
  interaction.shadingNormal = normal3f(0, 1, 0);
  interaction.dpdu = vec3f(1, 0, 0);
  interaction.dpdv = vec3f(0, 1, 0);
  interaction.shadingDpdu = vec3f(1, 0, 0);
  interaction.shadingDpdv = vec3f(0, 0, -1);
  interaction.wo = vec3f(0, 0, 1);
  interaction.faceIndex = 9;

  MaterialEvalContext ctx(interaction);
  CheckApprox(ctx.p.xyz.x, 1, static_cast<Float>(1e-6), "context copies position");
  CheckApprox(ctx.uv[0], static_cast<Float>(0.25), static_cast<Float>(1e-6), "context copies uv");
  CheckApprox(dot(convert_to_vec3(ctx.n), vec3f(0, 0, 1)), 1, static_cast<Float>(1e-6), "context keeps geometric normal");
  CheckApprox(dot(convert_to_vec3(ctx.ns), vec3f(0, 1, 0)), 1, static_cast<Float>(1e-6), "context uses shading normal");
  Check(ctx.faceIndex == 9, "context copies face index");

  Material material = Material::Diffuse(SpectrumTexture::Constant(static_cast<Float>(0.5)));
  SampledWavelengths lambda = TestWavelengths();
  UniversalTextureEvaluator evaluator;
  ScratchBuffer scratch(1024);
  BSDF bsdf = material.GetBSDF(evaluator, BaseContext(), lambda, scratch);
  std::optional<BSDFSample> sample = bsdf.Sample_f(
    vec3f(0, 0, 1),
    0,
    point2f(static_cast<Float>(0.25), static_cast<Float>(0.75))
  );
  Check(sample.has_value(), "scratch closure works before reset");
  Check(scratch.Used() > 0, "scratch buffer records closure allocation");
  scratch.Reset();
  Check(scratch.Used() == 0, "scratch reset releases closure storage after sample completion");

  BSDF second = material.GetBSDF(evaluator, BaseContext(), lambda, scratch);
  Check(second.Sample_f(
    vec3f(0, 0, 1),
    0,
    point2f(static_cast<Float>(0.5), static_cast<Float>(0.25))
  ).has_value(), "scratch buffer can allocate a fresh closure after reset");
}

} // namespace

void* operator new(std::size_t size) {
  if (g_countAllocations) {
    ++g_allocationCount;
  }
  if (void* pointer = std::malloc(size)) {
    return pointer;
  }
  throw std::bad_alloc();
}

void* operator new[](std::size_t size) {
  if (g_countAllocations) {
    ++g_allocationCount;
  }
  if (void* pointer = std::malloc(size)) {
    return pointer;
  }
  throw std::bad_alloc();
}

void operator delete(void* pointer) noexcept {
  std::free(pointer);
}

void operator delete[](void* pointer) noexcept {
  std::free(pointer);
}

void operator delete(void* pointer, std::size_t) noexcept {
  std::free(pointer);
}

void operator delete[](void* pointer, std::size_t) noexcept {
  std::free(pointer);
}

int main() {
  TestDiffuseMaterialMatchesDirectBxDF();
  TestReflectanceClamping();
  TestNoHeapAllocationPerDiffuseHit();
  TestAlphaAndBumpInvariants();
  TestInterfaceMaterialAndTable();
  TestSurfaceInteractionContextAndScratchReset();
  std::cout << "PR13 Material tests passed" << std::endl;
  return 0;
}
