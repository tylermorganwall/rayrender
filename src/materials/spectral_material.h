#ifndef RAYRENDER_MATERIALS_SPECTRAL_MATERIAL_H
#define RAYRENDER_MATERIALS_SPECTRAL_MATERIAL_H

#include "spectral_texture.h"

#include "../base/base.h"
#include "../render/spectral_bsdf.h"
#include "../render/spectral_scene.h"

#include <cstddef>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace rayrender {
namespace materials {

enum class MaterialType {
  Diffuse,
  Conductor,
  Dielectric,
  Interface
};

const char* MaterialTypeName(MaterialType type);

struct MaterialEvalContext : public TextureEvalContext {
  vec3f wo;
  normal3f n;
  normal3f ns;
  vec3f dpdus;
  vec3f dpdvs;
  normal3f dndus;
  normal3f dndvs;
  const render::ResolvedDielectricInterface* dielectric = nullptr;

  MaterialEvalContext() = default;
  explicit MaterialEvalContext(const render::SurfaceInteraction& interaction);
};

MaterialEvalContext MaterialEvalContextFromSurfaceInteraction(
  const render::SurfaceInteraction& interaction
);

struct MaterialTextureRequirements {
  bool spectrumTextures = false;
  bool floatTextures = false;
  bool alphaTexture = false;
  bool bumpTexture = false;

  bool Any() const;
};

struct MaterialAlphaResult {
  Float alpha = 1;
  bool accepted = true;
};

struct BumpMapResult {
  normal3f geometricNormal;
  normal3f shadingNormal;
  vec3f dpdus;
  vec3f dpdvs;
  bool hasBump = false;
};

BumpMapResult DefaultBumpMapResult(const MaterialEvalContext& ctx);
normal3f NormalFromShadingTangents(
  const vec3f& dpdus,
  const vec3f& dpdvs,
  normal3f referenceNormal
);

template <typename TextureEvaluator>
BumpMapResult EvaluateBumpMap(
  const TextureEvaluator& texEval,
  const FloatTexture& displacement,
  const MaterialEvalContext& ctx
) {
  TextureEvalContext shiftedCtx = ctx;
  Float du = static_cast<Float>(0.5) * (std::abs(ctx.dudx) + std::abs(ctx.dudy));
  if (du == 0) {
    du = static_cast<Float>(0.0005);
  }
  shiftedCtx.p = ctx.p + du * ctx.dpdus;
  shiftedCtx.uv = point2f(ctx.uv[0] + du, ctx.uv[1]);
  Float uDisplace = texEval(displacement, shiftedCtx);

  Float dv = static_cast<Float>(0.5) * (std::abs(ctx.dvdx) + std::abs(ctx.dvdy));
  if (dv == 0) {
    dv = static_cast<Float>(0.0005);
  }
  shiftedCtx = ctx;
  shiftedCtx.p = ctx.p + dv * ctx.dpdvs;
  shiftedCtx.uv = point2f(ctx.uv[0], ctx.uv[1] + dv);
  Float vDisplace = texEval(displacement, shiftedCtx);

  Float displace = texEval(displacement, ctx);

  BumpMapResult result = DefaultBumpMapResult(ctx);
  result.dpdus =
    ctx.dpdus + ((uDisplace - displace) / du) * convert_to_vec3(ctx.ns) +
    displace * convert_to_vec3(ctx.dndus);
  result.dpdvs =
    ctx.dpdvs + ((vDisplace - displace) / dv) * convert_to_vec3(ctx.ns) +
    displace * convert_to_vec3(ctx.dndvs);
  result.shadingNormal = NormalFromShadingTangents(result.dpdus, result.dpdvs, ctx.ns);
  result.hasBump = true;
  return result;
}

class DiffuseMaterial {
public:
  using BxDF = render::DiffuseBxDF;

  DiffuseMaterial() = default;
  explicit DiffuseMaterial(SpectrumTexture reflectance);
  DiffuseMaterial(
    SpectrumTexture reflectance,
    std::optional<FloatTexture> alpha,
    std::optional<FloatTexture> bump
  );

  template <typename TextureEvaluator>
  BxDF GetBxDF(
    const TextureEvaluator& texEval,
    const MaterialEvalContext& ctx,
    base::SampledWavelengths& lambda
  ) const {
    return BxDF(base::Clamp(texEval(reflectance_, ctx, lambda), 0, 1));
  }

  template <typename TextureEvaluator>
  MaterialAlphaResult EvaluateAlpha(
    const TextureEvaluator& texEval,
    const MaterialEvalContext& ctx,
    Float alphaSample = 0
  ) const {
    if (!alpha_) {
      return {};
    }
    Float alpha = ClampUnit(texEval(*alpha_, ctx));
    return MaterialAlphaResult{alpha, alphaSample < alpha};
  }

  template <typename TextureEvaluator>
  BumpMapResult EvaluateBump(
    const TextureEvaluator& texEval,
    const MaterialEvalContext& ctx
  ) const {
    if (!bump_) {
      return DefaultBumpMapResult(ctx);
    }
    return EvaluateBumpMap(texEval, *bump_, ctx);
  }

  bool HasAlpha() const;
  bool HasBump() const;
  const SpectrumTexture& Reflectance() const;
  const std::optional<FloatTexture>& AlphaTexture() const;
  const std::optional<FloatTexture>& BumpTexture() const;
  MaterialTextureRequirements TextureRequirements() const;

private:
  static Float ClampUnit(Float value);

  SpectrumTexture reflectance_;
  std::optional<FloatTexture> alpha_;
  std::optional<FloatTexture> bump_;
};

class ConductorMaterial {
public:
  using BxDF = render::ConductorBxDF;

  static ConductorMaterial FromEtaK(
    SpectrumTexture eta,
    SpectrumTexture k,
    FloatTexture roughness = FloatTexture::Constant(0),
    bool remapRoughness = true
  );
  static ConductorMaterial FromEtaK(
    SpectrumTexture eta,
    SpectrumTexture k,
    FloatTexture uRoughness,
    FloatTexture vRoughness,
    bool remapRoughness,
    std::optional<FloatTexture> alpha = std::nullopt,
    std::optional<FloatTexture> bump = std::nullopt
  );
  static ConductorMaterial FromReflectance(
    SpectrumTexture reflectance,
    FloatTexture roughness = FloatTexture::Constant(0),
    bool remapRoughness = true
  );
  static ConductorMaterial FromReflectance(
    SpectrumTexture reflectance,
    FloatTexture uRoughness,
    FloatTexture vRoughness,
    bool remapRoughness,
    std::optional<FloatTexture> alpha = std::nullopt,
    std::optional<FloatTexture> bump = std::nullopt
  );

  template <typename TextureEvaluator>
  BxDF GetBxDF(
    const TextureEvaluator& texEval,
    const MaterialEvalContext& ctx,
    base::SampledWavelengths& lambda
  ) const {
    Float uRoughness = texEval(uRoughness_, ctx);
    Float vRoughness = texEval(vRoughness_, ctx);
    if (remapRoughness_) {
      uRoughness = render::TrowbridgeReitzDistribution::RoughnessToAlpha(uRoughness);
      vRoughness = render::TrowbridgeReitzDistribution::RoughnessToAlpha(vRoughness);
    }

    base::SampledSpectrum eta;
    base::SampledSpectrum k;
    if (eta_ && k_) {
      eta = texEval(*eta_, ctx, lambda);
      k = texEval(*k_, ctx, lambda);
    } else {
      base::SampledSpectrum r =
        base::Clamp(texEval(*reflectance_, ctx, lambda), 0, static_cast<Float>(0.9999));
      eta = base::SampledSpectrum(1);
      k = static_cast<Float>(2) * base::Sqrt(r) /
          base::SafeSqrt(base::ClampZero(base::SampledSpectrum(1) - r));
    }

    return BxDF(render::TrowbridgeReitzDistribution(uRoughness, vRoughness), eta, k);
  }

  template <typename TextureEvaluator>
  MaterialAlphaResult EvaluateAlpha(
    const TextureEvaluator& texEval,
    const MaterialEvalContext& ctx,
    Float alphaSample = 0
  ) const {
    if (!alpha_) {
      return {};
    }
    Float alpha = ClampUnit(texEval(*alpha_, ctx));
    return MaterialAlphaResult{alpha, alphaSample < alpha};
  }

  template <typename TextureEvaluator>
  BumpMapResult EvaluateBump(
    const TextureEvaluator& texEval,
    const MaterialEvalContext& ctx
  ) const {
    if (!bump_) {
      return DefaultBumpMapResult(ctx);
    }
    return EvaluateBumpMap(texEval, *bump_, ctx);
  }

  bool UsesEtaK() const;
  bool UsesReflectance() const;
  bool HasAlpha() const;
  bool HasBump() const;
  bool RemapRoughness() const;
  MaterialTextureRequirements TextureRequirements() const;

private:
  ConductorMaterial(
    std::optional<SpectrumTexture> eta,
    std::optional<SpectrumTexture> k,
    std::optional<SpectrumTexture> reflectance,
    FloatTexture uRoughness,
    FloatTexture vRoughness,
    bool remapRoughness,
    std::optional<FloatTexture> alpha,
    std::optional<FloatTexture> bump
  );

  static Float ClampUnit(Float value);

  std::optional<SpectrumTexture> eta_;
  std::optional<SpectrumTexture> k_;
  std::optional<SpectrumTexture> reflectance_;
  FloatTexture uRoughness_ = FloatTexture::Constant(0);
  FloatTexture vRoughness_ = FloatTexture::Constant(0);
  bool remapRoughness_ = true;
  std::optional<FloatTexture> alpha_;
  std::optional<FloatTexture> bump_;
};

class DielectricMaterial {
public:
  using BxDF = render::DielectricBxDF;

  static DielectricMaterial Smooth();
  static DielectricMaterial FromRoughness(
    FloatTexture roughness = FloatTexture::Constant(0),
    bool remapRoughness = true
  );
  static DielectricMaterial FromRoughness(
    FloatTexture uRoughness,
    FloatTexture vRoughness,
    bool remapRoughness,
    std::optional<FloatTexture> alpha = std::nullopt,
    std::optional<FloatTexture> bump = std::nullopt
  );

  template <typename TextureEvaluator>
  BxDF GetBxDF(
    const TextureEvaluator& texEval,
    const MaterialEvalContext& ctx,
    base::SampledWavelengths& lambda
  ) const {
    if (ctx.dielectric == nullptr) {
      throw std::runtime_error("DielectricMaterial requires a resolved dielectric interface");
    }
    if (!ctx.dielectric->ratioIsConstant) {
      lambda.TerminateSecondary();
    }
    Float uRoughness = texEval(uRoughness_, ctx);
    Float vRoughness = texEval(vRoughness_, ctx);
    if (remapRoughness_) {
      uRoughness = render::TrowbridgeReitzDistribution::RoughnessToAlpha(uRoughness);
      vRoughness = render::TrowbridgeReitzDistribution::RoughnessToAlpha(vRoughness);
    }
    render::TrowbridgeReitzDistribution distribution(uRoughness, vRoughness);
    if (!distribution.EffectivelySmooth()) {
      throw std::runtime_error("rough spectral dielectric materials are deferred to PR 18");
    }
    Float eta = ctx.dielectric->Eta();
    if (!(eta > 0) || !std::isfinite(eta)) {
      throw std::runtime_error("resolved dielectric eta must be positive finite");
    }
    return BxDF(eta, distribution);
  }

  template <typename TextureEvaluator>
  MaterialAlphaResult EvaluateAlpha(
    const TextureEvaluator& texEval,
    const MaterialEvalContext& ctx,
    Float alphaSample = 0
  ) const {
    if (!alpha_) {
      return {};
    }
    Float alpha = ClampUnit(texEval(*alpha_, ctx));
    return MaterialAlphaResult{alpha, alphaSample < alpha};
  }

  template <typename TextureEvaluator>
  BumpMapResult EvaluateBump(
    const TextureEvaluator& texEval,
    const MaterialEvalContext& ctx
  ) const {
    if (!bump_) {
      return DefaultBumpMapResult(ctx);
    }
    return EvaluateBumpMap(texEval, *bump_, ctx);
  }

  bool HasAlpha() const;
  bool HasBump() const;
  bool RemapRoughness() const;
  MaterialTextureRequirements TextureRequirements() const;

private:
  DielectricMaterial(
    FloatTexture uRoughness,
    FloatTexture vRoughness,
    bool remapRoughness,
    std::optional<FloatTexture> alpha,
    std::optional<FloatTexture> bump
  );

  static Float ClampUnit(Float value);

  FloatTexture uRoughness_ = FloatTexture::Constant(0);
  FloatTexture vRoughness_ = FloatTexture::Constant(0);
  bool remapRoughness_ = true;
  std::optional<FloatTexture> alpha_;
  std::optional<FloatTexture> bump_;
};

class InterfaceMaterial {
public:
  using BxDF = render::NullBxDF;

  template <typename TextureEvaluator>
  BxDF GetBxDF(
    const TextureEvaluator&,
    const MaterialEvalContext&,
    base::SampledWavelengths&
  ) const {
    return BxDF();
  }

  template <typename TextureEvaluator>
  MaterialAlphaResult EvaluateAlpha(
    const TextureEvaluator&,
    const MaterialEvalContext&,
    Float = 0
  ) const {
    return {};
  }

  template <typename TextureEvaluator>
  BumpMapResult EvaluateBump(
    const TextureEvaluator&,
    const MaterialEvalContext& ctx
  ) const {
    return DefaultBumpMapResult(ctx);
  }

  bool HasAlpha() const;
  bool HasBump() const;
  MaterialTextureRequirements TextureRequirements() const;
};

class Material {
public:
  Material() = default;

  static Material Diffuse(SpectrumTexture reflectance);
  static Material Diffuse(DiffuseMaterial material);
  static Material Conductor(ConductorMaterial material);
  static Material Dielectric(DielectricMaterial material = DielectricMaterial::Smooth());
  static Material Interface();

  bool IsValid() const;
  MaterialType Type() const;
  MaterialTextureRequirements TextureRequirements() const;
  bool CanEvaluateTextures(const UniversalTextureEvaluator& evaluator) const;

  template <typename TextureEvaluator>
  render::BSDF GetBSDF(
    const TextureEvaluator& texEval,
    MaterialEvalContext ctx,
    base::SampledWavelengths& lambda,
    base::ScratchBuffer& scratch
  ) const {
    return std::visit(
      [&](const auto& material) -> render::BSDF {
        using ConcreteMaterial = std::decay_t<decltype(material)>;
        if constexpr (std::is_same<ConcreteMaterial, std::monostate>::value) {
          return render::BSDF();
        } else {
          BumpMapResult bump = material.EvaluateBump(texEval, ctx);
          if (bump.hasBump) {
            ctx.ns = bump.shadingNormal;
            ctx.dpdus = bump.dpdus;
            ctx.dpdvs = bump.dpdvs;
            ctx.dpdu = bump.dpdus;
            ctx.dpdv = bump.dpdvs;
          }

          using ConcreteBxDF = typename ConcreteMaterial::BxDF;
          ConcreteBxDF value = material.GetBxDF(texEval, ctx, lambda);
          ConcreteBxDF* bxdf = scratch.Create<ConcreteBxDF>(std::move(value));
          return render::BSDF(ctx.n, ctx.ns, ctx.dpdus, render::BxDF(bxdf));
        }
      },
      material_
    );
  }

  template <typename TextureEvaluator>
  MaterialAlphaResult EvaluateAlpha(
    const TextureEvaluator& texEval,
    const MaterialEvalContext& ctx,
    Float alphaSample = 0
  ) const {
    return std::visit(
      [&](const auto& material) -> MaterialAlphaResult {
        using ConcreteMaterial = std::decay_t<decltype(material)>;
        if constexpr (std::is_same<ConcreteMaterial, std::monostate>::value) {
          return {};
        } else {
          return material.EvaluateAlpha(texEval, ctx, alphaSample);
        }
      },
      material_
    );
  }

  template <typename TextureEvaluator>
  BumpMapResult EvaluateBump(
    const TextureEvaluator& texEval,
    const MaterialEvalContext& ctx
  ) const {
    return std::visit(
      [&](const auto& material) -> BumpMapResult {
        using ConcreteMaterial = std::decay_t<decltype(material)>;
        if constexpr (std::is_same<ConcreteMaterial, std::monostate>::value) {
          return DefaultBumpMapResult(ctx);
        } else {
          return material.EvaluateBump(texEval, ctx);
        }
      },
      material_
    );
  }

private:
  explicit Material(DiffuseMaterial material);
  explicit Material(ConductorMaterial material);
  explicit Material(DielectricMaterial material);
  explicit Material(InterfaceMaterial material);

  using Variant = std::variant<
    std::monostate,
    DiffuseMaterial,
    ConductorMaterial,
    DielectricMaterial,
    InterfaceMaterial
  >;

  Variant material_;
};

class SpectralMaterialTable {
public:
  base::MaterialHandle Add(Material material);
  const Material& Get(base::MaterialHandle handle) const;
  Material& Get(base::MaterialHandle handle);
  std::size_t Size() const;

private:
  void Validate(base::MaterialHandle handle) const;

  std::vector<Material> materials_;
};

} // namespace materials
} // namespace rayrender

#endif
