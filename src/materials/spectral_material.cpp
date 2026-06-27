#include "spectral_material.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace rayrender {
namespace materials {
namespace {

Float Clamp(Float value, Float low, Float high) {
  return std::min(std::max(value, low), high);
}

bool NormalHasLength(normal3f n) {
  return n.squared_length() > 0 && std::isfinite(n.squared_length());
}

bool VectorHasLength(vec3f v) {
  return v.squared_length() > 0 && std::isfinite(v.squared_length());
}

vec3f UnitOrFallback(vec3f v, vec3f fallback) {
  Float length = v.length();
  if (length == 0 || !std::isfinite(length)) {
    return fallback;
  }
  return v / length;
}

normal3f UnitNormalOrFallback(normal3f n, normal3f fallback) {
  Float length = n.length();
  if (length == 0 || !std::isfinite(length)) {
    return fallback;
  }
  return n / length;
}

vec3f OrthogonalFallback(vec3f z) {
  if (std::abs(z.xyz.x) > std::abs(z.xyz.y)) {
    return UnitOrFallback(vec3f(-z.xyz.z, 0, z.xyz.x), vec3f(1, 0, 0));
  }
  return UnitOrFallback(vec3f(0, z.xyz.z, -z.xyz.y), vec3f(1, 0, 0));
}

} // namespace

const char* MaterialTypeName(MaterialType type) {
  switch (type) {
  case MaterialType::Diffuse:
    return "Diffuse";
  case MaterialType::Conductor:
    return "Conductor";
  case MaterialType::Interface:
    return "Interface";
  }
  return "Unknown";
}

MaterialEvalContext::MaterialEvalContext(const render::SurfaceInteraction& interaction)
  : TextureEvalContext() {
  p = interaction.p;
  uv = interaction.uv;
  dpdx = interaction.dpdx;
  dpdy = interaction.dpdy;
  dpdu = interaction.dpdu;
  dpdv = interaction.dpdv;
  dudx = interaction.dudx;
  dvdx = interaction.dvdx;
  dudy = interaction.dudy;
  dvdy = interaction.dvdy;
  faceIndex = interaction.faceIndex;

  wo = interaction.wo;
  n = UnitNormalOrFallback(interaction.n, normal3f(0, 0, 1));
  ns = UnitNormalOrFallback(
    NormalHasLength(interaction.shadingNormal) ? interaction.shadingNormal : interaction.n,
    n
  );
  dpdus = VectorHasLength(interaction.shadingDpdu) ? interaction.shadingDpdu : interaction.dpdu;
  dpdvs = VectorHasLength(interaction.shadingDpdv) ? interaction.shadingDpdv : interaction.dpdv;
  if (!VectorHasLength(dpdus)) {
    dpdus = OrthogonalFallback(convert_to_vec3(ns));
  }
  if (!VectorHasLength(dpdvs)) {
    dpdvs = cross(convert_to_vec3(ns), UnitOrFallback(dpdus, vec3f(1, 0, 0)));
  }
  dndus = NormalHasLength(interaction.shadingDndu) ? interaction.shadingDndu : interaction.dndu;
  dndvs = NormalHasLength(interaction.shadingDndv) ? interaction.shadingDndv : interaction.dndv;
}

MaterialEvalContext MaterialEvalContextFromSurfaceInteraction(
  const render::SurfaceInteraction& interaction
) {
  return MaterialEvalContext(interaction);
}

bool MaterialTextureRequirements::Any() const {
  return spectrumTextures || floatTextures || alphaTexture || bumpTexture;
}

BumpMapResult DefaultBumpMapResult(const MaterialEvalContext& ctx) {
  BumpMapResult result;
  result.geometricNormal = ctx.n;
  result.shadingNormal = ctx.ns;
  result.dpdus = ctx.dpdus;
  result.dpdvs = ctx.dpdvs;
  result.hasBump = false;
  return result;
}

normal3f NormalFromShadingTangents(
  const vec3f& dpdus,
  const vec3f& dpdvs,
  normal3f referenceNormal
) {
  vec3f n = cross(dpdus, dpdvs);
  if (!VectorHasLength(n)) {
    return UnitNormalOrFallback(referenceNormal, normal3f(0, 0, 1));
  }
  n = UnitOrFallback(n, convert_to_vec3(referenceNormal));
  if (dot(n, referenceNormal) < 0) {
    n = -n;
  }
  return convert_to_normal3(n);
}

DiffuseMaterial::DiffuseMaterial(SpectrumTexture reflectance)
  : reflectance_(std::move(reflectance)) {}

DiffuseMaterial::DiffuseMaterial(
  SpectrumTexture reflectance,
  std::optional<FloatTexture> alpha,
  std::optional<FloatTexture> bump
)
  : reflectance_(std::move(reflectance)),
    alpha_(std::move(alpha)),
    bump_(std::move(bump)) {}

bool DiffuseMaterial::HasAlpha() const {
  return alpha_.has_value();
}

bool DiffuseMaterial::HasBump() const {
  return bump_.has_value();
}

const SpectrumTexture& DiffuseMaterial::Reflectance() const {
  return reflectance_;
}

const std::optional<FloatTexture>& DiffuseMaterial::AlphaTexture() const {
  return alpha_;
}

const std::optional<FloatTexture>& DiffuseMaterial::BumpTexture() const {
  return bump_;
}

MaterialTextureRequirements DiffuseMaterial::TextureRequirements() const {
  MaterialTextureRequirements requirements;
  requirements.spectrumTextures = true;
  requirements.floatTextures = alpha_.has_value() || bump_.has_value();
  requirements.alphaTexture = alpha_.has_value();
  requirements.bumpTexture = bump_.has_value();
  return requirements;
}

Float DiffuseMaterial::ClampUnit(Float value) {
  return Clamp(value, 0, 1);
}

ConductorMaterial ConductorMaterial::FromEtaK(
  SpectrumTexture eta,
  SpectrumTexture k,
  FloatTexture roughness,
  bool remapRoughness
) {
  return FromEtaK(
    std::move(eta),
    std::move(k),
    roughness,
    std::move(roughness),
    remapRoughness
  );
}

ConductorMaterial ConductorMaterial::FromEtaK(
  SpectrumTexture eta,
  SpectrumTexture k,
  FloatTexture uRoughness,
  FloatTexture vRoughness,
  bool remapRoughness,
  std::optional<FloatTexture> alpha,
  std::optional<FloatTexture> bump
) {
  return ConductorMaterial(
    std::move(eta),
    std::move(k),
    std::nullopt,
    std::move(uRoughness),
    std::move(vRoughness),
    remapRoughness,
    std::move(alpha),
    std::move(bump)
  );
}

ConductorMaterial ConductorMaterial::FromReflectance(
  SpectrumTexture reflectance,
  FloatTexture roughness,
  bool remapRoughness
) {
  return FromReflectance(
    std::move(reflectance),
    roughness,
    std::move(roughness),
    remapRoughness
  );
}

ConductorMaterial ConductorMaterial::FromReflectance(
  SpectrumTexture reflectance,
  FloatTexture uRoughness,
  FloatTexture vRoughness,
  bool remapRoughness,
  std::optional<FloatTexture> alpha,
  std::optional<FloatTexture> bump
) {
  return ConductorMaterial(
    std::nullopt,
    std::nullopt,
    std::move(reflectance),
    std::move(uRoughness),
    std::move(vRoughness),
    remapRoughness,
    std::move(alpha),
    std::move(bump)
  );
}

ConductorMaterial::ConductorMaterial(
  std::optional<SpectrumTexture> eta,
  std::optional<SpectrumTexture> k,
  std::optional<SpectrumTexture> reflectance,
  FloatTexture uRoughness,
  FloatTexture vRoughness,
  bool remapRoughness,
  std::optional<FloatTexture> alpha,
  std::optional<FloatTexture> bump
)
  : eta_(std::move(eta)),
    k_(std::move(k)),
    reflectance_(std::move(reflectance)),
    uRoughness_(std::move(uRoughness)),
    vRoughness_(std::move(vRoughness)),
    remapRoughness_(remapRoughness),
    alpha_(std::move(alpha)),
    bump_(std::move(bump)) {
  bool hasEtaK = eta_.has_value() || k_.has_value();
  if (hasEtaK && !(eta_.has_value() && k_.has_value())) {
    throw std::invalid_argument("ConductorMaterial requires both eta and k spectra");
  }
  if (reflectance_.has_value() == (eta_.has_value() && k_.has_value())) {
    throw std::invalid_argument("ConductorMaterial requires either eta/k or reflectance");
  }
}

bool ConductorMaterial::UsesEtaK() const {
  return eta_.has_value() && k_.has_value();
}

bool ConductorMaterial::UsesReflectance() const {
  return reflectance_.has_value();
}

bool ConductorMaterial::HasAlpha() const {
  return alpha_.has_value();
}

bool ConductorMaterial::HasBump() const {
  return bump_.has_value();
}

bool ConductorMaterial::RemapRoughness() const {
  return remapRoughness_;
}

MaterialTextureRequirements ConductorMaterial::TextureRequirements() const {
  MaterialTextureRequirements requirements;
  requirements.spectrumTextures = true;
  requirements.floatTextures = true;
  requirements.alphaTexture = alpha_.has_value();
  requirements.bumpTexture = bump_.has_value();
  return requirements;
}

Float ConductorMaterial::ClampUnit(Float value) {
  return Clamp(value, 0, 1);
}

bool InterfaceMaterial::HasAlpha() const {
  return false;
}

bool InterfaceMaterial::HasBump() const {
  return false;
}

MaterialTextureRequirements InterfaceMaterial::TextureRequirements() const {
  return {};
}

Material Material::Diffuse(SpectrumTexture reflectance) {
  return Material(DiffuseMaterial(std::move(reflectance)));
}

Material Material::Diffuse(DiffuseMaterial material) {
  return Material(std::move(material));
}

Material Material::Conductor(ConductorMaterial material) {
  return Material(std::move(material));
}

Material Material::Interface() {
  return Material(InterfaceMaterial());
}

Material::Material(DiffuseMaterial material) : material_(std::move(material)) {}

Material::Material(ConductorMaterial material) : material_(std::move(material)) {}

Material::Material(InterfaceMaterial material) : material_(material) {}

bool Material::IsValid() const {
  return !std::holds_alternative<std::monostate>(material_);
}

MaterialType Material::Type() const {
  if (std::holds_alternative<DiffuseMaterial>(material_)) {
    return MaterialType::Diffuse;
  }
  if (std::holds_alternative<ConductorMaterial>(material_)) {
    return MaterialType::Conductor;
  }
  if (std::holds_alternative<InterfaceMaterial>(material_)) {
    return MaterialType::Interface;
  }
  throw std::runtime_error("Invalid spectral Material has no type");
}

MaterialTextureRequirements Material::TextureRequirements() const {
  return std::visit(
    [](const auto& material) -> MaterialTextureRequirements {
      using ConcreteMaterial = std::decay_t<decltype(material)>;
      if constexpr (std::is_same<ConcreteMaterial, std::monostate>::value) {
        return {};
      } else {
        return material.TextureRequirements();
      }
    },
    material_
  );
}

bool Material::CanEvaluateTextures(const UniversalTextureEvaluator&) const {
  return IsValid();
}

base::MaterialHandle SpectralMaterialTable::Add(Material material) {
  materials_.push_back(std::move(material));
  return base::MaterialHandle::FromIndex(
    static_cast<base::MaterialHandle::IndexType>(materials_.size() - 1),
    1
  );
}

const Material& SpectralMaterialTable::Get(base::MaterialHandle handle) const {
  Validate(handle);
  return materials_[handle.Index()];
}

Material& SpectralMaterialTable::Get(base::MaterialHandle handle) {
  Validate(handle);
  return materials_[handle.Index()];
}

std::size_t SpectralMaterialTable::Size() const {
  return materials_.size();
}

void SpectralMaterialTable::Validate(base::MaterialHandle handle) const {
  if (!handle.IsValid() || handle.Index() >= materials_.size()) {
    throw std::out_of_range("Material handle is invalid");
  }
  if (handle.Generation() != 1) {
    throw std::out_of_range("Material handle generation is stale");
  }
}

} // namespace materials
} // namespace rayrender
