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
  case MaterialType::DiffuseTransmission:
    return "DiffuseTransmission";
  case MaterialType::Conductor:
    return "Conductor";
  case MaterialType::Dielectric:
    return "Dielectric";
  case MaterialType::CoatedDiffuse:
    return "CoatedDiffuse";
  case MaterialType::CoatedConductor:
    return "CoatedConductor";
  case MaterialType::ThinDielectric:
    return "ThinDielectric";
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

DiffuseTransmissionMaterial::DiffuseTransmissionMaterial(
  SpectrumTexture reflectance,
  SpectrumTexture transmittance
)
  : reflectance_(std::move(reflectance)),
    transmittance_(std::move(transmittance)) {}

DiffuseTransmissionMaterial::DiffuseTransmissionMaterial(
  SpectrumTexture reflectance,
  SpectrumTexture transmittance,
  Float scale,
  std::optional<FloatTexture> alpha,
  std::optional<FloatTexture> bump
)
  : reflectance_(std::move(reflectance)),
    transmittance_(std::move(transmittance)),
    scale_(scale),
    alpha_(std::move(alpha)),
    bump_(std::move(bump)) {}

bool DiffuseTransmissionMaterial::HasAlpha() const {
  return alpha_.has_value();
}

bool DiffuseTransmissionMaterial::HasBump() const {
  return bump_.has_value();
}

MaterialTextureRequirements DiffuseTransmissionMaterial::TextureRequirements() const {
  MaterialTextureRequirements requirements;
  requirements.spectrumTextures = true;
  requirements.floatTextures = alpha_.has_value() || bump_.has_value();
  requirements.alphaTexture = alpha_.has_value();
  requirements.bumpTexture = bump_.has_value();
  return requirements;
}

Float DiffuseTransmissionMaterial::ClampUnit(Float value) {
  return Clamp(value, 0, 1);
}

ConductorMaterial ConductorMaterial::FromEtaK(
  SpectrumTexture eta,
  SpectrumTexture k,
  FloatTexture roughness,
  bool remapRoughness
) {
  FloatTexture vRoughness = roughness;
  return FromEtaK(
    std::move(eta),
    std::move(k),
    roughness,
    std::move(vRoughness),
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
  FloatTexture vRoughness = roughness;
  return FromReflectance(
    std::move(reflectance),
    roughness,
    std::move(vRoughness),
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

DielectricMaterial DielectricMaterial::Smooth() {
  return FromRoughness(FloatTexture::Constant(0), true);
}

DielectricMaterial DielectricMaterial::FromRoughness(
  FloatTexture roughness,
  bool remapRoughness
) {
  FloatTexture vRoughness = roughness;
  return FromRoughness(
    roughness,
    std::move(vRoughness),
    remapRoughness
  );
}

DielectricMaterial DielectricMaterial::FromRoughness(
  FloatTexture uRoughness,
  FloatTexture vRoughness,
  bool remapRoughness,
  std::optional<FloatTexture> alpha,
  std::optional<FloatTexture> bump
) {
  return DielectricMaterial(
    std::move(uRoughness),
    std::move(vRoughness),
    remapRoughness,
    std::move(alpha),
    std::move(bump)
  );
}

DielectricMaterial::DielectricMaterial(
  FloatTexture uRoughness,
  FloatTexture vRoughness,
  bool remapRoughness,
  std::optional<FloatTexture> alpha,
  std::optional<FloatTexture> bump
)
  : uRoughness_(std::move(uRoughness)),
    vRoughness_(std::move(vRoughness)),
    remapRoughness_(remapRoughness),
    alpha_(std::move(alpha)),
    bump_(std::move(bump)) {}

bool DielectricMaterial::HasAlpha() const {
  return alpha_.has_value();
}

bool DielectricMaterial::HasBump() const {
  return bump_.has_value();
}

bool DielectricMaterial::RemapRoughness() const {
  return remapRoughness_;
}

MaterialTextureRequirements DielectricMaterial::TextureRequirements() const {
  MaterialTextureRequirements requirements;
  requirements.floatTextures = true;
  requirements.alphaTexture = alpha_.has_value();
  requirements.bumpTexture = bump_.has_value();
  return requirements;
}

Float DielectricMaterial::ClampUnit(Float value) {
  return Clamp(value, 0, 1);
}

CoatedDiffuseMaterial::CoatedDiffuseMaterial(
  SpectrumTexture reflectance,
  FloatTexture roughness,
  FloatTexture thickness,
  SpectrumTexture albedo,
  FloatTexture g,
  render::EtaSpectrumHandle eta,
  bool remapRoughness,
  int maxDepth,
  int nSamples
) : CoatedDiffuseMaterial(
      std::move(reflectance),
      roughness,
      std::move(roughness),
      std::move(thickness),
      std::move(albedo),
      std::move(g),
      std::move(eta),
      remapRoughness,
      maxDepth,
      nSamples
    ) {}

CoatedDiffuseMaterial::CoatedDiffuseMaterial(
  SpectrumTexture reflectance,
  FloatTexture uRoughness,
  FloatTexture vRoughness,
  FloatTexture thickness,
  SpectrumTexture albedo,
  FloatTexture g,
  render::EtaSpectrumHandle eta,
  bool remapRoughness,
  int maxDepth,
  int nSamples,
  std::optional<FloatTexture> alpha,
  std::optional<FloatTexture> bump
)
  : reflectance_(std::move(reflectance)),
    uRoughness_(std::move(uRoughness)),
    vRoughness_(std::move(vRoughness)),
    thickness_(std::move(thickness)),
    albedo_(std::move(albedo)),
    g_(std::move(g)),
    eta_(std::move(eta)),
    remapRoughness_(remapRoughness),
    maxDepth_(std::max(1, maxDepth)),
    nSamples_(std::max(1, nSamples)),
    alpha_(std::move(alpha)),
    bump_(std::move(bump)) {
  if (!eta_) {
    throw std::invalid_argument("CoatedDiffuseMaterial requires an eta spectrum");
  }
}

bool CoatedDiffuseMaterial::HasAlpha() const {
  return alpha_.has_value();
}

bool CoatedDiffuseMaterial::HasBump() const {
  return bump_.has_value();
}

bool CoatedDiffuseMaterial::RemapRoughness() const {
  return remapRoughness_;
}

MaterialTextureRequirements CoatedDiffuseMaterial::TextureRequirements() const {
  MaterialTextureRequirements requirements;
  requirements.spectrumTextures = true;
  requirements.floatTextures = true;
  requirements.alphaTexture = alpha_.has_value();
  requirements.bumpTexture = bump_.has_value();
  return requirements;
}

Float CoatedDiffuseMaterial::ClampUnit(Float value) {
  return Clamp(value, 0, 1);
}

Float CoatedDiffuseMaterial::ClampSignedUnit(Float value) {
  return Clamp(value, -1, 1);
}

CoatedConductorMaterial CoatedConductorMaterial::FromEtaK(
  SpectrumTexture eta,
  SpectrumTexture k,
  FloatTexture roughness,
  bool remapRoughness
) {
  FloatTexture interfaceVRoughness = roughness;
  FloatTexture conductorURoughness = roughness;
  FloatTexture conductorVRoughness = roughness;
  return FromEtaK(
    std::move(eta),
    std::move(k),
    roughness,
    std::move(interfaceVRoughness),
    std::move(conductorURoughness),
    std::move(conductorVRoughness),
    FloatTexture::Constant(static_cast<Float>(0.01)),
    SpectrumTexture::Constant(0),
    FloatTexture::Constant(0),
    render::ConstantEtaSpectrum(static_cast<Float>(1.5)),
    remapRoughness
  );
}

CoatedConductorMaterial CoatedConductorMaterial::FromEtaK(
  SpectrumTexture eta,
  SpectrumTexture k,
  FloatTexture interfaceURoughness,
  FloatTexture interfaceVRoughness,
  FloatTexture conductorURoughness,
  FloatTexture conductorVRoughness,
  FloatTexture thickness,
  SpectrumTexture albedo,
  FloatTexture g,
  render::EtaSpectrumHandle interfaceEta,
  bool remapRoughness,
  int maxDepth,
  int nSamples,
  std::optional<FloatTexture> alpha,
  std::optional<FloatTexture> bump
) {
  return CoatedConductorMaterial(
    std::move(eta),
    std::move(k),
    std::nullopt,
    std::move(interfaceURoughness),
    std::move(interfaceVRoughness),
    std::move(conductorURoughness),
    std::move(conductorVRoughness),
    std::move(thickness),
    std::move(albedo),
    std::move(g),
    std::move(interfaceEta),
    remapRoughness,
    maxDepth,
    nSamples,
    std::move(alpha),
    std::move(bump)
  );
}

CoatedConductorMaterial CoatedConductorMaterial::FromReflectance(
  SpectrumTexture reflectance,
  FloatTexture roughness,
  bool remapRoughness
) {
  FloatTexture interfaceVRoughness = roughness;
  FloatTexture conductorURoughness = roughness;
  FloatTexture conductorVRoughness = roughness;
  return FromReflectance(
    std::move(reflectance),
    roughness,
    std::move(interfaceVRoughness),
    std::move(conductorURoughness),
    std::move(conductorVRoughness),
    FloatTexture::Constant(static_cast<Float>(0.01)),
    SpectrumTexture::Constant(0),
    FloatTexture::Constant(0),
    render::ConstantEtaSpectrum(static_cast<Float>(1.5)),
    remapRoughness
  );
}

CoatedConductorMaterial CoatedConductorMaterial::FromReflectance(
  SpectrumTexture reflectance,
  FloatTexture interfaceURoughness,
  FloatTexture interfaceVRoughness,
  FloatTexture conductorURoughness,
  FloatTexture conductorVRoughness,
  FloatTexture thickness,
  SpectrumTexture albedo,
  FloatTexture g,
  render::EtaSpectrumHandle interfaceEta,
  bool remapRoughness,
  int maxDepth,
  int nSamples,
  std::optional<FloatTexture> alpha,
  std::optional<FloatTexture> bump
) {
  return CoatedConductorMaterial(
    std::nullopt,
    std::nullopt,
    std::move(reflectance),
    std::move(interfaceURoughness),
    std::move(interfaceVRoughness),
    std::move(conductorURoughness),
    std::move(conductorVRoughness),
    std::move(thickness),
    std::move(albedo),
    std::move(g),
    std::move(interfaceEta),
    remapRoughness,
    maxDepth,
    nSamples,
    std::move(alpha),
    std::move(bump)
  );
}

CoatedConductorMaterial::CoatedConductorMaterial(
  std::optional<SpectrumTexture> conductorEta,
  std::optional<SpectrumTexture> conductorK,
  std::optional<SpectrumTexture> reflectance,
  FloatTexture interfaceURoughness,
  FloatTexture interfaceVRoughness,
  FloatTexture conductorURoughness,
  FloatTexture conductorVRoughness,
  FloatTexture thickness,
  SpectrumTexture albedo,
  FloatTexture g,
  render::EtaSpectrumHandle interfaceEta,
  bool remapRoughness,
  int maxDepth,
  int nSamples,
  std::optional<FloatTexture> alpha,
  std::optional<FloatTexture> bump
)
  : conductorEta_(std::move(conductorEta)),
    conductorK_(std::move(conductorK)),
    reflectance_(std::move(reflectance)),
    interfaceURoughness_(std::move(interfaceURoughness)),
    interfaceVRoughness_(std::move(interfaceVRoughness)),
    conductorURoughness_(std::move(conductorURoughness)),
    conductorVRoughness_(std::move(conductorVRoughness)),
    thickness_(std::move(thickness)),
    albedo_(std::move(albedo)),
    g_(std::move(g)),
    interfaceEta_(std::move(interfaceEta)),
    remapRoughness_(remapRoughness),
    maxDepth_(std::max(1, maxDepth)),
    nSamples_(std::max(1, nSamples)),
    alpha_(std::move(alpha)),
    bump_(std::move(bump)) {
  bool hasEtaK = conductorEta_.has_value() || conductorK_.has_value();
  if (hasEtaK && !(conductorEta_.has_value() && conductorK_.has_value())) {
    throw std::invalid_argument("CoatedConductorMaterial requires both eta and k spectra");
  }
  if (reflectance_.has_value() == (conductorEta_.has_value() && conductorK_.has_value())) {
    throw std::invalid_argument("CoatedConductorMaterial requires either eta/k or reflectance");
  }
  if (!interfaceEta_) {
    throw std::invalid_argument("CoatedConductorMaterial requires an interface eta spectrum");
  }
}

bool CoatedConductorMaterial::UsesEtaK() const {
  return conductorEta_.has_value() && conductorK_.has_value();
}

bool CoatedConductorMaterial::UsesReflectance() const {
  return reflectance_.has_value();
}

bool CoatedConductorMaterial::HasAlpha() const {
  return alpha_.has_value();
}

bool CoatedConductorMaterial::HasBump() const {
  return bump_.has_value();
}

bool CoatedConductorMaterial::RemapRoughness() const {
  return remapRoughness_;
}

MaterialTextureRequirements CoatedConductorMaterial::TextureRequirements() const {
  MaterialTextureRequirements requirements;
  requirements.spectrumTextures = true;
  requirements.floatTextures = true;
  requirements.alphaTexture = alpha_.has_value();
  requirements.bumpTexture = bump_.has_value();
  return requirements;
}

Float CoatedConductorMaterial::ClampUnit(Float value) {
  return Clamp(value, 0, 1);
}

Float CoatedConductorMaterial::ClampSignedUnit(Float value) {
  return Clamp(value, -1, 1);
}

ThinDielectricMaterial::ThinDielectricMaterial(render::EtaSpectrumHandle eta)
  : eta_(std::move(eta)) {
  if (!eta_) {
    throw std::invalid_argument("ThinDielectricMaterial requires an eta spectrum");
  }
}

bool ThinDielectricMaterial::HasAlpha() const {
  return false;
}

bool ThinDielectricMaterial::HasBump() const {
  return false;
}

MaterialTextureRequirements ThinDielectricMaterial::TextureRequirements() const {
  return {};
}

const render::EtaSpectrumHandle& ThinDielectricMaterial::Eta() const {
  return eta_;
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

Material Material::DiffuseTransmission(DiffuseTransmissionMaterial material) {
  return Material(std::move(material));
}

Material Material::Conductor(ConductorMaterial material) {
  return Material(std::move(material));
}

Material Material::Dielectric(DielectricMaterial material) {
  return Material(std::move(material));
}

Material Material::CoatedDiffuse(CoatedDiffuseMaterial material) {
  return Material(std::move(material));
}

Material Material::CoatedConductor(CoatedConductorMaterial material) {
  return Material(std::move(material));
}

Material Material::ThinDielectric(ThinDielectricMaterial material) {
  return Material(std::move(material));
}

Material Material::Interface() {
  return Material(InterfaceMaterial());
}

Material::Material(DiffuseMaterial material) : material_(std::move(material)) {}

Material::Material(DiffuseTransmissionMaterial material) : material_(std::move(material)) {}

Material::Material(ConductorMaterial material) : material_(std::move(material)) {}

Material::Material(DielectricMaterial material) : material_(std::move(material)) {}

Material::Material(CoatedDiffuseMaterial material) : material_(std::move(material)) {}

Material::Material(CoatedConductorMaterial material) : material_(std::move(material)) {}

Material::Material(ThinDielectricMaterial material) : material_(std::move(material)) {}

Material::Material(InterfaceMaterial material) : material_(material) {}

bool Material::IsValid() const {
  return !std::holds_alternative<std::monostate>(material_);
}

MaterialType Material::Type() const {
  if (std::holds_alternative<DiffuseMaterial>(material_)) {
    return MaterialType::Diffuse;
  }
  if (std::holds_alternative<DiffuseTransmissionMaterial>(material_)) {
    return MaterialType::DiffuseTransmission;
  }
  if (std::holds_alternative<ConductorMaterial>(material_)) {
    return MaterialType::Conductor;
  }
  if (std::holds_alternative<DielectricMaterial>(material_)) {
    return MaterialType::Dielectric;
  }
  if (std::holds_alternative<CoatedDiffuseMaterial>(material_)) {
    return MaterialType::CoatedDiffuse;
  }
  if (std::holds_alternative<CoatedConductorMaterial>(material_)) {
    return MaterialType::CoatedConductor;
  }
  if (std::holds_alternative<ThinDielectricMaterial>(material_)) {
    return MaterialType::ThinDielectric;
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
