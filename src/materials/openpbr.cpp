#include "openpbr.h"
#include "../volumes/path_diagnostics.h"
#include <glm/glm.hpp>
#include <algorithm>
#include <cassert>
#include <cstdint>

// Namespace isolation avoids collisions with rayrender's vec3 and math helpers.
// Assertions in worker-side reference code terminate only the affected ray.
#define OPENPBR_ASSERT(expr, message) do { if (!(expr)) throw PathFailure(PathFailureKind::Radiance, message); } while (false)
#define OPENPBR_ASSERT_UNREACHABLE(message) throw PathFailure(PathFailureKind::Radiance, message)
namespace openpbr_reference {
#include <openpbr/openpbr.h>
}

namespace {
namespace ref = openpbr_reference;
glm::vec3 vector(const vec3f& v) { return {v[0], v[1], v[2]}; }
glm::vec3 vector(const normal3f& v) { return {v[0], v[1], v[2]}; }
glm::vec3 vector(const point3f& v) { return {v[0], v[1], v[2]}; }
point3f color(const glm::vec3& v) { return {v.x, v.y, v.z}; }

// R is the authoring boundary; validate again when loading serialized descriptors
// so malformed values cannot reach reference-library assertions or LUT indices.
ref::OpenPBR_ResolvedInputs inputs(const Rcpp::List& parameters) {
  auto result = ref::openpbr_make_default_resolved_inputs();
  auto scalar = [&](const char* name, double low, double high, bool positive = false) {
    Rcpp::NumericVector value = Rcpp::as<Rcpp::NumericVector>(parameters[name]);
    if (value.size() != 1 || !std::isfinite(value[0]) || !std::isfinite(float(value[0])) ||
        value[0] < low || value[0] > high || (positive && !(float(value[0]) > 0)))
      throw std::runtime_error(std::string("Invalid OpenPBR parameter: ") + name);
    return float(value[0]);
  };
  auto rgb = [&](const char* name) {
    Rcpp::NumericVector value = Rcpp::as<Rcpp::NumericVector>(parameters[name]);
    if (value.size() != 3) throw std::runtime_error(std::string("OpenPBR RGB triple required: ") + name);
    glm::vec3 value3;
    for (int c = 0; c < 3; ++c) {
      if (!std::isfinite(value[c]) || value[c] < 0 || value[c] > 1)
        throw std::runtime_error(std::string("Invalid OpenPBR color: ") + name);
      value3[c] = float(value[c]);
    }
    return value3;
  };
#define OPENPBR_SCALAR(name, lo, hi) result.name = scalar(#name, lo, hi)
#define OPENPBR_COLOR(name) result.name = rgb(#name)
  OPENPBR_SCALAR(base_weight, 0, 1);
  OPENPBR_COLOR(base_color);
  OPENPBR_SCALAR(base_diffuse_roughness, 0, 1);
  OPENPBR_SCALAR(base_metalness, 0, 1);
  OPENPBR_SCALAR(specular_weight, 0, INFINITY);
  OPENPBR_COLOR(specular_color);
  OPENPBR_SCALAR(specular_roughness, 0, 1);
  result.specular_ior = scalar("specular_ior", 0, INFINITY, true);
  OPENPBR_SCALAR(specular_roughness_anisotropy, 0, 1);
  OPENPBR_SCALAR(transmission_weight, 0, 1);
  OPENPBR_COLOR(transmission_color);
  OPENPBR_SCALAR(transmission_depth, 0, INFINITY);
  OPENPBR_COLOR(transmission_scatter);
  OPENPBR_SCALAR(transmission_scatter_anisotropy, -1, 1);
  OPENPBR_SCALAR(transmission_dispersion_scale, 0, 1);
  result.transmission_dispersion_abbe_number = scalar("transmission_dispersion_abbe_number", 0, INFINITY, true);
  OPENPBR_SCALAR(subsurface_weight, 0, 1);
  OPENPBR_COLOR(subsurface_color);
  OPENPBR_SCALAR(subsurface_radius, 0, INFINITY);
  OPENPBR_COLOR(subsurface_radius_scale);
  OPENPBR_SCALAR(subsurface_scatter_anisotropy, -1, 1);
  OPENPBR_SCALAR(fuzz_weight, 0, 1);
  OPENPBR_COLOR(fuzz_color);
  OPENPBR_SCALAR(fuzz_roughness, 0, 1);
  OPENPBR_SCALAR(coat_weight, 0, 1);
  OPENPBR_COLOR(coat_color);
  OPENPBR_SCALAR(coat_roughness, 0, 1);
  OPENPBR_SCALAR(coat_roughness_anisotropy, 0, 1);
  result.coat_ior = scalar("coat_ior", 0, INFINITY, true);
  OPENPBR_SCALAR(coat_darkening, 0, 1);
  OPENPBR_SCALAR(thin_film_weight, 0, 1);
  OPENPBR_SCALAR(thin_film_thickness, 0, INFINITY);
  result.thin_film_ior = scalar("thin_film_ior", 0, INFINITY, true);
  OPENPBR_SCALAR(emission_luminance, 0, INFINITY);
  OPENPBR_COLOR(emission_color);
  OPENPBR_SCALAR(geometry_opacity, 0, 1);
#undef OPENPBR_SCALAR
#undef OPENPBR_COLOR
  SEXP thin = parameters["geometry_thin_walled"];
  if (TYPEOF(thin) != LGLSXP || Rf_xlength(thin) != 1 || LOGICAL(thin)[0] == NA_LOGICAL)
    throw std::runtime_error("OpenPBR geometry_thin_walled must be TRUE or FALSE.");
  result.geometry_thin_walled = LOGICAL(thin)[0];
  return result;
}

// Use the UV tangent when available, with Gram-Schmidt orthogonalization. A
// degenerate tangent uses the reference library's deterministic normal-only frame.
ref::OpenPBR_Basis basis(glm::vec3 n, glm::vec3 tangent) {
  n = glm::normalize(n);
  tangent -= n * glm::dot(n, tangent);
  return glm::dot(tangent, tangent) > 1e-12f
      ? ref::openpbr_make_basis(n, glm::normalize(tangent), 1.f)
      : ref::openpbr_make_basis(n);
}
}

struct OpenPBRMaterial::Impl {
  ref::OpenPBR_ResolvedInputs parameters;
  std::shared_ptr<texture> base;
  std::shared_ptr<roughness_texture> roughness;
  bool has_roughness;
  point2f texture_repeat{1, 1};
  // Optional world-space geometry overrides. Zero means use the hit geometry.
  glm::vec3 normal{0}, tangent{0}, coat_normal{0}, coat_tangent{0};
};
struct OpenPBRInteraction::Impl {
  ref::OpenPBR_PreparedBsdf bsdf;
  glm::vec3 geometric_normal, shading_normal, view;
  glm::vec3 transmission_radiance_scale{1};
  Float eta_squared = 1;
  bool Supports(const glm::vec3& direction) const {
    return (glm::dot(view, geometric_normal) * glm::dot(direction, geometric_normal) < 0) ==
           (glm::dot(view, shading_normal) * glm::dot(direction, shading_normal) < 0);
  }
};

OpenPBRMaterial::OpenPBRMaterial(const Rcpp::List& parameters, std::shared_ptr<texture> base,
                               std::shared_ptr<roughness_texture> roughness, bool has_roughness,
                               point2f texture_repeat)
    : dielectric(point3f(1), 1.5f, point3f(0), 0), impl(new Impl) {
  impl->parameters = inputs(parameters);
  impl->base = std::move(base);
  impl->roughness = std::move(roughness);
  impl->has_roughness = has_roughness;
  impl->texture_repeat = texture_repeat;
  ref_idx = impl->parameters.specular_ior;
  double p = Rcpp::as<double>(parameters["priority"]);
  if (!std::isfinite(p) || p < 0 || p > INT_MAX || p != std::floor(p))
    throw std::runtime_error("OpenPBR priority must be a nonnegative integer.");
  priority = size_t(p);
  for (auto item : {std::make_pair("geometry_normal", &impl->normal),
                    std::make_pair("geometry_tangent", &impl->tangent),
                    std::make_pair("geometry_coat_normal", &impl->coat_normal),
                    std::make_pair("geometry_coat_tangent", &impl->coat_tangent)}) {
    if (!parameters.containsElementNamed(item.first) || Rf_isNull(parameters[item.first])) continue;
    Rcpp::NumericVector v = parameters[item.first];
    if (v.size() != 3) throw std::runtime_error("OpenPBR geometry vectors require three components.");
    for (int c = 0; c < 3; ++c) {
      if (!std::isfinite(v[c]) || !std::isfinite(float(v[c])))
        throw std::runtime_error("OpenPBR geometry vectors must be finite.");
      (*item.second)[c] = v[c];
    }
    if (!(glm::dot(*item.second, *item.second) > 0))
      throw std::runtime_error("OpenPBR geometry vectors must be nonzero.");
  }
}
OpenPBRMaterial::~OpenPBRMaterial() = default;
size_t OpenPBRMaterial::GetSize() { return sizeof(*this) + sizeof(Impl); }
bool OpenPBRMaterial::is_dielectric() const {
  const auto& p = impl->parameters;
  return !p.geometry_thin_walled && p.base_metalness < 1 &&
         (p.transmission_weight > 0 || p.subsurface_weight > 0);
}
point3f OpenPBRMaterial::get_albedo(const hit_record& h) const {
  point3f value = impl->base->value(h.u, h.v, h.p);
  for (int c = 0; c < 3; ++c)
    value[c] = std::isfinite(value[c]) ? std::clamp(value[c], Float(0), Float(1)) : 0;
  return value;
}
Float OpenPBRMaterial::EmissionEstimate() const {
  const auto& p = impl->parameters;
  return p.emission_luminance * (p.emission_color.x + p.emission_color.y + p.emission_color.z) / 3;
}

OpenPBRInteraction OpenPBRMaterial::Prepare(const Ray& ray, const hit_record& h, Float exterior_ior) const {
  auto p = impl->parameters;
  p.base_color = vector(get_albedo(h));
  if (impl->has_roughness)
    p.specular_roughness = impl->roughness->raw_value(
      h.u * impl->texture_repeat[0],
      h.v * impl->texture_repeat[1])[0];
  const auto geometric = vector(h.geometric_normal.squared_length() > 0 ? h.geometric_normal : h.normal);
  auto n = vector(h.physical_shading_normal.squared_length() > 0 ? h.physical_shading_normal :
                  h.has_bump ? h.bump_normal : h.normal);
  if (glm::dot(impl->normal, impl->normal) > 0) n = impl->normal;
  if (glm::dot(n, geometric) < 0) n = -n;
  auto t = glm::dot(impl->tangent, impl->tangent) > 0 ? impl->tangent : vector(h.dpdu);
  p.geometry_basis = basis(n, t);
  auto cn = glm::dot(impl->coat_normal, impl->coat_normal) > 0 ? impl->coat_normal : n;
  if (glm::dot(cn, geometric) < 0) cn = -cn;
  auto ct = glm::dot(impl->coat_tangent, impl->coat_tangent) > 0 ? impl->coat_tangent : t;
  p.geometry_coat_basis = basis(cn, ct);
  auto result = std::make_unique<OpenPBRInteraction::Impl>();
  result->geometric_normal = glm::normalize(geometric);
  result->shading_normal = p.geometry_basis.n;
  result->view = -glm::normalize(vector(ray.d));
  // Fixed RGB wavelengths are the reference RGB mode, not spectral transport.
  // Unit throughput changes lobe selection efficiency, not the estimator.
  result->bsdf = ref::openpbr_prepare(p, glm::vec3(1), ref::OpenPBR_BaseRgbWavelengths_nm,
                                    exterior_ior, result->view);
  if (!p.geometry_thin_walled) {
    // The pinned reference's microfacet BTDF includes the refraction Jacobian
    // in both f*cos and PDF, but not the additional radiance-mode 1/eta^2.
    // Read its prepared RGB ratios (including weight/coat/dispersion changes)
    // rather than independently reconstructing and risking a different IOR.
    // This is the sole dependency on a non-public field of the pinned version.
    const auto eta = result->bsdf.fuzz_lobe.coating_lobe.base_lobe.specular_lobe.eta_t_over_eta_i;
    result->transmission_radiance_scale = glm::vec3(1) / (eta * eta);
    result->eta_squared = (eta.x * eta.x + eta.y * eta.y + eta.z * eta.z) / 3;
  }
  return OpenPBRInteraction(std::move(result));
}
point3f OpenPBRMaterial::emitted(const Ray& ray, const hit_record& h, Float, Float,
                              const point3f&, bool& invisible) {
  invisible = false;
  if (!(EmissionEstimate() > 0)) return point3f(0);
  // Resolve the exterior on either side, including a light seen from inside an
  // overlapping dielectric. Remove only the most recent placement on exit.
  const dielectric* exterior = nullptr;
  if (ray.pri_stack) {
    size_t remove = ray.pri_stack->size();
    if (is_dielectric() && dot(ray.d, h.geometric_normal) > 0)
      for (size_t i = remove; i > 0; --i)
        if ((*ray.pri_stack)[i - 1] == this) { remove = i - 1; break; }
    for (size_t i = 0; i < ray.pri_stack->size(); ++i) {
      auto* d = (*ray.pri_stack)[i];
      if (i != remove && (!exterior || d->priority <= exterior->priority)) exterior = d;
    }
  }
  return Prepare(ray, h, exterior ? exterior->ref_idx : 1).Emission();
}
bool OpenPBRMaterial::scatter(const Ray&, const hit_record&, scatter_record&, random_gen&) {
  throw std::runtime_error("OpenPBR requires the NEE integrator.");
}
bool OpenPBRMaterial::scatter(const Ray&, const hit_record&, scatter_record&, Sampler*) {
  throw std::runtime_error("OpenPBR requires the NEE integrator.");
}

OpenPBRInteraction::OpenPBRInteraction(std::unique_ptr<Impl> state) : impl(std::move(state)) {}
OpenPBRInteraction::~OpenPBRInteraction() = default;
OpenPBRInteraction::OpenPBRInteraction(OpenPBRInteraction&&) noexcept = default;
point3f OpenPBRInteraction::Evaluate(const vec3f& direction) const {
  if (!impl->Supports(vector(direction))) return point3f(0);
  auto value = ref::openpbr_get_sum_of_diffuse_specular(ref::openpbr_eval(impl->bsdf, vector(direction)));
  if (glm::dot(impl->view, impl->geometric_normal) * glm::dot(vector(direction), impl->geometric_normal) < 0)
    value *= impl->transmission_radiance_scale;
  return color(value);
}
Float OpenPBRInteraction::Pdf(const vec3f& direction) const {
  return impl->Supports(vector(direction)) ? ref::openpbr_pdf(impl->bsdf, vector(direction)) : 0;
}
point3f OpenPBRInteraction::Emission() const { return color(impl->bsdf.emission); }
OpenPBRSample OpenPBRInteraction::Sample(Float branch, Float u, Float v) const {
  glm::vec3 direction;
  ref::OpenPBR_DiffuseSpecular weight;
  float density;
  std::uint32_t type;
  ref::openpbr_sample(impl->bsdf, glm::vec3(branch, u, v), direction, weight, density, type);
  if (!(density > 0)) return {};
  OpenPBRSample result;
  result.direction = vec3f(direction.x, direction.y, direction.z);
  result.weight = color(ref::openpbr_get_sum_of_diffuse_specular(weight));
  result.pdf = density;
  result.specular = type & ref::OpenPBR_BsdfLobeTypeSpecular;
  result.transmission = type & ref::OpenPBR_BsdfLobeTypeTransmission;
  // A perturbed frame must not turn a reflection into a geometric transmission
  // (or vice versa). Reject without resampling so the PDF stays unchanged.
  bool crossed = glm::dot(impl->view, impl->geometric_normal) * glm::dot(direction, impl->geometric_normal) < 0;
  if (crossed != result.transmission || !impl->Supports(direction)) return {};
  if (crossed) {
    result.weight *= color(impl->transmission_radiance_scale);
    result.eta_squared = impl->eta_squared;
  }
  return result;
}

OpenPBRVolume OpenPBRInterior(const Rcpp::List& parameters) {
  auto p = inputs(parameters);
  ref::OpenPBR_PreparedBsdf prepared;
  ref::OpenPBR_VolumeDerivedProps derived;
  ref::openpbr_prepare_volume(p, derived, prepared, true);
  OpenPBRVolume result;
  result.absorption = color(prepared.volume.extinction_coefficient * (glm::vec3(1) - prepared.volume.albedo));
  result.scattering = color(prepared.volume.extinction_coefficient * prepared.volume.albedo);
  result.anisotropy = prepared.volume.anisotropy;
  result.ior = p.specular_ior;
  return result;
}
