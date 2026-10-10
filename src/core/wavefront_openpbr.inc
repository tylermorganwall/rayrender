R"MSL(
  // The reference has already been preprocessed with its MSL backend selected.
  // Only these adapters know the packed layout and rayrender's radiance convention.
  struct PbrState {
    wf_openpbr::OpenPBR_PreparedBsdf prepared;
    float3 geometric, normal, view, transmission_scale;
    float eta_squared;
  };
  wf_openpbr::OpenPBR_Basis pbr_basis(float3 n, float3 tangent) {
    n = normalize(n);
    tangent -= n * dot(n, tangent);
    return dot(tangent, tangent) > 1e-12f
               ? wf_openpbr::openpbr_make_basis(n, normalize(tangent), 1.f)
               : wf_openpbr::openpbr_make_basis(n);
  }
  PbrState prepare_pbr(Material m, float3 rho, float3 view, float3 normal, float3 geometric,
                       float3 tangent, float exterior, float roughness, device const float4 *data) {
    auto p = wf_openpbr::openpbr_make_default_resolved_inputs();
    uint index = m.data;
)MSL"
#define WF_SCALAR(name) "p." #name " = data[index++].x;\n"
#define WF_COLOR(name) "p." #name " = data[index++].xyz;\n"
#include "wavefront_openpbr_fields.inc"
#undef WF_SCALAR
#undef WF_COLOR
    R"MSL(
  float3 override_n = data[index++].xyz, override_t = data[index++].xyz;
  float3 coat_n = data[index++].xyz, coat_t = data[index++].xyz;
  if (dot(override_n, override_n) > 0)
    normal = override_n;
  if (dot(normal, geometric) < 0)
    normal = -normal;
  if (dot(override_t, override_t) > 0)
    tangent = override_t;
  if (dot(coat_n, coat_n) == 0)
    coat_n = normal;
  if (dot(coat_n, geometric) < 0)
    coat_n = -coat_n;
  if (dot(coat_t, coat_t) == 0)
    coat_t = tangent;
  p.geometry_basis = pbr_basis(normal, tangent);
  p.geometry_coat_basis = pbr_basis(coat_n, coat_t);
  p.base_color = clamp(rho, 0.f, 1.f);
  if (roughness >= 0)
    p.specular_roughness = roughness;
  PbrState result;
  result.view = view;
  result.geometric = geometric;
  result.normal = p.geometry_basis.n;
  result.prepared = wf_openpbr::openpbr_prepare(
      p, float3(1), wf_openpbr::OpenPBR_BaseRgbWavelengths_nm, exterior, view);
  result.transmission_scale = 1;
  result.eta_squared = 1;
  if (!p.geometry_thin_walled) {
    float3 eta = result.prepared.fuzz_lobe.coating_lobe.base_lobe.specular_lobe.eta_t_over_eta_i;
    result.transmission_scale = 1 / (eta * eta);
    result.eta_squared = dot(eta * eta, float3(1.f / 3));
  }
  return result;
  }
  bool pbr_crossed(thread const PbrState &b, float3 wi) {
    return dot(b.view, b.geometric) * dot(wi, b.geometric) < 0;
  }
  bool pbr_supports(thread const PbrState &b, float3 wi) {
    return pbr_crossed(b, wi) == (dot(b.view, b.normal) * dot(wi, b.normal) < 0);
  }
  BsdfValue pbr_eval(thread const PbrState &b, float3 wi) {
    if (!pbr_supports(b, wi))
      return {0, 0};
    float3 value =
        wf_openpbr::openpbr_get_sum_of_diffuse_specular(wf_openpbr::openpbr_eval(b.prepared, wi));
    if (pbr_crossed(b, wi))
      value *= b.transmission_scale;
    return {value, wf_openpbr::openpbr_pdf(b.prepared, wi)};
  }
  BsdfSample pbr_sample(thread const PbrState &b, float3 u) {
    BsdfSample out = {0, 0, 0, 1, false};
    wf_openpbr::OpenPBR_DiffuseSpecular weight;
    uint type;
    wf_openpbr::openpbr_sample(b.prepared, u, out.direction, weight, out.pdf, type);
    bool transmission = type & wf_openpbr::OpenPBR_BsdfLobeTypeTransmission;
    if (!(out.pdf > 0) || transmission != pbr_crossed(b, out.direction) ||
        !pbr_supports(b, out.direction)) {
      out.pdf = 0;
      return out;
    }
    out.weight = wf_openpbr::openpbr_get_sum_of_diffuse_specular(weight);
    out.delta = type & wf_openpbr::OpenPBR_BsdfLobeTypeSpecular;
    if (transmission) {
      out.weight *= b.transmission_scale;
      out.eta_squared = b.eta_squared;
    }
    return out;
  }
    )MSL"
