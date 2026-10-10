#include "material.h"
#include <algorithm>

namespace {
normalmap::Vector vector(const vec3f& v) { return {v[0],v[1],v[2]}; }
normalmap::Vector vector(const normal3f& v) { return {v[0],v[1],v[2]}; }
normalmap::Model interaction(const Ray& ray, const hit_record& h) {
  const normal3f base = h.base_shading_normal.squared_length() > 0
      ? h.base_shading_normal : h.geometric_normal;
  const auto outgoing = -vector(ray.direction());
  // Both frames must see the same side of the surface. Reject a smooth-normal
  // silhouette mismatch rather than flipping the shading frame through the mesh.
  if (dot(vector(h.geometric_normal), outgoing) * dot(vector(base), outgoing) <= 0)
    return normalmap::Model(normalmap::Vector(0), normalmap::Vector(0), outgoing);
  normal3f shading = h.physical_shading_normal;
  if (!(shading.squared_length() > 0)) {
    shading = h.has_bump ? h.bump_normal : base;
  }
  // Conservation and the projected-area cosine are measured in the smooth
  // base frame. Only the actual bump tilt enters the microgeometry model;
  // using a triangle's flat normal here makes triangulation visible in lighting.
  return normalmap::Model(vector(base), vector(shading), outgoing);
}
class normalmap_pdf final : public pdf {
public:
  normalmap_pdf(normalmap::Model model, const Ray& ray, const hit_record& h)
      : model(std::move(model)), boundary(dot(h.geometric_normal, ray.direction()) > 0
                                         ? -h.geometric_normal : h.geometric_normal) {}
  Float value(const vec3f& d, random_gen&, Float) override {
    return dot(boundary, d) > 0 ? model.pdf(vector(d)) : 0;
  }
  Float value(const vec3f& d, Sampler*, Float) override {
    return dot(boundary, d) > 0 ? model.pdf(vector(d)) : 0;
  }
  vec3f generate(random_gen& rng, bool& diffuse, Float) override {
    const double u = rng.unif_rand(), v = rng.unif_rand();
    const double facet = rng.unif_rand(), escape = rng.unif_rand();
    diffuse = true;
    auto d = model.sample(u,v,facet,escape);
    const vec3f direction(d[0],d[1],d[2]);
    return dot(boundary, direction) > 0 ? direction : vec3f(0);
  }
  vec3f generate(Sampler* sampler, bool& diffuse, Float) override {
    auto uv = sampler->Get2D();
    const double facet = sampler->Get1D(), escape = sampler->Get1D();
    diffuse = true;
    auto d = model.sample(uv[0],uv[1],facet,escape);
    const vec3f direction(d[0],d[1],d[2]);
    return dot(boundary, direction) > 0 ? direction : vec3f(0);
  }
private:
  const normalmap::Model model;
  // Geometric clipping adds null events, without resampling or renormalizing
  // the PDF. This preserves the energy bound and prevents wrong-side reflection.
  const normal3f boundary;
};
}

point3f diffuse_material::get_albedo(const hit_record& h) const {
  point3f a = albedo->value(h);
  // Physical reflectance domain, including HDR images and procedural textures.
  // Invalid/nonfinite channels absorb. R validates constant input colors too.
  for (int i=0;i<3;++i) a[i] = std::isfinite(a[i]) ? std::clamp(a[i],Float(0),Float(1)) : 0;
  return a;
}
point3f diffuse_material::f(const Ray& ray, const hit_record& h, const vec3f& direction) const {
  if (dot(h.geometric_normal, ray.direction()) * dot(h.geometric_normal, direction) >= 0)
    return point3f(0);
  auto model = interaction(ray,h);
  auto wi = normalmap::normalized(vector(direction));
  const double cosine = std::max(0.0,dot(model.geometric(),wi));
  point3f result = get_albedo(h);
  for (int channel=0; channel<3; ++channel)
    result[channel] = model.eval_raw(wi,child,result[channel])*cosine;
  return result;
}
bool diffuse_material::scatter(const Ray& ray, const hit_record& h, scatter_record& s, random_gen&) {
  auto model = interaction(ray,h);
  if (!model.is_valid()) return false;
  s.is_specular = false;
  s.attenuation = get_albedo(h);
  s.pdf_ptr = new normalmap_pdf(std::move(model), ray, h);
  return true;
}
bool diffuse_material::scatter(const Ray& ray, const hit_record& h, scatter_record& s, Sampler*) {
  auto model = interaction(ray,h);
  if (!model.is_valid()) return false;
  s.is_specular = false;
  s.attenuation = get_albedo(h);
  s.pdf_ptr = new normalmap_pdf(std::move(model), ray, h);
  return true;
}

void SetPhysicalBump(hit_record &h, const material *mat, const bump_texture *bump, const Ray &ray) {
  h.ComputeDifferentials(ray);
  h.base_shading_normal = h.geometric_normal;
  h.physical_shading_normal =
      mat && mat->physical_normal_mapping() ? h.geometric_normal : normal3f(0);
  if (!bump)
    return;
  const TextureFootprint footprint{h.dudx, h.dvdx, h.dudy, h.dvdy, h.has_differentials};
  TextureEvalContext context;
  if (bump->height_texture) {
    context = TextureEvalContext::FromHit(h);
    context.face_index = -1; // Analytic primitives have no Ptex source face.
  }
  h.bump_normal = bump->perturb(h.u, h.v, h.p, unit_vector(h.geometric_normal), h.dpdu, h.dpdv,
                                h.dndu, h.dndv, footprint, bump->height_texture ? &context : nullptr);
  h.has_bump = true;
  if (mat && mat->physical_normal_mapping())
    h.physical_shading_normal = h.bump_normal;
}
