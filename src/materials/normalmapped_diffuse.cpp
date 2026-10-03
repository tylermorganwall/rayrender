#include "material.h"
#include <algorithm>

namespace {
normalmap::Vector vector(const vec3f& v) { return {v[0],v[1],v[2]}; }
normalmap::Vector vector(const normal3f& v) { return {v[0],v[1],v[2]}; }
normalmap::Model interaction(const Ray& ray, const hit_record& h) {
  normal3f shading = h.physical_shading_normal;
  if (!(shading.squared_length() > 0)) {
    // Analytic primitives have no interpolated/view-repaired smooth normals.
    shading = h.has_bump ? h.bump_normal : h.geometric_normal;
  }
  return normalmap::Model(vector(h.geometric_normal), vector(shading), -vector(ray.direction()));
}
class normalmap_pdf final : public pdf {
public:
  explicit normalmap_pdf(normalmap::Model model) : model(std::move(model)) {}
  Float value(const vec3f& d, random_gen&, Float) override { return model.pdf(vector(d)); }
  Float value(const vec3f& d, Sampler*, Float) override { return model.pdf(vector(d)); }
  vec3f generate(random_gen& rng, bool& diffuse, Float) override {
    const double u = rng.unif_rand(), v = rng.unif_rand();
    const double facet = rng.unif_rand(), escape = rng.unif_rand();
    diffuse = true;
    auto d = model.sample(u,v,facet,escape);
    return vec3f(d[0],d[1],d[2]);
  }
  vec3f generate(Sampler* sampler, bool& diffuse, Float) override {
    auto uv = sampler->Get2D();
    const double facet = sampler->Get1D(), escape = sampler->Get1D();
    diffuse = true;
    auto d = model.sample(uv[0],uv[1],facet,escape);
    return vec3f(d[0],d[1],d[2]);
  }
private:
  const normalmap::Model model;
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
  s.pdf_ptr = new normalmap_pdf(std::move(model));
  return true;
}
bool diffuse_material::scatter(const Ray& ray, const hit_record& h, scatter_record& s, Sampler*) {
  auto model = interaction(ray,h);
  if (!model.is_valid()) return false;
  s.is_specular = false;
  s.attenuation = get_albedo(h);
  s.pdf_ptr = new normalmap_pdf(std::move(model));
  return true;
}

void SetPhysicalBump(hit_record &h, const material *mat, const bump_texture *bump, const Ray &ray) {
  h.ComputeDifferentials(ray);
  h.physical_shading_normal =
      mat && mat->physical_normal_mapping() ? h.geometric_normal : normal3f(0);
  if (!bump)
    return;
  const TextureFootprint footprint{h.dudx, h.dvdx, h.dudy, h.dvdy, h.has_differentials};
  h.bump_normal = bump->perturb(h.u, h.v, h.p, unit_vector(h.geometric_normal), h.dpdu, h.dpdv,
                                h.dndu, h.dndv, footprint);
  h.has_bump = true;
  if (mat && mat->physical_normal_mapping())
    h.physical_shading_normal = h.bump_normal;
}
