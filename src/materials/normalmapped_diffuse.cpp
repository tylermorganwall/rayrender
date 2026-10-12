#include "material.h"
#include <algorithm>

namespace {
normalmap::Vector vector(const vec3f& v) { return {v[0],v[1],v[2]}; }
normalmap::Vector vector(const normal3f& v) { return {v[0],v[1],v[2]}; }
}
normalmap::Model NormalMapInteraction(const Ray& ray, const hit_record& h) {
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
namespace {
class normalmap_pdf final : public ScatteringDistribution {
public:
  template<class Next>
  normalmap_pdf(normalmap::Model model, const Ray& ray, const hit_record& h,
                normalmap::DiffuseChild child, point3f color, Next next)
      : scattering(std::move(model)), boundary(dot(h.geometric_normal, ray.direction()) > 0
                                         ? -h.geometric_normal : h.geometric_normal),
        child(child), color(color) {
    scattering.Prepare([&](const normalmap::Vector& view, const normalmap::Vector& normal) {
      const double u = next(), v = next();
      const auto wi = ChildSample(normal, u, v);
      point3f weight;
      for (int c=0; c<3; ++c) weight[c] = child.eval(view, wi, normal, color[c]) / normalmap::inv_pi;
      return normalmap::Scattering::ChildSample{wi, weight};
    });
  }
  point3f Evaluate(const vec3f& d) const override {
    if (dot(boundary, d) <= 0) return point3f(0);
    return scattering.Evaluate(normalmap::normalized(vector(d)),
        [&](const normalmap::Vector& view, const normalmap::Vector& wi, const normalmap::Vector& normal) {
          point3f value;
          for (int c=0; c<3; ++c)
            value[c] = child.eval(view, wi, normal, color[c]) * std::max(0.0, dot(normal, wi));
          return value;
        });
  }
  Float value(const vec3f& d, random_gen&, Float) override {
    return Value(d);
  }
  Float value(const vec3f& d, Sampler*, Float) override {
    return Value(d);
  }
  vec3f generate(random_gen& rng, bool& diffuse, Float) override {
    const double u = rng.unif_rand(), v = rng.unif_rand();
    const double facet = rng.unif_rand(), escape = rng.unif_rand();
    diffuse = true;
    auto d = scattering.Sample(facet, escape, [&](const normalmap::Vector&, const normalmap::Vector& n) {
      return ChildSample(n, u, v);
    });
    const vec3f direction(d[0],d[1],d[2]);
    return dot(boundary, direction) > 0 ? direction : vec3f(0);
  }
  vec3f generate(Sampler* sampler, bool& diffuse, Float) override {
    auto uv = sampler->Get2D();
    const double facet = sampler->Get1D(), escape = sampler->Get1D();
    diffuse = true;
    auto d = scattering.Sample(facet, escape, [&](const normalmap::Vector&, const normalmap::Vector& n) {
      return ChildSample(n, uv[0], uv[1]);
    });
    const vec3f direction(d[0],d[1],d[2]);
    return dot(boundary, direction) > 0 ? direction : vec3f(0);
  }
private:
  static normalmap::Vector ChildSample(const normalmap::Vector& n, double u, double v) {
    onb basis;
    basis.build_from_w_normalized(normal3f(n[0],n[1],n[2]));
    const double radius = std::sqrt(u), phi = 2*M_PI*v;
    return vector(basis.local(vec3f(radius*std::cos(phi), radius*std::sin(phi),
                                   std::sqrt(std::max(0.0,1-u)))));
  }
  Float Value(const vec3f& d) const {
    if (dot(boundary, d) <= 0) return 0;
    return scattering.Pdf(normalmap::normalized(vector(d)),
        [](const normalmap::Vector&, const normalmap::Vector& wi, const normalmap::Vector& n) {
          return std::max(0.0, dot(n, wi)) * normalmap::inv_pi;
        });
  }
  normalmap::Scattering scattering;
  // Geometric clipping adds null events, without resampling or renormalizing
  // the PDF. This preserves the energy bound and prevents wrong-side reflection.
  const normal3f boundary;
  const normalmap::DiffuseChild child;
  const point3f color;
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
  auto model = NormalMapInteraction(ray,h);
  auto wi = normalmap::normalized(vector(direction));
  const double cosine = std::max(0.0,dot(model.geometric(),wi));
  point3f result = get_albedo(h);
  for (int channel=0; channel<3; ++channel)
    result[channel] = model.eval_raw(wi,child,result[channel])*cosine;
  return result;
}
bool diffuse_material::scatter(const Ray& ray, const hit_record& h, scatter_record& s, random_gen& rng) {
  auto model = NormalMapInteraction(ray,h);
  if (!model.is_valid()) return false;
  s.is_specular = false;
  s.attenuation = get_albedo(h);
  auto closure = new normalmap_pdf(std::move(model), ray, h, child, s.attenuation,
                                 [&] { return rng.unif_rand(); });
  s.pdf_ptr = closure;
  s.distribution = closure;
  return true;
}
bool diffuse_material::scatter(const Ray& ray, const hit_record& h, scatter_record& s, Sampler* sampler) {
  auto model = NormalMapInteraction(ray,h);
  if (!model.is_valid()) return false;
  s.is_specular = false;
  s.attenuation = get_albedo(h);
  auto closure = new normalmap_pdf(std::move(model), ray, h, child, s.attenuation,
                                 [&] { return sampler->Get1D(); });
  s.pdf_ptr = closure;
  s.distribution = closure;
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
