#include "material.h"
#include <algorithm>

namespace {
using normalmap::Vector;
Vector vector(const vec3f& v) { return {v[0], v[1], v[2]}; }
vec3f direction(const Vector& v) { return {Float(v[0]), Float(v[1]), Float(v[2])}; }

// The primary facet carries the glossy mixture; the secondary facet is a
// perfect mirror. Both sampling and density include every reachable facet path.
class mapped_glossy_pdf final : public ScatteringDistribution {
public:
  template<class Next>
  mapped_glossy_pdf(normalmap::Model model, const Ray& ray, const hit_record& hit,
                   MicrofacetDistribution* distribution, const glossy& material, Next next)
      : scattering(std::move(model)),
        boundary(dot(hit.geometric_normal, ray.direction()) > 0
                     ? -hit.geometric_normal : hit.geometric_normal),
        distribution(distribution), alphas(distribution->Resolve(hit)), material(material),
        color(material.get_albedo(hit)) {
    scattering.Prepare([&](const Vector& view, const Vector& normal) {
      auto child = Child(view, normal);
      bool diffuse = true;
      const Float branch = next(), u = next(), v = next();
      const vec3f wi = child.Generate(u, v, branch, diffuse);
      const Float density = child.Value(wi);
      return normalmap::Scattering::ChildSample{vector(wi), density > 0
          ? ChildValue(view, vector(wi), normal) / density : point3f(0)};
    });
  }
  point3f Evaluate(const vec3f& d) const override {
    if (dot(boundary, d) <= 0) return point3f(0);
    return scattering.Evaluate(normalmap::normalized(vector(d)),
        [&](const Vector& view, const Vector& wi, const Vector& normal) {
          return ChildValue(view, wi, normal);
        });
  }
  Float value(const vec3f& d, random_gen&, Float) override { return Value(d); }
  Float value(const vec3f& d, Sampler*, Float) override { return Value(d); }
  vec3f generate(random_gen& rng, bool& diffuse, Float) override {
    const double choice = rng.unif_rand(), escape = rng.unif_rand(), branch = rng.unif_rand();
    const Float u = rng.unif_rand(), v = rng.unif_rand();
    return Sample(u, v, choice, escape, branch, diffuse);
  }
  vec3f generate(Sampler* sampler, bool& diffuse, Float) override {
    const double choice = sampler->Get1D(), escape = sampler->Get1D(), branch = sampler->Get1D();
    const vec2f uv = sampler->Get2D();
    return Sample(uv[0], uv[1], choice, escape, branch, diffuse);
  }
private:
  point3f ChildValue(const Vector& view, const Vector& wi, const Vector& normal) const {
    return material.Evaluate(view, wi, normal, color, alphas) *
           Float(std::max(0.0, dot(normal, wi)));
  }
  glossy_pdf Child(const Vector& view, const Vector& normal) const {
    glossy_pdf child(normal3f(normal[0], normal[1], normal[2]), -direction(view),
                     distribution, 0, 0);
    child.alphas = alphas;
    return child;
  }
  Float Value(const vec3f& d) const {
    if (dot(boundary, d) <= 0) return 0;
    return scattering.Pdf(normalmap::normalized(vector(d)), [&](const Vector& view, const Vector& wi, const Vector& n) {
      if (dot(n, view) <= 0 || dot(n, wi) <= 0) return Float(0);
      return Child(view, n).Value(direction(wi));
    });
  }
  vec3f Sample(Float u, Float v, double choice, double escape, double branch, bool& diffuse) const {
    const auto sampled = scattering.Sample(choice, escape,
        [&](const Vector& view, const Vector& normal) {
          if (dot(normal, view) <= 0) return Vector(0);
          return vector(Child(view, normal).Generate(u, v, Float(branch), diffuse));
        });
    const vec3f result = direction(sampled);
    // Geometric clipping adds null events, never retries or renormalizes.
    return dot(boundary, result) > 0 ? result : vec3f(0);
  }
  normalmap::Scattering scattering;
  const normal3f boundary;
  MicrofacetDistribution* distribution;
  const point2f alphas;
  const glossy& material;
  const point3f color;
};
}

bool glossy::scatter(const Ray& ray, const hit_record& hit, scatter_record& s, random_gen& rng) {
  auto model = NormalMapInteraction(ray, hit);
  if (!model.is_valid()) return false;
  s.is_specular = false;
  s.attenuation = get_albedo(hit);
  auto closure = new mapped_glossy_pdf(std::move(model), ray, hit, distribution, *this,
                                     [&] { return rng.unif_rand(); });
  s.pdf_ptr = closure;
  s.distribution = closure;
  return true;
}
bool glossy::scatter(const Ray& ray, const hit_record& hit, scatter_record& s, Sampler* sampler) {
  auto model = NormalMapInteraction(ray, hit);
  if (!model.is_valid()) return false;
  s.is_specular = false;
  s.attenuation = get_albedo(hit);
  auto closure = new mapped_glossy_pdf(std::move(model), ray, hit, distribution, *this,
                                     [&] { return sampler->Get1D(); });
  s.pdf_ptr = closure;
  s.distribution = closure;
  return true;
}

point3f glossy::Evaluate(const Vector& view, const Vector& incoming, const Vector& normal,
                         point3f color, point2f alphas) const {
  onb frame;
  frame.build_from_w_normalized(normal3f(normal[0], normal[1], normal[2]));
  const vec3f wo = unit_vector(frame.world_to_local(direction(view)));
  const vec3f wi = unit_vector(frame.world_to_local(direction(incoming)));
  // A child only reflects on the front of its primary facet. The normal-map
  // model, rather than an independent normal flip, handles back-facing bumps.
  if (!(wo[2] > 0 && wi[2] > 0)) return point3f(0);
  const vec3f half = unit_vector(wo + wi);
  const Float cosine = dot(wo, half);
  if (!(cosine > 0)) return point3f(0);
  auto pow5 = [](Float v) { return v*v*v*v*v; };
  const point3f diffuse = Float(28 / (23*M_PI)) * Rd * color * (point3f(1)-Rs) *
      (1-pow5(1-.5f*wi[2])) * (1-pow5(1-.5f*wo[2]));
  const point3f specular = distribution->D(half, alphas) /
      (4*cosine*std::max(wo[2], wi[2])) * SchlickFresnel(dot(wi, half));
  // Raw Ashikhmin--Shirley/Fresnel-blend BRDF. The adapter applies the surface
  // cosine exactly once, after the normal-map projected-area transport.
  return diffuse + specular;
}

point3f glossy::f(const Ray& ray, const hit_record& hit, const vec3f& wi) const {
  if (dot(hit.geometric_normal, ray.direction()) * dot(hit.geometric_normal, wi) >= 0)
    return point3f(0);
  const auto model = NormalMapInteraction(ray, hit);
  const auto incoming = normalmap::normalized(vector(wi));
  const point3f color = get_albedo(hit);
  const point2f alphas = distribution->Resolve(hit);
  auto value = model.evaluate(incoming,
      [&](const Vector& view, const Vector& light, const Vector& n) {
        return Evaluate(view, light, n, color, alphas);
      });
  value *= std::max(0.0, dot(model.geometric(), incoming));
  return value;
}

point3f glossy::SchlickFresnel(Float cosine) const {
  const Float x = 1-cosine;
  return Rs + x*x*x*x*x * (point3f(1)-Rs);
}
point3f glossy::get_albedo(const hit_record& hit) const {
  point3f color = albedo->value(hit);
  for (int c = 0; c < 3; ++c)
    color[c] = std::isfinite(color[c]) ? std::clamp(color[c], Float(0), Float(1)) : 0;
  return color;
}
size_t glossy::GetSize() { return sizeof(*this); }
