#ifndef RAYRENDER_NORMALMAP_H
#define RAYRENDER_NORMALMAP_H

#include "../math/vec3.h"
#include "../math/point3.h"
#include <algorithm>
#include <cmath>
#include <array>

class WavefrontSceneCompiler;
namespace normalmap {
using Vector = vec3<double>;
constexpr double inv_pi = 0.31830988618379067154;

// Pure mathematics: unit directions point away from the surface. Double
// intermediates protect projected-area ratios near grazing; no mutable state.
bool valid(const Vector& v);
Vector normalized(const Vector& v);

// EON raw child BRDF, including single-scattering reflectance rho.
// sigma is in radians; roughness = clamp(sigma / (pi/2), 0, 1).
// No cosine or normal-map terminator correction belongs in the child.
struct DiffuseChild {
  friend class ::WavefrontSceneCompiler;
  explicit DiffuseChild(double sigma = 0);
  double eval(const Vector& wo, const Vector& wi, const Vector& normal,
              double rho = 1) const;
  double directional_albedo(double cosine, double rho = 1) const;
private:
  double missing_energy(double cosine) const;
  double a, b, average_loss;
};

class Model {
  friend class Scattering;
public:
  Model(Vector geometric, Vector shading, Vector outgoing);
  double eval_raw(const Vector& wi, const DiffuseChild& child, double rho = 1) const;
  // General reflective child: evaluate a raw BRDF in the primary facet's
  // frame. RGB and scalar children share the same projected-area transport.
  template<class Evaluate>
  auto evaluate(const Vector& direction, const Evaluate& child) const
      -> decltype(child(Vector(0), Vector(0), Vector(0))) {
    using Value = decltype(child(wo, wo, p));
    if (!usable || !valid(direction)) return Value(0);
    const Vector wi = normalized(direction);
    const double cosine = dot(g, wi);
    if (!(cosine > 0)) return Value(0);
    if (identity) return child(wo, wi, g);
    const double cp = std::max(0.0, dot(p, wi));
    auto value = child(wo, wi, p);
    value *= lambda * cp;
    if (dot(t, wi) > 0) {
      const Vector reflected = mirror(wi);
      auto term = child(wo, reflected, p);
      term *= lambda * std::max(0.0, dot(p, reflected)) * (1-masking(reflected));
      value += term;
    }
    if (dot(t, wo) > 0) {
      auto term = child(mirror(wo), wi, p);
      term *= (1-lambda) * cp;
      value += term;
    }
    const double area = cp + s * std::max(0.0, dot(t, wi));
    value *= std::min(1/cosine, c/area);
    return value;
  }
  double pdf(const Vector& wi) const;
  template<class Density>
  double pdf(const Vector& direction, const Density& child) const {
    if (!usable || !valid(direction)) return 0;
    const Vector wi = normalized(direction);
    if (!(dot(g, wi) > 0)) return 0;
    if (identity) return child(wo, wi, g);
    double density = lambda * child(wo, wi, p) * masking(wi);
    if (dot(t, wi) > 0) {
      const Vector reflected = mirror(wi);
      density += lambda * child(wo, reflected, p) * (1-masking(reflected));
    }
    if (dot(t, wo) > 0) density += (1-lambda) * child(mirror(wo), wi, p);
    return density;
  }
  // Remap the facet coordinate to a uniform child-lobe coordinate. The child
  // callback captures its independent directional coordinates. Mirror events
  // have unit solid-angle Jacobian; rejected paths remain explicit null events.
  template<class Sample>
  Vector sample(double choice, double escape, const Sample& child) const {
    if (!usable) return Vector(0);
    const bool primary = identity || choice < lambda;
    if (!identity) choice = primary ? choice/lambda : (choice-lambda)/(1-lambda);
    Vector wi = child(primary ? wo : mirror(wo), p, choice);
    if (!valid(wi)) return Vector(0);
    if (!identity && primary && escape >= masking(wi)) wi = mirror(wi);
    return dot(g, wi) > 0 ? normalized(wi) : Vector(0);
  }
  // Four independent uniforms: cosine disk, first facet, escape decision.
  // A zero vector is the explicit null event, never a resampling request.
  Vector sample(double u, double v, double facet, double escape) const;
  bool is_valid() const { return usable; }
  bool is_identity() const { return identity; }
  const Vector& geometric() const { return g; }
  double masking(const Vector& w) const;
private:
  Vector mirror(const Vector& w) const;
  Vector g{0}, p{0}, t{0}, wo{0}, tangent{0}, bitangent{0};
  double c = 1, s = 0, lambda = 1;
  bool usable = false, identity = false;
};

// A positive stochastic expansion of the two-facet transport through 16 hits.
// Each stored lobe connects an internal primary-facet vertex to the exterior,
// either directly or via one mirror hit. Forced internal continuations carry
// their non-escape probabilities instead of discarding the whole prefix.
// Sampling a prefix once per scene vertex makes Evaluate and Pdf deterministic
// conditional on it: ordinary light/BSDF MIS then integrates that same closure.
class Scattering {
public:
  static constexpr int order = 16;
  explicit Scattering(Model model) : model(std::move(model)) {}
  struct ChildSample { Vector direction; point3f weight; };

  template<class Sample>
  void Prepare(const Sample& sample) {
    count = 0;
    total = 0;
    if (!model.usable) return;
    for (int branch = 0; branch < 2; ++branch) {
      if (branch && model.identity) break;
      Vector view = branch ? model.mirror(model.wo) : model.wo;
      point3f weight(branch ? 1-model.lambda : model.lambda);
      int remaining = model.identity ? 1 : order-branch;
      while (remaining > 0) {
        // Channel maximum gives every nonzero RGB contribution proposal support.
        // This mass only selects a lobe; it is not the PDF of the internal walk.
        const double mass = std::max({double(weight[0]), double(weight[1]), double(weight[2])});
        if (!(mass > 0) || !std::isfinite(mass)) break;
        lobes[count++] = {view, weight, mass, remaining > 1};
        total += mass;
        if (remaining <= 2) break;
        const auto child = sample(view, model.p);
        if (!valid(child.direction)) break;
        const Vector reflected = model.mirror(child.direction);
        // A failed escape from both primary and mirror reaches the next primary.
        weight *= child.weight * Float((1-model.masking(child.direction)) *
                                        (1-model.masking(reflected)));
        view = -reflected;
        remaining -= 2;
      }
    }
  }

  template<class Evaluator>
  point3f Evaluate(const Vector& wi, const Evaluator& child) const {
    if (!valid(wi) || dot(model.g, wi) <= 0) return point3f(0);
    const Vector reflected = model.mirror(wi);
    point3f value(0);
    for (int i=0; i<count; ++i) {
      const auto& lobe = lobes[i];
      point3f term = child(lobe.view, wi, model.p);
      if (lobe.mirror_exit && dot(model.t, wi) > 0)
        term += child(lobe.view, reflected, model.p) * Float(1-model.masking(reflected));
      value += lobe.weight * term;
    }
    return value * Float(model.masking(wi));
  }

  template<class Density>
  double Pdf(const Vector& wi, const Density& child) const {
    if (!(total > 0) || !valid(wi) || dot(model.g, wi) <= 0) return 0;
    const Vector reflected = model.mirror(wi);
    double density = 0;
    for (int i=0; i<count; ++i) {
      const auto& lobe = lobes[i];
      double term = child(lobe.view, wi, model.p);
      if (lobe.mirror_exit) {
        term *= model.masking(wi);
        if (dot(model.t, wi) > 0)
          term += child(lobe.view, reflected, model.p) * (1-model.masking(reflected));
      }
      density += lobe.mass * term;
    }
    return density / total;
  }

  template<class Draw>
  Vector Sample(double component, double escape, const Draw& child) const {
    if (!(total > 0)) return Vector(0);
    double target = component * total;
    int selected = 0;
    while (selected+1 < count && target >= lobes[selected].mass)
      target -= lobes[selected++].mass;
    const auto& lobe = lobes[selected];
    Vector wi = child(lobe.view, model.p);
    if (!valid(wi)) return Vector(0);
    if (lobe.mirror_exit && escape >= model.masking(wi)) wi = model.mirror(wi);
    return dot(model.g, wi) > 0 ? normalized(wi) : Vector(0);
  }
private:
  struct Lobe { Vector view; point3f weight; double mass; bool mirror_exit; };
  Model model;
  // Eight primary vertices from each possible initial facet; order counts both.
  std::array<Lobe, order> lobes;
  int count = 0;
  double total = 0;
};

// Metric dual-tangent construction around a raw smooth normal. bu/bv use the
// existing texture's per-texel slopes; bv increases down the image, hence +bv.
Vector perturb(Vector normal, Vector dpdu, Vector dpdv, double bu, double bv);
}
#endif
