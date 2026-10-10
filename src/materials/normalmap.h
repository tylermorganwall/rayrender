#ifndef RAYRENDER_NORMALMAP_H
#define RAYRENDER_NORMALMAP_H

#include "../math/vec3.h"
#include <algorithm>
#include <cmath>

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
public:
  Model(Vector geometric, Vector shading, Vector outgoing);
  double eval_raw(const Vector& wi, const DiffuseChild& child, double rho = 1) const;
  double pdf(const Vector& wi) const;
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

// Metric dual-tangent construction around a raw smooth normal. bu/bv use the
// existing texture's per-texel slopes; bv increases down the image, hence +bv.
Vector perturb(Vector normal, Vector dpdu, Vector dpdv, double bu, double bv);
}
#endif
