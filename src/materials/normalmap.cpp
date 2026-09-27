#include "normalmap.h"

namespace normalmap {
bool valid(const Vector& v) {
  return std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]) &&
         v.squared_length() > 0;
}
Vector normalized(const Vector& v) {
  return valid(v) ? v / v.length() : Vector(0);
}

// Portsmouth, Kutz & Hill, JCGT 14(1), 2025 (revised 2026), Eqs.10–19.
// https://jcgt.org/published/0014/01/06/
// Exact analytic albedo, not the optional polynomial fit. Factor the missing
// energy as B*(C1-G/pi) to avoid cancellation as roughness approaches zero.
namespace {
constexpr double fon_c1 = .5 - 2*inv_pi/3;
constexpr double fon_c2 = 2./3 - 28*inv_pi/15;
}
DiffuseChild::DiffuseChild(double sigma) {
  const double roughness = std::clamp(sigma * (2*inv_pi), 0.0, 1.0);
  a = 1 / (1 + fon_c1*roughness);
  b = roughness*a;
  average_loss = b*(fon_c1-fon_c2);
}
double DiffuseChild::missing_energy(double cosine) const {
  const double mu = std::clamp(cosine, 0.0, 1.0);
  const double sine = std::sqrt(std::max(0.0,1-mu*mu));
  // (1-sin^3(theta))/cos(theta) = cos(theta)*(1+s+s^2)/(1+s).
  // This rationalization is finite at exact grazing, with no epsilon bias.
  const double g = sine*(std::acos(mu)-sine*mu) + (2./3)*
                   (sine*mu*(1+sine+sine*sine)/(1+sine)-sine);
  return b*std::max(0.0,fon_c1-inv_pi*g);
}
double DiffuseChild::directional_albedo(double cosine, double rho) const {
  const double loss = missing_energy(cosine);
  const double multiple = rho*rho*(1-average_loss)/(1-rho*average_loss);
  return rho*(1-loss) + multiple*loss;
}
double DiffuseChild::eval(const Vector& wo, const Vector& wi, const Vector& n,
                          double rho) const {
  const double co = dot(n, wo), ci = dot(n, wi);
  if (!(co > 0 && ci > 0)) return 0;
  if (b == 0) return rho*inv_pi;
  const double tangent_dot = dot(wo-co*n, wi-ci*n);
  const double angular = tangent_dot > 0 ? tangent_dot/std::max(ci,co) : tangent_dot;
  const double single = rho*(a+b*angular);
  // Eq.19 is nonlinear in reflectance: apply it independently per RGB channel.
  // average_loss is positive for b>0. Divide one loss first to remain stable
  // for arbitrarily small positive roughness, without a discontinuous cutoff.
  const double multiple = rho*rho*(1-average_loss)/(1-rho*average_loss) *
                          (missing_energy(co)/average_loss)*missing_energy(ci);
  return inv_pi*(single+multiple);
}

Model::Model(Vector geometric, Vector shading, Vector outgoing) {
  if (!valid(geometric) || !valid(shading) || !valid(outgoing)) return;
  g = normalized(geometric);
  p = normalized(shading);
  wo = normalized(outgoing);
  // Orient the entire fixed frame to the geometric side, never p alone to wo.
  if (dot(g, wo) < 0) { g = -g; p = -p; }
  if (!(dot(g, wo) > 0)) return;
  c = std::min(1.0, dot(g,p));
  Vector tilt = p - c*g;
  s = tilt.length();
  // Author default: invalid tilt uses the geometric child. The near-identity
  // threshold corresponds to a vector error <=1e-7, covered by continuity tests.
  identity = c <= 0 || s <= 1e-7;
  if (identity) { p = g; c = 1; s = 0; }
  else {
    t = -tilt / s;
    const double ap = std::max(0.0, dot(p,wo));
    const double at = s * std::max(0.0, dot(t,wo));
    lambda = ap / (ap + at);
  }
  Vector axis = std::abs(p[0]) > .9 ? Vector(0,1,0) : Vector(1,0,0);
  tangent = normalized(cross(axis,p));
  bitangent = cross(p,tangent);
  usable = true;
}
Vector Model::mirror(const Vector& w) const { return w - 2*dot(w,t)*t; }
double Model::masking(const Vector& w) const {
  if (identity) return dot(g,w) > 0 ? 1 : 0;
  const double area = std::max(0.0,dot(p,w)) + s*std::max(0.0,dot(t,w));
  return area > 0 ? std::min(1.0, std::max(0.0,dot(g,w))*c/area) : 0;
}
double Model::eval_raw(const Vector& direction, const DiffuseChild& child, double rho) const {
  if (!usable || !valid(direction)) return 0;
  const Vector wi = normalized(direction);
  const double cosine = dot(g,wi);
  if (!(cosine > 0)) return 0;
  if (identity) return child.eval(wo,wi,g,rho);
  // Eq.23: primary, primary→mirror, mirror→primary. Each child is raw;
  // the primary-facet cosine and macro masking are explicit here.
  const double cp = std::max(0.0,dot(p,wi));
  double value = lambda * child.eval(wo,wi,p,rho) * cp;
  if (dot(t,wi) > 0) {
    const Vector reflected = mirror(wi);
    value += lambda * child.eval(wo,reflected,p,rho) *
             std::max(0.0,dot(p,reflected)) * (1-masking(reflected));
  }
  if (dot(t,wo) > 0)
    value += (1-lambda) * child.eval(mirror(wo),wi,p,rho) * cp;
  // G(wi)/cosine analytically cancels the grazing cosine in Eq.13.
  // Evaluating this ratio directly avoids 0/0 and lost significant digits.
  const double area = cp + s*std::max(0.0,dot(t,wi));
  return value * std::min(1/cosine, c/area);
}
double Model::pdf(const Vector& direction) const {
  if (!usable || !valid(direction)) return 0;
  const Vector wi = normalized(direction);
  if (!(dot(g,wi) > 0)) return 0;
  const double cp = std::max(0.0,dot(p,wi))*inv_pi;
  if (identity) return cp;
  double density = lambda * cp * masking(wi);
  if (dot(t,wi) > 0) {
    const Vector reflected = mirror(wi);
    density += lambda * std::max(0.0,dot(p,reflected))*inv_pi *
               (1-masking(reflected));
  }
  if (dot(t,wo) > 0) density += (1-lambda)*cp;
  return density;
}
Vector Model::sample(double u, double v, double facet, double escape) const {
  if (!usable) return Vector(0);
  const double radius = std::sqrt(u), phi = 2 / inv_pi * v;
  Vector wi = radius*std::cos(phi)*tangent + radius*std::sin(phi)*bitangent +
              std::sqrt(std::max(0.0,1-u))*p;
  if (!identity && facet < lambda && escape >= masking(wi)) wi = mirror(wi);
  return dot(g,wi) > 0 ? normalized(wi) : Vector(0);
}
Vector perturb(Vector n, Vector u, Vector v, double bu, double bv) {
  n = normalized(n);
  u -= dot(u,n)*n;
  v -= dot(v,n)*n;
  const double uu = dot(u,u), uv = dot(u,v), vv = dot(v,v);
  const double determinant = uu*vv - uv*uv;
  if (!(determinant > 1e-12*uu*vv) || !std::isfinite(determinant)) return n;
  return normalized(n - bu*(vv*u - uv*v)/determinant +
                         bv*(uu*v - uv*u)/determinant);
}
}
