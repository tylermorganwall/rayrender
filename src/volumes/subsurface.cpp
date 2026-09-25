#include "subsurface.h"
#include "../materials/material.h"
#include "boundary.h"
#include <limits>

namespace {
double unit_uniform(double u) {
  return std::clamp(u, 0.0, std::nextafter(1.0, 0.0));
}
double log_sum(double a, double b) {
  if (a == -INFINITY) return b;
  if (b == -INFINITY) return a;
  double m = std::max(a, b);
  return m + std::log1p(std::exp(std::min(a, b) - m));
}
double log_exponential(double rate, double t, bool collision) {
  if (rate == 0) return collision ? -INFINITY : 0;
  return (collision ? std::log(rate) : 0) - rate * t;
}
}

point3f SubsurfaceCollisionPoint(const point3f &p, const Ray &ray, double t,
                                const hit_record &endpoint, bool &adjusted) {
  double distance2 = 0, error2 = 0;
  for (int a = 0; a < 3; ++a) {
    double delta = double(p[a]) - endpoint.p[a];
    double error = endpoint.pError[a] + gamma(3) *
        (std::abs(double(ray.o[a])) + std::abs(double(ray.d[a]) * t));
    distance2 += delta * delta;
    error2 += error * error;
  }
  adjusted = Float(t) >= endpoint.t || distance2 <= error2;
  if (adjusted) return OffsetMediumOrigin(endpoint, -ray.d);

  // A grazing collision can round onto the boundary plane while remaining far
  // from the ray's endpoint along that plane. A Euclidean endpoint-distance
  // check misses this case, leaving a subsequent inward ray at t=0 with the
  // medium already active. Correct only plane crossings within position error,
  // retaining the sampled distance/density and the collision's tangent position.
  normal3f n = endpoint.geometric_normal.squared_length() > 0
                   ? endpoint.geometric_normal : endpoint.normal;
  double plane_distance = 0, plane_error = 0, direction_dot = 0, normal_squared = 0;
  hit_record local = endpoint;
  local.p = p;
  for (int a = 0; a < 3; ++a) {
    double error = endpoint.pError[a] + gamma(3) *
        (std::abs(double(ray.o[a])) + std::abs(double(ray.d[a]) * t));
    local.pError[a] = Float(error);
    plane_distance += (double(p[a]) - endpoint.p[a]) * n[a];
    plane_error += error * std::abs(double(n[a]));
    direction_dot += double(ray.d[a]) * n[a];
    normal_squared += double(n[a]) * n[a];
  }
  if (normal_squared > 0 && direction_dot != 0 &&
      plane_distance * direction_dot >= 0 && std::abs(plane_distance) <= plane_error) {
    for (int a = 0; a < 3; ++a)
      local.p[a] = Float(double(p[a]) - plane_distance * n[a] / normal_squared);
    adjusted = true;
    return OffsetMediumOrigin(local, -ray.d);
  }
  return p;
}

double SubsurfaceProposal::Pole(double albedo) {
  if (!(albedo > 0 && albedo < 1)) return 0;
  // Solve a*atanh(k)/k=1. A bounded approximate pole remains a valid guide;
  // limiting only the proposal avoids singular poles for strongly absorbing RGB.
  double lo = 0, hi = .95;
  for (int i = 0; i < 48; ++i) {
    double k = (lo + hi) / 2;
    if (albedo * std::atanh(k) < k) lo = k;
    else hi = k;
  }
  return (lo + hi) / 2;
}
SubsurfaceProposal SubsurfaceProposal::Ordinary(const Medium &m) {
  SubsurfaceProposal p;
  p.phase = HGPhaseFunction(m.g);
  for (int c = 0; c < 3; ++c) {
    p.extinction[c] = double(m.sigma_a[c]) + m.sigma_s[c];
    p.pole[c] = m.subsurface_pole[c];
  }
  return p;
}
double SubsurfaceProposal::GuidePdf(int c, const vec3f &wi) const {
  double k = pole[c];
  if (k < 1e-7) return 1 / (4 * M_PI);
  double mu = std::clamp(double(dot(axis, wi)), -1.0, 1.0);
  return k / (4 * M_PI * std::atanh(k) * (1 - k * mu));
}
double SubsurfaceProposal::DirectionPdf(int c, const vec3f &wi) const {
  double physical = phase.p(wo, wi);
  return guided ? .5 * (physical + GuidePdf(c, wi)) : physical;
}
vec3f SubsurfaceProposal::SampleDirection(int hero, double strategy, double u, double v) const {
  if (!guided || strategy < .5) return phase.Sample(wo, u, v).wi;
  double k = pole[hero], mu = 2 * u - 1;
  if (k >= 1e-7)
    mu = (-std::expm1(std::log1p(k) - 2 * u * std::atanh(k))) / k;
  mu = std::clamp(mu, -1.0, 1.0);
  double s = std::sqrt(std::max(0.0, 1 - mu * mu));
  onb basis;
  basis.build_from_w(axis);
  return unit_vector(basis.local(s * std::cos(2 * M_PI * v), s * std::sin(2 * M_PI * v), mu));
}
double SubsurfaceProposal::GuideRate(int c, const vec3f &wi) const {
  return extinction[c] * (1 - pole[c] * std::clamp(double(dot(axis, wi)), -1.0, 1.0));
}
double SubsurfaceProposal::LogDistancePdf(int c, const vec3f &wi, double t, bool collision) const {
  double ordinary = log_exponential(extinction[c], t, collision);
  if (!guided) return ordinary;
  double pg = .5 * GuidePdf(c, wi) / DirectionPdf(c, wi);
  return log_sum(std::log1p(-pg) + ordinary,
                 std::log(pg) + log_exponential(GuideRate(c, wi), t, collision));
}
double SubsurfaceProposal::SampleDistance(int hero, const vec3f &wi, double strategy, double u) const {
  double rate = extinction[hero];
  if (guided && strategy < .5 * GuidePdf(hero, wi) / DirectionPdf(hero, wi))
    rate = GuideRate(hero, wi);
  return rate == 0 ? INFINITY : -std::log1p(-unit_uniform(u)) / rate;
}

SubsurfaceBoundaryBSDF::SubsurfaceBoundaryBSDF(const vec3f &incoming,
    const normal3f &normal, double e, double roughness) : eta(e), alpha(roughness * roughness) {
  normal3f n = dot(incoming, normal) < 0 ? normal : -normal;
  frame.build_from_w_normalized(n);
  wo = -unit_vector(frame.world_to_local(incoming));
}
double SubsurfaceBoundaryBSDF::D(const vec3f &m) const {
  if (m[2] <= 0) return 0;
  double a2 = alpha * alpha;
  double q = double(m[0]) * m[0] + double(m[1]) * m[1] + a2 * double(m[2]) * m[2];
  return a2 / (M_PI * q * q);
}
double SubsurfaceBoundaryBSDF::G1(const vec3f &w) const {
  double z2 = double(w[2]) * w[2];
  if (z2 == 0) return 0;
  double tan2 = (double(w[0]) * w[0] + double(w[1]) * w[1]) / z2;
  return 2 / (1 + std::sqrt(1 + alpha * alpha * tan2));
}
double SubsurfaceBoundaryBSDF::EvaluateLocal(const vec3f &wi, bool density) const {
  if (IsSpecular() || wo[2] <= 0 || wi[2] == 0) return 0;
  bool reflection = wi[2] > 0;
  vec3f sum = reflection ? wo + wi : wo + Float(eta) * wi;
  if (sum.squared_length() == 0) return 0;
  vec3f m = unit_vector(sum);
  if (m[2] < 0) m = -m;
  double om = dot(wo, m), im = dot(wi, m);
  if (om <= 0 || (reflection ? im <= 0 : im >= 0)) return 0;
  double fresnel = FrDielectric(om, 1 / eta);
  if (reflection) {
    return density ? D(m) * m[2] * fresnel / (4 * om)
                   : fresnel * D(m) * G1(wo) * G1(wi) / (4 * wo[2]);
  }
  double denom = im + om / eta;
  if (denom == 0) return 0;
  return density ? (1 - fresnel) * D(m) * m[2] * std::abs(im) / (denom * denom)
                 : (1 - fresnel) * D(m) * G1(wo) * G1(wi) * std::abs(im * om) /
                       (wo[2] * denom * denom * eta * eta);
}
double SubsurfaceBoundaryBSDF::Evaluate(const vec3f &wi) const {
  return EvaluateLocal(unit_vector(frame.world_to_local(wi)), false);
}
double SubsurfaceBoundaryBSDF::Pdf(const vec3f &wi) const {
  return EvaluateLocal(unit_vector(frame.world_to_local(wi)), true);
}
SubsurfaceBoundarySample SubsurfaceBoundaryBSDF::Sample(double branch, double u, double v) const {
  SubsurfaceBoundarySample s;
  s.specular = IsSpecular();
  vec3f m(0, 0, 1);
  if (!s.specular) {
    double uu = unit_uniform(u);
    double tan2 = alpha * alpha * uu / (1 - uu);
    double z = 1 / std::sqrt(1 + tan2), xy = std::sqrt(std::max(0.0, 1 - z * z));
    m = vec3f(xy * std::cos(2 * M_PI * v), xy * std::sin(2 * M_PI * v), z);
  }
  if (dot(wo, m) <= 0) return s; // NDF's invisible facets are zero-weight samples.
  double fresnel = FrDielectric(dot(wo, m), 1 / eta);
  vec3f wi;
  if (branch < fresnel) {
    wi = Reflect(wo, convert_to_normal3(m));
    if (wi[2] <= 0) return s;
  } else {
    if (eta == 1) wi = -wo;
    else if (!refract(-wo, m, 1 / eta, wi)) return s;
    if (wi[2] >= 0) return s;
    s.transmission = true;
  }
  s.wi = unit_vector(frame.local_to_world(wi));
  if (s.specular) {
    s.pdf = s.transmission ? 1 - fresnel : fresnel;
    s.weight = s.transmission ? 1 / (eta * eta) : 1;
  } else {
    s.pdf = Pdf(s.wi);
    s.weight = s.pdf > 0 ? Evaluate(s.wi) / s.pdf : 0;
  }
  return s;
}
