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

Float SubsurfaceFlightLimit(const Ray &ray, double distance) {
  if (!std::isfinite(distance) || distance >= MaxT) return MaxT;
  // The physical distance stays double precision. Round the traversal limit
  // upward and inspect a small neighborhood beyond it so a near-endpoint hit
  // still supplies the normal/error bounds used by SubsurfaceCollisionPoint.
  // This is a query guard, never an extra flight distance or a resampled event.
  double magnitude = std::max({1.0, std::abs(double(ray.o[0])),
      std::abs(double(ray.o[1])), std::abs(double(ray.o[2])), distance});
  double guarded = distance + 32 * std::numeric_limits<Float>::epsilon() * magnitude;
  if (guarded >= MaxT) return MaxT;
  return std::nextafter(Float(guarded), std::numeric_limits<Float>::infinity());
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

// The radial density 2*pi*r*Sr/A is a 1:3 mixture of exponentials with
// scales d and 3d. Sampling the full tail avoids a radius-dependent cutoff bias.
double NormalizedDiffusionProfile::SampleRadius(int c, double component, double u) const {
  return -(component < .25 ? radius[c] : 3 * radius[c]) * std::log1p(-unit_uniform(u));
}
double NormalizedDiffusionProfile::LogAreaPdf(int c, double r) const {
  if (!(r > 0)) return INFINITY; // integrable 1/r singularity; callers handle its limit
  double x = r / radius[c];
  return -x / 3 + std::log1p(std::exp(-2 * x / 3)) -
         std::log(8 * M_PI) - std::log(radius[c]) - std::log(r);
}
double NormalizedDiffusionProfile::Cdf(int c, double r) const {
  return -.25 * std::expm1(-r / radius[c]) - .75 * std::expm1(-r / (3 * radius[c]));
}

std::optional<DiffusionSurfaceSample> SampleDiffusionSurface(
    const VolumeScene &scene, const MediumEntry &body, const point3f &entry,
    const normal3f &entry_normal, Float time, Sampler *sampler, random_gen &rng,
    const std::atomic<bool> *cancel) {
  const Medium &medium = *body.boundary->medium;
  NormalizedDiffusionProfile profile{medium.diffusion_color, medium.diffusion_radius};
  if (*std::max_element(profile.color.begin(), profile.color.end()) == 0) return {};
  onb frame;
  frame.build_from_w(entry_normal);
  const double axis_probability[] = {.25, .25, .5};
  double axis_uniform = sampler->Get1D();
  int axis = axis_uniform < .5 ? 2 : (axis_uniform < .75 ? 0 : 1);
  int channel = std::min(2, int(sampler->Get1D() * 3));
  double component = sampler->Get1D(), radial_uniform = sampler->Get1D();
  // A sampler may return exactly zero. Use an open variate for the singular
  // profile, without changing the density on any interval of positive measure.
  radial_uniform = std::max(radial_uniform, std::numeric_limits<double>::epsilon());
  double radius = profile.SampleRadius(channel, component, radial_uniform);
  double azimuth = 2 * M_PI * sampler->Get1D();
  vec3f direction = frame[axis];
  // The body bounds delimit the whole chord, including transformed instances.
  // Other winning dielectrics can cut this region but cannot enlarge it.
  aabb bounds;
  if (!body.boundary->bounding_box(time, time, bounds) || !scene.boundary_bvh) return {};
  Transform placement = body.medium_to_world * Inverse(body.boundary->medium_to_world);
  bounds = placement(bounds);
  double farthest = 0;
  for (int corner = 0; corner < 8; ++corner) {
    point3f p = bounds.Corner(corner);
    farthest = std::max(farthest, std::hypot(double(p[0]) - entry[0],
        double(p[1]) - entry[1], double(p[2]) - entry[2]));
  }
  // A projected offset beyond every corner cannot reach this object. Test in
  // double before converting to positions, so even enormous valid radii miss
  // cleanly without generating infinite ray coordinates. This is a bounds
  // rejection, not a truncation of the diffusion profile.
  if (radius > farthest * (1 + 8 * std::numeric_limits<Float>::epsilon())) return {};
  point3f center = entry + Float(radius * std::cos(azimuth)) * frame[(axis + 1) % 3] +
                          Float(radius * std::sin(azimuth)) * frame[(axis + 2) % 3];
  double lower = INFINITY, upper = -INFINITY, magnitude = 1;
  for (int corner = 0; corner < 8; ++corner) {
    point3f p = bounds.Corner(corner);
    double projected = 0;
    for (int a = 0; a < 3; ++a) {
      projected += (double(p[a]) - center[a]) * direction[a];
      magnitude = std::max(magnitude, std::abs(double(p[a])));
    }
    lower = std::min(lower, projected);
    upper = std::max(upper, projected);
  }
  // Roundoff padding ensures the first and last endpoints lie outside the body.
  double padding = 64 * std::numeric_limits<Float>::epsilon() * magnitude;
  Ray probe(center + Float(lower - padding) * direction, direction, time);
  probe.segment_absorption = true;
  Float end = Float(upper - lower + 2 * padding);
  if (!bounds.hit(probe, 0, end, rng)) return {};
  Float query_min = 0;
  std::vector<DiffusionSurfaceSample> candidates;
  while (!(cancel && cancel->load(std::memory_order_relaxed))) {
    hit_record h;
    if (!scene.boundary_bvh->hit(probe, query_min, end, h, rng)) break;
    // An invalid primitive hit must not turn this ordered walk into an
    // unbounded candidate list. Report the path failure instead of retrying it.
    if (!std::isfinite(h.OrderedDistance()) || h.OrderedDistance() < probe.medium_t_min)
      throw PathFailure(PathFailureKind::PositionPrecision,
                        "Non-finite or non-advancing diffusion boundary intersection.");
    // Classify both sides geometrically, rather than replaying just this hit.
    // A coincident glass/liquid face can have two crossings at the same t; a
    // closest-hit query reports only one. These offset memberships resolve the
    // complete priority interface and also handle contacts at mesh seams.
    Ray before_ray(OffsetMediumOrigin(h, -direction), -direction, time);
    Ray after_ray(OffsetMediumOrigin(h, direction), direction, time);
    auto previous = scene.InitialState(before_ray, cancel);
    auto state = scene.InitialState(after_ray, cancel);
    const auto *before_entry = previous.Active(), *after_entry = state.Active();
    bool before = before_entry && before_entry->boundary_id == body.boundary_id;
    bool after = after_entry && after_entry->boundary_id == body.boundary_id;
    if (before != after) {
      DiffusionSurfaceSample candidate;
      candidate.hit = h;
      candidate.outward = h.geometric_normal;
      if ((dot(direction, candidate.outward) > 0) != before)
        candidate.outward = -candidate.outward;
      candidate.outside = before ? state : previous;
      // Two faces separated by less than their position error can both report
      // this same effective transition after side classification. Count it once:
      // duplicated exits would multiply the BSSRDF's energy. Distinct resolved
      // intervals still retain every crossing, including opposite orientations.
      bool duplicate = false;
      if (!candidates.empty()) {
        const auto &last = candidates.back();
        double tolerance = h.pError.length() + last.hit.pError.length() +
                           8 * std::numeric_limits<Float>::epsilon() * magnitude;
        duplicate = ((dot(direction, last.outward) > 0) == before) &&
                    (h.p - last.hit.p).length() <= tolerance;
      }
      if (!duplicate) candidates.push_back(std::move(candidate));
    }
    // Keep the exact probe line: normal offsets can skip thin contact regions.
    probe.medium_t_min = std::nextafter(h.OrderedDistance(), INFINITY);
    query_min = Float(probe.medium_t_min);
    if (double(query_min) > probe.medium_t_min)
      query_min = std::nextafter(query_min, Float(-INFINITY));
  }
  if (candidates.empty() || (cancel && cancel->load(std::memory_order_relaxed))) return {};
  // Select crossings in proportion to their actual RGB reflectance profile.
  // Uniform selection wastes half the normal-axis samples on the opposite side
  // of a thick solid. The compensating weights then multiply on repeated
  // glass/liquid visits, producing severe bright outliers. Importance selection
  // retains every crossing and the full profile tail without that extra variance.
  std::vector<double> log_scores(candidates.size(), -INFINITY);
  double log_total = -INFINITY;
  for (size_t i = 0; i < candidates.size(); ++i) {
    const point3f &p = candidates[i].hit.p;
    double r = std::max(std::hypot(double(p[0]) - entry[0], double(p[1]) - entry[1],
                                   double(p[2]) - entry[2]), std::numeric_limits<double>::min());
    for (int c = 0; c < 3; ++c)
      if (profile.color[c] > 0)
        log_scores[i] = log_sum(log_scores[i], std::log(profile.color[c]) + profile.LogAreaPdf(c, r));
    log_total = log_sum(log_total, log_scores[i]);
  }
  if (!std::isfinite(log_total)) return {};
  double select = unit_uniform(sampler->Get1D()), cumulative = 0;
  // If the accumulated CDF rounds just below one, fall back to a positive mass.
  size_t chosen = size_t(std::max_element(log_scores.begin(), log_scores.end()) - log_scores.begin());
  for (size_t i = 0; i < candidates.size(); ++i) {
    cumulative += std::exp(log_scores[i] - log_total);
    if (select < cumulative) { chosen = i; break; }
  }
  double log_selection_pdf = log_scores[chosen] - log_total;
  auto result = std::move(candidates[chosen]);
  result.candidates = candidates.size();

  // Combine axis/channel proposals in area measure, and divide by the chosen
  // crossing's conditional probability. The MIS weights use projected densities
  // without crossing probabilities: they still sum to one at every surface point,
  // so this remains unbiased even though each axis has a different candidate set.
  double coordinates[3];
  for (int a = 0; a < 3; ++a) {
    coordinates[a] = 0;
    for (int b = 0; b < 3; ++b)
      coordinates[a] += (double(result.hit.p[b]) - entry[b]) * frame[a][b];
  }
  double distance = std::hypot(coordinates[0], coordinates[1], coordinates[2]);
  // At exact coincidence the limiting normal-axis density cancels the same
  // 1/r singularity in Sr. Such rounded samples use a common positive distance.
  double minimum_distance = std::numeric_limits<double>::min();
  distance = std::max(distance, minimum_distance);
  double log_pdf = -INFINITY;
  for (int a = 0; a < 3; ++a) {
    double projected = std::hypot(coordinates[(a + 1) % 3], coordinates[(a + 2) % 3]);
    projected = std::max(projected, minimum_distance);
    double cosine = std::abs(double(dot(result.outward, frame[a])));
    if (!(cosine > 0)) continue;
    for (int c = 0; c < 3; ++c)
      log_pdf = log_sum(log_pdf, profile.LogAreaPdf(c, projected) +
                       std::log(axis_probability[a] * cosine / 3));
  }
  if (!std::isfinite(log_pdf)) return {};
  for (int c = 0; c < 3; ++c)
    result.weight[c] = profile.color[c] *
        std::exp(profile.LogAreaPdf(c, distance) - log_pdf - log_selection_pdf);
  return result;
}

namespace {
// Evaluate in double precision, including index-matched interfaces and the
// critical angle. eta is interior/exterior; directions lie on the exterior side.
double diffusion_transmission(double cosine, double eta) {
  if (eta == 1) return 1;
  if (!(cosine > 0)) return 0;
  double sin2_transmitted = (1 - cosine * cosine) / (eta * eta);
  if (sin2_transmitted >= 1) return 0;
  double transmitted = std::sqrt(std::max(0.0, 1 - sin2_transmitted));
  double parallel = (eta * cosine - transmitted) / (eta * cosine + transmitted);
  double perpendicular = (cosine - eta * transmitted) / (cosine + eta * transmitted);
  return 1 - .5 * (parallel * parallel + perpendicular * perpendicular);
}
double diffusion_fresnel_normalization(double eta) {
  if (eta == 1) return 1;
  // Per-thread cache: material IOR pairs are constant throughout a render.
  // Integrate 2*mu*(1-F) dmu only above the critical angle. Substituting
  // mu^2 = critical^2 + (1-critical^2)*t^2 removes the critical-angle root.
  thread_local std::vector<std::pair<double, double>> cache;
  for (const auto &entry : cache) if (entry.first == eta) return entry.second;
  double critical2 = std::max(0.0, 1 - eta * eta), integral = 0;
  constexpr int steps = 256;
  for (int i = 0; i <= steps; ++i) {
    double t = double(i) / steps;
    double mu = std::sqrt(critical2 + (1 - critical2) * t * t);
    double f = 2 * (1 - critical2) * t * diffusion_transmission(mu, eta);
    integral += (i == 0 || i == steps ? 1 : (i % 2 ? 4 : 2)) * f;
  }
  integral /= 3 * steps;
  cache.emplace_back(eta, integral);
  return integral;
}
// Importance distribution for the external exit lobe. In t coordinates,
// mu^2 = critical^2 + (1-critical^2)*t^2, the critical-angle square root is
// smooth. Tabulate g(t) = 2*(1-critical^2)*t*T(mu) as a piecewise-linear density.
// Inverting its piecewise-quadratic CDF and evaluating that SAME density keeps
// the estimator unbiased for the unchanged analytic Fresnel lobe, even when
// the table approximates it. The table affects variance, not the material.
struct DiffusionExitDistribution {
  static constexpr int intervals = 1024;
  double eta, critical2, support, mass;
  std::array<double, intervals + 1> density{}, cdf{};

  explicit DiffusionExitDistribution(double ratio)
      : eta(ratio), critical2(std::max(0.0, 1 - ratio * ratio)),
        support(1 - critical2), mass(0) {
    for (int i = 0; i <= intervals; ++i) {
      double t = double(i) / intervals;
      double mu = std::sqrt(critical2 + support * t * t);
      density[i] = 2 * support * t * diffusion_transmission(mu, eta);
      if (i > 0) cdf[i] = cdf[i - 1] + (density[i - 1] + density[i]) / (2 * intervals);
    }
    mass = cdf.back();
  }

  double SampleCosine(double u) const {
    // Reverse the variate so the index-matched case agrees with the existing
    // cosine sampler. Avoid sampling exactly at the zero-density critical edge.
    double target = (1 - unit_uniform(u)) * mass;
    int bin = std::clamp(int(std::upper_bound(cdf.begin(), cdf.end(), target) - cdf.begin()) - 1,
                         0, intervals - 1);
    double area = std::max(0.0, (target - cdf[bin]) * intervals);
    double a = density[bin], delta = density[bin + 1] - a;
    // Rationalized quadratic inverse is stable for both flat and sloped bins.
    double denominator = a + std::sqrt(std::max(0.0, a * a + 2 * delta * area));
    double fraction = denominator > 0 ? 2 * area / denominator : 0;
    double t = (bin + std::clamp(fraction, 0.0, 1.0)) / intervals;
    double mu = std::sqrt(critical2 + support * t * t);
    // Float directions must remain on the transmitting side after rounding.
    Float first_transmitting = std::nextafter(Float(std::sqrt(critical2)), Float(1));
    return std::clamp(mu, double(first_transmitting), 1.0);
  }

  double Pdf(double mu) const {
    if (!(mu > 0) || mu * mu <= critical2) return 0;
    double t = std::sqrt(std::max(0.0, (mu * mu - critical2) / support));
    double coordinate = std::min(t, 1.0) * intervals;
    int bin = std::min(int(coordinate), intervals - 1);
    double fraction = coordinate - bin;
    double g = density[bin] + fraction * (density[bin + 1] - density[bin]);
    // p(omega) = p(t) * dt/dmu / (2*pi); azimuth is uniform.
    return g * mu / (2 * M_PI * support * t * mass);
  }
};

const DiffusionExitDistribution& diffusion_exit_distribution(double eta) {
  // Independent immutable tables per worker, keyed only by the IOR ratio.
  // Heap ownership keeps returned references stable as new ratios are cached.
  thread_local std::vector<std::unique_ptr<DiffusionExitDistribution>> cache;
  for (const auto& entry : cache) if (entry->eta == eta) return *entry;
  cache.emplace_back(std::make_unique<DiffusionExitDistribution>(eta));
  return *cache.back();
}

}
DiffusionExitBSDF::DiffusionExitBSDF(const normal3f &normal, double inside, double outside)
    : eta(inside / outside), normalization(diffusion_fresnel_normalization(eta)) {
  frame.build_from_w(normal);
}
DiffusionExitTable ExportDiffusionExitTable(double eta) {
  const auto &source = diffusion_exit_distribution(eta);
  DiffusionExitTable result;
  result.eta = eta;
  result.critical2 = source.critical2;
  result.support = source.support;
  result.mass = source.mass;
  result.normalization = diffusion_fresnel_normalization(eta);
  for (size_t i = 0; i < result.samples.size(); ++i)
    result.samples[i] = {float(source.density[i]), float(source.cdf[i])};
  return result;
}
double DiffusionExitBSDF::Pdf(const vec3f &wi) const {
  double cosine = std::clamp(double(dot(wi, frame.w())), 0.0, 1.0);
  return eta == 1 ? cosine / M_PI : diffusion_exit_distribution(eta).Pdf(cosine);
}
double DiffusionExitBSDF::Evaluate(const vec3f &wi) const {
  double cosine = std::clamp(double(dot(wi, frame.w())), 0.0, 1.0);
  return cosine / M_PI * diffusion_transmission(cosine, eta) / normalization * eta * eta;
}
vec3f DiffusionExitBSDF::Sample(double u, double v) const {
  double cosine = eta == 1 ? std::sqrt(1 - unit_uniform(u)) :
                            diffusion_exit_distribution(eta).SampleCosine(u);
  double radius = std::sqrt(std::max(0.0, 1 - cosine * cosine));
  double phi = 2 * M_PI * unit_uniform(v);
  return frame.local(Float(radius * std::cos(phi)), Float(radius * std::sin(phi)),
                     Float(cosine));
}
