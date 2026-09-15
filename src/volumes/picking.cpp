#include "picking.h"
#include "../materials/material.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {
using Depth = std::array<double, 3>;

double Opacity(const Depth &depth) {
  // expm1 also resolves optically thin segments without subtracting from one.
  return -(std::expm1(-depth[0]) + std::expm1(-depth[1]) + std::expm1(-depth[2])) / 3;
}

point3f At(const Ray &ray, double t) {
  return point3f(double(ray.o[0]) + t * ray.d[0],
                 double(ray.o[1]) + t * ray.d[1],
                 double(ray.o[2]) + t * ray.d[2]);
}

double NextKnot(const DensityIndexRay &ray, double t, double end) {
  for (int a = 0; a < 3; ++a) {
    double d = ray.direction[a];
    if (d == 0)
      continue;
    double x = std::fma(t, d, ray.origin[a]);
    double plane = d > 0 ? std::floor(x) + 1 : std::ceil(x) - 1;
    double next = (plane - ray.origin[a]) / d;
    if (next <= t) // Roundoff at a knot must not stall traversal.
      next = (plane + (d > 0 ? 1 : -1) - ray.origin[a]) / d;
    if (next > t)
      end = std::min(end, next);
  }
  return end;
}

// Within one interpolation cell, trilinear density along a ray is cubic.
// Two-point Gauss quadrature integrates it exactly, up to field lookup roundoff.
double DensityIntegral(const Medium &medium, const Ray &ray, double a, double b) {
  double half = (b - a) / 2, mid = a + half;
  constexpr double node = 0.57735026918962576451;
  return half * (double(medium.Density(At(ray, mid - node * half))) +
                 double(medium.Density(At(ray, mid + node * half))));
}

// Accumulate a single homogeneous/interpolation interval, bisecting only the
// interval containing the requested opacity. RGB extinction is kept throughout.
std::optional<double> IntegrateInterval(const Medium *medium, const Ray &local,
                                        const point3f &glass, double a, double b,
                                        Depth &depth, double opacity) {
  if (!(b > a))
    return {};
  point3f sigma = medium ? medium->sigma_a + medium->sigma_s : point3f(0);
  auto integrated = [&](double end) {
    double density = medium ? (medium->IsHomogeneous() ? end - a
                                                      : DensityIntegral(*medium, local, a, end))
                            : 0;
    Depth result = depth;
    for (int c = 0; c < 3; ++c)
      result[c] += double(glass[c]) * (end - a) + double(sigma[c]) * density;
    return result;
  };
  Depth total = integrated(b);
  if (Opacity(total) < opacity) {
    depth = total;
    return {};
  }
  double lo = a, hi = b;
  for (int i = 0; i < 40; ++i) {
    double mid = lo + (hi - lo) / 2;
    if (Opacity(integrated(mid)) >= opacity)
      hi = mid;
    else
      lo = mid;
  }
  return hi;
}

std::optional<double> IntegrateSegment(const Ray &ray, double distance,
                                       const VolumePathState &state, Depth &depth,
                                       double opacity, const std::function<bool()> &cancel) {
  const dielectric *active_glass = nullptr;
  for (auto *glass : state.glass)
    if (!active_glass || glass->priority < active_glass->priority)
      active_glass = glass;
  point3f glass = active_glass ? active_glass->attenuation : point3f(0);
  const auto *entry = state.Active();
  if (!entry)
    return IntegrateInterval(nullptr, ray, glass, 0, distance, depth, opacity);
  const Medium &medium = *entry->boundary->medium;
  Transform inverse = Inverse(entry->medium_to_world);
  Ray local(inverse(ray.o), inverse(ray.d), ray.time());
  if (medium.IsHomogeneous())
    return IntegrateInterval(&medium, local, glass, 0, distance, depth, opacity);
  DensityIndexRay index_ray = medium.DensityRay(local);
  auto iterator = medium.SampleRay(local, distance);
  double previous = 0;
  while (auto segment = iterator.Next()) {
    if (cancel && cancel())
      return {};
    if (auto hit = IntegrateInterval(nullptr, local, glass, previous, segment->t_min,
                                     depth, opacity))
      return hit;
    bool empty = segment->sigma_maj[0] == 0 && segment->sigma_maj[1] == 0 &&
                 segment->sigma_maj[2] == 0;
    double t = segment->t_min;
    while (t < segment->t_max) {
      if (cancel && cancel())
        return {};
      double end = empty ? segment->t_max : NextKnot(index_ray, t, segment->t_max);
      if (auto hit = IntegrateInterval(empty ? nullptr : &medium, local, glass, t, end,
                                       depth, opacity))
        return hit;
      t = end;
    }
    previous = segment->t_max;
  }
  return IntegrateInterval(nullptr, local, glass, previous, distance, depth, opacity);
}
} // namespace

std::optional<RayPick> PickRay(const Ray &input, hitable *world, const VolumeScene *scene,
                              double opacity, const std::function<bool()> &cancel) {
  if (!std::isfinite(opacity) || opacity <= 0 || opacity >= 1)
    throw std::invalid_argument("Picking opacity must be strictly between zero and one.");
  if (!world || !std::isfinite(input.d.length()) || input.d.length() == 0 ||
      (cancel && cancel()))
    return {};
  // Reconstruct the ray so its cached intersection data matches the unit
  // direction. All following integration intervals are world-space distances.
  Ray ray(input.o, unit_vector(input.d), input.time());
  VolumePathState state = scene ? scene->InitialState(ray, nullptr, cancel) : VolumePathState{};
  random_gen rng(0); // Only geometry/alpha intersections use this local RNG.
  Depth depth{};
  while (!(cancel && cancel())) {
    if (scene)
      state.SetRay(ray);
    hit_record h;
    bool hit = world->hit(ray, 0, MaxT, h, rng);
    double distance = hit ? h.t : MaxT;
    if (scene) {
      auto t = IntegrateSegment(ray, distance, state, depth, opacity, cancel);
      if (cancel && cancel())
        return {};
      if (t)
        return RayPick{At(ray, *t), true};
    }
    if (!hit)
      return {};
    if (h.infinite_area_hit || (h.shape && h.shape->GetName() == "EnvironmentLight"))
      return RayPick{h.p, false, true};
    if (h.medium_boundary && !h.medium_boundary->keep_surface) {
      state.Cross(h, ray.d);
      ray = Ray(OffsetMediumOrigin(h, ray.d), ray.d, ray.time());
      continue;
    }
    if (h.alpha_miss) {
      ray = Ray(OffsetMediumOrigin(h, ray.d), ray.d, ray.time());
      continue;
    }
    return RayPick{h.p, false};
  }
  return {};
}
