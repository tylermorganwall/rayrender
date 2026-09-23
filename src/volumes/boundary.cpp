#include "boundary.h"
#include "../core/bvh.h"
#include "../hitables/box.h"
#include "../hitables/ellipsoid.h"
#include "../hitables/mesh3d.h"
#include "../hitables/plymesh.h"
#include "../hitables/raymesh.h"
#include "../hitables/sphere.h"
#include "../hitables/trimesh.h"
#include "../materials/material.h"
#include <map>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace {
TriangleMesh *boundary_mesh(hitable *geometry) {
  if (auto *g = dynamic_cast<trimesh *>(geometry))
    return g->mesh.get();
  if (auto *g = dynamic_cast<plymesh *>(geometry))
    return g->mesh.get();
  if (auto *g = dynamic_cast<mesh3d *>(geometry))
    return g->mesh.get();
  if (auto *g = dynamic_cast<raymesh *>(geometry))
    return g->mesh.get();
  return nullptr;
}
double signed_mesh_volume(const TriangleMesh &mesh) {
  double sum = 0;
  const point3f origin = mesh.p[0];
  for (size_t f = 0; f < mesh.nTriangles; ++f) {
    double p[3][3];
    for (int v = 0; v < 3; ++v)
      for (int a = 0; a < 3; ++a)
        p[v][a] = double(mesh.p[mesh.vertexIndices[3 * f + v]][a]) - origin[a];
    sum += p[0][0] * (p[1][1] * p[2][2] - p[1][2] * p[2][1]) +
           p[0][1] * (p[1][2] * p[2][0] - p[1][0] * p[2][2]) +
           p[0][2] * (p[1][0] * p[2][1] - p[1][1] * p[2][0]);
  }
  return sum / 6;
}
struct MediumBoxInterval {
  point3f origin;
  vec3f direction;
  double t_near = -INFINITY, t_far = INFINITY;
  int near_axis = -1, far_axis = -1;
};
// A contact with no resolvable interior interval is not a boundary crossing.
bool box_interval(const box &geometry, const Ray &ray, MediumBoxInterval &interval) {
  interval.origin = (*geometry.WorldToObject)(ray.o);
  interval.direction = (*geometry.WorldToObject)(ray.d);
  for (int axis = 0; axis < 3; ++axis) {
    if (interval.direction[axis] == 0) {
      if (interval.origin[axis] <= geometry.pmin[axis] ||
          interval.origin[axis] >= geometry.pmax[axis])
        return false;
      continue;
    }
    double a = (double(geometry.pmin[axis]) - interval.origin[axis]) / interval.direction[axis];
    double b = (double(geometry.pmax[axis]) - interval.origin[axis]) / interval.direction[axis];
    if (a > b)
      std::swap(a, b);
    if (a > interval.t_near) {
      interval.t_near = a;
      interval.near_axis = axis;
    }
    if (b < interval.t_far) {
      interval.t_far = b;
      interval.far_axis = axis;
    }
  }
  double tolerance = 8 * std::numeric_limits<Float>::epsilon() *
                     std::max({1.0, std::abs(interval.t_near), std::abs(interval.t_far)});
  return interval.t_far - interval.t_near > tolerance;
}
bool has_interior_interval(const hitable *geometry, const Ray &ray) {
  const auto *b = dynamic_cast<const box *>(geometry);
  MediumBoxInterval interval;
  return !b || box_interval(*b, ray, interval);
}
// Independent Float rectangle tests can both reject an entry at a box edge,
// leaving only the exit hit. Invisible containers need a single consistent slab
// interval; their intersections do not require surface texture or bump evaluation.
bool hit_invisible_box(const box &geometry, const Ray &ray, Float lo, Float hi, hit_record &h) {
  MediumBoxInterval interval;
  if (!box_interval(geometry, ray, interval))
    return false;
  bool entering = interval.t_near >= lo;
  double t = entering ? interval.t_near : interval.t_far;
  int axis = entering ? interval.near_axis : interval.far_axis;
  if (axis < 0 || t < lo || t > hi)
    return false;
  const auto &o = interval.origin;
  const auto &d = interval.direction;
  normal3f normal(0);
  normal[axis] = (d[axis] > 0 ? 1 : -1) * (entering ? -1 : 1);
  point3f p;
  vec3f error;
  for (int a = 0; a < 3; ++a) {
    p[a] = Float(double(o[a]) + double(d[a]) * t);
    error[a] = Float(gamma(3) * (std::abs(double(o[a])) + std::abs(double(d[a]) * t)));
  }
  p[axis] = normal[axis] < 0 ? geometry.pmin[axis] : geometry.pmax[axis];
  h.t = Float(t);
  h.p = (*geometry.ObjectToWorld)(p, error, &h.pError);
  h.normal = h.geometric_normal = unit_vector((*geometry.ObjectToWorld)(normal));
  h.bump_normal = h.normal;
  h.dpdu = h.dpdv = vec3f(0);
  h.u = h.v = 0;
  h.mat_ptr = geometry.mat_ptr.get();
  h.shape = &geometry;
  h.alpha_miss = h.has_bump = h.infinite_area_hit = false;
  return true;
}
} // namespace

point3f OffsetMediumOrigin(const hit_record &h, const vec3f &direction) {
  normal3f n = h.geometric_normal.squared_length() > 0 ? h.geometric_normal : h.normal;
  point3f origin = OffsetRayOrigin(h.p, h.pError, n, direction);
  // At a box edge or corner, offset every incident face. Moving off only the
  // selected rectangle leaves a second face at t=0, causing a duplicate crossing.
  if (h.medium_boundary)
    if (auto *box_geometry = dynamic_cast<const box *>(h.medium_boundary->geometry.get())) {
      Transform object_to_world = h.MediumToWorld();
      if (h.medium_boundary->medium)
        object_to_world = object_to_world * Inverse(h.medium_boundary->medium->medium_to_object);
      Transform inverse = Inverse(object_to_world);
      point3f p = inverse(h.p), q = inverse(origin);
      vec3f error(0);
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
          error[a] += std::abs(inverse.GetMatrix().m[a][b]) * h.pError[b];
      bool entering = dot(direction, n) < 0, changed = false;
      for (int a = 0; a < 3; ++a) {
        Float margin = 8 * std::numeric_limits<Float>::epsilon() *
                           std::max({Float(1), std::abs(p[a]), std::abs(box_geometry->pmin[a]),
                                     std::abs(box_geometry->pmax[a])}) +
                       error[a];
        if (std::abs(p[a] - box_geometry->pmin[a]) <= margin) {
          q[a] = box_geometry->pmin[a] + (entering ? margin : -margin);
          changed = true;
        } else if (std::abs(p[a] - box_geometry->pmax[a]) <= margin) {
          q[a] = box_geometry->pmax[a] + (entering ? -margin : margin);
          changed = true;
        }
      }
      if (changed)
        origin = object_to_world(q);
    }
  for (int a = 0; a < 3; ++a)
    if (origin[a] == h.p[a] && n[a] != 0) {
      Float side = (dot(direction, n) < 0 ? -1 : 1) * n[a];
      // A face on a coordinate-zero plane has zero normal error. One subnormal
      // ULP is lost in subsequent intersection arithmetic, particularly at a
      // concave edge. Use the point's scale, without a fixed world-unit epsilon.
      Float scale = std::max({std::abs(h.p[0]), std::abs(h.p[1]), std::abs(h.p[2]),
                              h.pError.length()});
      origin[a] += side * (8 * std::numeric_limits<Float>::epsilon() * scale);
      origin[a] = std::nextafter(origin[a], side > 0 ? Float(INFINITY) : Float(-INFINITY));
    }
  return origin;
}

MediumBoundary::MediumBoundary(std::shared_ptr<hitable> geom, std::shared_ptr<const Medium> m,
                               const Transform &object_to_world, bool surface, uint64_t id)
    : geometry(std::move(geom)), medium(std::move(m)),
      medium_to_world(object_to_world * (medium ? medium->medium_to_object : Transform())),
      keep_surface(surface), boundary_id(id) {
  if (boundary_id == 0)
    throw std::invalid_argument("Medium boundaries require a nonzero scene-assigned ID.");
  ValidateMediumTransform(medium_to_world);
  if (auto *mesh = boundary_mesh(geometry.get()))
    orientation_sign = signed_mesh_volume(*mesh) < 0 ? -1 : 1;
  else if ((dynamic_cast<sphere *>(geometry.get()) || dynamic_cast<ellipsoid *>(geometry.get())) &&
           geometry->reverseOrientation)
    orientation_sign = -1;
}
void MediumBoundary::Annotate(hit_record &h) const {
  h.medium_boundary = this;
  h.medium_to_world = medium_to_world;
  h.boundary_id = boundary_id;
  if (h.geometric_normal.squared_length() == 0)
    h.geometric_normal = h.normal;
  h.geometric_normal = unit_vector(h.geometric_normal) * orientation_sign;
}
const bool MediumBoundary::hit(const Ray &r, Float lo, Float hi, hit_record &h,
                               random_gen &rng) const {
  if (!r.segment_absorption && !keep_surface)
    return false;
  const auto *b = dynamic_cast<const box *>(geometry.get());
  if (b && !keep_surface) {
    if (!hit_invisible_box(*b, r, lo, hi, h))
      return false;
  } else if (!has_interior_interval(geometry.get(), r) || !geometry->hit(r, lo, hi, h, rng))
    return false;
  if (h.OrderedDistance() < r.medium_t_min) {
    // Analytic primitives report Float distances. Round their next query up;
    // triangles apply the precise limit within their intersection predicate.
    Float next = Float(r.medium_t_min);
    if (double(next) < r.medium_t_min) next = std::nextafter(next, Float(INFINITY));
    if (next <= lo || next > hi) return false;
    return hit(r, next, hi, h, rng);
  }
  Annotate(h);
  return true;
}
const bool MediumBoundary::hit(const Ray &r, Float lo, Float hi, hit_record &h,
                               Sampler *sampler) const {
  if (!r.segment_absorption && !keep_surface)
    return false;
  const auto *b = dynamic_cast<const box *>(geometry.get());
  if (b && !keep_surface) {
    if (!hit_invisible_box(*b, r, lo, hi, h))
      return false;
  } else if (!has_interior_interval(geometry.get(), r) || !geometry->hit(r, lo, hi, h, sampler))
    return false;
  if (h.OrderedDistance() < r.medium_t_min) {
    Float next = Float(r.medium_t_min);
    if (double(next) < r.medium_t_min) next = std::nextafter(next, Float(INFINITY));
    if (next <= lo || next > hi) return false;
    return hit(r, next, hi, h, sampler);
  }
  Annotate(h);
  return true;
}
double PriorityInterface::Eta() const {
  return (after ? double(after->ref_idx) : 1) / (before ? double(before->ref_idx) : 1);
}
const dielectric *VolumePathState::ActiveDielectric() const {
  const dielectric *active = nullptr;
  for (const auto *d : glass)
    if (!active || d->priority <= active->priority) active = d;
  return active;
}
bool VolumePathState::ContainsSubsurface() const {
  for (const auto &entry : media)
    if (entry.boundary->medium->subsurface) return true;
  return false;
}
const MediumEntry *VolumePathState::Active() const {
  if (media.empty()) return nullptr;
  // Explicit media retain their nesting semantics. SSS instead follows the
  // same winning dielectric as refraction and Beer-Lambert attenuation.
  if (!media.back().boundary->medium->subsurface) return &media.back();
  const dielectric *winner = ActiveDielectric();
  for (auto i = media.rbegin(); i != media.rend(); ++i) {
    if (!i->boundary->medium->subsurface) return &*i;
    if (i->surface == winner) return &*i;
  }
  return nullptr;
}
PriorityInterface VolumePathState::Interface(const hit_record &h, const vec3f &direction) const {
  const auto *d = static_cast<const dielectric *>(h.mat_ptr);
  PriorityInterface result;
  result.entering = dot(direction, h.geometric_normal) < 0;
  result.before = ActiveDielectric();
  // Remove only the last placement of this material on exit. Instances may
  // share a material pointer, so removing every occurrence loses containment.
  size_t removed = glass.size();
  if (!result.entering)
    for (size_t i = glass.size(); i > 0; --i)
      if (glass[i - 1] == d) { removed = i - 1; break; }
  for (size_t i = 0; i < glass.size(); ++i)
    if (i != removed && (!result.after || glass[i]->priority <= result.after->priority))
      result.after = glass[i];
  if (result.entering && (!result.after || d->priority <= result.after->priority)) result.after = d;
  return result;
}
void VolumePathState::CrossDielectric(const hit_record &h, const vec3f &direction) {
  auto *d = static_cast<dielectric *>(h.mat_ptr);
  auto interface = Interface(h, direction);
  if (h.medium_boundary && h.medium_boundary->medium &&
      h.medium_boundary->medium->subsurface && interface.Hidden() &&
      interface.before && interface.before != d) {
    bool contained = std::any_of(media.begin(), media.end(), [&](const MediumEntry &entry) {
      return entry.boundary == h.medium_boundary && entry.boundary_id == h.boundary_id;
    });
    // A grazing origin offset can cross a tiny wedge of a losing solid without
    // an intervening resolved hit. Its membership is a set, not a nesting count:
    // assign the known outgoing side idempotently. This is safe only while a
    // different dielectric wins on BOTH sides, so no active segment, IOR, or
    // extinction is changed. Real SSS interfaces retain strict crossing checks.
    if (contained == interface.entering) return;
  }
  const auto *previous = Active();
  uint64_t previous_id = previous ? previous->boundary_id : 0;
  if (dot(direction, h.geometric_normal) < 0) glass.push_back(d);
  else {
    auto i = std::find(glass.rbegin(), glass.rend(), d);
    if (i != glass.rend()) glass.erase(std::next(i).base());
  }
  Cross(h, direction);
  auto *next = Active();
  if (next && next->boundary_id != previous_id &&
      (!h.medium_boundary || next->boundary_id != h.boundary_id)) next->guide_valid = false;
}
void VolumePathState::Cross(const hit_record &h, const vec3f &direction) {
  if (!h.medium_boundary || !h.medium_boundary->medium)
    return;
  bool entering = dot(direction, h.geometric_normal) < 0;
  if (entering) {
    for (const auto &entry : media)
      if (entry.boundary == h.medium_boundary && entry.boundary_id == h.boundary_id) {
        std::ostringstream message;
        message << "Repeated entry into " << h.medium_boundary->geometry->GetName() << " at ("
                << h.p[0] << ", " << h.p[1] << ", " << h.p[2] << "), t=" << h.t
                << ". Check mesh orientation, self intersections, and nesting.";
        throw std::runtime_error(message.str());
      }
    if (!media.empty()) media.back().guide_valid = false;
    const auto *surface = h.mat_ptr && h.mat_ptr->is_dielectric()
                              ? static_cast<const dielectric *>(h.mat_ptr) : nullptr;
    media.push_back({h.medium_boundary, h.boundary_id, h.MediumToWorld(),
                     convert_to_vec3(h.geometric_normal), true, surface});
  } else {
    auto found = std::find_if(media.begin(), media.end(), [&](const MediumEntry &entry) {
      return entry.boundary == h.medium_boundary && entry.boundary_id == h.boundary_id;
    });
    bool valid = found != media.end();
    if (valid && std::next(found) != media.end()) {
      valid = h.medium_boundary->medium->subsurface;
      for (auto i = std::next(found); valid && i != media.end(); ++i)
        valid = i->boundary->medium->subsurface;
    }
    if (!valid) {
      std::ostringstream message;
      message << "Non-nested or inconsistently oriented medium boundaries encountered at ("
              << h.p[0] << ", " << h.p[1] << ", " << h.p[2] << ") exiting "
              << h.medium_boundary->geometry->GetName() << ". Active boundary: "
              << (media.empty() ? "vacuum" : media.back().boundary->geometry->GetName())
              << ". Use disjoint or nested closed volumes.";
      throw std::runtime_error(message.str());
    }
    media.erase(found);
    if (!media.empty()) media.back().guide_valid = false;
  }
}
uint64_t VolumeScene::ReserveBoundaryIds(uint64_t count) {
  if (count > std::numeric_limits<uint64_t>::max() - boundary_count)
    throw std::overflow_error("Too many boundary placements for 64-bit scene IDs.");
  uint64_t offset = boundary_count;
  boundary_count += count;
  return offset;
}
void VolumeScene::Finish(Float t0, Float t1) {
  if (!boundaries.objects.empty())
    boundary_bvh = std::make_shared<BVHAggregate>(boundaries.objects, t0, t1, 1, true);
}
VolumePathState VolumeScene::InitialState(const Ray &ray, const std::atomic<bool> *cancel,
                                          const std::function<bool()> &poll) const {
  VolumePathState state;
  if (!boundary_bvh)
    return state;
  vec3f direction = unit_vector(vec3f(1, 0.317f, 0.129f));
  Ray probe(ray.o, direction, ray.time());
  probe.segment_absorption = true;
  Float lower = 0;
  random_gen rng(0);
  std::vector<hit_record> crossings;
  bool first_crossing = true;
  // Trace outward from the exact camera origin, then replay crossings inward
  // from known vacuum. This avoids cancellation when constructing a distant
  // starting point whose line must terminate exactly at the camera.
  auto cancelled = [&] {
    return (cancel && cancel->load(std::memory_order_relaxed)) || (poll && poll());
  };
  while (!cancelled()) {
    hit_record h;
    if (!boundary_bvh->hit(probe, lower, MaxT, h, rng))
      break;
    // A ray starting exactly on a boundary still processes its t=0 hit.
    // Initialize the side immediately before that crossing, using the actual
    // ray direction rather than the arbitrary containment probe direction.
    // This also applies when a surface spawn rounds onto a nearby boundary.
    bool at_origin = first_crossing && h.t == 0;
    bool ray_outgoing = dot(ray.d, h.geometric_normal) > 0;
    bool probe_outgoing = dot(direction, h.geometric_normal) > 0;
    if (!at_origin || ray_outgoing == probe_outgoing)
      crossings.push_back(h);
    first_crossing = false;
    // Normal offsets change the probe's line and can jump over the exit of
    // a thin wedge at a pointed mesh corner. Advance only the query limit.
    probe.medium_t_min = std::nextafter(h.OrderedDistance(), INFINITY);
    lower = Float(probe.medium_t_min);
    if (double(lower) > probe.medium_t_min)
      lower = std::nextafter(lower, Float(-INFINITY));
  }
  for (auto i = crossings.rbegin(); i != crossings.rend(); ++i) {
    if (cancelled())
      break;
    const hit_record &h = *i;
    if (h.medium_boundary && h.medium_boundary->keep_surface && h.mat_ptr &&
        h.mat_ptr->is_dielectric()) {
      state.CrossDielectric(h, -direction);
    } else state.Cross(h, -direction);
  }
  // Containment probes are not physical entries and cannot supply a guide plane.
  for (auto &entry : state.media) entry.guide_valid = false;
  return state;
}
bool ValidateMediumBoundary(hitable *geometry, bool required) {
  if (dynamic_cast<sphere *>(geometry) || dynamic_cast<box *>(geometry) ||
      dynamic_cast<ellipsoid *>(geometry))
    return true;
  TriangleMesh *mesh = boundary_mesh(geometry);
  auto invalid = [&] {
    if (required)
      throw std::runtime_error("Medium boundary must be a supported closed shape or a watertight, "
                               "consistently oriented triangle mesh without alpha cutouts.");
    return false;
  };
  if (!mesh || !mesh->nTriangles)
    return invalid();
  for (const auto &a : mesh->alpha_textures)
    if (a)
      return invalid();
  // Weld equal positions so OBJ normal/UV seams do not create false open edges.
  std::map<std::tuple<Float, Float, Float>, size_t> positions;
  std::vector<size_t> ids(mesh->nVertices);
  for (size_t i = 0; i < mesh->nVertices; ++i) {
    const auto &p = mesh->p[i];
    for (int a = 0; a < 3; ++a)
      if (!std::isfinite(p[a]))
        return invalid();
    ids[i] = positions.emplace(std::make_tuple(p[0], p[1], p[2]), positions.size()).first->second;
  }
  std::map<std::pair<size_t, size_t>, std::pair<int, int>> edges;
  for (size_t f = 0; f < mesh->nTriangles; ++f)
    for (int j = 0; j < 3; ++j) {
      size_t a = ids[mesh->vertexIndices[3 * f + j]],
             b = ids[mesh->vertexIndices[3 * f + (j + 1) % 3]];
      if (a == b)
        return invalid();
      auto &count = edges[std::minmax(a, b)];
      ++count.first;
      count.second += a < b ? 1 : -1;
    }
  for (const auto &e : edges)
    if (e.second.first != 2 || e.second.second != 0)
      return invalid();
  double volume = signed_mesh_volume(*mesh);
  if (!std::isfinite(volume) || volume == 0)
    return invalid();
  return true;
}

void ValidateMediumTransform(const Transform &transform) {
  const auto &m = transform.GetMatrix();
  const auto &inv = transform.GetInverseMatrix();
  for (int r = 0; r < 4; ++r)
    for (int c = 0; c < 4; ++c)
      if (!std::isfinite(m.m[r][c]) || !std::isfinite(inv.m[r][c]))
        throw std::runtime_error("Medium boundary transforms must be finite and invertible.");
  if (m.m[3][0] != 0 || m.m[3][1] != 0 || m.m[3][2] != 0 || m.m[3][3] != 1)
    throw std::runtime_error("Medium boundary transforms must be affine.");
}

Rcpp::List VolumeScene::Statistics() const {
  size_t bytes = 0;
  for (const auto &entry : *medium_cache)
    if (entry.second)
      bytes += entry.second->MemoryBytes();
  auto event_quantile_upper = [&](double fraction) {
    uint64_t paths = statistics.subsurface_event_paths.load(), cumulative = 0;
    if (!paths) return 0.0;
    for (size_t i = 0; i < statistics.subsurface_event_histogram.size(); ++i) {
      cumulative += statistics.subsurface_event_histogram[i].load();
      if (double(cumulative) >= fraction * double(paths)) return std::ldexp(1.0, i);
    }
    return double(statistics.max_subsurface_events.load());
  };
  return Rcpp::List::create(
      Rcpp::Named("paths") = double(statistics.paths.load()),
      Rcpp::Named("majorant_segments") = double(statistics.segments.load()),
      Rcpp::Named("camera_null_events") = double(statistics.null_events.load()),
      Rcpp::Named("scattering_events") = double(statistics.scattering_events.load()),
      Rcpp::Named("shadow_candidates") = double(statistics.shadow_candidates.load()),
      Rcpp::Named("subsurface_events") = double(statistics.subsurface_events.load()),
      Rcpp::Named("subsurface_boundaries") = double(statistics.subsurface_boundaries.load()),
      Rcpp::Named("guide_eligible") = double(statistics.guide_eligible.load()),
      Rcpp::Named("guide_fallback") = double(statistics.guide_fallback.load()),
      Rcpp::Named("subsurface_intersections") = double(statistics.subsurface_intersections.load()),
      Rcpp::Named("max_subsurface_events") = double(statistics.max_subsurface_events.load()),
      Rcpp::Named("subsurface_event_paths") = double(statistics.subsurface_event_paths.load()),
      Rcpp::Named("total_subsurface_events") = double(statistics.total_subsurface_events.load()),
      Rcpp::Named("subsurface_events_p95_upper") = event_quantile_upper(.95),
      Rcpp::Named("subsurface_events_p99_upper") = event_quantile_upper(.99),
      Rcpp::Named("rounded_subsurface_flights") = double(statistics.rounded_subsurface_flights.load()),
      Rcpp::Named("medium_bytes") = double(bytes));
}
