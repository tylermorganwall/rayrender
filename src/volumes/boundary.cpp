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
// Face-by-face box intersections can report a contact at an exterior edge as
// an exit. A contact with no resolvable interior interval is not a crossing.
bool has_interior_interval(const hitable *geometry, const Ray &ray) {
  const auto *b = dynamic_cast<const box *>(geometry);
  if (!b)
    return true;
  point3f o = (*b->WorldToObject)(ray.o);
  vec3f d = (*b->WorldToObject)(ray.d);
  double near = -INFINITY, far = INFINITY;
  for (int axis = 0; axis < 3; ++axis) {
    if (d[axis] == 0) {
      if (o[axis] <= b->pmin[axis] || o[axis] >= b->pmax[axis])
        return false;
      continue;
    }
    double a = (double(b->pmin[axis]) - o[axis]) / d[axis];
    double z = (double(b->pmax[axis]) - o[axis]) / d[axis];
    if (a > z)
      std::swap(a, z);
    near = std::max(near, a);
    far = std::min(far, z);
  }
  double tolerance =
      8 * std::numeric_limits<Float>::epsilon() * std::max({1.0, std::abs(near), std::abs(far)});
  return far - near > tolerance;
}
} // namespace

point3f OffsetMediumOrigin(const hit_record &h, const vec3f &direction) {
  normal3f n = h.geometric_normal.squared_length() > 0 ? h.geometric_normal : h.normal;
  point3f origin = OffsetRayOrigin(h.p, h.pError, n, direction);
  // At a box edge or corner, offset every incident face. Moving off only the
  // selected rectangle leaves a second face at t=0, causing a duplicate crossing.
  if (h.medium_boundary)
    if (auto *box_geometry = dynamic_cast<const box *>(h.medium_boundary->geometry.get())) {
      Transform object_to_world = h.medium_to_world;
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
  if (!has_interior_interval(geometry.get(), r) || !geometry->hit(r, lo, hi, h, rng))
    return false;
  Annotate(h);
  return true;
}
const bool MediumBoundary::hit(const Ray &r, Float lo, Float hi, hit_record &h,
                               Sampler *sampler) const {
  if (!r.segment_absorption && !keep_surface)
    return false;
  if (!has_interior_interval(geometry.get(), r) || !geometry->hit(r, lo, hi, h, sampler))
    return false;
  Annotate(h);
  return true;
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
    media.push_back({h.medium_boundary, h.boundary_id, h.medium_to_world});
  } else {
    if (media.empty() || media.back().boundary != h.medium_boundary ||
        media.back().boundary_id != h.boundary_id) {
      std::ostringstream message;
      message << "Non-nested or inconsistently oriented medium boundaries encountered at ("
              << h.p[0] << ", " << h.p[1] << ", " << h.p[2] << ") exiting "
              << h.medium_boundary->geometry->GetName() << ". Active boundary: "
              << (media.empty() ? "vacuum" : media.back().boundary->geometry->GetName())
              << ". Use disjoint or nested closed volumes.";
      throw std::runtime_error(message.str());
    }
    media.pop_back();
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
VolumePathState VolumeScene::InitialState(const Ray &ray, const std::atomic<bool> *cancel) const {
  VolumePathState state;
  if (!boundary_bvh)
    return state;
  vec3f direction = unit_vector(vec3f(1, 0.317f, 0.129f));
  Ray probe(ray.o, direction, ray.time());
  probe.segment_absorption = true;
  random_gen rng(0);
  std::vector<hit_record> crossings;
  bool first_crossing = true;
  // Trace outward from the exact camera origin, then replay crossings inward
  // from known vacuum. This avoids cancellation when constructing a distant
  // starting point whose line must terminate exactly at the camera.
  while (!(cancel && cancel->load(std::memory_order_relaxed))) {
    hit_record h;
    if (!boundary_bvh->hit(probe, 0, MaxT, h, rng))
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
    point3f origin = OffsetMediumOrigin(h, direction);
    probe = Ray(origin, direction, ray.time());
    probe.segment_absorption = true;
  }
  for (auto i = crossings.rbegin(); i != crossings.rend(); ++i) {
    if (cancel && cancel->load(std::memory_order_relaxed))
      break;
    const hit_record &h = *i;
    state.Cross(h, -direction);
    if (h.medium_boundary && h.medium_boundary->keep_surface && h.mat_ptr &&
        h.mat_ptr->is_dielectric()) {
      auto *d = static_cast<dielectric *>(h.mat_ptr);
      if (dot(-direction, h.geometric_normal) < 0)
        state.glass.push_back(d);
      else {
        auto j = std::find(state.glass.rbegin(), state.glass.rend(), d);
        if (j != state.glass.rend())
          state.glass.erase(std::next(j).base());
      }
    }
  }
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
  return Rcpp::List::create(
      Rcpp::Named("paths") = double(statistics.paths.load()),
      Rcpp::Named("majorant_segments") = double(statistics.segments.load()),
      Rcpp::Named("camera_null_events") = double(statistics.null_events.load()),
      Rcpp::Named("scattering_events") = double(statistics.scattering_events.load()),
      Rcpp::Named("shadow_candidates") = double(statistics.shadow_candidates.load()),
      Rcpp::Named("medium_bytes") = double(bytes));
}
