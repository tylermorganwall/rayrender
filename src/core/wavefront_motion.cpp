#include "wavefront_motion.h"
#include <cstring>
#include <set>

namespace {
Transform unpack(const WFVector *rows) {
  Float m[4][4] = {{0}};
  m[3][3] = 1;
  for (int i = 0; i < 3; ++i) {
    m[i][0] = rows[i].x;
    m[i][1] = rows[i].y;
    m[i][2] = rows[i].z;
    m[i][3] = rows[i].w;
  }
  return Transform(m);
}
void pack(const Matrix4x4 &m, WFVector *rows) {
  for (int i = 0; i < 3; ++i)
    rows[i] = {float(m.m[i][0]), float(m.m[i][1]), float(m.m[i][2]), float(m.m[i][3])};
}
void compose(const Transform &delta, WFVector *forward, WFVector *inverse) {
  Transform result = delta * unpack(forward);
  pack(result.GetMatrix(), forward);
  pack(result.GetInverseMatrix(), inverse);
}
void point(const Transform &t, WFVector &v) {
  auto p = t(point3f(v.x, v.y, v.z));
  v.x = p[0];
  v.y = p[1];
  v.z = p[2];
}
void normal(const Transform &t, WFVector &v) {
  auto n = t(normal3f(v.x, v.y, v.z));
  v.x = n[0];
  v.y = n[1];
  v.z = n[2];
}
void extend(WFVector &lo, WFVector &hi, const WFVector &p) {
  lo.x = std::min(lo.x, p.x);
  lo.y = std::min(lo.y, p.y);
  lo.z = std::min(lo.z, p.z);
  hi.x = std::max(hi.x, p.x);
  hi.y = std::max(hi.y, p.y);
  hi.z = std::max(hi.z, p.z);
}
} // namespace

void WavefrontMotion::EnableGPU() {
  // Keep only geometry needed by the host light-selection tree. Boundary
  // bounds can be conservatively transformed once per body/owner, rather than
  // revisiting every triangle. Metal already owns all immutable references.
  std::set<std::pair<uint32_t, uint32_t>> groups;
  auto collect = [&](const auto &records) {
    for (const auto &r : records)
      if (r.value.boundary)
        groups.emplace(r.value.boundary, r.owner);
  };
  collect(triangles);
  collect(quadrics);
  collect(instances);
  boundary_owners.assign(groups.begin(), groups.end());
  host_triangles.clear();
  host_quadrics.clear();
  for (const auto &r : triangles)
    if (r.value.light != UINT32_MAX)
      host_triangles.push_back(r);
  for (const auto &r : quadrics)
    if (r.value.light != UINT32_MAX)
      host_quadrics.push_back(r);
  gpu_geometry = true;
}

bool WavefrontMotion::Update(WavefrontScene &scene, float time) {
  if (time == last_time)
    return false;
  std::vector<Transform> transforms(nodes.size() + 1);
  for (size_t i = 0; i < nodes.size(); ++i) {
    const auto &n = nodes[i];
    Transform value;
    n.animation->Interpolate(time, &value);
    transforms[i + 1] = transforms[n.parent] * n.before * value * n.reference_inverse;
  }
  if (whole_scene_frame)
    pack(transforms[1].GetInverseMatrix(), inverse_frame);
  poses.resize(transforms.size());
  for (size_t i = 0; i < transforms.size(); ++i) {
    pack(transforms[i].GetMatrix(), poses[i].forward);
    pack(transforms[i].GetInverseMatrix(), poses[i].inverse);
  }
  scene.boundaries = reference_bounds;
  scene.boundary_lo = reference_lo;
  scene.boundary_hi = reference_hi;
  auto bounds = [&](uint32_t body, const WFVector &p) {
    if (!body)
      return;
    extend(scene.boundaries[body].lo, scene.boundaries[body].hi, p);
    extend(scene.boundary_lo, scene.boundary_hi, p);
  };
  if (gpu_geometry)
    for (const auto &group : boundary_owners) {
      const auto &reference = reference_bounds[group.first];
      for (int corner = 0; corner < 8; ++corner) {
        WFVector p{corner & 1 ? reference.hi.x : reference.lo.x, corner & 2 ? reference.hi.y : reference.lo.y,
                   corner & 4 ? reference.hi.z : reference.lo.z, 0};
        point(transforms[group.second], p);
        // Device positions use float arithmetic. Expand host membership bounds
        // past both rounding paths, including cancellation around the origin.
        float scale = 1;
        for (const auto &row : poses[group.second].forward)
          scale = std::max(scale, std::abs(row.w) +
                                      std::abs(row.x) * std::max(std::abs(reference.lo.x), std::abs(reference.hi.x)) +
                                      std::abs(row.y) * std::max(std::abs(reference.lo.y), std::abs(reference.hi.y)) +
                                      std::abs(row.z) * std::max(std::abs(reference.lo.z), std::abs(reference.hi.z)));
        float pad = 64 * std::numeric_limits<float>::epsilon() * scale;
        bounds(group.first, {p.x - pad, p.y - pad, p.z - pad, 0});
        bounds(group.first, {p.x + pad, p.y + pad, p.z + pad, 0});
      }
    }
  for (const auto &r : gpu_geometry ? host_triangles : triangles) {
    auto value = r.value;
    const auto &t = transforms[r.owner];
    for (int k = 0; k < 3; ++k) {
      point(t, value.p[k]);
      normal(t, value.n[k]);
      bounds(value.boundary, value.p[k]);
    }
    normal(t, value.geometric);
    scene.triangles[r.index] = value;
  }
  for (const auto &r : gpu_geometry ? host_quadrics : quadrics) {
    auto value = r.value;
    compose(transforms[r.owner], value.forward, value.inverse);
    auto m = unpack(value.forward).GetMatrix();
    vec3f x(m.m[0][0], m.m[1][0], m.m[2][0]), y(m.m[0][1], m.m[1][1], m.m[2][1]), z(m.m[0][2], m.m[1][2], m.m[2][2]);
    value.measure.x = std::abs(dot(x, cross(y, z)));
    Float radius = x.length();
    value.measure.z = std::abs(y.length() - radius) < 1e-6 * radius && std::abs(z.length() - radius) < 1e-6 * radius &&
                              std::abs(dot(x, y)) + std::abs(dot(x, z)) + std::abs(dot(y, z)) < 1e-6 * radius * radius
                          ? radius
                          : 0;
    std::array<float, 3> l, h;
    for (int a = 0; a < 3; ++a) {
      const auto &row = value.forward[a];
      // Match the compiler's conservative bounds: forward/inverse transforms
      // round separately on upload, including after a shutter pose update.
      double extent = std::sqrt(double(row.x) * row.x + double(row.y) * row.y + double(row.z) * row.z);
      double padding = 8 * std::numeric_limits<float>::epsilon() * (std::abs(row.w) + extent);
      l[a] = std::nextafter(float(double(row.w) - extent - padding), -INFINITY);
      h[a] = std::nextafter(float(double(row.w) + extent + padding), INFINITY);
    }
    WFVector lo{l[0], l[1], l[2], 0}, hi{h[0], h[1], h[2], 0};
    if (value.csg_count || value.curve_offset)
      WavefrontCSGBounds(value, lo, hi);
    scene.quadric_bounds[2 * r.index] = lo;
    scene.quadric_bounds[2 * r.index + 1] = hi;
    bounds(value.boundary, lo);
    bounds(value.boundary, hi);
    scene.quadrics[r.index] = value;
  }
  if (!gpu_geometry)
    for (const auto &r : surfaces) {
      auto value = r.value;
      const auto &t = transforms[r.owner];
      compose(t, value.forward, value.inverse);
      compose(t, value.bump_forward, value.bump_inverse);
      value.bump_transformed = 1;
      scene.surfaces[r.index] = value;
    }
  for (const auto &r : instances) {
    auto value = r.value;
    const auto &t = transforms[r.owner];
    compose(t, value.forward, value.inverse);
    scene.instances[r.index] = value;
    if (value.boundary) {
      const auto &body = reference_bounds[value.boundary];
      for (int c = 0; c < 8; ++c) {
        WFVector p = {c & 1 ? body.hi.x : body.lo.x, c & 2 ? body.hi.y : body.lo.y, c & 4 ? body.hi.z : body.lo.z, 0};
        point(t, p);
        bounds(value.boundary, p);
      }
    }
  }
  for (const auto &r : media) {
    auto value = r.value;
    auto transform = unpack(value.inverse) * Inverse(transforms[r.owner]);
    pack(transform.GetMatrix(), value.inverse);
    scene.media[r.index] = value;
  }
  last_time = time;
  return true;
}
