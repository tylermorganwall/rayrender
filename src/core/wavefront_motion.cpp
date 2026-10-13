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
void extend(WFVector &lo, WFVector &hi, const WFVector &p) {
  lo.x = std::min(lo.x, p.x);
  lo.y = std::min(lo.y, p.y);
  lo.z = std::min(lo.z, p.z);
  hi.x = std::max(hi.x, p.x);
  hi.y = std::max(hi.y, p.y);
  hi.z = std::max(hi.z, p.z);
}
} // namespace

std::vector<unsigned char> WavefrontMotion::Export(WavefrontScene &scene) {
  struct Header {
    uint32_t count, owners[5], nodes, bounds;
  } header{};
  struct Animation {
    WFVector endpoints[2][3], translation[2], rotation[2], scale[2][3];
    WFVector before[3], reference_inverse[3];
    uint32_t parent;
    float start, end;
    uint32_t padding = 0;
  };
  static_assert(sizeof(Header) == 32 && sizeof(Animation) == 368, "motion ABI mismatch");
  std::vector<unsigned char> bytes(sizeof(header));
  auto append = [&](const void *data, size_t size) {
    size_t offset = (bytes.size() + 15) & ~size_t(15);
    if (offset + size > UINT32_MAX)
      throw std::runtime_error("Metal motion metadata exceeds 32-bit capacity");
    bytes.resize(offset + size);
    if (size)
      std::memcpy(bytes.data() + offset, data, size);
    return uint32_t(offset);
  };
  auto counts = sizes(scene);
  for (size_t k = 0; k < 5; ++k) {
    owners[k].resize(counts[k]);
    header.owners[k] = append(owners[k].data(), owners[k].size() * sizeof(uint32_t));
  }
  std::vector<Animation> animations;
  for (const auto &node : nodes) {
    Animation a{};
    a.parent = node.parent;
    a.start = node.animation->StartTime();
    a.end = node.animation->EndTime();
    pack(node.before.GetMatrix(), a.before);
    pack(node.reference_inverse.GetMatrix(), a.reference_inverse);
    Quaternion rotations[2];
    for (int k = 0; k < 2; ++k) {
      const auto &t = k ? node.animation->EndTransform() : node.animation->StartTransform();
      pack(t.GetMatrix(), a.endpoints[k]);
      vec3f translation;
      Matrix4x4 scale;
      AnimatedTransform::Decompose(t.GetMatrix(), &translation, &rotations[k], &scale);
      a.translation[k] = {float(translation[0]), float(translation[1]), float(translation[2]), 0};
      pack(scale, a.scale[k]);
    }
    if (dot(rotations[0], rotations[1]) < 0)
      rotations[1] = -rotations[1];
    for (int k = 0; k < 2; ++k)
      a.rotation[k] = {float(rotations[k].v[0]), float(rotations[k].v[1]), float(rotations[k].v[2]),
                       float(rotations[k].w)};
    animations.push_back(a);
  }
  header.count = nodes.size();
  header.nodes = append(animations.data(), animations.size() * sizeof(Animation));

  // Bounds in each owner's reference frame reject unrelated animation groups
  // before querying the shared immutable hardware BVH. Prototype bounds are
  // computed once, irrespective of the number of placements.
  std::vector<aabb> groups(nodes.size() + 1), meshes(scene.meshes.size());
  for (size_t i = 0; i < scene.meshes.size(); ++i) {
    const auto &mesh = scene.meshes[i];
    for (uint32_t j = 0; j < mesh.count; ++j)
      for (const auto &p : scene.triangles[mesh.start + j].p)
        meshes[i] = surrounding_box(meshes[i], aabb(point3f(p.x, p.y, p.z)));
  }
  size_t loose = scene.meshes.empty() ? scene.triangles.size() : scene.loose_triangles;
  for (size_t i = 0; i < loose; ++i)
    for (const auto &p : scene.triangles[i].p)
      groups[owners[0][i]] = surrounding_box(groups[owners[0][i]], aabb(point3f(p.x, p.y, p.z)));
  for (size_t i = 0; i < scene.quadrics.size(); ++i) {
    const auto &lo = scene.quadric_bounds[2 * i], &hi = scene.quadric_bounds[2 * i + 1];
    groups[owners[1][i]] =
        surrounding_box(groups[owners[1][i]], aabb(point3f(lo.x, lo.y, lo.z), point3f(hi.x, hi.y, hi.z)));
  }
  for (size_t i = 0; i < scene.instances.size(); ++i) {
    const auto &instance = scene.instances[i];
    if (instance.triangle_count)
      groups[owners[3][i]] = surrounding_box(groups[owners[3][i]], unpack(instance.forward)(meshes[instance.mesh]));
  }
  std::vector<WFVector> boxes;
  for (const auto &box : groups) {
    auto lo = box.min(), hi = box.max();
    if (lo[0] > hi[0]) {
      boxes.push_back({INFINITY, INFINITY, INFINITY, 0});
      boxes.push_back({-INFINITY, -INFINITY, -INFINITY, 0});
      continue;
    }
    Float pad = 64 * std::numeric_limits<float>::epsilon() *
                std::max(Float(1), std::max(MaxComponent(Abs(lo)), MaxComponent(Abs(hi))));
    boxes.push_back({float(lo[0] - pad), float(lo[1] - pad), float(lo[2] - pad), 0});
    boxes.push_back({float(hi[0] + pad), float(hi[1] + pad), float(hi[2] + pad), 0});
  }
  header.bounds = append(boxes.data(), boxes.size() * sizeof(WFVector));
  std::memcpy(bytes.data(), &header, sizeof(header));

  // Membership and diffusion bounds enclose the full native motion, including
  // nested rotation. They are only broad-phase bounds; membership still tests
  // the exact geometry at the path's time.
  std::set<std::pair<uint32_t, uint32_t>> bodies;
  auto collect = [&](const auto &values, size_t kind) {
    for (size_t i = 0; i < values.size(); ++i)
      if (values[i].boundary && owners[kind][i])
        bodies.emplace(values[i].boundary, owners[kind][i]);
  };
  collect(scene.triangles, 0);
  collect(scene.quadrics, 1);
  collect(scene.instances, 3);
  const auto reference_bounds = scene.boundaries;
  for (const auto &body : bodies) {
    const auto &b = reference_bounds[body.first];
    aabb box(point3f(b.lo.x, b.lo.y, b.lo.z), point3f(b.hi.x, b.hi.y, b.hi.z));
    for (uint32_t owner = body.second; owner; owner = nodes[owner - 1].parent) {
      const auto &node = nodes[owner - 1];
      box = node.before(node.animation->MotionBounds(node.reference_inverse(box)));
    }
    auto lo = box.min(), hi = box.max();
    float pad = 64 * std::numeric_limits<float>::epsilon();
    for (int k = 0; k < 3; ++k)
      pad =
          std::max(pad, 64 * std::numeric_limits<float>::epsilon() * float(std::max(std::abs(lo[k]), std::abs(hi[k]))));
    WFVector l{float(lo[0]) - pad, float(lo[1]) - pad, float(lo[2]) - pad, 0};
    WFVector h{float(hi[0]) + pad, float(hi[1]) + pad, float(hi[2]) + pad, 0};
    extend(scene.boundaries[body.first].lo, scene.boundaries[body.first].hi, l);
    extend(scene.boundaries[body.first].lo, scene.boundaries[body.first].hi, h);
    extend(scene.boundary_lo, scene.boundary_hi, l);
    extend(scene.boundary_lo, scene.boundary_hi, h);
  }
  return bytes;
}
