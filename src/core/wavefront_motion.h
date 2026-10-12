#ifndef RAYRENDER_WAVEFRONT_MOTION_H
#define RAYRENDER_WAVEFRONT_MOTION_H
#include "../math/animatedtransform.h"
#include "wavefront.h"
#include <array>
#include <limits>

// Moving geometry uses coherent shutter samples: one exact instant per image
// sample, shared by camera, visibility, lights and media. This is Monte Carlo
// integration over time, not a finite set of baked motion poses. It trades
// spatially correlated temporal noise for exact native transform interpolation
// (including nested rotation/shear) on existing devices. Independent placements
// refit only changing structures; a shared whole-scene motion transforms rays
// into one immutable reference BVH instead.
class WavefrontMotion {
public:
  struct Node {
    uint32_t parent;
    const AnimatedTransform *animation;
    Transform before, reference_inverse;
  };
  std::vector<Node> nodes;
  // Collection order: triangles, quadrics, surfaces, instances, media.
  using Sizes = std::array<size_t, 5>;
  std::array<std::vector<uint32_t>, 5> owners;
  static Sizes sizes(const WavefrontScene &s) {
    return {s.triangles.size(), s.quadrics.size(), s.surfaces.size(), s.instances.size(), s.media.size()};
  }
  uint32_t Add(uint32_t parent, const AnimatedTransform &animation, const Transform &before) {
    // animate_objects() attaches an identical animation to each row. Coalesce
    // consecutive equal placements in constant time, preserving distinct times
    // and parents; do not introduce a quadratic search for independent motion.
    if (!nodes.empty()) {
      const auto &previous = nodes.back();
      if (previous.parent == parent && previous.before == before && previous.animation->SameMotion(animation))
        return nodes.size();
    }
    Transform reference;
    animation.Interpolate(0, &reference);
    nodes.push_back({parent, &animation, before, Inverse(before * reference)});
    return nodes.size();
  }
  void Claim(Sizes begin, Sizes end, uint32_t owner) {
    for (size_t k = 0; k < 5; ++k) {
      owners[k].resize(end[k]);
      for (size_t i = begin[k]; i < end[k]; ++i)
        if (!owners[k][i])
          owners[k][i] = owner;
    }
  }
  void PrefixInstances(size_t prefix) { owners[3].insert(owners[3].begin(), prefix, 0); }
  template <class T> struct Record {
    uint32_t index, owner;
    T value;
  };
  template <class T>
  static std::vector<Record<T>> Save(const std::vector<T> &values, const std::vector<uint32_t> &ids) {
    std::vector<Record<T>> result;
    for (size_t i = 0; i < ids.size(); ++i)
      if (ids[i])
        result.push_back({uint32_t(i), ids[i], values[i]});
    return result;
  }
  bool whole_scene_frame = false;
  WFVector inverse_frame[3] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}};
  void Seal(WavefrontScene &s) {
    gpu_geometry = false;
    last_time = std::numeric_limits<float>::quiet_NaN();
    host_triangles.clear();
    host_quadrics.clear();
    boundary_owners.clear();
    triangles = Save(s.triangles, owners[0]);
    quadrics = Save(s.quadrics, owners[1]);
    surfaces = Save(s.surfaces, owners[2]);
    instances = Save(s.instances, owners[3]);
    media = Save(s.media, owners[4]);
    reference_bounds = s.boundaries;
    reference_lo = s.boundary_lo;
    reference_hi = s.boundary_hi;
    // A single transform owning every rendered primitive can move rays instead
    // of the acceleration structure. Prototype vertices are owned by placements.
    auto owns = [&](size_t kind, size_t count) {
      for (size_t i = 0; i < count; ++i)
        if (i >= owners[kind].size() || owners[kind][i] != 1)
          return false;
      return true;
    };
    whole_scene_frame = nodes.size() == 1 && owns(0, s.meshes.empty() ? s.triangles.size() : s.loose_triangles) &&
                        owns(1, s.quadrics.size());
    for (size_t i = 0; i < s.instances.size() && whole_scene_frame; ++i)
      if (s.instances[i].triangle_count && (i >= owners[3].size() || owners[3][i] != 1))
        whole_scene_frame = false;
    s.texture_features &= ~32u;
    if (whole_scene_frame)
      s.texture_features |= 32; // Specialize reference-frame traversal.
  }
  bool Update(WavefrontScene &, float time);

  // Per-owner affine deltas are small host metadata. Metal applies them to
  // immutable geometry records, then refits the affected acceleration trees.
  struct Pose {
    WFVector forward[3], inverse[3];
  };
  std::vector<Pose> poses;
  const auto &TriangleReferences() const { return triangles; }
  const auto &QuadricReferences() const { return quadrics; }
  const auto &SurfaceReferences() const { return surfaces; }
  const auto &InstanceReferences() const { return instances; }
  void EnableGPU();
  bool gpu_geometry = false;

private:
  std::vector<Record<WFTriangle>> triangles;
  std::vector<Record<WFQuadric>> quadrics;
  std::vector<Record<WFTriangle>> host_triangles;
  std::vector<Record<WFQuadric>> host_quadrics;
  std::vector<Record<WFSurface>> surfaces;
  std::vector<Record<WFMeshInstance>> instances;
  std::vector<Record<WFMedium>> media;
  std::vector<WFBoundary> reference_bounds;
  std::vector<std::pair<uint32_t, uint32_t>> boundary_owners;
  WFVector reference_lo, reference_hi;
  float last_time = std::numeric_limits<float>::quiet_NaN();
};
#endif
