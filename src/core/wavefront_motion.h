#ifndef RAYRENDER_WAVEFRONT_MOTION_H
#define RAYRENDER_WAVEFRONT_MOTION_H
#include "../math/animatedtransform.h"
#include "wavefront.h"
#include <array>
#include <limits>
#include <functional>
#include <unordered_map>

// Immutable geometry is traversed in each animated owner's reference frame.
// Camera paths retain independent shutter times through all transport stages.
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
    if (parent > nodes.size() || nodes.size() >= UINT32_MAX)
      throw std::runtime_error("Invalid animation parent or motion owner capacity exceeded.");
    // Preserve the cheap adjacent-repeat path, then reuse earlier motion IDs
    // across the whole scene. BVH traversal can interleave identical motions;
    // assigning each run a new owner makes every ray revisit the same motion.
    if (!nodes.empty()) {
      const auto &previous = nodes.back();
      if (previous.parent == parent && previous.before == before && previous.animation->SameMotion(animation))
        return nodes.size();
    }
    const size_t hash = MotionHash(parent, animation, before);
    const auto candidates = motion_lookup.equal_range(hash);
    for (auto entry = candidates.first; entry != candidates.second; ++entry) {
      const auto &previous = nodes[entry->second - 1];
      // Hash collisions must never merge distinct transforms, parents or times.
      if (previous.parent == parent && previous.before == before && previous.animation->SameMotion(animation))
        return entry->second;
    }
    Transform reference;
    animation.Interpolate(0, &reference);
    nodes.push_back({parent, &animation, before, Inverse(before * reference)});
    const uint32_t owner = nodes.size();
    motion_lookup.emplace(hash, owner);
    return owner;
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
  void Seal(WavefrontScene &scene) {
    for (const auto &collection : owners)
      for (uint32_t owner : collection)
        if (owner > nodes.size()) throw std::runtime_error("Invalid motion owner reference.");
    scene.texture_features |= 1024;
    packed = Export(scene);
  }
  // Native counterpart for the optional batched CPU atmosphere transport.
  Transform At(uint32_t owner, Float time) const {
    Transform result;
    while (owner) {
      const auto &node = nodes[owner - 1];
      Transform value;
      node.animation->Interpolate(time, &value);
      result = node.before * value * node.reference_inverse * result;
      owner = node.parent;
    }
    return result;
  }
  std::vector<unsigned char> packed;
  // Validate each submitted exposure, including later frames reusing storage.
  void ValidateTimes(double begin, double end) const;
  // Refresh only per-frame clock fields, never keys or immutable geometry.
  void UpdateClocks();
  void CopyClockState(void *destination) const;
  // Recompute only changed mesh boxes, then union placements in their original
  // traversal frames. Authored transforms and packed owner IDs never change.
  void RefitBounds(WavefrontScene &, const std::vector<WFMeshRange> &);

private:
  // Store IDs rather than keys or pointers into nodes: vector growth keeps IDs
  // valid, and the retained node already owns the transform comparison data.
  std::unordered_multimap<size_t, uint32_t> motion_lookup;
  std::vector<uint32_t> packed_transform_offsets;
  std::vector<aabb> mesh_bounds, loose_groups, loose_boundaries;
  std::vector<std::pair<uint32_t, uint32_t>> boundary_owners;
  std::vector<WFVector> ReferenceBounds(WavefrontScene &, const std::vector<WFMeshRange> &, bool initialize);
  void SweepBoundaries(WavefrontScene &) const;

  static size_t MotionHash(uint32_t parent, const AnimatedTransform &animation, const Transform &before) {
    size_t hash = std::hash<uint32_t>{}(parent);
    auto combine = [&](Float value) {
      // Numeric hashing gives +0 and -0 the same key, matching exact equality.
      hash ^= std::hash<Float>{}(value) + 0x9e3779b9u + (hash << 6) + (hash >> 2);
    };
    combine(animation.StartTime());
    combine(animation.EndTime());
    hash ^= std::hash<const void *>{}(animation.Track().get());
    hash ^= std::hash<double>{}(animation.TimeOrigin());
    for (const Transform *transform : {&before, &animation.StartTransform(), &animation.EndTransform()}) {
      // Transform equality includes the stored inverse as well as the forward
      // matrix; do not approximate motion or reconstruct a different inverse.
      for (const Matrix4x4 *matrix : {&transform->GetMatrix(), &transform->GetInverseMatrix()})
        for (int row = 0; row < 4; ++row)
          for (int column = 0; column < 4; ++column)
            combine(matrix->m[row][column]);
    }
    return hash;
  }

  std::vector<unsigned char> Export(WavefrontScene &);
};
#endif
