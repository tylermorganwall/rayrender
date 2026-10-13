#ifndef RAYRENDER_WAVEFRONT_MOTION_H
#define RAYRENDER_WAVEFRONT_MOTION_H
#include "../math/animatedtransform.h"
#include "wavefront.h"
#include <array>
#include <limits>

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
  void Seal(WavefrontScene &scene) {
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

private:
  std::vector<unsigned char> Export(WavefrontScene &);
};
#endif
