#ifndef RAYRENDER_MODEL_DEFORMATION_H
#define RAYRENDER_MODEL_DEFORMATION_H

#include "keyframe.h"

namespace animation {
struct SkinInfluence { size_t joint; double weight; };
struct MorphTarget {
  std::vector<vec3f> positions, normals, tangents;
};
// Immutable mesh-local source arrays. Every evaluation starts here, never at
// the previous pose. Influences retain every JOINTS_n/WEIGHTS_n set.
struct DeformationMesh {
  std::vector<point3f> positions;
  std::vector<normal3f> normals;
  std::vector<vec3f> tangents;
  std::vector<double> handedness;
  std::vector<int> indices;
  std::vector<MorphTarget> targets;
  std::vector<std::vector<SkinInfluence>> influences;
  bool recalculate_normals = false;
  void Validate(size_t joints) const;
};
struct ModelSkin {
  std::vector<size_t> joints;
  std::vector<Transform> inverse_bind;
};
struct DeformedMesh {
  std::vector<point3f> positions;
  std::vector<normal3f> normals;
  std::vector<vec3f> tangents;
  std::vector<double> handedness;
};
// Palette maps source mesh coordinates back into the posed mesh node frame:
// inverse(mesh_global) * joint_global * inverse_bind. Outer placement is absent.
// Morph first, then linear blend skinning. Normals use the inverse transpose of
// the blended linear map; tangents use that map and Gram-Schmidt against normals.
DeformedMesh DeformMesh(const DeformationMesh &, const std::vector<double> &weights,
                       const std::vector<Transform> &palette);
} // namespace animation
#endif
