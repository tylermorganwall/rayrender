#include "wavefront_deformation.h"
#include "../hitables/mesh3d.h"
#include <cstring>

bool UpdateWavefrontTriangle(const triangle &tri, const Transform &placement, WFTriangle &out) {
  const auto &mesh = *tri.mesh;
  const auto a = mesh.p[tri.v[0]], b = mesh.p[tri.v[1]], c = mesh.p[tri.v[2]];
  normal3f geometric = convert_to_normal3(cross(b - a, c - a));
  const bool nondegenerate = geometric.squared_length() > 0;
  if (tri.reverseOrientation ^ tri.transformSwapsHandedness) geometric = -geometric;
  const auto world_normal = nondegenerate ? unit_vector(placement(geometric)) : normal3f(0, 0, 0);
  out.geometric.x = world_normal[0];
  out.geometric.y = world_normal[1];
  out.geometric.z = world_normal[2];
  bool smooth = mesh.has_normals && tri.n;
  for (int k = 0; k < 3 && smooth; ++k)
    smooth = tri.n[k] >= 0 && tri.n[k] < mesh.nNormals;
  out.flags = (out.flags & ~8u) | (smooth ? 8u : 0u);
  const auto &inverse = placement.GetInverseMatrix();
  for (int k = 0; k < 3; ++k) {
    const auto p = placement(mesh.p[tri.v[k]]);
    normal3f n = geometric;
    if (smooth) {
      n = mesh.n[tri.n[k]];
      if (tri.reverseOrientation) n = -n;
    }
    // Keep corner lengths under nonuniform scale: the CPU transforms the
    // interpolated normal, so independently normalizing corners changes it.
    normal3f transformed = n;
    if (!placement.IsIdentity())
      for (int axis = 0; axis < 3; ++axis)
        transformed[axis] = inverse.m[0][axis] * n[0] + inverse.m[1][axis] * n[1] +
                            inverse.m[2][axis] * n[2];
    else if (n.squared_length() > 0)
      transformed = placement(n);
    out.p[k].x = p[0]; out.p[k].y = p[1]; out.p[k].z = p[2];
    out.n[k].x = transformed[0]; out.n[k].y = transformed[1]; out.n[k].z = transformed[2];
  }
  return nondegenerate;
}

std::vector<WFMeshRange> WavefrontDeformation::Update(WavefrontScene &scene) {
  std::vector<WFMeshRange> changed;
  for (auto &mesh : meshes) {
    if (mesh.revision == mesh.source->DeformationRevision()) continue;
    for (size_t i = 0; i < mesh.faces.size(); ++i) {
      auto &out = scene.triangles[mesh.range.start + i];
      const auto &face = *mesh.faces[i];
      UpdateWavefrontTriangle(face, mesh.placement, out);
      out.geometric.x *= mesh.orientation;
      out.geometric.y *= mesh.orientation;
      out.geometric.z *= mesh.orientation;
      if (out.flags & 128) {
        uint32_t offset;
        std::memcpy(&offset, &out.geometric.w, sizeof(offset));
        scene.texels[offset] = {float(face.mesh->alpha_v[face.n[0]]),
                                float(face.mesh->alpha_v[face.n[1]]),
                                float(face.mesh->alpha_v[face.n[2]]), 0};
      }
    }
    mesh.revision = mesh.source->DeformationRevision();
    changed.push_back(mesh.range);
  }
  return changed;
}
