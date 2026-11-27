#include "model_deformation.h"
#include <stdexcept>

namespace animation {
namespace {
template<class T> bool finite_vector(const T &v) {
  return std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]);
}
void require(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}
}
void DeformationMesh::Validate(size_t joints) const {
  const size_t count = positions.size();
  require(count > 0 && indices.size() % 3 == 0, "Deformation requires triangle geometry.");
  for (const auto &v : positions) require(finite_vector(v), "Nonfinite rest position.");
  require(normals.empty() || normals.size() == count, "Normal count differs from positions.");
  require(tangents.empty() || (tangents.size() == count && handedness.size() == count && !normals.empty()),
          "Tangents require matching normals and handedness.");
  for (const auto &v : normals) require(finite_vector(v) && v.squared_length() > 0, "Invalid rest normal.");
  for (const auto &v : tangents) require(finite_vector(v) && v.squared_length() > 0, "Invalid rest tangent.");
  for (double h : handedness) require(h == -1 || h == 1, "Tangent handedness must be -1 or 1.");
  for (int index : indices) require(index >= 0 && size_t(index) < count, "Deformation index out of range.");
  for (const auto &target : targets) {
    for (const auto *values : {&target.positions, &target.normals, &target.tangents}) {
      require(values->empty() || values->size() == count, "Morph target count differs from positions.");
      for (const auto &v : *values) require(finite_vector(v), "Nonfinite morph delta.");
    }
    require(target.normals.empty() || !normals.empty(), "Morph normals require source normals.");
    require(target.tangents.empty() || !tangents.empty(), "Morph tangents require source tangents.");
  }
  require(influences.empty() || (joints && influences.size() == count), "Invalid skin influence count.");
  for (const auto &vertex : influences) {
    double sum = 0;
    for (const auto &influence : vertex) {
      require(influence.joint < joints && std::isfinite(influence.weight) && influence.weight >= 0,
              "Invalid joint index or skin weight.");
      sum += influence.weight;
    }
    require(std::isfinite(sum) && sum > 0, "Skin vertex has zero total weight.");
  }
}

DeformedMesh DeformMesh(const DeformationMesh &mesh, const std::vector<double> &weights,
                       const std::vector<Transform> &palette) {
  require(weights.size() == mesh.targets.size(), "Morph weight count differs from targets.");
  for (double weight : weights) require(std::isfinite(weight), "Nonfinite evaluated morph weight.");
  DeformedMesh result{mesh.positions, mesh.normals, mesh.tangents, mesh.handedness};
  for (size_t i = 0; i < result.positions.size(); ++i) {
    vec3f normal, tangent;
    if (!mesh.normals.empty()) normal = vec3f(mesh.normals[i][0], mesh.normals[i][1], mesh.normals[i][2]);
    if (!mesh.tangents.empty()) tangent = mesh.tangents[i];
    for (size_t target = 0; target < weights.size(); ++target) {
      const auto &delta = mesh.targets[target];
      if (!delta.positions.empty()) result.positions[i] += Float(weights[target]) * delta.positions[i];
      if (!delta.normals.empty()) normal += Float(weights[target]) * delta.normals[i];
      if (!delta.tangents.empty()) tangent += Float(weights[target]) * delta.tangents[i];
    }
    double determinant = 1;
    if (!mesh.influences.empty()) {
      // Accumulate the affine map in double precision and normalize all weights
      // together, including additional influence sets and quantized weights.
      double matrix[3][4] = {}, total = 0;
      for (const auto &influence : mesh.influences[i]) {
        require(influence.joint < palette.size(), "Skin palette is incomplete.");
        const auto &joint = palette[influence.joint].GetMatrix();
        total += influence.weight;
        for (int r = 0; r < 3; ++r) for (int c = 0; c < 4; ++c)
          matrix[r][c] += influence.weight * joint.m[r][c];
      }
      require(total > 0, "Skin vertex has zero total weight.");
      for (auto &row : matrix) for (double &value : row) value /= total;
      const auto &m = matrix;
      determinant = m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1]) -
        m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0]) +
        m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
      require(std::isfinite(determinant) && std::abs(determinant) > 1e-12,
              "Evaluated skin blend is singular; normals cannot be transformed exactly.");
      Float affine[4][4] = {{0}}; affine[3][3] = 1;
      for (int r = 0; r < 3; ++r) for (int c = 0; c < 4; ++c) affine[r][c] = matrix[r][c];
      Transform skin(affine);
      result.positions[i] = skin(result.positions[i]);
      if (!mesh.normals.empty()) {
        auto n = skin(normal3f(normal[0], normal[1], normal[2]));
        normal = vec3f(n[0], n[1], n[2]);
      }
      if (!mesh.tangents.empty()) tangent = skin(tangent);
    }
    require(finite_vector(result.positions[i]), "Nonfinite deformed position.");
    if (!mesh.normals.empty() && !mesh.recalculate_normals) {
      require(finite_vector(normal) && normal.squared_length() > 1e-20, "Degenerate deformed normal.");
      normal = unit_vector(normal);
      result.normals[i] = normal3f(normal[0], normal[1], normal[2]);
    }
    if (!mesh.tangents.empty()) {
      if (!mesh.recalculate_normals) tangent -= dot(tangent, normal) * normal;
      require(finite_vector(tangent) && tangent.squared_length() > 1e-20, "Degenerate deformed tangent.");
      result.tangents[i] = unit_vector(tangent);
      result.handedness[i] *= determinant < 0 ? -1 : 1;
    }
  }
  if (mesh.recalculate_normals) {
    // Rebuild from this pose only. Corner-angle weighting avoids making shading
    // depend on the areas of skinny neighboring triangles. Index splits retain
    // hard edges; we deliberately do not weld coincident imported vertices.
    std::vector<vec3f> sums(result.positions.size(), vec3f(0,0,0));
    for (size_t face = 0; face < mesh.indices.size(); face += 3) {
      for (int corner = 0; corner < 3; ++corner) {
        const int vertex = mesh.indices[face + corner];
        const auto &origin = result.positions[vertex];
        const auto a = result.positions[mesh.indices[face + (corner+1)%3]] - origin;
        const auto b = result.positions[mesh.indices[face + (corner+2)%3]] - origin;
        const auto area = cross(a,b);
        const double length = area.length();
        // Collapsed faces contribute nothing; an undefined vertex normal fails
        // before publication rather than retaining a previous frame's normal.
        if (!(length > 0)) continue;
        const double angle = std::atan2(length, double(dot(a,b)));
        sums[vertex] += Float(angle/length) * area;
      }
    }
    result.normals.resize(result.positions.size());
    for (size_t i = 0; i < sums.size(); ++i) {
      require(finite_vector(sums[i]) && sums[i].squared_length() > 1e-20,
              "Cannot recalculate a normal without a nondegenerate incident surface.");
      const auto normal = unit_vector(sums[i]);
      result.normals[i] = normal3f(normal[0],normal[1],normal[2]);
      if (!result.tangents.empty()) {
        auto tangent = result.tangents[i] - dot(result.tangents[i],normal)*normal;
        require(finite_vector(tangent) && tangent.squared_length() > 1e-20,
                "Degenerate tangent after normal recalculation.");
        result.tangents[i] = unit_vector(tangent);
      }
    }
  }
  return result;
}
} // namespace animation
