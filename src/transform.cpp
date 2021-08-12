#include "transform.h"

#include "mathinline.h"
#include "aabb.h"

Transform Translate(const vec3f &delta) {
  Matrix4x4 m(1, 0, 0, delta.x(), 0, 1, 0, delta.y(), 0, 0, 1, delta.z(), 0, 0, 0,
              1);
  Matrix4x4 minv(1, 0, 0, -delta.x(), 0, 1, 0, -delta.y(), 0, 0, 1, -delta.z(), 0,
                 0, 0, 1);
  return Transform(m, minv);
}

Transform Scale(Float x, Float y, Float z) {
  Matrix4x4 m(x, 0, 0, 0, 0, y, 0, 0, 0, 0, z, 0, 0, 0, 0, 1);
  Matrix4x4 minv(1 / x, 0, 0, 0, 0, 1 / y, 0, 0, 0, 0, 1 / z, 0, 0, 0, 0, 1);
  return Transform(m, minv);
}

Transform RotateX(Float theta) {
  Float sinTheta = std::sin(Radians(theta));
  Float cosTheta = std::cos(Radians(theta));
  Matrix4x4 m(1, 0, 0, 0, 0, cosTheta, -sinTheta, 0, 0, sinTheta, cosTheta, 0,
              0, 0, 0, 1);
  return Transform(m, Transpose(m));
}

Transform RotateY(Float theta) {
  Float sinTheta = std::sin(Radians(theta));
  Float cosTheta = std::cos(Radians(theta));
  Matrix4x4 m(cosTheta, 0, sinTheta, 0, 0, 1, 0, 0, -sinTheta, 0, cosTheta, 0,
              0, 0, 0, 1);
  return Transform(m, Transpose(m));
}

Transform RotateZ(Float theta) {
  Float sinTheta = std::sin(Radians(theta));
  Float cosTheta = std::cos(Radians(theta));
  Matrix4x4 m(cosTheta, -sinTheta, 0, 0, sinTheta, cosTheta, 0, 0, 0, 0, 1, 0,
              0, 0, 0, 1);
  return Transform(m, Transpose(m));
}

Transform Rotate(Float theta, const vec3f &axis) {
  vec3f a = unit_vector(axis);
  Float sinTheta = std::sin(Radians(theta));
  Float cosTheta = std::cos(Radians(theta));
  Matrix4x4 m;
  // Compute rotation of first basis vector
  m.m[0][0] = a.x() * a.x() + (1 - a.x() * a.x()) * cosTheta;
  m.m[0][1] = a.x() * a.y() * (1 - cosTheta) - a.z() * sinTheta;
  m.m[0][2] = a.x() * a.z() * (1 - cosTheta) + a.y() * sinTheta;
  m.m[0][3] = 0;
  
  // Compute rotations of second and third basis vectors
  m.m[1][0] = a.x() * a.y() * (1 - cosTheta) + a.z() * sinTheta;
  m.m[1][1] = a.y() * a.y() + (1 - a.y() * a.y()) * cosTheta;
  m.m[1][2] = a.y() * a.z() * (1 - cosTheta) - a.x() * sinTheta;
  m.m[1][3] = 0;
  
  m.m[2][0] = a.x() * a.z() * (1 - cosTheta) - a.y() * sinTheta;
  m.m[2][1] = a.y() * a.z() * (1 - cosTheta) + a.x() * sinTheta;
  m.m[2][2] = a.z() * a.z() + (1 - a.z() * a.z()) * cosTheta;
  m.m[2][3] = 0;
  return Transform(m, Transpose(m));
}

Transform LookAt(const point3f &pos, const point3f &look, const vec3f &up) {
  Matrix4x4 cameraToWorld;
  // Initialize fourth column of viewing matrix
  cameraToWorld.m[0][3] = pos.x();
  cameraToWorld.m[1][3] = pos.y();
  cameraToWorld.m[2][3] = pos.z();
  cameraToWorld.m[3][3] = 1;
  
  // Initialize first three columns of viewing matrix
  vec3f dir = unit_vector(look - pos);
  if (cross(unit_vector(up), dir).length() == 0) {
    throw std::runtime_error("\"up\" vector and viewing direction passed to LookAt are pointing in the same direction.  Using the identity transformation.");
    return Transform();
  }
  vec3f right = unit_vector(cross(unit_vector(up), dir));
  vec3f newUp = cross(dir, right);
  cameraToWorld.m[0][0] = right.x();
  cameraToWorld.m[1][0] = right.y();
  cameraToWorld.m[2][0] = right.z();
  cameraToWorld.m[3][0] = 0.;
  cameraToWorld.m[0][1] = newUp.x();
  cameraToWorld.m[1][1] = newUp.y();
  cameraToWorld.m[2][1] = newUp.z();
  cameraToWorld.m[3][1] = 0.;
  cameraToWorld.m[0][2] = dir.x();
  cameraToWorld.m[1][2] = dir.y();
  cameraToWorld.m[2][2] = dir.z();
  cameraToWorld.m[3][2] = 0.;
  return Transform(Inverse(cameraToWorld), cameraToWorld);
}

aabb Transform::operator()(const aabb &b) const {
  const Transform &M = *this;
  vec3f pMin = b.min();
  vec3f pMax = b.max();
  
  aabb ret(M(pMin));
  ret = surrounding_box(ret, M(point3f(pMax.x(), pMin.y(), pMin.z())));
  ret = surrounding_box(ret, M(point3f(pMin.x(), pMax.y(), pMin.z())));
  ret = surrounding_box(ret, M(point3f(pMin.x(), pMin.y(), pMax.z())));
  ret = surrounding_box(ret, M(point3f(pMin.x(), pMax.y(), pMax.z())));
  ret = surrounding_box(ret, M(point3f(pMax.x(), pMax.y(), pMin.z())));
  ret = surrounding_box(ret, M(point3f(pMax.x(), pMin.y(), pMax.z())));
  ret = surrounding_box(ret, M(point3f(pMax.x(), pMax.y(), pMax.z())));
  return ret;
}

Transform Transform::operator*(const Transform &t2) const {
  return Transform(Matrix4x4::Mul(m, t2.m), Matrix4x4::Mul(t2.mInv, mInv));
}

bool Transform::SwapsHandedness() const {
  Float det = m.m[0][0] * (m.m[1][1] * m.m[2][2] - m.m[1][2] * m.m[2][1]) -
    m.m[0][1] * (m.m[1][0] * m.m[2][2] - m.m[1][2] * m.m[2][0]) +
    m.m[0][2] * (m.m[1][0] * m.m[2][1] - m.m[1][1] * m.m[2][0]);
  return det < 0;
}


Transform Orthographic(Float zNear, Float zFar) {
  return Scale(1, 1, 1 / (zFar - zNear)) * Translate(vec3f(0, 0, -zNear));
}

Transform Perspective(Float fov, Float n, Float f) {
  // Perform projective divide for perspective projection
  Matrix4x4 persp(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, f / (f - n), -f * n / (f - n),
                  0, 0, 1, 0);
  
  // Scale canonical perspective view to specified field of view
  Float invTanAng = 1 / std::tan(Radians(fov) / 2);
  return Scale(invTanAng, invTanAng, 1) * Transform(persp);
}