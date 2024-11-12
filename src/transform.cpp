#include "transform.h"

#include "mathinline.h"
#include "aabb.h"
#include "hitable.h"

Transform::Transform(const Float mat[4][4]) {
  m = Matrix4x4(mat[0][0], mat[0][1], mat[0][2], mat[0][3],
                mat[1][0], mat[1][1], mat[1][2], mat[1][3],
                mat[2][0], mat[2][1], mat[2][2], mat[2][3],
                mat[3][0], mat[3][1], mat[3][2], mat[3][3]);
  mInv = Inverse(m);
}

Transform::Transform(const Rcpp::NumericMatrix& mat) {
  m = Matrix4x4(mat(0,0), mat(0,1), mat(0,2), mat(0,3),
                mat(1,0), mat(1,1), mat(1,2), mat(1,3),
                mat(2,0), mat(2,1), mat(2,2), mat(2,3),
                mat(3,0), mat(3,1), mat(3,2), mat(3,3));
  mInv = Inverse(m);
}

bool Transform::operator==(const Transform &t) const {
  return t.m == m && t.mInv == mInv;
}
bool Transform::operator!=(const Transform &t) const {
  return t.m != m || t.mInv != mInv;
}
bool Transform::operator<(const Transform &t2) const {
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      if (m.m[i][j] > t2.m.m[i][j]) return false;
    }
  }
  return true;
}

Transform Translate(const vec3f &delta) {
  Matrix4x4 m(1, 0, 0, delta.x(), 
              0, 1, 0, delta.y(), 
              0, 0, 1, delta.z(), 
              0, 0, 0, 1);
  Matrix4x4 minv(1, 0, 0, -delta.x(), 
                 0, 1, 0, -delta.y(), 
                 0, 0, 1, -delta.z(), 
                 0, 0, 0, 1);
  return Transform(m, minv);
}

Transform Scale(Float x, Float y, Float z) {
  Matrix4x4 m(x, 0, 0, 0, 
              0, y, 0, 0, 
              0, 0, z, 0, 
              0, 0, 0, 1);
  Matrix4x4 minv(1 / x, 0, 0, 0, 
                 0, 1 / y, 0, 0, 
                 0, 0, 1 / z, 0, 
                 0, 0, 0,     1);
  return Transform(m, minv);
}

Transform RotateX(Float theta) {
  Float sinTheta = std::sin(Radians(theta));
  Float cosTheta = std::cos(Radians(theta));
  Matrix4x4 m(1,        0,         0, 0, 
              0, cosTheta, -sinTheta, 0, 
              0, sinTheta,  cosTheta, 0,
              0,        0,         0, 1);
  return Transform(m, Transpose(m));
}

Transform RotateY(Float theta) {
  Float sinTheta = std::sin(Radians(theta));
  Float cosTheta = std::cos(Radians(theta));
  Matrix4x4 m(cosTheta,  0, sinTheta, 0, 
              0,         1,        0, 0, 
              -sinTheta, 0, cosTheta, 0,
              0,         0,        0, 1);
  return Transform(m, Transpose(m));
}

Transform RotateZ(Float theta) {
  Float sinTheta = std::sin(Radians(theta));
  Float cosTheta = std::cos(Radians(theta));
  Matrix4x4 m(cosTheta, -sinTheta, 0, 0, 
              sinTheta,  cosTheta, 0, 0, 
                     0,         0, 1, 0,
                     0,         0, 0, 1);
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


Transform Orthographic(Float zNear, Float zFar) {
  return Scale(1, 1, 1 / (zFar - zNear)) * Translate(vec3f(0, 0, -zNear));
}

Transform Perspective(Float fov, Float n, Float f) {
  // Perform projective divide for perspective projection
  Matrix4x4 persp(1, 0, 0, 0, 
                  0, 1, 0, 0, 
                  0, 0, f / (f - n), -f * n / (f - n),
                  0, 0, 1, 0);
  
  // Scale canonical perspective view to specified field of view
  Float invTanAng = 1 / std::tan(Radians(fov) / 2);
  return Scale(invTanAng, invTanAng, 1) * Transform(persp);
}


bool Transform::SwapsHandedness() const {
  Float det = m.m[0][0] * (m.m[1][1] * m.m[2][2] - m.m[1][2] * m.m[2][1]) -
    m.m[0][1] * (m.m[1][0] * m.m[2][2] - m.m[1][2] * m.m[2][0]) +
    m.m[0][2] * (m.m[1][0] * m.m[2][1] - m.m[1][1] * m.m[2][0]);
  return det < 0;
}


bool Transform::IsIdentity() const {
  return (m.m[0][0] == 1.f && m.m[0][1] == 0.f &&
          m.m[0][2] == 0.f && m.m[0][3] == 0.f &&
          m.m[1][0] == 0.f && m.m[1][1] == 1.f &&
          m.m[1][2] == 0.f && m.m[1][3] == 0.f &&
          m.m[2][0] == 0.f && m.m[2][1] == 0.f &&
          m.m[2][2] == 1.f && m.m[2][3] == 0.f &&
          m.m[3][0] == 0.f && m.m[3][1] == 0.f &&
          m.m[3][2] == 0.f && m.m[3][3] == 1.f);
}
const Matrix4x4& Transform::GetMatrix() const { return m; }
const Matrix4x4& Transform::GetInverseMatrix() const { return mInv; }
bool Transform::HasScale() const {
  Float la2 = (*this)(vec3f(1, 0, 0)).squared_length();
  Float lb2 = (*this)(vec3f(0, 1, 0)).squared_length();
  Float lc2 = (*this)(vec3f(0, 0, 1)).squared_length();
#define NOT_ONE(x) ((x) < .999f || (x) > 1.001f)
  return (NOT_ONE(la2) || NOT_ONE(lb2) || NOT_ONE(lc2));
#undef NOT_ONE
}

// Transform Inline Functions

aabb Transform::operator()(const aabb &b) const {
  const Transform &M = *this;
  point3f pMin = b.min();
  point3f pMax = b.max();
  
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


// ray Transform::operator()(const ray &r) const {
//   point3f o = (*this)(r.origin());
//   vec3f d = (*this)(r.direction());
//   Float tMax = r.tMax;
//   return ray(o, d, r.pri_stack, r.time(), tMax);
// }


ray Transform::operator()(const ray &r) const {
  vec3f oError;
  point3f o = (*this)(r.origin(), &oError);
  vec3f d = (*this)(r.direction());
  // Offset ray origin to edge of error bounds and compute _tMax_
  Float lengthSquared = d.squared_length();
  Float tMax = r.tMax;
  if (lengthSquared > 0) {
    Float dt = dot(Abs(d), oError) / lengthSquared;
    o += d * dt;
    tMax -= dt;
  }
  return ray(o, d, r.pri_stack, r.time(), tMax);
}


 ray Transform::operator()(const ray &r, vec3f *oError,
                               vec3f *dError) const {
  point3f o = (*this)(r.origin(), oError);
  vec3f d = (*this)(r.direction(), dError);
  Float tMax = r.tMax;
  Float lengthSquared = d.squared_length();
  if (lengthSquared > 0) {
    Float dt = dot(Abs(d), *oError) / lengthSquared;
    o += d * dt;
    //        tMax -= dt;
  }
  return ray(o, d, r.pri_stack, r.time(), tMax);
}

 ray Transform::operator()(const ray &r, const vec3f &oErrorIn,
                               const vec3f &dErrorIn, vec3f *oErrorOut,
                               vec3f *dErrorOut) const {
  point3f o = (*this)(r.origin(), oErrorIn, oErrorOut);
  vec3f d = (*this)(r.direction(), dErrorIn, dErrorOut);
  Float tMax = r.tMax;
  Float lengthSquared = d.squared_length();
  if (lengthSquared > 0) {
    Float dt = dot(Abs(d), *oErrorOut) / lengthSquared;
    o += d * dt;
    //        tMax -= dt;
  }
  return ray(o, d, r.pri_stack, r.time(), tMax);
}


hit_record Transform::operator()(const hit_record &r) const {
  hit_record hr;
  hr.p = (*this)(r.p, r.pError, &hr.pError);
  hr.normal = (*this)(r.normal);
  hr.bump_normal = (*this)(r.bump_normal);
  hr.dpdu = (*this)(r.dpdu);
  hr.dpdv = (*this)(r.dpdv);
  hr.mat_ptr = r.mat_ptr;
  hr.has_bump = r.has_bump;
#ifdef DEBUGBVH
  hr.bvh_nodes = r.bvh_nodes;
#endif
  hr.u = r.u;
  hr.v = r.v;
  hr.t = r.t;
  hr.shape = r.shape;
  hr.alpha_miss = r.alpha_miss;
  
  
  
  //Need to transform wo if used
  return(hr);
}

hit_record Transform::operator()(hit_record &r) const {
  hit_record hr;
  hr.p = (*this)(r.p, r.pError, &hr.pError);

  hr.normal = (*this)(r.normal);
  hr.bump_normal = (*this)(r.bump_normal);
  hr.dpdu = (*this)(r.dpdu);
  hr.dpdv = (*this)(r.dpdv);
  hr.mat_ptr = r.mat_ptr;
  hr.has_bump = r.has_bump;
#ifdef DEBUGBVH
  hr.bvh_nodes = r.bvh_nodes;
#endif
  hr.u = r.u;
  hr.v = r.v;
  hr.t = r.t;
  hr.shape = r.shape;
  hr.alpha_miss = r.alpha_miss;
  
  
  //Need to transform wo if used
  return(hr);
}

// inline rayDifferential Transform::operator()(const rayDifferential &r) const {
//   ray tr = (*this)(ray(r));
//   rayDifferential ret(tr.o, tr.d, tr.tMax, tr.time, tr.medium);
//   ret.hasDifferentials = r.hasDifferentials;
//   ret.rxOrigin = (*this)(r.rxOrigin);
//   ret.ryOrigin = (*this)(r.ryOrigin);
//   ret.rxDirection = (*this)(r.rxDirection);
//   ret.ryDirection = (*this)(r.ryDirection);
//   return ret;
// }

template <typename T>
point3<T> Transform::operator()(const point3<T> &p,
                                vec3<T> *pError) const {
  T x = p.x(), y = p.y(), z = p.z();
  // Compute transformed coordinates from point _pt_
  T xp = (m.m[0][0] * x + m.m[0][1] * y) + (m.m[0][2] * z + m.m[0][3]);
  T yp = (m.m[1][0] * x + m.m[1][1] * y) + (m.m[1][2] * z + m.m[1][3]);
  T zp = (m.m[2][0] * x + m.m[2][1] * y) + (m.m[2][2] * z + m.m[2][3]);
  T wp = (m.m[3][0] * x + m.m[3][1] * y) + (m.m[3][2] * z + m.m[3][3]);

  // Compute absolute error for transformed point
  T xAbsSum = (std::fabs(m.m[0][0] * x) + std::fabs(m.m[0][1] * y) +
    std::fabs(m.m[0][2] * z) + std::fabs(m.m[0][3]));
  T yAbsSum = (std::fabs(m.m[1][0] * x) + std::fabs(m.m[1][1] * y) +
    std::fabs(m.m[1][2] * z) + std::fabs(m.m[1][3]));
  T zAbsSum = (std::fabs(m.m[2][0] * x) + std::fabs(m.m[2][1] * y) +
    std::fabs(m.m[2][2] * z) + std::fabs(m.m[2][3]));
  *pError = gamma(3) * vec3<T>(xAbsSum, yAbsSum, zAbsSum);
  // CHECK_NE(wp, 0);
  if (wp == 1) {
    return point3<T>(xp, yp, zp);
  } else {
    return point3<T>(xp, yp, zp) / wp;
  }
}

template <typename T>
point3<T> Transform::operator()(const point3<T> &pt,
                                const vec3<T> &ptError,
                                vec3<T> *absError) const {
  T x = pt.x(), y = pt.y(), z = pt.z();
  T xp = (m.m[0][0] * x + m.m[0][1] * y) + (m.m[0][2] * z + m.m[0][3]);
  T yp = (m.m[1][0] * x + m.m[1][1] * y) + (m.m[1][2] * z + m.m[1][3]);
  T zp = (m.m[2][0] * x + m.m[2][1] * y) + (m.m[2][2] * z + m.m[2][3]);
  T wp = (m.m[3][0] * x + m.m[3][1] * y) + (m.m[3][2] * z + m.m[3][3]);
  absError->e[0] =
    (gamma(3) + (T)1) *
    (std::fabs(m.m[0][0]) * ptError.x()  + 
     std::fabs(m.m[0][1]) * ptError.y()  +
     std::fabs(m.m[0][2]) * ptError.z()) +
    gamma(3) * (std::fabs(m.m[0][0] * x) + 
                std::fabs(m.m[0][1] * y) +
                std::fabs(m.m[0][2] * z) + std::fabs(m.m[0][3]));
  absError->e[1] =
    (gamma(3) + (T)1) *
    (std::fabs(m.m[1][0]) * ptError.x()  + 
     std::fabs(m.m[1][1]) * ptError.y()  +
     std::fabs(m.m[1][2]) * ptError.z()) +
    gamma(3) * (std::fabs(m.m[1][0] * x) + 
                std::fabs(m.m[1][1] * y) +
                std::fabs(m.m[1][2] * z) + std::fabs(m.m[1][3]));
  absError->e[2] =
    (gamma(3) + (T)1) *
    (std::fabs(m.m[2][0]) * ptError.x()  + 
     std::fabs(m.m[2][1]) * ptError.y()  +
     std::fabs(m.m[2][2]) * ptError.z()) +
    gamma(3) * (std::fabs(m.m[2][0] * x) + 
                std::fabs(m.m[2][1] * y) +
                std::fabs(m.m[2][2] * z) + std::fabs(m.m[2][3]));
  // CHECK_NE(wp, 0);
  if (wp == 1.) {
    return point3<T>(xp, yp, zp);
  } else {
    return point3<T>(xp, yp, zp) / wp;
  }
}

template <typename T>
vec3<T> Transform::operator()(const vec3<T> &v,
                              vec3<T> *absError) const {
  T x = v.x(), y = v.y(), z = v.z();
  absError->e[0] =
    gamma(3) * (std::fabs(m.m[0][0] * v.x()) + 
                std::fabs(m.m[0][1] * v.y()) +
                std::fabs(m.m[0][2] * v.z()));
  absError->e[1] =
    gamma(3) * (std::fabs(m.m[1][0] * v.x()) + 
                std::fabs(m.m[1][1] * v.y()) +
                std::fabs(m.m[1][2] * v.z()));
  absError->e[2] =
    gamma(3) * (std::fabs(m.m[2][0] * v.x()) + 
                std::fabs(m.m[2][1] * v.y()) +
                std::fabs(m.m[2][2] * v.z()));
  return vec3<T>(m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z,
                 m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z,
                 m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z);
}

template <typename T>
vec3<T> Transform::operator()(const vec3<T> &v,
                              const vec3<T> &vError,
                              vec3<T> *absError) const {
  T x = v.x(), y = v.y(), z = v.z();
  absError->e[0] =
    (gamma(3) + (T)1) *
    (std::fabs(m.m[0][0]) * vError.x() + 
     std::fabs(m.m[0][1]) * vError.y() +
     std::fabs(m.m[0][2]) * vError.z()) +
    gamma(3) * (std::fabs(m.m[0][0] * v.x()) + 
                std::fabs(m.m[0][1] * v.y()) +
                std::fabs(m.m[0][2] * v.z()));
  absError->e[1] =
    (gamma(3) + (T)1) *
    (std::fabs(m.m[1][0]) * vError.x() + 
     std::fabs(m.m[1][1]) * vError.y() +
     std::fabs(m.m[1][2]) * vError.z()) +
    gamma(3) * (std::fabs(m.m[1][0] * v.x()) + 
                std::fabs(m.m[1][1] * v.y()) +
                std::fabs(m.m[1][2] * v.z()));
  absError->e[2] =
    (gamma(3) + (T)1) *
    (std::fabs(m.m[2][0]) * vError.x() + 
     std::fabs(m.m[2][1]) * vError.y() +
     std::fabs(m.m[2][2]) * vError.z()) +
    gamma(3) * (std::fabs(m.m[2][0] * v.x()) + 
                std::fabs(m.m[2][1] * v.y()) +
                std::fabs(m.m[2][2] * v.z()));
  return vec3<T>(m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z,
                 m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z,
                 m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z);
}

vec3f Transform::w() {
  return(vec3f(m.m[0][2],
               m.m[1][2],
               m.m[2][2]));
}
vec3f Transform::u() {
  return(vec3f(m.m[0][0],
               m.m[1][0],
               m.m[2][0]));
}
vec3f Transform::v() {
  return(vec3f(m.m[0][1],
               m.m[1][1],
               m.m[2][1]));
}

// Unit Tests
#include <testthat.h>

context("Translate function works correctly") {
  test_that("Translate applies correct translation") {
      vec3f delta(1.0f, 2.0f, 3.0f);
      Transform translate = Translate(delta);

      // Expected translation matrix
      Matrix4x4 expected(
          1, 0, 0, 1.0f,
          0, 1, 0, 2.0f,
          0, 0, 1, 3.0f,
          0, 0, 0, 1);

      expect_true(translate.GetMatrix() == expected);

      // Test inverse
      Matrix4x4 expectedInv(
          1, 0, 0, -1.0f,
          0, 1, 0, -2.0f,
          0, 0, 1, -3.0f,
          0, 0, 0, 1);

      expect_true(translate.GetInverseMatrix() == expectedInv);

      // Test transforming a point
      point3f p(0.0f, 0.0f, 0.0f);
      point3f pTransformed = translate(p);
      point3f expectedPoint(1.0f, 2.0f, 3.0f);

      // Compare individual components
      expect_true(pTransformed.x() == Approx(expectedPoint.x()));
      expect_true(pTransformed.y() == Approx(expectedPoint.y()));
      expect_true(pTransformed.z() == Approx(expectedPoint.z()));
  }
};

// Test Scale function
context("Scale function works correctly") {
  test_that("Scale applies correct scaling") {
      Float sx = 2.0f, sy = 3.0f, sz = 4.0f;
      Transform scale = Scale(sx, sy, sz);

      // Expected scaling matrix
      Matrix4x4 expected(
          2.0f, 0,    0,    0,
          0,    3.0f, 0,    0,
          0,    0,    4.0f, 0,
          0,    0,    0,    1);

      expect_true(scale.GetMatrix() == expected);

      // Test inverse
      Matrix4x4 expectedInv(
          1.0f / 2.0f, 0,            0,            0,
          0,           1.0f / 3.0f,  0,            0,
          0,           0,            1.0f / 4.0f,  0,
          0,           0,            0,            1);

      expect_true(scale.GetInverseMatrix() == expectedInv);

      // Test transforming a point
      point3f p(1.0f, 1.0f, 1.0f);
      point3f pTransformed = scale(p);
      point3f expectedPoint(2.0f, 3.0f, 4.0f);

      expect_true(pTransformed.x() == Approx(expectedPoint.x()));
      expect_true(pTransformed.y() == Approx(expectedPoint.y()));
      expect_true(pTransformed.z() == Approx(expectedPoint.z()));
  }
}

// Test RotateX function
context("RotateX function works correctly") {
  test_that("RotateX rotates point around X-axis") {
      Float theta = 90.0f; // degrees
      Transform rotateX = RotateX(theta);

      // Expected rotation matrix (90 degrees around X-axis)
      Float rad = Radians(theta);
      Float cosTheta = std::cos(rad);
      Float sinTheta = std::sin(rad);

      Matrix4x4 expected(
          1, 0,         0,        0,
          0, cosTheta, -sinTheta, 0,
          0, sinTheta,  cosTheta, 0,
          0, 0,         0,        1);

      expect_true(rotateX.GetMatrix() == expected);

      // Test transforming a point
      point3f p(0.0f, 1.0f, 0.0f);
      point3f pTransformed = rotateX(p);
      point3f expectedPoint(0.0f, 0.0f, 1.0f);

      expect_true(pTransformed.x() == Approx(expectedPoint.x()));
      expect_true(pTransformed.y() == Approx(expectedPoint.y()));
      expect_true(pTransformed.z() == Approx(expectedPoint.z()));
  }
}

// Test RotateY function
context("RotateY function works correctly") {
  test_that("RotateY rotates point around Y-axis") {
      Float theta = 90.0f; // degrees
      Transform rotateY = RotateY(theta);

      // Expected rotation matrix (90 degrees around Y-axis)
      Float rad = Radians(theta);
      Float cosTheta = std::cos(rad);
      Float sinTheta = std::sin(rad);

      Matrix4x4 expected(
          cosTheta,  0, sinTheta, 0,
          0,         1, 0,        0,
          -sinTheta, 0, cosTheta, 0,
          0,         0, 0,        1);

      expect_true(rotateY.GetMatrix() == expected);

      // Test transforming a point
      point3f p(1.0f, 0.0f, 0.0f);
      point3f pTransformed = rotateY(p);
      point3f expectedPoint(0.0f, 0.0f, -1.0f);

      expect_true(pTransformed.x() == Approx(expectedPoint.x()));
      expect_true(pTransformed.y() == Approx(expectedPoint.y()));
      expect_true(pTransformed.z() == Approx(expectedPoint.z()));
  }
}

// Test RotateZ function
context("RotateZ function works correctly") {
  test_that("RotateZ rotates point around Z-axis") {
      Float theta = 90.0f; // degrees
      Transform rotateZ = RotateZ(theta);

      // Expected rotation matrix (90 degrees around Z-axis)
      Float rad = Radians(theta);
      Float cosTheta = std::cos(rad);
      Float sinTheta = std::sin(rad);

      Matrix4x4 expected(
          cosTheta, -sinTheta, 0, 0,
          sinTheta,  cosTheta, 0, 0,
          0,         0,        1, 0,
          0,         0,        0, 1);

      expect_true(rotateZ.GetMatrix() == expected);

      // Test transforming a point
      point3f p(1.0f, 0.0f, 0.0f);
      point3f pTransformed = rotateZ(p);
      point3f expectedPoint(0.0f, 1.0f, 0.0f);

      expect_true(pTransformed.x() == Approx(expectedPoint.x()));
      expect_true(pTransformed.y() == Approx(expectedPoint.y()));
      expect_true(pTransformed.z() == Approx(expectedPoint.z()));
  }
}

// Test Rotate function with arbitrary axis
context("Rotate function works correctly with arbitrary axis") {
  test_that("Rotate rotates point around arbitrary axis") {
      Float theta = 90.0f; // degrees
      vec3f axis(0.0f, 0.0f, 1.0f); // Z-axis
      Transform rotate = Rotate(theta, axis);

      // Should be equivalent to RotateZ(90)
      Transform rotateZ = RotateZ(theta);

      // Compare matrices
      expect_true(rotate.GetMatrix() == rotateZ.GetMatrix());

      // Test transforming a point
      point3f p(1.0f, 0.0f, 0.0f);
      point3f pTransformed = rotate(p);
      point3f expectedPoint(0.0f, 1.0f, 0.0f);

      expect_true(pTransformed.x() == Approx(expectedPoint.x()));
      expect_true(pTransformed.y() == Approx(expectedPoint.y()));
      expect_true(pTransformed.z() == Approx(expectedPoint.z()));
  }
}

// Test LookAt function
context("LookAt function works correctly") {
  test_that("LookAt creates correct camera transform") {
      point3f pos(0.0f, 0.0f, 0.0f);
      point3f look(0.0f, 0.0f, -1.0f);
      vec3f up(0.0f, 1.0f, 0.0f);

      Transform lookAt = LookAt(pos, look, up);

      // Test transforming a point in front of the camera
      point3f p(0.0f, 0.0f, -5.0f);
      point3f pTransformed = lookAt(p);
      point3f expectedPoint(0.0f, 0.0f, -5.0f);

      expect_true(pTransformed.x() == Approx(expectedPoint.x()));
      expect_true(pTransformed.y() == Approx(expectedPoint.y()));
      expect_true(pTransformed.z() == Approx(expectedPoint.z()));
  }
}

// Test Orthographic function
context("Orthographic function works correctly") {
  test_that("Orthographic projection transforms point correctly") {
      Float zNear = 0.0f;
      Float zFar = 1.0f;

      Transform ortho = Orthographic(zNear, zFar);

      // Test transforming a point
      point3f p(1.0f, 2.0f, 3.0f);
      point3f pTransformed = ortho(p);

      // Since ortho scales Z by 1 / (zFar - zNear)
      Float scaleZ = 1 / (zFar - zNear);
      point3f expectedPoint(1.0f, 2.0f, 3.0f * scaleZ - zNear * scaleZ);

      expect_true(pTransformed.x() == Approx(expectedPoint.x()));
      expect_true(pTransformed.y() == Approx(expectedPoint.y()));
      expect_true(pTransformed.z() == Approx(expectedPoint.z()));
  }
}

// Test SwapsHandedness function
context("SwapsHandedness function works correctly") {
  test_that("SwapsHandedness detects handedness swap correctly") {
      Transform scale = Scale(-1.0f, 1.0f, 1.0f);
      expect_true(scale.SwapsHandedness() == true);

      Transform identity;
      expect_true(identity.SwapsHandedness() == false);
  }
}

// Test IsIdentity function
context("IsIdentity function works correctly") {
  test_that("IsIdentity detects identity transform correctly") {
      Transform identity;
      expect_true(identity.IsIdentity() == true);

      Transform translate = Translate(vec3f(1.0f, 2.0f, 3.0f));
      expect_true(translate.IsIdentity() == false);
  }
}

// Test HasScale function
context("HasScale function works correctly") {
  test_that("HasScale detects scaling correctly") {
      Transform scale = Scale(2.0f, 2.0f, 2.0f);
      expect_true(scale.HasScale() == true);

      Transform translate = Translate(vec3f(1.0f, 2.0f, 3.0f));
      expect_true(translate.HasScale() == false);
  }
}

// Test operator* (composition of transforms)
context("Transform composition works correctly") {
  test_that("Composition of transforms applies correctly") {
      Transform translate = Translate(vec3f(1.0f, 0.0f, 0.0f));
      Transform scale = Scale(2.0f, 2.0f, 2.0f);

      Transform composed = translate * scale;

      // Test transforming a point
      point3f p(1.0f, 1.0f, 1.0f);
      point3f pTransformed = composed(p);

      // Expected result: first scales to (2,2,2), then translates to (3,2,2)
      point3f expectedPoint(3.0f, 2.0f, 2.0f);

      expect_true(pTransformed.x() == Approx(expectedPoint.x()));
      expect_true(pTransformed.y() == Approx(expectedPoint.y()));
      expect_true(pTransformed.z() == Approx(expectedPoint.z()));
  }
}