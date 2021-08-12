#ifndef TRANSFORMH
#define TRANSFORMH

#include "vec3.h"
#include "point3.h"
#include "matrix.h"
#include "normal.h"
#include "aabb.h"

class Transform {
public:
Transform() { }
  Transform(const Float mat[4][4]) {
    m = Matrix4x4(mat[0][0], mat[0][1], mat[0][2], mat[0][3],
                  mat[1][0], mat[1][1], mat[1][2], mat[1][3],
                  mat[2][0], mat[2][1], mat[2][2], mat[2][3],
                  mat[3][0], mat[3][1], mat[3][2], mat[3][3]);
    mInv = Inverse(m);
  }
  Transform(const Matrix4x4 &m) : m(m), mInv(Inverse(m)) { }
  Transform(const Matrix4x4 &m, const Matrix4x4 &mInv) 
    : m(m), mInv(mInv) {
  }
  friend Transform Inverse(const Transform &t) {
    return Transform(t.mInv, t.m);
  }
  friend Transform Transpose(const Transform &t) {
    return Transform(Transpose(t.m), Transpose(t.mInv));
  }
  bool operator==(const Transform &t) const {
    return t.m == m && t.mInv == mInv;
  }
  bool operator!=(const Transform &t) const {
    return t.m != m || t.mInv != mInv;
  }
  bool operator<(const Transform &t2) const {
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) {
        if (m.m[i][j] < t2.m.m[i][j]) return true;
        if (m.m[i][j] > t2.m.m[i][j]) return false;
      }
      return false;
  }
  bool IsIdentity() const {
    return (m.m[0][0] == 1.f && m.m[0][1] == 0.f &&
            m.m[0][2] == 0.f && m.m[0][3] == 0.f &&
            m.m[1][0] == 0.f && m.m[1][1] == 1.f &&
            m.m[1][2] == 0.f && m.m[1][3] == 0.f &&
            m.m[2][0] == 0.f && m.m[2][1] == 0.f &&
            m.m[2][2] == 1.f && m.m[2][3] == 0.f &&
            m.m[3][0] == 0.f && m.m[3][1] == 0.f &&
            m.m[3][2] == 0.f && m.m[3][3] == 1.f);
  }
  const Matrix4x4 &GetMatrix() const { return m; }
  const Matrix4x4 &GetInverseMatrix() const { return mInv; }
  bool HasScale() const {
    Float la2 = (*this)(vec3f(1, 0, 0)).squared_length();
    Float lb2 = (*this)(vec3f(0, 1, 0)).squared_length();
    Float lc2 = (*this)(vec3f(0, 0, 1)).squared_length();
#define NOT_ONE(x) ((x) < .999f || (x) > 1.001f)
    return (NOT_ONE(la2) || NOT_ONE(lb2) || NOT_ONE(lc2));
#undef NOT_ONE
  }
  template <typename T> inline point3<T> operator()(const point3<T> &p) const;
  template <typename T> inline vec3<T> operator()(const vec3<T> &v) const;
  template <typename T> inline normal3<T> operator()(const normal3<T> &) const;
  template <typename T> inline void operator()(const normal3<T> &, normal3<T> *nt) const;
  inline ray operator()(const ray &r) const;
  // inline rayDifferential operator()(const rayDifferential &r) const;
  aabb operator()(const aabb &b) const;
  Transform operator*(const Transform &t2) const;
  bool SwapsHandedness() const;
  // SurfaceInteraction operator()(const SurfaceInteraction &si) const;
  template <typename T> inline point3<T> operator()(const point3<T> &pt, vec3<T> *absError) const;
  template <typename T> inline point3<T>
  operator()(const point3<T> &p, const vec3<T> &pError,
           vec3<T> *pTransError) const;
  template <typename T> inline vec3<T>
  operator()(const vec3<T> &v, vec3<T> *vTransError) const;
  template <typename T> inline vec3<T>
  operator()(const vec3<T> &v, const vec3<T> &vError,
           vec3<T> *vTransError) const;
  inline ray operator()(const ray &r, vec3f *oError,
                      vec3f *dError) const;
  inline ray operator()(const ray &r, const vec3f &oErrorIn,
                      const vec3f &dErrorIn, vec3f *oErrorOut,
                      vec3f *dErrorOut) const;
  
private:
  Matrix4x4 m, mInv;
  friend class AnimatedTransform;
  friend struct Quaternion;
  
};


Transform Translate(const vec3f &delta);
Transform Scale(Float x, Float y, Float z);
Transform RotateX(Float theta);
Transform RotateY(Float theta);
Transform RotateZ(Float theta);
Transform Rotate(Float theta, const vec3f &axis);
Transform LookAt(const point3f &pos, const point3f &look, const vec3f &up);
Transform Orthographic(Float znear, Float zfar);
Transform Perspective(Float fov, Float znear, Float zfar);
// bool SolveLinearSystem2x2(const Float A[2][2], const Float B[2], Float *x0,
//                           Float *x1);

// Transform Inline Functions
template <typename T>
inline point3<T> Transform::operator()(const point3<T> &p) const {
  T x = p.x(), y = p.y(), z = p.z();
  T xp = m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z + m.m[0][3];
  T yp = m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z + m.m[1][3];
  T zp = m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z + m.m[2][3];
  T wp = m.m[3][0] * x + m.m[3][1] * y + m.m[3][2] * z + m.m[3][3];
  // CHECK_NE(wp, 0);
  if (wp == 1) {
    return point3<T>(xp, yp, zp);
  } else {
    return point3<T>(xp, yp, zp) / wp;
  }
}

template <typename T>
inline vec3<T> Transform::operator()(const vec3<T> &v) const {
  T x = v.x(), y = v.y(), z = v.z();
  return vec3<T>(m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z,
                 m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z,
                 m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z);
}

template <typename T>
inline normal3<T> Transform::operator()(const normal3<T> &n) const {
  T x = n.x(), y = n.y(), z = n.z();
  return normal3<T>(mInv.m[0][0] * x + mInv.m[1][0] * y + mInv.m[2][0] * z,
                    mInv.m[0][1] * x + mInv.m[1][1] * y + mInv.m[2][1] * z,
                    mInv.m[0][2] * x + mInv.m[1][2] * y + mInv.m[2][2] * z);
}

inline ray Transform::operator()(const ray &r) const {
  vec3f oError;
  vec3f o = (*this)(r.origin(), &oError);
  vec3f d = (*this)(r.direction());
  // Offset ray origin to edge of error bounds and compute _tMax_
  Float lengthSquared = d.squared_length();
  Float tMax = r.tMax;
  if (lengthSquared > 0) {
    Float dt = dot(Abs(d), oError) / lengthSquared;
    o += d * dt;
    tMax -= dt;
  }
  return ray(o, d, r.pri_stack, tMax, r.time());
}

// inline rayDifferential Transform::operator()(const rayDifferential &r) const {
//   ray tr = (*this)(ray(r));
//   RayDifferential ret(tr.o, tr.d, tr.tMax, tr.time, tr.medium);
//   ret.hasDifferentials = r.hasDifferentials;
//   ret.rxOrigin = (*this)(r.rxOrigin);
//   ret.ryOrigin = (*this)(r.ryOrigin);
//   ret.rxDirection = (*this)(r.rxDirection);
//   ret.ryDirection = (*this)(r.ryDirection);
//   return ret;
// }

template <typename T>
inline point3<T> Transform::operator()(const point3<T> &p,
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
  *pError = tgamma(3) * vec3<T>(xAbsSum, yAbsSum, zAbsSum);
  // CHECK_NE(wp, 0);
  if (wp == 1) {
    return point3<T>(xp, yp, zp);
  } else {
    return point3<T>(xp, yp, zp) / wp;
  }
}

template <typename T>
inline point3<T> Transform::operator()(const point3<T> &pt,
                                     const vec3<T> &ptError,
                                     vec3<T> *absError) const {
  T x = pt.x(), y = pt.y(), z = pt.z();
  T xp = (m.m[0][0] * x + m.m[0][1] * y) + (m.m[0][2] * z + m.m[0][3]);
  T yp = (m.m[1][0] * x + m.m[1][1] * y) + (m.m[1][2] * z + m.m[1][3]);
  T zp = (m.m[2][0] * x + m.m[2][1] * y) + (m.m[2][2] * z + m.m[2][3]);
  T wp = (m.m[3][0] * x + m.m[3][1] * y) + (m.m[3][2] * z + m.m[3][3]);
  absError->e[0] =
    (tgamma(3) + (T)1) *
    (std::fabs(m.m[0][0]) * ptError.x() + std::fabs(m.m[0][1]) * ptError.y() +
    std::fabs(m.m[0][2]) * ptError.z()) +
    tgamma(3) * (std::fabs(m.m[0][0] * x) + std::fabs(m.m[0][1] * y) +
    std::fabs(m.m[0][2] * z) + std::fabs(m.m[0][3]));
  absError->e[1] =
    (tgamma(3) + (T)1) *
    (std::fabs(m.m[1][0]) * ptError.x() + std::fabs(m.m[1][1]) * ptError.y() +
    std::fabs(m.m[1][2]) * ptError.z()) +
    tgamma(3) * (std::fabs(m.m[1][0] * x) + std::fabs(m.m[1][1] * y) +
    std::fabs(m.m[1][2] * z) + std::fabs(m.m[1][3]));
  absError->e[2] =
    (tgamma(3) + (T)1) *
    (std::fabs(m.m[2][0]) * ptError.x() + std::fabs(m.m[2][1]) * ptError.y() +
    std::fabs(m.m[2][2]) * ptError.z()) +
    tgamma(3) * (std::fabs(m.m[2][0] * x) + std::fabs(m.m[2][1] * y) +
    std::fabs(m.m[2][2] * z) + std::fabs(m.m[2][3]));
  // CHECK_NE(wp, 0);
  if (wp == 1.) {
    return point3<T>(xp, yp, zp);
  } else {
    return point3<T>(xp, yp, zp) / wp;
  }
}

template <typename T>
inline vec3<T> Transform::operator()(const vec3<T> &v,
                                      vec3<T> *absError) const {
  T x = v.x(), y = v.y(), z = v.z();
  absError->e[0] =
    tgamma(3) * (std::fabs(m.m[0][0] * v.x()) + std::fabs(m.m[0][1] * v.y()) +
    std::fabs(m.m[0][2] * v.z()));
  absError->e[1] =
    tgamma(3) * (std::fabs(m.m[1][0] * v.x()) + std::fabs(m.m[1][1] * v.y()) +
    std::fabs(m.m[1][2] * v.z()));
  absError->e[2] =
    tgamma(3) * (std::fabs(m.m[2][0] * v.x()) + std::fabs(m.m[2][1] * v.y()) +
    std::fabs(m.m[2][2] * v.z()));
  return vec3<T>(m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z,
                    m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z,
                    m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z);
}

template <typename T>
inline vec3<T> Transform::operator()(const vec3<T> &v,
                                      const vec3<T> &vError,
                                      vec3<T> *absError) const {
  T x = v.x(), y = v.y(), z = v.z();
  absError->e[0] =
    (tgamma(3) + (T)1) *
    (std::fabs(m.m[0][0]) * vError.x() + std::fabs(m.m[0][1]) * vError.y() +
    std::fabs(m.m[0][2]) * vError.z()) +
    tgamma(3) * (std::fabs(m.m[0][0] * v.x()) + std::fabs(m.m[0][1] * v.y()) +
    std::fabs(m.m[0][2] * v.z()));
  absError->e[1] =
    (tgamma(3) + (T)1) *
    (std::fabs(m.m[1][0]) * vError.x() + std::fabs(m.m[1][1]) * vError.y() +
    std::fabs(m.m[1][2]) * vError.z()) +
    tgamma(3) * (std::fabs(m.m[1][0] * v.x()) + std::fabs(m.m[1][1] * v.y()) +
    std::fabs(m.m[1][2] * v.z()));
  absError->e[2] =
    (tgamma(3) + (T)1) *
    (std::fabs(m.m[2][0]) * vError.x() + std::fabs(m.m[2][1]) * vError.y() +
    std::fabs(m.m[2][2]) * vError.z()) +
    tgamma(3) * (std::fabs(m.m[2][0] * v.x()) + std::fabs(m.m[2][1] * v.y()) +
    std::fabs(m.m[2][2] * v.z()));
  return vec3<T>(m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z,
                 m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z,
                 m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z);
}

inline ray Transform::operator()(const ray &r, vec3f *oError,
                               vec3f *dError) const {
  vec3f o = (*this)(r.origin(), oError);
  vec3f d = (*this)(r.direction(), dError);
  Float tMax = r.tMax;
  Float lengthSquared = d.squared_length();
  if (lengthSquared > 0) {
    Float dt = dot(Abs(d), *oError) / lengthSquared;
    o += d * dt;
    //        tMax -= dt;
  }
  return ray(o, d, r.pri_stack, tMax, r.time());
}

inline ray Transform::operator()(const ray &r, const vec3f &oErrorIn,
                               const vec3f &dErrorIn, vec3f *oErrorOut,
                               vec3f *dErrorOut) const {
  vec3f o = (*this)(r.origin(), oErrorIn, oErrorOut);
  vec3f d = (*this)(r.direction(), dErrorIn, dErrorOut);
  Float tMax = r.tMax;
  Float lengthSquared = d.squared_length();
  if (lengthSquared > 0) {
    Float dt = dot(Abs(d), *oErrorOut) / lengthSquared;
    o += d * dt;
    //        tMax -= dt;
  }
  return ray(o, d, r.pri_stack, tMax, r.time());
}

// AnimatedTransform Declarations
// class AnimatedTransform {
// public:
//   // AnimatedTransform Public Methods
//   AnimatedTransform(const Transform *startTransform, Float startTime,
//                     const Transform *endTransform, Float endTime);
//   static void Decompose(const Matrix4x4 &m, vec3f *T, Quaternion *R,
//                         Matrix4x4 *S);
//   void Interpolate(Float time, Transform *t) const;
//   ray operator()(const ray &r) const;
//   rayDifferential operator()(const rayDifferential &r) const;
//   point3f operator()(Float time, const point3f &p) const;
//   vec3f operator()(Float time, const vec3f &v) const;
//   bool HasScale() const {
//     return startTransform->HasScale() || endTransform->HasScale();
//   }
//   Bounds3f MotionBounds(const Bounds3f &b) const;
//   Bounds3f BoundPointMotion(const point3f &p) const;
//   
// private:
//   // AnimatedTransform Private Data
//   const Transform *startTransform, *endTransform;
//   const Float startTime, endTime;
//   const bool actuallyAnimated;
//   vec3f T[2];
//   Quaternion R[2];
//   Matrix4x4 S[2];
//   bool hasRotation;
//   struct DerivativeTerm {
//     DerivativeTerm() {}
//     DerivativeTerm(Float c, Float x, Float y, Float z)
//       : kc(c), kx(x), ky(y), kz(z) {}
//     Float kc, kx, ky, kz;
//     Float Eval(const point3f &p) const {
//       return kc + kx * p.x + ky * p.y + kz * p.z;
//     }
//   };
//   DerivativeTerm c1[3], c2[3], c3[3], c4[3], c5[3];
// };

#endif