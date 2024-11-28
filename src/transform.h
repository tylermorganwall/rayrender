#ifndef TRANSFORMH
#define TRANSFORMH

#include "matrix.h"
#include "aabb.h"
#include "point3.h"
#include "vec3.h"
#include "normal.h"

struct hit_record;

class Transform {
public:
Transform() { }
  Transform(const Float mat[4][4]);
  Transform(const Rcpp::NumericMatrix& mat);
  
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
  bool operator==(const Transform &t) const;
  bool operator!=(const Transform &t) const;
  bool operator<(const Transform &t2) const;
  bool IsIdentity() const;
  const Matrix4x4& GetMatrix() const;
  const Matrix4x4& GetInverseMatrix() const;
  bool HasScale() const;
  bool SwapsHandedness() const;
  vec3f w();
  vec3f u();
  vec3f v();
  
  //Transformations
  aabb operator()(const aabb &b) const;
  Transform operator*(const Transform &t2) const;
  
  template <typename T> point3<T> operator()(const point3<T> &p) const {
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
  
  template <typename T> vec3<T> operator()(const vec3<T> &v) const {
    T x = v.x(), y = v.y(), z = v.z();
    return vec3<T>(m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z,
                   m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z,
                   m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z);
  }
  
  #ifndef RAYSIMDVEC
  normal3f operator()(const normal3f &n) const {
    Float x = n.x(), y = n.y(), z = n.z();
    return normal3f(mInv.m[0][0] * x + mInv.m[1][0] * y + mInv.m[2][0] * z,
                    mInv.m[0][1] * x + mInv.m[1][1] * y + mInv.m[2][1] * z,
                    mInv.m[0][2] * x + mInv.m[1][2] * y + mInv.m[2][2] * z);
  }
  #endif

  #ifdef RAYSIMDVEC

  template <>
  point3<float> operator()(const point3<float> &p) const {
    // Load point into FVec4 (including w = 1.0f)
    FVec4 p_v = p.e; // [x, y, z, 1.0f]
    p_v[3] = 1.0;

    // Compute transformed point using SIMD dot products
    float xp = simd_dot(m.m[0], p_v);
    float yp = simd_dot(m.m[1], p_v);
    float zp = simd_dot(m.m[2], p_v);
    float wp = simd_dot(m.m[3], p_v);

    // Return the transformed point
    if (wp == 1.0f) {
        return point3<float>(xp, yp, zp);
    } else {
        float inv_wp = 1.0f / wp;
        return point3<float>(xp * inv_wp, yp * inv_wp, zp * inv_wp);
    }
  }

  template<>
  vec3<float> operator()(const vec3<float> &v) const {
    // Load vector into FVec4 (w component is 0.0f)
    FVec4 v_v = v.e; // [x, y, z, 0.0f]

    // Zero out the w components of the matrix rows
    FVec4 x_row = m.m[0]; x_row.xyzw[3] = 0.0f;
    FVec4 y_row = m.m[1]; y_row.xyzw[3] = 0.0f;
    FVec4 z_row = m.m[2]; z_row.xyzw[3] = 0.0f;

    // Compute transformed vector components
    float x = simd_dot(x_row, v_v);
    float y = simd_dot(y_row, v_v);
    float z = simd_dot(z_row, v_v);

    return vec3<float>(x, y, z);
  }

  
  normal3f operator()(const normal3f &n) const {
    // Load normal into FVec4 (w component is 0.0f)
    FVec4 n_v = n.e; // [x, y, z, 0.0f]

    // Transpose the inverse matrix (only the upper-left 3x3 part)
    FVec4 col0 = simd_set(mInv.m[0][0], mInv.m[1][0], mInv.m[2][0], 0.0f);
    FVec4 col1 = simd_set(mInv.m[0][1], mInv.m[1][1], mInv.m[2][1], 0.0f);
    FVec4 col2 = simd_set(mInv.m[0][2], mInv.m[1][2], mInv.m[2][2], 0.0f);

    // Compute transformed normal components
    float x = simd_dot(col0, n_v);
    float y = simd_dot(col1, n_v);
    float z = simd_dot(col2, n_v);

    return normal3f(x, y, z);
  }

  #endif
  
  // template <typename T> void operator()(const normal3<T> &, normal3<T> *nt) const;
  ray operator()(const ray &r) const;
  
  hit_record operator()(const hit_record &r) const;
  hit_record operator()(hit_record &r) const;
  
  // inline rayDifferential operator()(const rayDifferential &r) const;
  // SurfaceInteraction operator()(const SurfaceInteraction &si) const;
  template <typename T> point3<T> operator()(const point3<T> &pt, vec3<T> *absError) const;
  template <typename T> point3<T> operator()(const point3<T> &p, const vec3<T> &pError,
           vec3<T> *pTransError) const;
  template <typename T> vec3<T> operator()(const vec3<T> &v, vec3<T> *vTransError) const;
  template <typename T> vec3<T> operator()(const vec3<T> &v, const vec3<T> &vError,
           vec3<T> *vTransError) const;

  ray operator()(const ray &r, vec3f *oError,vec3f *dError) const;
  ray operator()(const ray &r, const vec3f &oErrorIn,
                 const vec3f &dErrorIn, vec3f *oErrorOut,
                 vec3f *dErrorOut) const;
  
  
  friend std::ostream& operator<<(std::ostream& o, Transform const& t) {
    return o << t.m;
  }
  
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

#endif