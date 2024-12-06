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
  #else 
  normal3f operator()(const normal3f &n) const;
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