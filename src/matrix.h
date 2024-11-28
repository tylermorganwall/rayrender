#ifndef MATRIXH
#define MATRIXH

#include "vectypes.h"
#include "simd.h"
#include <cstring>

struct alignas(16) Matrix4x4 {
  Matrix4x4() {
    m[0][0] = m[1][1] = m[2][2] = m[3][3] = 1.f;
    m[0][1] = m[0][2] = m[0][3] = m[1][0] =
    m[1][2] = m[1][3] = m[2][0] = m[2][1] = m[2][3] =
    m[3][0] = m[3][1] = m[3][2] = 0.f;
  }
  Matrix4x4(Float mat[4][4]);
#ifdef RAYSIMDVEC
  Matrix4x4(FVec4 mat[4]) {
    m[0] = mat[0];
    m[1] = mat[1];
    m[2] = mat[2];
    m[3] = mat[3];
  }
#endif
  Matrix4x4(Float t00, Float t01, Float t02, Float t03,
            Float t10, Float t11, Float t12, Float t13,
            Float t20, Float t21, Float t22, Float t23,
            Float t30, Float t31, Float t32, Float t33);
  bool operator==(const Matrix4x4 &m2) const {
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        if (m[i][j] != m2.m[i][j]) {
          return false;
        }
      }
    }
    return true;
  }
  bool operator!=(const Matrix4x4 &m2) const {
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        if (m[i][j] != m2.m[i][j]) {
          return true;
        }
      }
    }
    return false;
  }
  friend Float Det(const Matrix4x4 &m1) {
    return(
        m1.m[0][3] * m1.m[1][2] * m1.m[2][1] * m1.m[3][0] - m1.m[0][2] * m1.m[1][3] * m1.m[2][1] * m1.m[3][0] -
        m1.m[0][3] * m1.m[1][1] * m1.m[2][2] * m1.m[3][0] + m1.m[0][1] * m1.m[1][3] * m1.m[2][2] * m1.m[3][0] +
        m1.m[0][2] * m1.m[1][1] * m1.m[2][3] * m1.m[3][0] - m1.m[0][1] * m1.m[1][2] * m1.m[2][3] * m1.m[3][0] -
        m1.m[0][3] * m1.m[1][2] * m1.m[2][0] * m1.m[3][1] + m1.m[0][2] * m1.m[1][3] * m1.m[2][0] * m1.m[3][1] +
        m1.m[0][3] * m1.m[1][0] * m1.m[2][2] * m1.m[3][1] - m1.m[0][0] * m1.m[1][3] * m1.m[2][2] * m1.m[3][1] -
        m1.m[0][2] * m1.m[1][0] * m1.m[2][3] * m1.m[3][1] + m1.m[0][0] * m1.m[1][2] * m1.m[2][3] * m1.m[3][1] +
        m1.m[0][3] * m1.m[1][1] * m1.m[2][0] * m1.m[3][2] - m1.m[0][1] * m1.m[1][3] * m1.m[2][0] * m1.m[3][2] -
        m1.m[0][3] * m1.m[1][0] * m1.m[2][1] * m1.m[3][2] + m1.m[0][0] * m1.m[1][3] * m1.m[2][1] * m1.m[3][2] +
        m1.m[0][1] * m1.m[1][0] * m1.m[2][3] * m1.m[3][2] - m1.m[0][0] * m1.m[1][1] * m1.m[2][3] * m1.m[3][2] -
        m1.m[0][2] * m1.m[1][1] * m1.m[2][0] * m1.m[3][3] + m1.m[0][1] * m1.m[1][2] * m1.m[2][0] * m1.m[3][3] +
        m1.m[0][2] * m1.m[1][0] * m1.m[2][1] * m1.m[3][3] - m1.m[0][0] * m1.m[1][2] * m1.m[2][1] * m1.m[3][3] -
        m1.m[0][1] * m1.m[1][0] * m1.m[2][2] * m1.m[3][3] + m1.m[0][0] * m1.m[1][1] * m1.m[2][2] * m1.m[3][3]
    );
  }
  friend std::ostream& operator<<(std::ostream& o, Matrix4x4 const& m1)  {
    return (o << m1.m[0][0] << " " << m1.m[0][1] << " " << m1.m[0][2] << " " << m1.m[0][3] << "\n" <<
      m1.m[1][0] << " " << m1.m[1][1] << " " << m1.m[1][2] << " " << m1.m[1][3] << "\n" <<
      m1.m[2][0] << " " << m1.m[2][1] << " " << m1.m[2][2] << " " << m1.m[2][3] << "\n" << 
      m1.m[3][0] << " " << m1.m[3][1] << " " << m1.m[3][2] << " " << m1.m[3][3] << "\n");
  }
  friend Matrix4x4 Transpose(const Matrix4x4 &);

  static Matrix4x4 Mul(const Matrix4x4 &m1, const Matrix4x4 &m2) {
    Matrix4x4 r;
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
#ifdef RAYSIMDVEC
        FVec4 col{m2.m[0][j], m2.m[1][j], m2.m[2][j], m2.m[3][j]};
        r.m[i][j] = simd_dot(m1.m[i], col);
#else
        r.m[i][j] = m1.m[i][0] * m2.m[0][j] + 
                    m1.m[i][1] * m2.m[1][j] + 
                    m1.m[i][2] * m2.m[2][j] + 
                    m1.m[i][3] * m2.m[3][j];
#endif
      }
    }
    return r;
  }
  
  friend Matrix4x4 Inverse(const Matrix4x4 &);
  
  // Float m[4][4];
#ifdef RAYSIMDVEC
  FVec4 m[4];
#else
  Float m[4][4];
#endif
};

#endif