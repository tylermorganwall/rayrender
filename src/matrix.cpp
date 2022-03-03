#include "matrix.h"

Matrix4x4::Matrix4x4(Float mat[4][4]) { memcpy(m, mat, 16 * sizeof(Float)); }

Matrix4x4::Matrix4x4(Float t00, Float t01, Float t02, Float t03, Float t10,
                     Float t11, Float t12, Float t13, Float t20, Float t21,
                     Float t22, Float t23, Float t30, Float t31, Float t32,
                     Float t33) {
  m[0][0] = t00;
  m[0][1] = t01;
  m[0][2] = t02;
  m[0][3] = t03;
  m[1][0] = t10;
  m[1][1] = t11;
  m[1][2] = t12;
  m[1][3] = t13;
  m[2][0] = t20;
  m[2][1] = t21;
  m[2][2] = t22;
  m[2][3] = t23;
  m[3][0] = t30;
  m[3][1] = t31;
  m[3][2] = t32;
  m[3][3] = t33;
}

#include "RcppThread.h"

Matrix4x4 Inverse(const Matrix4x4 &m) {
  int indxc[4], indxr[4];
  int ipiv[4] = {0, 0, 0, 0};
  Float minv[4][4];
  memcpy(minv, m.m, 4 * 4 * sizeof(Float));
  for (int i = 0; i < 4; i++) {
    int irow = 0, icol = 0;
    Float big = 0.f;
    // Choose pivot
    for (int j = 0; j < 4; j++) {
      if (ipiv[j] != 1) {
        for (int k = 0; k < 4; k++) {
          if (ipiv[k] == 0) {
            if (std::fabs(minv[j][k]) >= big) {
              big = Float(std::fabs(minv[j][k]));
              irow = j;
              icol = k;
            }
          } else if (ipiv[k] > 1) {
            // RcppThread::Rcout << "Singular Matrix:\n" << m << "\n";
            throw std::runtime_error("Singular matrix in MatrixInvert");          
          }
        }
      }
    }
    ++ipiv[icol];
    // Swap rows _irow_ and _icol_ for pivot
    if (irow != icol) {
      for (int k = 0; k < 4; ++k) std::swap(minv[irow][k], minv[icol][k]);
    }
    indxr[i] = irow;
    indxc[i] = icol;
    if (minv[icol][icol] == 0.f) {
      // RcppThread::Rcout << "Singular Matrix:\n" << m << "\n";
      throw std::runtime_error("Singular matrix in MatrixInvert");
    }
    
    // Set $m[icol][icol]$ to one by scaling row _icol_ appropriately
    Float pivinv = 1. / minv[icol][icol];
    minv[icol][icol] = 1.;
    for (int j = 0; j < 4; j++) {
      minv[icol][j] *= pivinv;
    }
    
    // Subtract this row from others to zero out their columns
    for (int j = 0; j < 4; j++) {
      if (j != icol) {
        Float save = minv[j][icol];
        minv[j][icol] = 0;
        for (int k = 0; k < 4; k++) minv[j][k] -= minv[icol][k] * save;
      }
    }
  }
  // Swap columns to reflect permutation
  for (int j = 3; j >= 0; j--) {
    if (indxr[j] != indxc[j]) {
      for (int k = 0; k < 4; k++) {
        std::swap(minv[k][indxr[j]], minv[k][indxc[j]]);
      }
    }
  }
  return Matrix4x4(minv);
}

Matrix4x4 Transpose(const Matrix4x4 &m) {
  return Matrix4x4(m.m[0][0], m.m[1][0], m.m[2][0], m.m[3][0], 
                   m.m[0][1], m.m[1][1], m.m[2][1], m.m[3][1], 
                   m.m[0][2], m.m[1][2], m.m[2][2], m.m[3][2], 
                   m.m[0][3], m.m[1][3], m.m[2][3], m.m[3][3]);
}
