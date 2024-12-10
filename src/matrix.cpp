#include "matrix.h"
#include "simd.h"
#include <testthat.h>

static_assert(std::is_trivially_copyable<FVec4>::value,
              "FVec4 must be trivially copyable to use memcpy safely.");

Matrix4x4::Matrix4x4(Float mat[4][4]) { memcpy((Float*)m, mat, 16 * sizeof(Float)); }

Matrix4x4::Matrix4x4(Float t00, Float t01, Float t02, Float t03, 
                     Float t10, Float t11, Float t12, Float t13, 
                     Float t20, Float t21, Float t22, Float t23, 
                     Float t30, Float t31, Float t32, Float t33) {
#ifdef RAYSIMDVEC
  m[0] = FVec4(t00, t01, t02, t03);
  m[1] = FVec4(t10, t11, t12, t13);
  m[2] = FVec4(t20, t21, t22, t23);
  m[3] = FVec4(t30, t31, t32, t33);
#else
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
#endif
}

#include "RcppThread.h"

// #ifndef RAYSIMDVEC

// Matrix4x4 Inverse(const Matrix4x4 &m) {
//     int indxc[4], indxr[4];
//     int ipiv[4] = {0, 0, 0, 0};
//     FVec4 minv[4];
//     memcpy(minv, m.m, 4 * sizeof(FVec4));

//     for (int i = 0; i < 4; i++) {
//         int irow = 0, icol = 0;
//         Float big = 0.f;
//         // Pivot selection
//         for (int j = 0; j < 4; j++) {
//             if (ipiv[j] != 1) {
//                 for (int k = 0; k < 4; k++) {
//                     if (ipiv[k] == 0) {
//                         Float absval = std::fabsf(minv[j][k]);

//                         if (absval >= big) {
//                             big = absval;
//                             irow = j;
//                             icol = k;
//                         }
//                     } else if (ipiv[k] > 1) [[unlikely]] {
//                         throw std::runtime_error("Singular matrix in MatrixInvert");
//                     }
//                 }
//             }
//         }
//         ++ipiv[icol];
//         // Row swapping
//         if (irow != icol) {
//             std::swap(minv[irow], minv[icol]);
//         }
//         indxr[i] = irow;
//         indxc[i] = icol;
//         if (minv[icol][icol] == 0.f) [[unlikely]] {
//             throw std::runtime_error("Singular matrix in MatrixInvert");
//         }

//         // Scaling pivot row
//         Float pivinv = 1.0f / minv[icol][icol];
//         minv[icol][icol] = 1.0f;
//         // Scale the entire row
//         minv[icol] = simd_mul(minv[icol], simd_set1(pivinv));

//         // Eliminating other rows
//         for (int j = 0; j < 4; j++) {
//             if (j != icol) {
//                 Float save = minv[j][icol];
//                 minv[j][icol] = 0.0f;
//                 FVec4 scaled_row = simd_mul(minv[icol], simd_set1(save));
//                 minv[j] = simd_sub(minv[j], scaled_row);
//             }
//         }
//     }
//     // Column swapping to reflect permutation
//     for (int j = 3; j >= 0; j--) {
//         if (indxr[j] != indxc[j]) {
//             for (int k = 0; k < 4; k++) {
//                 std::swap(minv[k][indxr[j]], minv[k][indxc[j]]);
//             }
//         }
//     }
//     return Matrix4x4(minv);
// }
// #else
Matrix4x4 Inverse(const Matrix4x4 &m) {
  int indxc[4], indxr[4];
  int ipiv[4] = {0, 0, 0, 0};
  Float minv[4][4];
  std::memcpy(&minv, (Float*)&m.m, 4 * 4 * sizeof(Float));
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
          } else if (ipiv[k] > 1) [[unlikely]] {
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
    if (minv[icol][icol] == 0.f)[[unlikely]] {
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
// #endif

Matrix4x4 Transpose(const Matrix4x4 &m) {
  return Matrix4x4(m.m[0][0], m.m[1][0], m.m[2][0], m.m[3][0], 
                   m.m[0][1], m.m[1][1], m.m[2][1], m.m[3][1], 
                   m.m[0][2], m.m[1][2], m.m[2][2], m.m[3][2], 
                   m.m[0][3], m.m[1][3], m.m[2][3], m.m[3][3]);
}

#include <testthat.h>
#include <Rcpp.h>

context("Inverse computes matrix inverse correctly") {
    test_that("[Inverse of identity matrix]") {
        // Create identity matrix
        Matrix4x4 identity(
            1.0f, 0.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 0.0f, 1.0f
        );

        // Compute inverse
        Matrix4x4 inv_identity = Inverse(identity);

        // Check that inverse is the identity matrix
        for (int i = 0; i < 4; ++i) {
            #ifdef RAYSIMDVEC
            for (int j = 0; j < 4; ++j) {
                expect_true(std::fabs(inv_identity.m[i][j] - identity.m[i][j]) < 1e-6f);
            }
            #else
            for (int j = 0; j < 4; ++j) {
                expect_true(std::fabs(inv_identity.m[i][j] - identity.m[i][j]) < 1e-6f);
            }
            #endif
        }
    }

    test_that("[Inverse of a known matrix]") {
        // Define a known invertible matrix
        Matrix4x4 mat(
            4.0f, 7.0f, 2.0f, 3.0f,
            0.0f, 5.0f, 0.0f, 1.0f,
            0.0f, 0.0f, 6.0f, 0.0f,
            0.0f, 0.0f, 0.0f, 1.0f
        );

        // Compute inverse
        Matrix4x4 inv_mat = Inverse(mat);

        // Multiply original matrix by its inverse
        Matrix4x4 product = mat.Mul(mat,inv_mat);

        // Create identity matrix
        Matrix4x4 identity(
            1.0f, 0.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 0.0f, 1.0f
        );

        // Check that product is close to identity matrix
        for (int i = 0; i < 4; ++i) {
            #ifdef RAYSIMDVEC
            for (int j = 0; j < 4; ++j) {
                expect_true(std::fabs(product.m[i][j] - identity.m[i][j]) < 1e-4f);
            }
            #else
            for (int j = 0; j < 4; ++j) {
                expect_true(std::fabs(product.m[i][j] - identity.m[i][j]) < 1e-4f);
            }
            #endif
        }
    }

    test_that("[Inverse throws exception for singular matrix]") {
        // Define a singular matrix (determinant is zero)
        Matrix4x4 singular_mat(
            1.0f, 2.0f, 3.0f, 4.0f,
            2.0f, 4.0f, 6.0f, 8.0f,
            3.0f, 6.0f, 9.0f, 12.0f,
            4.0f, 8.0f, 12.0f, 16.0f
        );

        // Attempt to compute inverse and expect an exception
        expect_error_as(Inverse(singular_mat), std::runtime_error);
    }

    test_that("[Inverse of random matrices]") {
        // Define multiple random invertible matrices
        std::vector<Matrix4x4> matrices = {
            Matrix4x4(
                2.0f, 1.0f, 0.0f, 0.0f,
                1.0f, 2.0f, 1.0f, 0.0f,
                0.0f, 1.0f, 2.0f, 1.0f,
                0.0f, 0.0f, 1.0f, 2.0f
            ),
            Matrix4x4(
                5.0f, -2.0f, 2.0f, 7.0f,
                1.0f, 0.0f, 0.0f, 3.0f,
                -3.0f, 1.0f, 5.0f, 0.0f,
                3.0f, -1.0f, -9.0f, 4.0f
            )
            // Add more matrices if desired
        };

        // Create identity matrix
        Matrix4x4 identity(
            1.0f, 0.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 0.0f, 1.0f
        );

        for (const auto& mat : matrices) {
            // Compute inverse
            Matrix4x4 inv_mat = Inverse(mat);

            // Multiply original matrix by its inverse
            Matrix4x4 product = mat.Mul(mat,inv_mat);

            // Check that product is close to identity matrix
            for (int i = 0; i < 4; ++i) {
                #ifdef RAYSIMDVEC
                for (int j = 0; j < 4; ++j) {
                    expect_true(std::fabs(product.m[i][j] - identity.m[i][j]) < 1e-4f);
                }
                #else
                for (int j = 0; j < 4; ++j) {
                    expect_true(std::fabs(product.m[i][j] - identity.m[i][j]) < 1e-4f);
                }
                #endif
            }
        }
    }
}
