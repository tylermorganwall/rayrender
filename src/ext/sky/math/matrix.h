/**
 * Copyright (c) 2016 Eric Bruneton
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holders nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef MATH_MATRIX_H_
#define MATH_MATRIX_H_

#include "../../../ext/sky/math/scalar.h"
#include "../../../ext/sky/math/vector.h"

namespace dimensional {

template<typename T>
struct Matrix3 {
  T m00, m01, m02;
  T m10, m11, m12;
  T m20, m21, m22;
  inline Matrix3() {}
  constexpr inline Matrix3(T m00, T m01, T m02,
                           T m10, T m11, T m12,
                           T m20, T m21, T m22)
      : m00(m00), m01(m01), m02(m02),
        m10(m10), m11(m11), m12(m12),
        m20(m20), m21(m21), m22(m22) {}
};

template<int U1, int U2, int U3, int U4, int U5,
    int V1, int V2, int V3, int V4, int V5>
Vector3<Scalar<U1 + V1, U2 + V2, U3 + V3, U4 + V4, U5 + V5>> operator*(
    const Matrix3<Scalar<U1, U2, U3, U4, U5>>& lhs,
    const Vector3<Scalar<V1, V2, V3, V4, V5>>& rhs) {
  return Vector3<Scalar<U1 + V1, U2 + V2, U3 + V3, U4 + V4, U5 + V5>>(
      lhs.m00 * rhs.x + lhs.m01 * rhs.y + lhs.m02 * rhs.z,
      lhs.m10 * rhs.x + lhs.m11 * rhs.y + lhs.m12 * rhs.z,
      lhs.m20 * rhs.x + lhs.m21 * rhs.y + lhs.m22 * rhs.z);
}

template<int U1, int U2, int U3, int U4, int U5,
    int V1, int V2, int V3, int V4, int V5>
Matrix3<Scalar<U1 + V1, U2 + V2, U3 + V3, U4 + V4, U5 + V5>> operator*(
    const Matrix3<Scalar<U1, U2, U3, U4, U5>>& lhs,
    const Matrix3<Scalar<V1, V2, V3, V4, V5>>& rhs) {
  return Matrix3<Scalar<U1 + V1, U2 + V2, U3 + V3, U4 + V4, U5 + V5>>(
      lhs.m00 * rhs.m00 + lhs.m01 * rhs.m10 + lhs.m02 * rhs.m20,
      lhs.m00 * rhs.m01 + lhs.m01 * rhs.m11 + lhs.m02 * rhs.m21,
      lhs.m00 * rhs.m02 + lhs.m01 * rhs.m12 + lhs.m02 * rhs.m22,
      lhs.m10 * rhs.m00 + lhs.m11 * rhs.m10 + lhs.m12 * rhs.m20,
      lhs.m10 * rhs.m01 + lhs.m11 * rhs.m11 + lhs.m12 * rhs.m21,
      lhs.m10 * rhs.m02 + lhs.m11 * rhs.m12 + lhs.m12 * rhs.m22,
      lhs.m20 * rhs.m00 + lhs.m21 * rhs.m10 + lhs.m22 * rhs.m20,
      lhs.m20 * rhs.m01 + lhs.m21 * rhs.m11 + lhs.m22 * rhs.m21,
      lhs.m20 * rhs.m02 + lhs.m21 * rhs.m12 + lhs.m22 * rhs.m22);
}

template<int U1, int U2, int U3, int U4, int U5>
Matrix3<Scalar<-U1, -U2, -U3, -U4, -U5>> inverse(
    const Matrix3<Scalar<U1, U2, U3, U4, U5>>& mat) {
  auto c00 = mat.m11 * mat.m22 - mat.m12 * mat.m21;
  auto c01 = mat.m02 * mat.m21 - mat.m01 * mat.m22;
  auto c02 = mat.m01 * mat.m12 - mat.m02 * mat.m11;
  auto c10 = mat.m12 * mat.m20 - mat.m10 * mat.m22;
  auto c11 = mat.m00 * mat.m22 - mat.m02 * mat.m20;
  auto c12 = mat.m02 * mat.m10 - mat.m00 * mat.m12;
  auto c20 = mat.m10 * mat.m21 - mat.m11 * mat.m20;
  auto c21 = mat.m01 * mat.m20 - mat.m00 * mat.m21;
  auto c22 = mat.m00 * mat.m11 - mat.m01 * mat.m10;
  auto det = mat.m00 * c00 + mat.m01 * c10 + mat.m02 * c20;
  return Matrix3<Scalar<-U1, -U2, -U3, -U4, -U5>>(
      c00 / det, c01 / det, c02 / det,
      c10 / det, c11 / det, c12 / det,
      c20 / det, c21 / det, c22 / det);
}

}  // namespace dimensional

#endif  // MATH_MATRIX_H_
