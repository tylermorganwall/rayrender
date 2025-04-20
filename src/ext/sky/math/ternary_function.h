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

#ifndef MATH_TERNARY_FUNCTION_H_
#define MATH_TERNARY_FUNCTION_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "../../../ext/sky/math/vector.h"

namespace dimensional {

// A function from [0:1]x[0:1]x[0:1] to values of type T represented by its
// values at NXxNYxNZ uniformly distributed samples (i + 0.5) / NX,
// (j + 0.5) / NY, (k + 0.5) / NZ and trilinearly interpolated in between.
template<unsigned int NX, unsigned int NY, unsigned int NZ, class T>
class TernaryFunction {
 public:
  TernaryFunction() { value_.reset(new T[NX * NY * NZ]); }

  explicit TernaryFunction(const T& constant_value) {
    value_.reset(new T[NX * NY * NZ]);
    for (unsigned int i = 0; i < NX * NY * NZ; ++i) {
      value_[i] = constant_value;
    }
  }

  TernaryFunction(const TernaryFunction& rhs) = delete;
  TernaryFunction& operator=(const TernaryFunction& rhs) = delete;

  TernaryFunction& operator+=(const TernaryFunction& rhs) {
    for (unsigned int i = 0; i < NX * NY * NZ; ++i) {
      value_[i] += rhs.value_[i];
    }
    return *this;
  }

  inline unsigned int size_x() const { return NX; }
  inline unsigned int size_y() const { return NY; }
  inline unsigned int size_z() const { return NZ; }

  virtual inline const T& Get(int i, int j, int k) const {
    assert(i >= 0 && i < static_cast<int>(NX) &&
           j >= 0 && j < static_cast<int>(NY) &&
           k >= 0 && k < static_cast<int>(NZ));
    return value_[i + j * NX + k * NX * NY];
  }

  void Set(const T& constant_value) {
    for (unsigned int i = 0; i < NX * NY * NZ; ++i) {
      value_[i] = constant_value;
    }
  }

  void Set(const TernaryFunction& rhs) {
    for (unsigned int i = 0; i < NX * NY * NZ; ++i) {
      value_[i] = rhs.value_[i];
    }
  }

  inline void Set(int i, int j, int k, T value) {
    assert(i >= 0 && i < static_cast<int>(NX) &&
           j >= 0 && j < static_cast<int>(NY) &&
           k >= 0 && k < static_cast<int>(NZ));
    value_[i + j * NX + k * NX * NY] = value;
  }

  const T operator()(double x, double y, double z) const {
    double u = x * NX - 0.5;
    double v = y * NY - 0.5;
    double w = z * NZ - 0.5;
    int i = std::floor(u);
    int j = std::floor(v);
    int k = std::floor(w);
    u -= i;
    v -= j;
    w -= k;
    int i0 = std::max(0, std::min(static_cast<int>(NX) - 1, i));
    int i1 = std::max(0, std::min(static_cast<int>(NX) - 1, i + 1));
    int j0 = std::max(0, std::min(static_cast<int>(NY) - 1, j));
    int j1 = std::max(0, std::min(static_cast<int>(NY) - 1, j + 1));
    int k0 = std::max(0, std::min(static_cast<int>(NZ) - 1, k));
    int k1 = std::max(0, std::min(static_cast<int>(NZ) - 1, k + 1));
    return Get(i0, j0, k0) * ((1.0 - u) * (1.0 - v) * (1.0 - w)) +
        Get(i1, j0, k0) * (u * (1.0 - v) * (1.0 - w)) +
        Get(i0, j1, k0) * ((1.0 - u) * v * (1.0 - w)) +
        Get(i1, j1, k0) * (u * v * (1.0 - w)) +
        Get(i0, j0, k1) * ((1.0 - u) * (1.0 - v) * w) +
        Get(i1, j0, k1) * (u * (1.0 - v) * w) +
        Get(i0, j1, k1) * ((1.0 - u) * v * w) +
        Get(i1, j1, k1) * (u * v * w);
  }

  void Load(const std::string& filename) {
    std::ifstream file(filename, std::ifstream::binary | std::ifstream::in);
    file.read(reinterpret_cast<char*>(value_.get()), NX * NY * NZ * sizeof(T));
    file.close();
  }

  void Save(const std::string& filename) const {
    std::ofstream file(filename, std::ofstream::binary);
    file.write(
        reinterpret_cast<const char*>(value_.get()), NX * NY * NZ * sizeof(T));
    file.close();
  }

 protected:
  std::unique_ptr<T[]> value_;
};

template<unsigned int NX, unsigned int NY, unsigned int NZ, typename T>
T texture(const TernaryFunction<NX, NY, NZ, T>& table, const vec3& uvw) {
  return table(uvw.x(), uvw.y(), uvw.z());
}

}  // namespace dimensional

#endif  // MATH_TERNARY_FUNCTION_H_
