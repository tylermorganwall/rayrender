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

#ifndef MATH_FUNCTION_H_
#define MATH_FUNCTION_H_

#include <algorithm>
#include <cassert>
#include <cmath>

namespace dimensional {

// A function from [0:1] to values of type T represented by its values at N
// uniformly distributed samples (i + 0.5) / N, and linearly interpolated in
// between.
template<unsigned int N, typename T>
class Function {
 public:
  Function() {}

  explicit Function(const T& constant_value) {
    for (unsigned int i = 0; i < N; ++i) {
      value_[i] = constant_value;
    }
  }

  Function(const Function& rhs) = delete;
  Function& operator=(const Function& rhs) = delete;

  inline unsigned int size() const { return N; }

  inline const T& operator[](int index) const {
    assert(index >= 0 && index < static_cast<int>(N));
    return value_[index];
  }

  inline T& operator[](int index) {
    assert(index >= 0 && index < static_cast<int>(N));
    return value_[index];
  }

  const T operator()(double x) const {
    double u = x * N - 0.5;
    int i = std::floor(u);
    u -= i;
    int i0 = std::max(0, std::min(static_cast<int>(N) - 1, i));
    int i1 = std::max(0, std::min(static_cast<int>(N) - 1, i + 1));
    return value_[i0] * (1.0 - u) + value_[i1] * u;
  }

 private:
  T value_[N];
};

}  // namespace dimensional

#endif  // MATH_FUNCTION_H_
