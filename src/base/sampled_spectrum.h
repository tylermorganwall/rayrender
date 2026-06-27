#ifndef RAYRENDER_BASE_SAMPLED_SPECTRUM_H
#define RAYRENDER_BASE_SAMPLED_SPECTRUM_H

#include "../math/float.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <initializer_list>
#include <sstream>
#include <stdexcept>
#include <string>

namespace rayrender {
namespace base {

inline constexpr int NSpectrumSamples = 4;
inline constexpr Float LambdaMin = 360.f;
inline constexpr Float LambdaMax = 830.f;
inline constexpr Float CIEYIntegral = static_cast<Float>(106.856895);

inline Float Sqr(Float value) {
  return value * value;
}

inline Float Lerp(Float t, Float a, Float b) {
  return (static_cast<Float>(1) - t) * a + t * b;
}

inline Float VisibleWavelengthsPDF(Float lambda) {
  if (lambda < static_cast<Float>(360) || lambda > static_cast<Float>(830)) {
    return 0;
  }
  Float c = std::cosh(static_cast<Float>(0.0072) * (lambda - static_cast<Float>(538)));
  return static_cast<Float>(0.0039398042) / Sqr(c);
}

inline Float SampleVisibleWavelengths(Float u) {
  return static_cast<Float>(538) -
         static_cast<Float>(138.888889) *
           std::atanh(static_cast<Float>(0.85691062) - static_cast<Float>(1.82750197) * u);
}

class SampledSpectrum {
public:
  SampledSpectrum() = default;

  explicit SampledSpectrum(Float c) {
    values_.fill(c);
  }

  explicit SampledSpectrum(const std::array<Float, NSpectrumSamples>& values) : values_(values) {}

  SampledSpectrum(std::initializer_list<Float> values) {
    if (values.size() != NSpectrumSamples) {
      throw std::invalid_argument("SampledSpectrum initializer must have NSpectrumSamples values");
    }
    std::copy(values.begin(), values.end(), values_.begin());
  }

  Float& operator[](int index) {
    assert(index >= 0 && index < NSpectrumSamples);
    return values_[index];
  }

  Float operator[](int index) const {
    assert(index >= 0 && index < NSpectrumSamples);
    return values_[index];
  }

  explicit operator bool() const {
    for (Float value : values_) {
      if (value != 0) {
        return true;
      }
    }
    return false;
  }

  SampledSpectrum& operator+=(const SampledSpectrum& spectrum) {
    for (int i = 0; i < NSpectrumSamples; ++i) {
      values_[i] += spectrum.values_[i];
    }
    return *this;
  }

  SampledSpectrum& operator-=(const SampledSpectrum& spectrum) {
    for (int i = 0; i < NSpectrumSamples; ++i) {
      values_[i] -= spectrum.values_[i];
    }
    return *this;
  }

  SampledSpectrum& operator*=(const SampledSpectrum& spectrum) {
    for (int i = 0; i < NSpectrumSamples; ++i) {
      values_[i] *= spectrum.values_[i];
    }
    return *this;
  }

  SampledSpectrum& operator*=(Float value) {
    assert(!std::isnan(value));
    for (Float& component : values_) {
      component *= value;
    }
    return *this;
  }

  SampledSpectrum& operator/=(const SampledSpectrum& spectrum) {
    for (int i = 0; i < NSpectrumSamples; ++i) {
      assert(spectrum.values_[i] != 0);
      values_[i] /= spectrum.values_[i];
    }
    return *this;
  }

  SampledSpectrum& operator/=(Float value) {
    assert(value != 0);
    assert(!std::isnan(value));
    for (Float& component : values_) {
      component /= value;
    }
    return *this;
  }

  SampledSpectrum operator-() const {
    SampledSpectrum result;
    for (int i = 0; i < NSpectrumSamples; ++i) {
      result.values_[i] = -values_[i];
    }
    return result;
  }

  bool Equals(const SampledSpectrum& spectrum) const {
    return values_ == spectrum.values_;
  }

  bool HasNaNs() const {
    for (Float value : values_) {
      if (std::isnan(value)) {
        return true;
      }
    }
    return false;
  }

  bool IsInf() const {
    for (Float value : values_) {
      if (std::isinf(value)) {
        return true;
      }
    }
    return false;
  }

  bool IsFinite() const {
    for (Float value : values_) {
      if (!std::isfinite(value)) {
        return false;
      }
    }
    return true;
  }

  bool IsPositive() const {
    return MinComponentValue() > 0;
  }

  Float MinComponentValue() const {
    Float value = values_[0];
    for (int i = 1; i < NSpectrumSamples; ++i) {
      value = std::min(value, values_[i]);
    }
    return value;
  }

  Float MaxComponentValue() const {
    Float value = values_[0];
    for (int i = 1; i < NSpectrumSamples; ++i) {
      value = std::max(value, values_[i]);
    }
    return value;
  }

  Float Average() const {
    Float sum = values_[0];
    for (int i = 1; i < NSpectrumSamples; ++i) {
      sum += values_[i];
    }
    return sum / NSpectrumSamples;
  }

  std::string ToString() const {
    std::ostringstream out;
    out << "[ ";
    for (int i = 0; i < NSpectrumSamples; ++i) {
      out << values_[i];
      if (i + 1 < NSpectrumSamples) {
        out << ", ";
      }
    }
    out << " ]";
    return out.str();
  }

private:
  std::array<Float, NSpectrumSamples> values_{};
};

inline SampledSpectrum operator+(SampledSpectrum lhs, const SampledSpectrum& rhs) {
  lhs += rhs;
  return lhs;
}

inline SampledSpectrum operator-(SampledSpectrum lhs, const SampledSpectrum& rhs) {
  lhs -= rhs;
  return lhs;
}

inline SampledSpectrum operator-(Float lhs, const SampledSpectrum& rhs) {
  assert(!std::isnan(lhs));
  SampledSpectrum result;
  for (int i = 0; i < NSpectrumSamples; ++i) {
    result[i] = lhs - rhs[i];
  }
  return result;
}

inline SampledSpectrum operator*(SampledSpectrum lhs, const SampledSpectrum& rhs) {
  lhs *= rhs;
  return lhs;
}

inline SampledSpectrum operator*(SampledSpectrum lhs, Float rhs) {
  lhs *= rhs;
  return lhs;
}

inline SampledSpectrum operator*(Float lhs, SampledSpectrum rhs) {
  rhs *= lhs;
  return rhs;
}

inline SampledSpectrum operator/(SampledSpectrum lhs, const SampledSpectrum& rhs) {
  lhs /= rhs;
  return lhs;
}

inline SampledSpectrum operator/(SampledSpectrum lhs, Float rhs) {
  lhs /= rhs;
  return lhs;
}

inline SampledSpectrum SafeDiv(SampledSpectrum numerator, SampledSpectrum denominator) {
  SampledSpectrum result;
  for (int i = 0; i < NSpectrumSamples; ++i) {
    result[i] = denominator[i] != 0 ? numerator[i] / denominator[i] : 0;
  }
  return result;
}

template <typename U, typename V>
inline SampledSpectrum Clamp(const SampledSpectrum& spectrum, U low, V high) {
  SampledSpectrum result;
  for (int i = 0; i < NSpectrumSamples; ++i) {
    result[i] = std::min(static_cast<Float>(high), std::max(static_cast<Float>(low), spectrum[i]));
  }
  assert(!result.HasNaNs());
  return result;
}

inline SampledSpectrum ClampZero(const SampledSpectrum& spectrum) {
  SampledSpectrum result;
  for (int i = 0; i < NSpectrumSamples; ++i) {
    result[i] = std::max(static_cast<Float>(0), spectrum[i]);
  }
  assert(!result.HasNaNs());
  return result;
}

inline SampledSpectrum Sqrt(const SampledSpectrum& spectrum) {
  SampledSpectrum result;
  for (int i = 0; i < NSpectrumSamples; ++i) {
    result[i] = std::sqrt(spectrum[i]);
  }
  assert(!result.HasNaNs());
  return result;
}

inline SampledSpectrum SafeSqrt(const SampledSpectrum& spectrum) {
  return Sqrt(ClampZero(spectrum));
}

inline SampledSpectrum Pow(const SampledSpectrum& spectrum, Float exponent) {
  SampledSpectrum result;
  for (int i = 0; i < NSpectrumSamples; ++i) {
    result[i] = std::pow(spectrum[i], exponent);
  }
  return result;
}

inline SampledSpectrum Exp(const SampledSpectrum& spectrum) {
  SampledSpectrum result;
  for (int i = 0; i < NSpectrumSamples; ++i) {
    result[i] = std::exp(spectrum[i]);
  }
  assert(!result.HasNaNs());
  return result;
}

class SampledWavelengths {
public:
  SampledWavelengths() = default;

  static SampledWavelengths SampleUniform(
    Float u,
    Float lambdaMin = LambdaMin,
    Float lambdaMax = LambdaMax
  ) {
    assert(std::isfinite(u));
    assert(u >= 0 && u <= 1);
    assert(lambdaMin < lambdaMax);

    SampledWavelengths wavelengths;
    wavelengths.lambda_[0] = Lerp(u, lambdaMin, lambdaMax);

    Float delta = (lambdaMax - lambdaMin) / NSpectrumSamples;
    for (int i = 1; i < NSpectrumSamples; ++i) {
      wavelengths.lambda_[i] = wavelengths.lambda_[i - 1] + delta;
      if (wavelengths.lambda_[i] > lambdaMax) {
        wavelengths.lambda_[i] = lambdaMin + (wavelengths.lambda_[i] - lambdaMax);
      }
    }

    for (int i = 0; i < NSpectrumSamples; ++i) {
      wavelengths.pdf_[i] = static_cast<Float>(1) / (lambdaMax - lambdaMin);
    }
    wavelengths.AssertInvariants(lambdaMin, lambdaMax);
    return wavelengths;
  }

  static SampledWavelengths SampleVisible(Float u) {
    assert(std::isfinite(u));
    assert(u >= 0 && u <= 1);

    SampledWavelengths wavelengths;
    for (int i = 0; i < NSpectrumSamples; ++i) {
      Float up = u + static_cast<Float>(i) / NSpectrumSamples;
      if (up > 1) {
        up -= 1;
      }

      wavelengths.lambda_[i] = SampleVisibleWavelengths(up);
      wavelengths.pdf_[i] = VisibleWavelengthsPDF(wavelengths.lambda_[i]);
    }
    wavelengths.AssertInvariants();
    return wavelengths;
  }

  Float operator[](int index) const {
    assert(index >= 0 && index < NSpectrumSamples);
    return lambda_[index];
  }

  Float& operator[](int index) {
    assert(index >= 0 && index < NSpectrumSamples);
    return lambda_[index];
  }

  Float PDF(int index) const {
    assert(index >= 0 && index < NSpectrumSamples);
    return pdf_[index];
  }

  SampledSpectrum PDF() const {
    return SampledSpectrum(pdf_);
  }

  void TerminateSecondary() {
    AssertInvariants();
    if (SecondaryTerminated()) {
      return;
    }
    for (int i = 1; i < NSpectrumSamples; ++i) {
      pdf_[i] = 0;
    }
    pdf_[0] /= NSpectrumSamples;
    AssertInvariants();
  }

  bool SecondaryTerminated() const {
    for (int i = 1; i < NSpectrumSamples; ++i) {
      if (pdf_[i] != 0) {
        return false;
      }
    }
    return true;
  }

  bool InvariantsHold(Float lambdaMin = LambdaMin, Float lambdaMax = LambdaMax) const {
    bool terminated = SecondaryTerminated();
    bool firstPDFPositive = pdf_[0] > 0 && std::isfinite(pdf_[0]);
    if (!firstPDFPositive) {
      return false;
    }
    for (int i = 0; i < NSpectrumSamples; ++i) {
      if (!std::isfinite(lambda_[i]) || lambda_[i] < lambdaMin - static_cast<Float>(1e-3) ||
          lambda_[i] > lambdaMax + static_cast<Float>(1e-3)) {
        return false;
      }
      if (!std::isfinite(pdf_[i]) || pdf_[i] < 0) {
        return false;
      }
      if (!terminated && pdf_[i] <= 0) {
        return false;
      }
      if (terminated && i > 0 && pdf_[i] != 0) {
        return false;
      }
    }
    return true;
  }

  std::string ToString() const {
    std::ostringstream out;
    out << "[ SampledWavelengths lambda: [";
    for (int i = 0; i < NSpectrumSamples; ++i) {
      out << ' ' << lambda_[i] << (i + 1 < NSpectrumSamples ? ',' : ' ');
    }
    out << "] pdf: [";
    for (int i = 0; i < NSpectrumSamples; ++i) {
      out << ' ' << pdf_[i] << (i + 1 < NSpectrumSamples ? ',' : ' ');
    }
    out << "] ]";
    return out.str();
  }

private:
  void AssertInvariants(Float lambdaMin = LambdaMin, Float lambdaMax = LambdaMax) const {
#ifndef NDEBUG
    assert(InvariantsHold(lambdaMin, lambdaMax));
#else
    (void)lambdaMin;
    (void)lambdaMax;
#endif
  }

  std::array<Float, NSpectrumSamples> lambda_{};
  std::array<Float, NSpectrumSamples> pdf_{};
};

} // namespace base
} // namespace rayrender

#endif
