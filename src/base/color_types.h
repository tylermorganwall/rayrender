#ifndef RAYRENDER_BASE_COLOR_TYPES_H
#define RAYRENDER_BASE_COLOR_TYPES_H

#include "../math/float.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace rayrender {
namespace base {

struct RGB {
  Float r = 0;
  Float g = 0;
  Float b = 0;

  RGB() = default;
  explicit RGB(Float value) : r(value), g(value), b(value) {}
  RGB(Float red, Float green, Float blue) : r(red), g(green), b(blue) {}

  Float& operator[](int index) { return (&r)[index]; }
  Float operator[](int index) const { return (&r)[index]; }

  bool HasNaNs() const { return std::isnan(r) || std::isnan(g) || std::isnan(b); }
  bool IsInf() const { return std::isinf(r) || std::isinf(g) || std::isinf(b); }
  bool IsFinite() const { return std::isfinite(r) && std::isfinite(g) && std::isfinite(b); }
  bool IsNonNegative() const { return r >= 0 && g >= 0 && b >= 0; }
  Float MaxComponentValue() const { return std::max(r, std::max(g, b)); }
  Float MinComponentValue() const { return std::min(r, std::min(g, b)); }
  Float Average() const { return (r + g + b) / static_cast<Float>(3); }
};

struct XYZ {
  Float x = 0;
  Float y = 0;
  Float z = 0;

  XYZ() = default;
  explicit XYZ(Float value) : x(value), y(value), z(value) {}
  XYZ(Float xx, Float yy, Float zz) : x(xx), y(yy), z(zz) {}

  Float& operator[](int index) { return (&x)[index]; }
  Float operator[](int index) const { return (&x)[index]; }

  bool HasNaNs() const { return std::isnan(x) || std::isnan(y) || std::isnan(z); }
  bool IsInf() const { return std::isinf(x) || std::isinf(y) || std::isinf(z); }
  bool IsFinite() const { return std::isfinite(x) && std::isfinite(y) && std::isfinite(z); }
};

inline RGB operator+(const RGB& lhs, const RGB& rhs) {
  return RGB(lhs.r + rhs.r, lhs.g + rhs.g, lhs.b + rhs.b);
}

inline RGB operator-(const RGB& lhs, const RGB& rhs) {
  return RGB(lhs.r - rhs.r, lhs.g - rhs.g, lhs.b - rhs.b);
}

inline RGB operator-(const RGB& value) {
  return RGB(-value.r, -value.g, -value.b);
}

inline RGB operator*(const RGB& lhs, const RGB& rhs) {
  return RGB(lhs.r * rhs.r, lhs.g * rhs.g, lhs.b * rhs.b);
}

inline RGB operator/(const RGB& lhs, const RGB& rhs) {
  return RGB(lhs.r / rhs.r, lhs.g / rhs.g, lhs.b / rhs.b);
}

inline RGB operator*(const RGB& lhs, Float rhs) {
  return RGB(lhs.r * rhs, lhs.g * rhs, lhs.b * rhs);
}

inline RGB operator*(Float lhs, const RGB& rhs) {
  return rhs * lhs;
}

inline RGB operator/(const RGB& lhs, Float rhs) {
  return RGB(lhs.r / rhs, lhs.g / rhs, lhs.b / rhs);
}

inline RGB& operator+=(RGB& lhs, const RGB& rhs) {
  lhs.r += rhs.r;
  lhs.g += rhs.g;
  lhs.b += rhs.b;
  return lhs;
}

inline RGB& operator-=(RGB& lhs, const RGB& rhs) {
  lhs.r -= rhs.r;
  lhs.g -= rhs.g;
  lhs.b -= rhs.b;
  return lhs;
}

inline RGB& operator*=(RGB& lhs, Float rhs) {
  lhs.r *= rhs;
  lhs.g *= rhs;
  lhs.b *= rhs;
  return lhs;
}

inline RGB& operator/=(RGB& lhs, Float rhs) {
  lhs.r /= rhs;
  lhs.g /= rhs;
  lhs.b /= rhs;
  return lhs;
}

inline XYZ operator+(const XYZ& lhs, const XYZ& rhs) {
  return XYZ(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
}

inline XYZ operator-(const XYZ& lhs, const XYZ& rhs) {
  return XYZ(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
}

inline XYZ operator*(const XYZ& lhs, Float rhs) {
  return XYZ(lhs.x * rhs, lhs.y * rhs, lhs.z * rhs);
}

inline XYZ operator*(Float lhs, const XYZ& rhs) {
  return rhs * lhs;
}

inline XYZ operator/(const XYZ& lhs, Float rhs) {
  return XYZ(lhs.x / rhs, lhs.y / rhs, lhs.z / rhs);
}

inline bool IsFinite(const RGB& value) {
  return value.IsFinite();
}

inline bool IsFinite(const XYZ& value) {
  return value.IsFinite();
}

inline Float SafeDivComponent(Float numerator, Float denominator, Float fallback = 0) {
  if (!std::isfinite(numerator) || !std::isfinite(denominator) || denominator == 0) {
    return fallback;
  }
  Float value = numerator / denominator;
  return std::isfinite(value) ? value : fallback;
}

inline RGB SafeDiv(const RGB& numerator, const RGB& denominator, Float fallback = 0) {
  return RGB(
    SafeDivComponent(numerator.r, denominator.r, fallback),
    SafeDivComponent(numerator.g, denominator.g, fallback),
    SafeDivComponent(numerator.b, denominator.b, fallback)
  );
}

inline RGB ClampZero(const RGB& value) {
  return RGB(
    std::max(static_cast<Float>(0), value.r),
    std::max(static_cast<Float>(0), value.g),
    std::max(static_cast<Float>(0), value.b)
  );
}

inline RGB Sqrt(const RGB& value) {
  RGB clamped = ClampZero(value);
  return RGB(std::sqrt(clamped.r), std::sqrt(clamped.g), std::sqrt(clamped.b));
}

inline RGB Exp(const RGB& value) {
  return RGB(std::exp(value.r), std::exp(value.g), std::exp(value.b));
}

class ColorSpaceMatrix3x3 {
public:
  ColorSpaceMatrix3x3() = default;

  explicit ColorSpaceMatrix3x3(const std::array<std::array<Float, 3>, 3>& rows) : m_(rows) {}

  static ColorSpaceMatrix3x3 Identity() {
    ColorSpaceMatrix3x3 result;
    result.m_[0][0] = 1;
    result.m_[1][1] = 1;
    result.m_[2][2] = 1;
    return result;
  }

  static ColorSpaceMatrix3x3 FromRows(
    std::array<Float, 3> row0,
    std::array<Float, 3> row1,
    std::array<Float, 3> row2
  ) {
    return ColorSpaceMatrix3x3({row0, row1, row2});
  }

  Float operator()(int row, int column) const { return m_[row][column]; }
  Float& operator()(int row, int column) { return m_[row][column]; }

  bool IsFinite() const {
    for (const auto& row : m_) {
      for (Float value : row) {
        if (!std::isfinite(value)) {
          return false;
        }
      }
    }
    return true;
  }

  std::array<Float, 3> Apply(const std::array<Float, 3>& value) const {
    return {
      m_[0][0] * value[0] + m_[0][1] * value[1] + m_[0][2] * value[2],
      m_[1][0] * value[0] + m_[1][1] * value[1] + m_[1][2] * value[2],
      m_[2][0] * value[0] + m_[2][1] * value[1] + m_[2][2] * value[2]
    };
  }

private:
  std::array<std::array<Float, 3>, 3> m_{{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
};

inline ColorSpaceMatrix3x3 operator*(const ColorSpaceMatrix3x3& lhs, const ColorSpaceMatrix3x3& rhs) {
  ColorSpaceMatrix3x3 result;
  for (int row = 0; row < 3; ++row) {
    for (int column = 0; column < 3; ++column) {
      Float sum = 0;
      for (int inner = 0; inner < 3; ++inner) {
        sum += lhs(row, inner) * rhs(inner, column);
      }
      result(row, column) = sum;
    }
  }
  return result;
}

inline XYZ TransformRGBToXYZ(const ColorSpaceMatrix3x3& matrix, const RGB& rgb) {
  std::array<Float, 3> out = matrix.Apply({rgb.r, rgb.g, rgb.b});
  return XYZ(out[0], out[1], out[2]);
}

inline RGB TransformXYZToRGB(const ColorSpaceMatrix3x3& matrix, const XYZ& xyz) {
  std::array<Float, 3> out = matrix.Apply({xyz.x, xyz.y, xyz.z});
  return RGB(out[0], out[1], out[2]);
}

enum class RGBTransferFunction {
  Linear,
  SRGB
};

struct RGBColorEncoding {
  RGBTransferFunction transfer = RGBTransferFunction::Linear;

  static RGBColorEncoding Linear() { return RGBColorEncoding{RGBTransferFunction::Linear}; }
  static RGBColorEncoding SRGB() { return RGBColorEncoding{RGBTransferFunction::SRGB}; }

  RGB Decode(const RGB& encoded) const {
    if (transfer == RGBTransferFunction::Linear) {
      return encoded;
    }
    return RGB(DecodeSRGB(encoded.r), DecodeSRGB(encoded.g), DecodeSRGB(encoded.b));
  }

  RGB Encode(const RGB& linear) const {
    if (transfer == RGBTransferFunction::Linear) {
      return linear;
    }
    return RGB(EncodeSRGB(linear.r), EncodeSRGB(linear.g), EncodeSRGB(linear.b));
  }

private:
  static Float DecodeSRGB(Float value) {
    if (value <= static_cast<Float>(0.04045)) {
      return value / static_cast<Float>(12.92);
    }
    return std::pow((value + static_cast<Float>(0.055)) / static_cast<Float>(1.055), static_cast<Float>(2.4));
  }

  static Float EncodeSRGB(Float value) {
    if (value <= static_cast<Float>(0.0031308)) {
      return value * static_cast<Float>(12.92);
    }
    return static_cast<Float>(1.055) * std::pow(value, static_cast<Float>(1.0 / 2.4)) - static_cast<Float>(0.055);
  }
};

struct RGBColorSpace {
  const char* name = "linear-identity";
  ColorSpaceMatrix3x3 rgbToXYZ = ColorSpaceMatrix3x3::Identity();
  ColorSpaceMatrix3x3 xyzToRGB = ColorSpaceMatrix3x3::Identity();
  RGBColorEncoding encoding = RGBColorEncoding::Linear();

  XYZ ToXYZ(const RGB& linearRGB) const {
    return TransformRGBToXYZ(rgbToXYZ, linearRGB);
  }

  RGB ToLinearRGB(const XYZ& xyz) const {
    return TransformXYZToRGB(xyzToRGB, xyz);
  }

  RGB Decode(const RGB& encodedRGB) const {
    return encoding.Decode(encodedRGB);
  }

  RGB Encode(const RGB& linearRGB) const {
    return encoding.Encode(linearRGB);
  }

  static RGBColorSpace SRGB() {
    RGBColorSpace colorSpace;
    colorSpace.name = "sRGB";
    colorSpace.encoding = RGBColorEncoding::SRGB();
    colorSpace.rgbToXYZ = ColorSpaceMatrix3x3::FromRows(
      {static_cast<Float>(0.4124564), static_cast<Float>(0.3575761), static_cast<Float>(0.1804375)},
      {static_cast<Float>(0.2126729), static_cast<Float>(0.7151522), static_cast<Float>(0.0721750)},
      {static_cast<Float>(0.0193339), static_cast<Float>(0.1191920), static_cast<Float>(0.9503041)}
    );
    colorSpace.xyzToRGB = ColorSpaceMatrix3x3::FromRows(
      {static_cast<Float>(3.2404542), static_cast<Float>(-1.5371385), static_cast<Float>(-0.4985314)},
      {static_cast<Float>(-0.9692660), static_cast<Float>(1.8760108), static_cast<Float>(0.0415560)},
      {static_cast<Float>(0.0556434), static_cast<Float>(-0.2040259), static_cast<Float>(1.0572252)}
    );
    return colorSpace;
  }
};

} // namespace base
} // namespace rayrender

#endif
