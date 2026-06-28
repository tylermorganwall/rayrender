#ifndef RAYRENDER_BASE_COLOR_TYPES_H
#define RAYRENDER_BASE_COLOR_TYPES_H

#include "../math/float.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>

namespace rayrender {
namespace base {

class DenselySampledSpectrum;
class RGBSigmoidPolynomial;
class RGBToSpectrumTable;

struct Chromaticity {
  Float x = 0;
  Float y = 0;

  Chromaticity() = default;
  Chromaticity(Float xx, Float yy) : x(xx), y(yy) {}
};

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

  Chromaticity xy() const {
    Float sum = x + y + z;
    if (sum == 0) {
      return Chromaticity(0, 0);
    }
    return Chromaticity(x / sum, y / sum);
  }

  static XYZ FromxyY(const Chromaticity& xy, Float luminance = 1) {
    if (xy.y == 0) {
      return XYZ(0, 0, 0);
    }
    return XYZ(
      xy.x * luminance / xy.y,
      luminance,
      (static_cast<Float>(1) - xy.x - xy.y) * luminance / xy.y
    );
  }
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

  static ColorSpaceMatrix3x3 FromColumns(
    std::array<Float, 3> column0,
    std::array<Float, 3> column1,
    std::array<Float, 3> column2
  ) {
    return ColorSpaceMatrix3x3::FromRows(
      {column0[0], column1[0], column2[0]},
      {column0[1], column1[1], column2[1]},
      {column0[2], column1[2], column2[2]}
    );
  }

  static ColorSpaceMatrix3x3 Diagonal(Float d0, Float d1, Float d2) {
    ColorSpaceMatrix3x3 result;
    result(0, 0) = d0;
    result(1, 1) = d1;
    result(2, 2) = d2;
    return result;
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

  Float Determinant() const {
    return m_[0][0] * (m_[1][1] * m_[2][2] - m_[1][2] * m_[2][1]) -
           m_[0][1] * (m_[1][0] * m_[2][2] - m_[1][2] * m_[2][0]) +
           m_[0][2] * (m_[1][0] * m_[2][1] - m_[1][1] * m_[2][0]);
  }

  ColorSpaceMatrix3x3 Inverse() const {
    Float determinant = Determinant();
    if (determinant == 0 || !std::isfinite(determinant)) {
      throw std::invalid_argument("Color space matrix is singular");
    }

    Float invDet = static_cast<Float>(1) / determinant;
    return ColorSpaceMatrix3x3::FromRows(
      {
        (m_[1][1] * m_[2][2] - m_[1][2] * m_[2][1]) * invDet,
        (m_[0][2] * m_[2][1] - m_[0][1] * m_[2][2]) * invDet,
        (m_[0][1] * m_[1][2] - m_[0][2] * m_[1][1]) * invDet
      },
      {
        (m_[1][2] * m_[2][0] - m_[1][0] * m_[2][2]) * invDet,
        (m_[0][0] * m_[2][2] - m_[0][2] * m_[2][0]) * invDet,
        (m_[0][2] * m_[1][0] - m_[0][0] * m_[1][2]) * invDet
      },
      {
        (m_[1][0] * m_[2][1] - m_[1][1] * m_[2][0]) * invDet,
        (m_[0][1] * m_[2][0] - m_[0][0] * m_[2][1]) * invDet,
        (m_[0][0] * m_[1][1] - m_[0][1] * m_[1][0]) * invDet
      }
    );
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

  static Float SRGBToLinear(Float value) {
    if (value <= static_cast<Float>(0.04045)) {
      return value * static_cast<Float>(1.0 / 12.92);
    }
    Float p = EvaluatePolynomial(
      value,
      static_cast<Float>(-0.0163933279112946),
      static_cast<Float>(-0.7386328024653209),
      static_cast<Float>(-11.199318357635072),
      static_cast<Float>(-47.46726633009393),
      static_cast<Float>(-36.04572663838034)
    );
    Float q = EvaluatePolynomial(
      value,
      static_cast<Float>(-0.004261480793199332),
      static_cast<Float>(-19.140923959601675),
      static_cast<Float>(-59.096406619244426),
      static_cast<Float>(-18.225745396846637),
      static_cast<Float>(1)
    );
    return p / q * value;
  }

  static Float LinearToSRGB(Float value) {
    if (value <= static_cast<Float>(0.0031308)) {
      return value * static_cast<Float>(12.92);
    }
    Float sqrtValue = std::sqrt(std::max(static_cast<Float>(0), value));
    Float p = EvaluatePolynomial(
      sqrtValue,
      static_cast<Float>(-0.0016829072605308378),
      static_cast<Float>(0.03453868659826638),
      static_cast<Float>(0.7642611304733891),
      static_cast<Float>(2.0041169284241644),
      static_cast<Float>(0.7551545191665577),
      static_cast<Float>(-0.016202083165206348)
    );
    Float q = EvaluatePolynomial(
      sqrtValue,
      static_cast<Float>(4.178892964897981e-7),
      static_cast<Float>(-0.00004375359692957097),
      static_cast<Float>(0.03467195408529984),
      static_cast<Float>(0.6085338522168684),
      static_cast<Float>(1.8970238036421054),
      static_cast<Float>(1)
    );
    return p / q * value;
  }

private:
  static Float EvaluatePolynomial(Float, Float coefficient) {
    return coefficient;
  }

  template <typename... Coefficients>
  static Float EvaluatePolynomial(Float t, Float coefficient, Coefficients... remaining) {
    return t * EvaluatePolynomial(t, static_cast<Float>(remaining)...) + coefficient;
  }

  static Float DecodeSRGB(Float value) {
    return SRGBToLinear(value);
  }

  static Float EncodeSRGB(Float value) {
    return LinearToSRGB(value);
  }
};

struct RGBColorSpace {
  const char* name = "linear-identity";
  Chromaticity r;
  Chromaticity g;
  Chromaticity b;
  Chromaticity w;
  ColorSpaceMatrix3x3 rgbToXYZ = ColorSpaceMatrix3x3::Identity();
  ColorSpaceMatrix3x3 xyzToRGB = ColorSpaceMatrix3x3::Identity();
  RGBColorEncoding encoding = RGBColorEncoding::Linear();
  std::shared_ptr<const DenselySampledSpectrum> illuminant;
  std::shared_ptr<const RGBToSpectrumTable> rgbToSpectrumTable;

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

  bool HasSpectralReconstruction() const {
    return static_cast<bool>(illuminant) && static_cast<bool>(rgbToSpectrumTable);
  }

  RGBSigmoidPolynomial ToRGBCoeffs(const RGB& rgb) const;

  static RGBColorSpace FromPrimaries(
    const char* name,
    Chromaticity red,
    Chromaticity green,
    Chromaticity blue,
    const XYZ& whiteXYZ,
    RGBColorEncoding encoding = RGBColorEncoding::Linear()
  ) {
    XYZ R = XYZ::FromxyY(red);
    XYZ G = XYZ::FromxyY(green);
    XYZ B = XYZ::FromxyY(blue);
    ColorSpaceMatrix3x3 rgb = ColorSpaceMatrix3x3::FromColumns(
      {R.x, R.y, R.z},
      {G.x, G.y, G.z},
      {B.x, B.y, B.z}
    );
    std::array<Float, 3> scale = rgb.Inverse().Apply({whiteXYZ.x, whiteXYZ.y, whiteXYZ.z});

    RGBColorSpace colorSpace;
    colorSpace.name = name;
    colorSpace.r = red;
    colorSpace.g = green;
    colorSpace.b = blue;
    colorSpace.w = whiteXYZ.xy();
    colorSpace.encoding = encoding;
    colorSpace.rgbToXYZ = rgb * ColorSpaceMatrix3x3::Diagonal(scale[0], scale[1], scale[2]);
    colorSpace.xyzToRGB = colorSpace.rgbToXYZ.Inverse();
    return colorSpace;
  }

  static RGBColorSpace FromPrimaries(
    const char* name,
    Chromaticity red,
    Chromaticity green,
    Chromaticity blue,
    Chromaticity white,
    RGBColorEncoding encoding = RGBColorEncoding::Linear()
  ) {
    return FromPrimaries(name, red, green, blue, XYZ::FromxyY(white), encoding);
  }

  static RGBColorSpace SRGB() {
    return FromPrimaries(
      "sRGB",
      Chromaticity(static_cast<Float>(0.64), static_cast<Float>(0.33)),
      Chromaticity(static_cast<Float>(0.30), static_cast<Float>(0.60)),
      Chromaticity(static_cast<Float>(0.15), static_cast<Float>(0.06)),
      Chromaticity(static_cast<Float>(0.3127), static_cast<Float>(0.3290)),
      RGBColorEncoding::SRGB()
    );
  }

  static RGBColorSpace DCIP3() {
    return FromPrimaries(
      "DCI-P3",
      Chromaticity(static_cast<Float>(0.680), static_cast<Float>(0.320)),
      Chromaticity(static_cast<Float>(0.265), static_cast<Float>(0.690)),
      Chromaticity(static_cast<Float>(0.150), static_cast<Float>(0.060)),
      Chromaticity(static_cast<Float>(0.3127), static_cast<Float>(0.3290))
    );
  }

  static RGBColorSpace Rec2020() {
    return FromPrimaries(
      "Rec.2020",
      Chromaticity(static_cast<Float>(0.708), static_cast<Float>(0.292)),
      Chromaticity(static_cast<Float>(0.170), static_cast<Float>(0.797)),
      Chromaticity(static_cast<Float>(0.131), static_cast<Float>(0.046)),
      Chromaticity(static_cast<Float>(0.3127), static_cast<Float>(0.3290))
    );
  }

  static RGBColorSpace ACES2065_1() {
    return FromPrimaries(
      "ACES2065-1",
      Chromaticity(static_cast<Float>(0.7347), static_cast<Float>(0.2653)),
      Chromaticity(static_cast<Float>(0.0), static_cast<Float>(1.0)),
      Chromaticity(static_cast<Float>(0.0001), static_cast<Float>(-0.0770)),
      Chromaticity(static_cast<Float>(0.32168), static_cast<Float>(0.33767))
    );
  }

  static std::string NormalizedName(std::string name) {
    std::string normalized;
    normalized.reserve(name.size());
    for (char c : name) {
      if (c == '_' || c == '.' || c == ' ') {
        continue;
      }
      normalized.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    }
    return normalized;
  }

  static std::string CanonicalName(const std::string& name) {
    std::string normalized = NormalizedName(name);
    if (normalized == "srgb") {
      return "sRGB";
    }
    if (normalized == "dci-p3" || normalized == "dcip3") {
      return "DCI-P3";
    }
    if (normalized == "rec2020") {
      return "Rec.2020";
    }
    if (normalized == "aces2065-1" || normalized == "aces20651") {
      return "ACES2065-1";
    }
    throw std::invalid_argument("Unknown RGB color space: " + name);
  }

  static bool IsKnownName(const std::string& name) {
    try {
      (void)CanonicalName(name);
    } catch (const std::invalid_argument&) {
      return false;
    }
    return true;
  }

  static RGBColorSpace Named(const std::string& name) {
    std::string canonical = CanonicalName(name);
    if (canonical == "sRGB") {
      return SRGB();
    }
    if (canonical == "DCI-P3") {
      return DCIP3();
    }
    if (canonical == "Rec.2020") {
      return Rec2020();
    }
    return ACES2065_1();
  }
};

} // namespace base
} // namespace rayrender

#endif
