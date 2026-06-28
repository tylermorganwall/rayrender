#ifndef RAYRENDER_BASE_RGB_SPECTRUM_H
#define RAYRENDER_BASE_RGB_SPECTRUM_H

#include "spectrum.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace rayrender {
namespace base {

inline constexpr int RGBToSpectrumResolution = 64;
inline constexpr int RGBToSpectrumMaxComponentCount = 3;
inline constexpr int RGBToSpectrumCoefficientChannels = 3;
inline constexpr std::size_t RGBToSpectrumCoefficientCount =
  static_cast<std::size_t>(RGBToSpectrumMaxComponentCount) *
  RGBToSpectrumResolution *
  RGBToSpectrumResolution *
  RGBToSpectrumResolution *
  RGBToSpectrumCoefficientChannels;
inline constexpr std::size_t RGBToSpectrumPayloadFloatCount =
  RGBToSpectrumResolution + RGBToSpectrumCoefficientCount;

struct RGBToSpectrumTableLoadStats {
  std::uint64_t lookups = 0;
  std::uint64_t cacheHits = 0;
  std::uint64_t fileLoads = 0;
};

class RGBSigmoidPolynomial {
public:
  RGBSigmoidPolynomial() = default;
  RGBSigmoidPolynomial(Float c0, Float c1, Float c2) : c0_(c0), c1_(c1), c2_(c2) {}

  Float operator()(Float lambda) const {
    Float value = c2_ + lambda * (c1_ + lambda * c0_);
    return Sigmoid(value);
  }

  Float MaxValue() const {
    Float result = std::max((*this)(LambdaMin), (*this)(LambdaMax));
    if (c0_ != 0) {
      Float lambda = -c1_ / (static_cast<Float>(2) * c0_);
      if (lambda >= LambdaMin && lambda <= LambdaMax) {
        result = std::max(result, (*this)(lambda));
      }
    }
    return result;
  }

  Float C0() const { return c0_; }
  Float C1() const { return c1_; }
  Float C2() const { return c2_; }

  std::string ToString() const {
    std::ostringstream out;
    out << "[ RGBSigmoidPolynomial c0: " << c0_
        << " c1: " << c1_ << " c2: " << c2_ << " ]";
    return out.str();
  }

private:
  static Float Sigmoid(Float value) {
    if (std::isinf(value)) {
      return value > 0 ? static_cast<Float>(1) : static_cast<Float>(0);
    }
    return static_cast<Float>(0.5) +
           value / (static_cast<Float>(2) * std::sqrt(static_cast<Float>(1) + Sqr(value)));
  }

  Float c0_ = 0;
  Float c1_ = 0;
  Float c2_ = 0;
};

namespace detail {

inline constexpr std::array<unsigned char, 16> RGBToSpectrumTableMagic = {
  'R', 'A', 'Y', 'R', 'G', 'B', 'S', 'P', 'E', 'C', 'v', '1', 0, 0, 0, 0
};
inline constexpr std::uint32_t RGBToSpectrumTableVersion = 1;
inline constexpr std::uint32_t RGBToSpectrumTableEndianMarker = 0x01020304u;
inline constexpr std::size_t RGBToSpectrumTableHeaderSize = 16 + 10 * 4 + 32;

inline std::uint32_t ReadLittleEndianUInt32(const std::vector<unsigned char>& bytes, std::size_t offset) {
  if (offset + 4 > bytes.size()) {
    throw std::runtime_error("RGB-to-spectrum table header is truncated");
  }
  return static_cast<std::uint32_t>(bytes[offset]) |
         (static_cast<std::uint32_t>(bytes[offset + 1]) << 8) |
         (static_cast<std::uint32_t>(bytes[offset + 2]) << 16) |
         (static_cast<std::uint32_t>(bytes[offset + 3]) << 24);
}

inline Float ReadLittleEndianFloat32(const std::vector<unsigned char>& bytes, std::size_t offset) {
  if (!std::numeric_limits<float>::is_iec559) {
    throw std::runtime_error("RGB-to-spectrum table loading requires IEEE-754 float support");
  }
  std::uint32_t bits = ReadLittleEndianUInt32(bytes, offset);
  float value = 0;
  std::memcpy(&value, &bits, sizeof(value));
  if (!std::isfinite(value) && !std::isinf(value)) {
    throw std::runtime_error("RGB-to-spectrum table contains an invalid float32 value");
  }
  return static_cast<Float>(value);
}

inline std::uint32_t Adler32(const std::vector<unsigned char>& bytes, std::size_t offset, std::size_t count) {
  constexpr std::uint32_t modulus = 65521;
  std::uint32_t a = 1;
  std::uint32_t b = 0;
  for (std::size_t i = 0; i < count; ++i) {
    a = (a + bytes[offset + i]) % modulus;
    b = (b + a) % modulus;
  }
  return (b << 16) | a;
}

inline std::string ReadFixedString(const std::vector<unsigned char>& bytes, std::size_t offset, std::size_t count) {
  if (offset + count > bytes.size()) {
    throw std::runtime_error("RGB-to-spectrum table color-space id is truncated");
  }
  std::string value;
  for (std::size_t i = 0; i < count && bytes[offset + i] != 0; ++i) {
    value.push_back(static_cast<char>(bytes[offset + i]));
  }
  return value;
}

template <typename Predicate>
inline int FindInterval(int size, Predicate predicate) {
  int first = 0;
  int length = size;
  while (length > 0) {
    int half = length >> 1;
    int middle = first + half;
    if (predicate(middle)) {
      first = middle + 1;
      length -= half + 1;
    } else {
      length = half;
    }
  }
  return std::clamp(first - 1, 0, size - 2);
}

} // namespace detail

class RGBToSpectrumTable {
public:
  RGBToSpectrumTable() = default;

  RGBToSpectrumTable(
    std::string colorSpaceId,
    std::array<Float, RGBToSpectrumResolution> zNodes,
    std::vector<Float> coefficients
  )
    : colorSpaceId_(std::move(colorSpaceId)),
      zNodes_(zNodes),
      coefficients_(std::move(coefficients)) {
    if (coefficients_.size() != RGBToSpectrumCoefficientCount) {
      throw std::invalid_argument("RGB-to-spectrum table coefficient count does not match pbrt layout");
    }
  }

  static RGBToSpectrumTable LoadFromFile(const std::string& path, const std::string& expectedColorSpaceId = "sRGB") {
    std::ifstream input(path.c_str(), std::ios::binary);
    if (!input) {
      throw std::runtime_error("Unable to open RGB-to-spectrum table file: " + path);
    }

    input.seekg(0, std::ios::end);
    std::streamoff streamSize = input.tellg();
    if (streamSize < static_cast<std::streamoff>(detail::RGBToSpectrumTableHeaderSize)) {
      throw std::runtime_error("RGB-to-spectrum table file is truncated: " + path);
    }
    input.seekg(0, std::ios::beg);

    std::vector<unsigned char> bytes(static_cast<std::size_t>(streamSize));
    input.read(reinterpret_cast<char*>(bytes.data()), static_cast<std::streamsize>(streamSize));
    if (!input) {
      throw std::runtime_error("Unable to read RGB-to-spectrum table file: " + path);
    }

    for (std::size_t i = 0; i < detail::RGBToSpectrumTableMagic.size(); ++i) {
      if (bytes[i] != detail::RGBToSpectrumTableMagic[i]) {
        throw std::runtime_error("RGB-to-spectrum table has invalid magic: " + path);
      }
    }

    std::size_t offset = detail::RGBToSpectrumTableMagic.size();
    std::uint32_t version = detail::ReadLittleEndianUInt32(bytes, offset);
    offset += 4;
    std::uint32_t endianMarker = detail::ReadLittleEndianUInt32(bytes, offset);
    offset += 4;
    std::uint32_t scalarBytes = detail::ReadLittleEndianUInt32(bytes, offset);
    offset += 4;
    std::uint32_t resolution = detail::ReadLittleEndianUInt32(bytes, offset);
    offset += 4;
    std::uint32_t maxComponentCases = detail::ReadLittleEndianUInt32(bytes, offset);
    offset += 4;
    std::uint32_t coefficientChannels = detail::ReadLittleEndianUInt32(bytes, offset);
    offset += 4;
    std::uint32_t scaleCount = detail::ReadLittleEndianUInt32(bytes, offset);
    offset += 4;
    std::uint32_t coefficientCount = detail::ReadLittleEndianUInt32(bytes, offset);
    offset += 4;
    std::uint32_t payloadBytes = detail::ReadLittleEndianUInt32(bytes, offset);
    offset += 4;
    std::uint32_t expectedChecksum = detail::ReadLittleEndianUInt32(bytes, offset);
    offset += 4;
    std::string colorSpaceId = detail::ReadFixedString(bytes, offset, 32);
    offset += 32;

    if (version != detail::RGBToSpectrumTableVersion) {
      throw std::runtime_error("RGB-to-spectrum table version mismatch in " + path);
    }
    if (endianMarker != detail::RGBToSpectrumTableEndianMarker) {
      throw std::runtime_error("RGB-to-spectrum table byte-order marker mismatch in " + path);
    }
    if (scalarBytes != 4) {
      throw std::runtime_error("RGB-to-spectrum table scalar precision is not float32 in " + path);
    }
    if (resolution != RGBToSpectrumResolution || scaleCount != RGBToSpectrumResolution) {
      throw std::runtime_error("RGB-to-spectrum table resolution mismatch in " + path);
    }
    if (maxComponentCases != RGBToSpectrumMaxComponentCount ||
        coefficientChannels != RGBToSpectrumCoefficientChannels ||
        coefficientCount != RGBToSpectrumCoefficientCount) {
      throw std::runtime_error("RGB-to-spectrum table dimensions mismatch in " + path);
    }
    if (!expectedColorSpaceId.empty() && colorSpaceId != expectedColorSpaceId) {
      throw std::runtime_error("RGB-to-spectrum table color-space id mismatch in " + path);
    }
    if (payloadBytes != RGBToSpectrumPayloadFloatCount * sizeof(float)) {
      throw std::runtime_error("RGB-to-spectrum table payload size mismatch in " + path);
    }
    if (bytes.size() != offset + payloadBytes) {
      throw std::runtime_error("RGB-to-spectrum table file size does not match header in " + path);
    }
    std::uint32_t actualChecksum = detail::Adler32(bytes, offset, payloadBytes);
    if (actualChecksum != expectedChecksum) {
      throw std::runtime_error("RGB-to-spectrum table checksum mismatch in " + path);
    }

    std::array<Float, RGBToSpectrumResolution> zNodes{};
    std::size_t payloadOffset = offset;
    for (int i = 0; i < RGBToSpectrumResolution; ++i) {
      zNodes[static_cast<std::size_t>(i)] = detail::ReadLittleEndianFloat32(bytes, payloadOffset);
      payloadOffset += sizeof(float);
    }

    std::vector<Float> coefficients(RGBToSpectrumCoefficientCount);
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
      coefficients[i] = detail::ReadLittleEndianFloat32(bytes, payloadOffset);
      payloadOffset += sizeof(float);
    }

    return RGBToSpectrumTable(std::move(colorSpaceId), zNodes, std::move(coefficients));
  }

  RGBSigmoidPolynomial operator()(RGB rgb) const {
    if (!rgb.IsFinite() || rgb.r < 0 || rgb.g < 0 || rgb.b < 0 ||
        rgb.r > 1 || rgb.g > 1 || rgb.b > 1) {
      throw std::invalid_argument("RGB-to-spectrum table lookup requires finite RGB values in [0, 1]");
    }

    if (rgb.r == rgb.g && rgb.g == rgb.b) {
      Float c2 = (rgb.r - static_cast<Float>(0.5)) /
                 std::sqrt(rgb.r * (static_cast<Float>(1) - rgb.r));
      return RGBSigmoidPolynomial(0, 0, c2);
    }

    int maxComponent =
      (rgb[0] > rgb[1]) ? ((rgb[0] > rgb[2]) ? 0 : 2) : ((rgb[1] > rgb[2]) ? 1 : 2);
    Float z = rgb[maxComponent];
    Float x = rgb[(maxComponent + 1) % 3] * (RGBToSpectrumResolution - 1) / z;
    Float y = rgb[(maxComponent + 2) % 3] * (RGBToSpectrumResolution - 1) / z;

    int xi = std::min(static_cast<int>(x), RGBToSpectrumResolution - 2);
    int yi = std::min(static_cast<int>(y), RGBToSpectrumResolution - 2);
    int zi = detail::FindInterval(
      RGBToSpectrumResolution,
      [&](int i) { return zNodes_[static_cast<std::size_t>(i)] < z; }
    );
    Float dx = x - xi;
    Float dy = y - yi;
    Float dz = (z - zNodes_[static_cast<std::size_t>(zi)]) /
               (zNodes_[static_cast<std::size_t>(zi + 1)] - zNodes_[static_cast<std::size_t>(zi)]);

    std::array<Float, 3> c{};
    for (int i = 0; i < 3; ++i) {
      auto co = [&](int ddx, int ddy, int ddz) {
        return Coefficient(maxComponent, zi + ddz, yi + ddy, xi + ddx, i);
      };
      c[static_cast<std::size_t>(i)] = Lerp(
        dz,
        Lerp(dy, Lerp(dx, co(0, 0, 0), co(1, 0, 0)), Lerp(dx, co(0, 1, 0), co(1, 1, 0))),
        Lerp(dy, Lerp(dx, co(0, 0, 1), co(1, 0, 1)), Lerp(dx, co(0, 1, 1), co(1, 1, 1)))
      );
    }

    return RGBSigmoidPolynomial(c[0], c[1], c[2]);
  }

  Float ScaleNode(int index) const {
    return zNodes_.at(static_cast<std::size_t>(index));
  }

  Float Coefficient(int maxComponent, int z, int y, int x, int channel) const {
    return coefficients_.at(CoefficientIndex(maxComponent, z, y, x, channel));
  }

  const std::string& ColorSpaceId() const { return colorSpaceId_; }

private:
  static std::size_t CoefficientIndex(int maxComponent, int z, int y, int x, int channel) {
    return (((static_cast<std::size_t>(maxComponent) * RGBToSpectrumResolution + z) *
               RGBToSpectrumResolution + y) *
              RGBToSpectrumResolution + x) *
             RGBToSpectrumCoefficientChannels +
           channel;
  }

  std::string colorSpaceId_ = "sRGB";
  std::array<Float, RGBToSpectrumResolution> zNodes_{};
  std::vector<Float> coefficients_;
};

namespace detail {

inline std::mutex& RGBToSpectrumTableCacheMutex() {
  static std::mutex mutex;
  return mutex;
}

inline std::map<std::string, std::shared_ptr<const RGBToSpectrumTable>>& RGBToSpectrumTableCache() {
  static std::map<std::string, std::shared_ptr<const RGBToSpectrumTable>> cache;
  return cache;
}

inline RGBToSpectrumTableLoadStats& MutableRGBToSpectrumTableLoadStats() {
  static RGBToSpectrumTableLoadStats stats;
  return stats;
}

inline std::string RGBToSpectrumTableCacheKey(
  const std::string& path,
  const std::string& expectedColorSpaceId
) {
  return path + "\n" + expectedColorSpaceId;
}

} // namespace detail

inline RGBToSpectrumTableLoadStats GetRGBToSpectrumTableLoadStats() {
  std::lock_guard<std::mutex> lock(detail::RGBToSpectrumTableCacheMutex());
  return detail::MutableRGBToSpectrumTableLoadStats();
}

inline void ResetRGBToSpectrumTableLoadStats() {
  std::lock_guard<std::mutex> lock(detail::RGBToSpectrumTableCacheMutex());
  detail::MutableRGBToSpectrumTableLoadStats() = RGBToSpectrumTableLoadStats();
}

inline std::shared_ptr<const RGBToSpectrumTable> LoadRGBToSpectrumTable(
  const std::string& path,
  const std::string& expectedColorSpaceId = "sRGB"
) {
  std::lock_guard<std::mutex> lock(detail::RGBToSpectrumTableCacheMutex());
  std::map<std::string, std::shared_ptr<const RGBToSpectrumTable>>& cache =
    detail::RGBToSpectrumTableCache();
  RGBToSpectrumTableLoadStats& stats = detail::MutableRGBToSpectrumTableLoadStats();
  std::string cacheKey = detail::RGBToSpectrumTableCacheKey(path, expectedColorSpaceId);

  ++stats.lookups;
  auto found = cache.find(cacheKey);
  if (found != cache.end()) {
    ++stats.cacheHits;
    return found->second;
  }

  ++stats.fileLoads;
  std::shared_ptr<const RGBToSpectrumTable> table =
    std::make_shared<RGBToSpectrumTable>(RGBToSpectrumTable::LoadFromFile(path, expectedColorSpaceId));
  cache[cacheKey] = table;
  return table;
}

inline std::string RGBToSpectrumTableColorSpaceId(const std::string& colorSpaceName) {
  std::string canonical = RGBColorSpace::CanonicalName(colorSpaceName);
  if (canonical == "Rec.2020") {
    return "Rec2020";
  }
  return canonical;
}

inline std::string RGBToSpectrumTableFilename(const std::string& colorSpaceName) {
  std::string canonical = RGBColorSpace::CanonicalName(colorSpaceName);
  if (canonical == "sRGB") {
    return "rgb-to-spectrum-srgb-v1.bin";
  }
  if (canonical == "DCI-P3") {
    return "rgb-to-spectrum-dci-p3-v1.bin";
  }
  if (canonical == "Rec.2020") {
    return "rgb-to-spectrum-rec2020-v1.bin";
  }
  return "rgb-to-spectrum-aces2065-1-v1.bin";
}

inline std::shared_ptr<const RGBToSpectrumTable> LoadRGBToSpectrumTableForColorSpace(
  const std::string& colorSpaceName,
  const std::string& preferredAssetDirectory = ""
) {
  std::string canonical = RGBColorSpace::CanonicalName(colorSpaceName);
  std::string assetDirectory = FindSpectralAssetDirectory(preferredAssetDirectory);
  std::string path = detail::JoinPath(assetDirectory, RGBToSpectrumTableFilename(canonical));
  if (!detail::FileExists(path)) {
    throw std::runtime_error(
      "RGB-to-spectrum table for color space " + canonical +
      " is not packaged: " + path
    );
  }
  return LoadRGBToSpectrumTable(path, RGBToSpectrumTableColorSpaceId(canonical));
}

inline std::shared_ptr<const RGBToSpectrumTable> LoadSRGBRGBToSpectrumTable(
  const std::string& preferredAssetDirectory = ""
) {
  return LoadRGBToSpectrumTableForColorSpace("sRGB", preferredAssetDirectory);
}

inline RGBSigmoidPolynomial RGBColorSpace::ToRGBCoeffs(const RGB& rgb) const {
  if (!rgbToSpectrumTable) {
    throw std::invalid_argument("RGBColorSpace is missing an RGB-to-spectrum table");
  }
  if (!rgb.IsFinite() || rgb.r < 0 || rgb.g < 0 || rgb.b < 0) {
    throw std::invalid_argument("RGBColorSpace spectrum reconstruction requires finite non-negative RGB values");
  }
  return (*rgbToSpectrumTable)(ClampZero(rgb));
}

inline RGBColorSpace LoadRGBColorSpace(
  const std::string& colorSpaceName,
  const std::string& preferredAssetDirectory
);

inline RGBColorSpace LoadSRGBColorSpace(const std::string& preferredAssetDirectory = "") {
  return LoadRGBColorSpace("sRGB", preferredAssetDirectory);
}

inline RGBColorSpace LoadRGBColorSpace(
  const std::string& colorSpaceName,
  const std::string& preferredAssetDirectory = ""
) {
  std::string assetDirectory = FindSpectralAssetDirectory(preferredAssetDirectory);
  NamedSpectrumRegistry registry = NamedSpectrumRegistry::LoadFromDirectory(assetDirectory);
  RGBColorSpace namedColorSpace = RGBColorSpace::Named(colorSpaceName);
  std::string canonical = RGBColorSpace::CanonicalName(colorSpaceName);
  const Spectrum& illuminantSpectrum =
    registry.GetOrThrow(canonical == "ACES2065-1" ? "illum-acesD60" : "stdillum-D65");
  const Spectrum& x = registry.GetOrThrow("cie-x");
  const Spectrum& y = registry.GetOrThrow("cie-y");
  const Spectrum& z = registry.GetOrThrow("cie-z");

  XYZ white = SpectrumToXYZ(illuminantSpectrum, x, y, z);
  RGBColorSpace colorSpace = RGBColorSpace::FromPrimaries(
    namedColorSpace.name,
    namedColorSpace.r,
    namedColorSpace.g,
    namedColorSpace.b,
    white,
    namedColorSpace.encoding
  );
  colorSpace.illuminant = std::make_shared<DenselySampledSpectrum>(
    DenselySampledSpectrum::SampleFunction(
      [&](Float lambda) { return illuminantSpectrum(lambda); },
      static_cast<int>(LambdaMin),
      static_cast<int>(LambdaMax)
    )
  );
  colorSpace.rgbToSpectrumTable = LoadRGBToSpectrumTableForColorSpace(canonical, assetDirectory);
  return colorSpace;
}

inline void ValidateFiniteNonNegativeRGB(const RGB& rgb, const char* context) {
  if (!rgb.IsFinite() || rgb.r < 0 || rgb.g < 0 || rgb.b < 0) {
    throw std::invalid_argument(std::string(context) + " requires finite non-negative RGB values");
  }
}

class RGBAlbedoSpectrum {
public:
  RGBAlbedoSpectrum() = default;

  RGBAlbedoSpectrum(const RGBColorSpace& colorSpace, const RGB& rgb) {
    ValidateFiniteNonNegativeRGB(rgb, "RGBAlbedoSpectrum");
    if (rgb.r > 1 || rgb.g > 1 || rgb.b > 1) {
      throw std::invalid_argument("RGBAlbedoSpectrum requires RGB values in [0, 1]");
    }
    rsp_ = colorSpace.ToRGBCoeffs(rgb);
  }

  Float operator()(Float lambda) const { return rsp_(lambda); }

  SampledSpectrum Sample(const SampledWavelengths& lambda) const {
    SampledSpectrum spectrum;
    for (int i = 0; i < NSpectrumSamples; ++i) {
      spectrum[i] = rsp_(lambda[i]);
    }
    return spectrum;
  }

  Float MaxValue() const { return rsp_.MaxValue(); }

  std::string ToString() const {
    std::ostringstream out;
    out << "[ RGBAlbedoSpectrum rsp: " << rsp_.ToString() << " ]";
    return out.str();
  }

private:
  RGBSigmoidPolynomial rsp_;
};

class RGBUnboundedSpectrum {
public:
  RGBUnboundedSpectrum() = default;

  RGBUnboundedSpectrum(const RGBColorSpace& colorSpace, const RGB& rgb) {
    ValidateFiniteNonNegativeRGB(rgb, "RGBUnboundedSpectrum");
    Float maxComponent = rgb.MaxComponentValue();
    scale_ = static_cast<Float>(2) * maxComponent;
    rsp_ = colorSpace.ToRGBCoeffs(scale_ != 0 ? rgb / scale_ : RGB(0, 0, 0));
  }

  Float operator()(Float lambda) const { return scale_ * rsp_(lambda); }

  SampledSpectrum Sample(const SampledWavelengths& lambda) const {
    SampledSpectrum spectrum;
    for (int i = 0; i < NSpectrumSamples; ++i) {
      spectrum[i] = scale_ * rsp_(lambda[i]);
    }
    return spectrum;
  }

  Float MaxValue() const { return scale_ * rsp_.MaxValue(); }
  Float Scale() const { return scale_; }

  std::string ToString() const {
    std::ostringstream out;
    out << "[ RGBUnboundedSpectrum scale: " << scale_
        << " rsp: " << rsp_.ToString() << " ]";
    return out.str();
  }

private:
  Float scale_ = 0;
  RGBSigmoidPolynomial rsp_;
};

class RGBIlluminantSpectrum {
public:
  RGBIlluminantSpectrum() = default;

  RGBIlluminantSpectrum(const RGBColorSpace& colorSpace, const RGB& rgb)
    : illuminant_(colorSpace.illuminant) {
    ValidateFiniteNonNegativeRGB(rgb, "RGBIlluminantSpectrum");
    if (!illuminant_) {
      throw std::invalid_argument("RGBIlluminantSpectrum requires an RGBColorSpace illuminant");
    }
    Float maxComponent = rgb.MaxComponentValue();
    scale_ = static_cast<Float>(2) * maxComponent;
    rsp_ = colorSpace.ToRGBCoeffs(scale_ != 0 ? rgb / scale_ : RGB(0, 0, 0));
  }

  Float operator()(Float lambda) const {
    if (!illuminant_) {
      return 0;
    }
    return scale_ * rsp_(lambda) * (*illuminant_)(lambda);
  }

  SampledSpectrum Sample(const SampledWavelengths& lambda) const {
    if (!illuminant_) {
      return SampledSpectrum(0);
    }
    SampledSpectrum spectrum;
    for (int i = 0; i < NSpectrumSamples; ++i) {
      spectrum[i] = scale_ * rsp_(lambda[i]) * (*illuminant_)(lambda[i]);
    }
    return spectrum;
  }

  Float MaxValue() const {
    if (!illuminant_) {
      return 0;
    }
    return scale_ * rsp_.MaxValue() * illuminant_->MaxValue();
  }

  const DenselySampledSpectrum* Illuminant() const { return illuminant_.get(); }
  Float Scale() const { return scale_; }

  std::string ToString() const {
    std::ostringstream out;
    out << "[ RGBIlluminantSpectrum scale: " << scale_
        << " rsp: " << rsp_.ToString() << " ]";
    return out.str();
  }

private:
  Float scale_ = 0;
  RGBSigmoidPolynomial rsp_;
  std::shared_ptr<const DenselySampledSpectrum> illuminant_;
};

template <typename SpectrumLike>
inline Float InnerProductSpectrumLike(const SpectrumLike& f, const Spectrum& g) {
  Float integral = 0;
  for (int lambda = static_cast<int>(LambdaMin); lambda <= static_cast<int>(LambdaMax); ++lambda) {
    Float lambdaValue = static_cast<Float>(lambda);
    integral += f(lambdaValue) * g(lambdaValue);
  }
  return integral;
}

template <typename SpectrumLike>
inline XYZ SpectrumLikeToXYZ(const SpectrumLike& spectrum, const Spectrum& x, const Spectrum& y, const Spectrum& z) {
  return XYZ(
    InnerProductSpectrumLike(spectrum, x) / CIEYIntegral,
    InnerProductSpectrumLike(spectrum, y) / CIEYIntegral,
    InnerProductSpectrumLike(spectrum, z) / CIEYIntegral
  );
}

inline XYZ SpectrumToXYZ(const RGBAlbedoSpectrum& spectrum, const Spectrum& x, const Spectrum& y, const Spectrum& z) {
  return SpectrumLikeToXYZ(spectrum, x, y, z);
}

inline XYZ SpectrumToXYZ(const RGBUnboundedSpectrum& spectrum, const Spectrum& x, const Spectrum& y, const Spectrum& z) {
  return SpectrumLikeToXYZ(spectrum, x, y, z);
}

inline XYZ SpectrumToXYZ(
  const RGBIlluminantSpectrum& spectrum,
  const Spectrum& x,
  const Spectrum& y,
  const Spectrum& z
) {
  return SpectrumLikeToXYZ(spectrum, x, y, z);
}

} // namespace base
} // namespace rayrender

#endif
