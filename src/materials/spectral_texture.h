#ifndef RAYRENDER_MATERIALS_SPECTRAL_TEXTURE_H
#define RAYRENDER_MATERIALS_SPECTRAL_TEXTURE_H

#include "../base/base.h"
#include "../math/perlin.h"
#include "../math/vectypes.h"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

struct hit_record;

namespace rayrender {
namespace materials {

enum class TextureValueType {
  Float,
  Spectrum
};

enum class TextureKind {
  Constant,
  Scale,
  Mix,
  Checker,
  Noise,
  Image
};

enum class TextureImageEncoding {
  Auto,
  Linear,
  SRGB
};

enum class TextureSemanticRole {
  Scalar,
  Albedo,
  Illuminant,
  Unbounded
};

enum class TextureWrapMode {
  Repeat,
  Clamp,
  Black
};

enum class TextureFilterMode {
  Nearest,
  Bilinear,
  EWA
};

enum class TextureScalarChannel {
  Luminance,
  Red,
  Green,
  Blue,
  Alpha
};

const char* TextureValueTypeName(TextureValueType valueType);
const char* TextureKindName(TextureKind kind);
const char* TextureImageEncodingName(TextureImageEncoding encoding);
const char* TextureSemanticRoleName(TextureSemanticRole role);
const char* TextureWrapModeName(TextureWrapMode wrap);
const char* TextureFilterModeName(TextureFilterMode filter);
const char* TextureScalarChannelName(TextureScalarChannel channel);

struct TextureRayDifferentials {
  vec3f dpdx;
  vec3f dpdy;
  Float dudx = 0;
  Float dvdx = 0;
  Float dudy = 0;
  Float dvdy = 0;
};

struct TextureEvalContext {
  point3f p;
  point2f uv;
  vec3f dpdx;
  vec3f dpdy;
  vec3f dpdu;
  vec3f dpdv;
  Float dudx = 0;
  Float dvdx = 0;
  Float dudy = 0;
  Float dvdy = 0;
  int faceIndex = -1;
};

TextureEvalContext TextureEvalContextFromHitRecord(const hit_record& hit);
TextureEvalContext TextureEvalContextFromHitRecord(
  const hit_record& hit,
  const TextureRayDifferentials& differentials
);

struct TextureFilterFootprint {
  Float dudx = 0;
  Float dvdx = 0;
  Float dudy = 0;
  Float dvdy = 0;
  Float width = 0;
};

TextureFilterFootprint TextureFilterFootprintFromContext(
  const TextureEvalContext& ctx,
  int imageWidth,
  int imageHeight
);

struct ImageCacheKey {
  std::string filename;
  TextureValueType valueType = TextureValueType::Spectrum;
  TextureImageEncoding encoding = TextureImageEncoding::Auto;
  std::string colorSpace = "sRGB";
  TextureSemanticRole role = TextureSemanticRole::Albedo;
  TextureWrapMode wrap = TextureWrapMode::Repeat;
  TextureFilterMode filter = TextureFilterMode::EWA;
  TextureScalarChannel scalarChannel = TextureScalarChannel::Luminance;
  Float maxAnisotropy = 8;

  ImageCacheKey() = default;
  ImageCacheKey(
    std::string filename,
    TextureValueType valueType,
    TextureImageEncoding encoding,
    std::string colorSpace,
    TextureSemanticRole role,
    TextureWrapMode wrap,
    TextureFilterMode filter,
    TextureScalarChannel scalarChannel = TextureScalarChannel::Luminance,
    Float maxAnisotropy = 8
  );

  static ImageCacheKey Color(
    std::string filename,
    TextureSemanticRole role,
    TextureImageEncoding encoding = TextureImageEncoding::Auto,
    std::string colorSpace = "sRGB",
    TextureWrapMode wrap = TextureWrapMode::Repeat,
    TextureFilterMode filter = TextureFilterMode::EWA,
    Float maxAnisotropy = 8
  );

  static ImageCacheKey Scalar(
    std::string filename,
    TextureScalarChannel channel = TextureScalarChannel::Luminance,
    TextureImageEncoding encoding = TextureImageEncoding::Linear,
    TextureWrapMode wrap = TextureWrapMode::Repeat,
    TextureFilterMode filter = TextureFilterMode::EWA,
    Float maxAnisotropy = 8
  );

  bool operator==(const ImageCacheKey& other) const;
  bool operator!=(const ImageCacheKey& other) const { return !(*this == other); }
  std::string ToString() const;
};

struct ImageCacheKeyHash {
  std::size_t operator()(const ImageCacheKey& key) const;
};

class ImageTextureData {
public:
  static std::shared_ptr<const ImageTextureData> FromInterleaved(
    ImageCacheKey key,
    int width,
    int height,
    int channels,
    std::vector<Float> values,
    const base::RGBColorSpace& inputColorSpace = base::RGBColorSpace::SRGB()
  );

  static std::shared_ptr<const ImageTextureData> LoadFromFile(
    const ImageCacheKey& key,
    const base::RGBColorSpace& inputColorSpace = base::RGBColorSpace::SRGB()
  );

  int Width() const { return width_; }
  int Height() const { return height_; }
  int Channels() const { return channels_; }
  TextureValueType ValueType() const { return key_.valueType; }
  const ImageCacheKey& Key() const { return key_; }
  bool WasDecodedToLinear() const { return decodedColorToLinear_; }

  base::RGB SampleRGB(const TextureEvalContext& ctx) const;
  Float SampleScalar(const TextureEvalContext& ctx) const;
  TextureFilterFootprint Footprint(const TextureEvalContext& ctx) const;

private:
  ImageTextureData(
    ImageCacheKey key,
    int width,
    int height,
    int channels,
    std::vector<Float> values,
    bool decodedColorToLinear
  );

  Float TexelChannel(int x, int y, int channel) const;
  Float TexelChannelWrapped(int x, int y, int channel) const;
  base::RGB TexelRGBWrapped(int x, int y) const;
  Float TexelScalarWrapped(int x, int y) const;
  base::RGB SampleRGBNearest(const TextureEvalContext& ctx) const;
  base::RGB SampleRGBBilinear(const TextureEvalContext& ctx) const;
  Float SampleScalarNearest(const TextureEvalContext& ctx) const;
  Float SampleScalarBilinear(const TextureEvalContext& ctx) const;
  int WrapIndex(int index, int size, bool& outside) const;

  ImageCacheKey key_;
  int width_ = 0;
  int height_ = 0;
  int channels_ = 0;
  std::vector<Float> texels_;
  bool decodedColorToLinear_ = false;
};

class SpectralImageCache {
public:
  using Loader = std::function<std::shared_ptr<const ImageTextureData>(const ImageCacheKey&)>;

  std::shared_ptr<const ImageTextureData> LookupOrLoad(
    const ImageCacheKey& key,
    const base::RGBColorSpace& inputColorSpace = base::RGBColorSpace::SRGB()
  );

  std::shared_ptr<const ImageTextureData> LookupOrCreate(
    const ImageCacheKey& key,
    const Loader& loader
  );

  void Store(std::shared_ptr<const ImageTextureData> image);
  std::size_t Size() const;

private:
  mutable std::mutex mutex_;
  std::unordered_map<ImageCacheKey, std::shared_ptr<const ImageTextureData>, ImageCacheKeyHash> images_;
};

class FloatTexture;
class SpectrumTexture;

using FloatTexturePtr = std::shared_ptr<const FloatTexture>;
using SpectrumTexturePtr = std::shared_ptr<const SpectrumTexture>;

struct ConstantFloatTexture {
  Float value = 0;
};

struct ScaleFloatTexture {
  FloatTexturePtr texture;
  Float scale = 1;
};

struct MixFloatTexture {
  FloatTexturePtr tex1;
  FloatTexturePtr tex2;
  std::variant<Float, FloatTexturePtr> amount = static_cast<Float>(0.5);
};

struct CheckerFloatTexture {
  FloatTexturePtr tex1;
  FloatTexturePtr tex2;
  int dimension = 2;
  Float uscale = 1;
  Float vscale = 1;
};

struct NoiseFloatTexture {
  Float value1 = 0;
  Float value2 = 1;
  Float scale = 1;
  Float phase = 0;
  Float intensity = 1;
  perlin noise;
};

struct ImageFloatTexture {
  std::shared_ptr<const ImageTextureData> image;
  Float scale = 1;
  bool invert = false;
};

class FloatTexture {
public:
  FloatTexture();

  static FloatTexture Constant(Float value);
  static FloatTexture Scale(FloatTexture texture, Float scale);
  static FloatTexture Mix(FloatTexture tex1, FloatTexture tex2, Float amount);
  static FloatTexture Mix(FloatTexture tex1, FloatTexture tex2, FloatTexture amount);
  static FloatTexture Checker(
    FloatTexture tex1,
    FloatTexture tex2,
    int dimension = 2,
    Float uscale = 1,
    Float vscale = 1
  );
  static FloatTexture Noise(
    Float value1 = 0,
    Float value2 = 1,
    Float scale = 1,
    Float phase = 0,
    Float intensity = 1
  );
  static FloatTexture Image(std::shared_ptr<const ImageTextureData> image, Float scale = 1, bool invert = false);

  Float Evaluate(const TextureEvalContext& ctx) const;
  TextureKind Kind() const;

private:
  using Variant = std::variant<
    ConstantFloatTexture,
    ScaleFloatTexture,
    MixFloatTexture,
    CheckerFloatTexture,
    NoiseFloatTexture,
    ImageFloatTexture
  >;

  explicit FloatTexture(Variant data);

  friend class UniversalTextureEvaluator;
  friend struct TextureDebugEvaluation;
  Variant data_;
};

struct ConstantSpectrumTexture {
  base::Spectrum spectrum;
};

struct ScaleSpectrumTexture {
  SpectrumTexturePtr texture;
  Float scale = 1;
};

struct MixSpectrumTexture {
  SpectrumTexturePtr tex1;
  SpectrumTexturePtr tex2;
  std::variant<Float, FloatTexturePtr> amount = static_cast<Float>(0.5);
};

struct CheckerSpectrumTexture {
  SpectrumTexturePtr tex1;
  SpectrumTexturePtr tex2;
  int dimension = 2;
  Float uscale = 1;
  Float vscale = 1;
};

struct NoiseSpectrumTexture {
  SpectrumTexturePtr tex1;
  SpectrumTexturePtr tex2;
  Float scale = 1;
  Float phase = 0;
  Float intensity = 1;
  perlin noise;
};

struct ImageSpectrumTexture {
  std::shared_ptr<const ImageTextureData> image;
  TextureSemanticRole role = TextureSemanticRole::Albedo;
  base::RGBColorSpace colorSpace = base::RGBColorSpace::SRGB();
  Float scale = 1;
};

class SpectrumTexture {
public:
  SpectrumTexture();

  static SpectrumTexture Constant(base::Spectrum spectrum);
  static SpectrumTexture Constant(Float value);
  static SpectrumTexture Scale(SpectrumTexture texture, Float scale);
  static SpectrumTexture Mix(SpectrumTexture tex1, SpectrumTexture tex2, Float amount);
  static SpectrumTexture Mix(SpectrumTexture tex1, SpectrumTexture tex2, FloatTexture amount);
  static SpectrumTexture Checker(
    SpectrumTexture tex1,
    SpectrumTexture tex2,
    int dimension = 2,
    Float uscale = 1,
    Float vscale = 1
  );
  static SpectrumTexture Noise(
    SpectrumTexture tex1,
    SpectrumTexture tex2,
    Float scale = 1,
    Float phase = 0,
    Float intensity = 1
  );
  static SpectrumTexture Image(
    std::shared_ptr<const ImageTextureData> image,
    TextureSemanticRole role,
    base::RGBColorSpace colorSpace,
    Float scale = 1
  );

  base::SampledSpectrum Evaluate(
    const TextureEvalContext& ctx,
    const base::SampledWavelengths& lambda
  ) const;
  TextureKind Kind() const;

private:
  using Variant = std::variant<
    ConstantSpectrumTexture,
    ScaleSpectrumTexture,
    MixSpectrumTexture,
    CheckerSpectrumTexture,
    NoiseSpectrumTexture,
    ImageSpectrumTexture
  >;

  explicit SpectrumTexture(Variant data);

  friend class UniversalTextureEvaluator;
  friend struct TextureDebugEvaluation;
  Variant data_;
};

class UniversalTextureEvaluator {
public:
  Float operator()(const FloatTexture& texture, const TextureEvalContext& ctx) const;

  base::SampledSpectrum operator()(
    const SpectrumTexture& texture,
    const TextureEvalContext& ctx,
    const base::SampledWavelengths& lambda
  ) const;
};

class TextureHandleTable {
public:
  base::TextureHandle Add(FloatTexture texture);
  base::TextureHandle Add(SpectrumTexture texture);

  TextureValueType ValueType(base::TextureHandle handle) const;
  const FloatTexture& GetFloat(base::TextureHandle handle) const;
  const SpectrumTexture& GetSpectrum(base::TextureHandle handle) const;
  std::size_t Size() const { return entries_.size(); }

private:
  struct Entry {
    TextureValueType valueType = TextureValueType::Float;
    base::TextureHandle::GenerationType generation = 1;
    std::variant<FloatTexture, SpectrumTexture> texture;
  };

  const Entry& GetEntry(base::TextureHandle handle) const;

  std::vector<Entry> entries_;
};

struct TextureDebugEvaluation {
  TextureValueType valueType = TextureValueType::Float;
  TextureKind kind = TextureKind::Constant;
  Float scalar = 0;
  base::RGB linearRGB;
  bool hasLinearRGB = false;
  bool decodedColorToLinear = false;
  base::SampledSpectrum spectrum;
  TextureFilterFootprint footprint;
};

TextureDebugEvaluation DebugEvaluateFloatTexture(
  const FloatTexture& texture,
  const TextureEvalContext& ctx
);

TextureDebugEvaluation DebugEvaluateSpectrumTexture(
  const SpectrumTexture& texture,
  const TextureEvalContext& ctx,
  const base::SampledWavelengths& lambda
);

} // namespace materials
} // namespace rayrender

#endif
