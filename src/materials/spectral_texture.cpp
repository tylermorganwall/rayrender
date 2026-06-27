#include "../materials/spectral_texture.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace rayrender {
namespace materials {
namespace {

template <typename T>
void HashCombine(std::size_t& seed, const T& value) {
  seed ^= std::hash<T>{}(value) + 0x9e3779b97f4a7c15ull + (seed << 6) + (seed >> 2);
}

std::string LowercaseExtension(const std::string& filename) {
  std::size_t dot = filename.find_last_of('.');
  if (dot == std::string::npos) {
    return "";
  }
  std::string extension = filename.substr(dot);
  std::transform(extension.begin(), extension.end(), extension.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return extension;
}

bool IsHighDynamicRangeFile(const std::string& filename) {
  std::string extension = LowercaseExtension(filename);
  return extension == ".hdr" || extension == ".exr";
}

TextureImageEncoding ResolveEncoding(const ImageCacheKey& key) {
  if (key.encoding != TextureImageEncoding::Auto) {
    return key.encoding;
  }
  if (key.valueType == TextureValueType::Float) {
    return TextureImageEncoding::Linear;
  }
  return IsHighDynamicRangeFile(key.filename) ? TextureImageEncoding::Linear : TextureImageEncoding::SRGB;
}

base::RGB DecodeRGB(const base::RGB& rgb, TextureImageEncoding encoding) {
  if (encoding == TextureImageEncoding::SRGB) {
    return base::RGBColorEncoding::SRGB().Decode(rgb);
  }
  return rgb;
}

bool NearlyNeutral(const base::RGB& rgb) {
  constexpr Float tolerance = static_cast<Float>(1e-6);
  return std::fabs(rgb.r - rgb.g) <= tolerance && std::fabs(rgb.r - rgb.b) <= tolerance;
}

bool IsOdd(int value) {
  return std::abs(value % 2) == 1;
}

bool UseSecondCheckerTexel(const TextureEvalContext& ctx, int dimension, Float uscale, Float vscale) {
  if (dimension == 2) {
    int u = static_cast<int>(std::floor(ctx.uv[0] * uscale));
    int v = static_cast<int>(std::floor(ctx.uv[1] * vscale));
    return IsOdd(u + v);
  }

  int x = static_cast<int>(std::floor(ctx.p.xyz.x * uscale));
  int y = static_cast<int>(std::floor(ctx.p.xyz.y * vscale));
  int z = static_cast<int>(std::floor(ctx.p.xyz.z * vscale));
  return IsOdd(x + y + z);
}

Float NoiseWeight(const perlin& noise, const TextureEvalContext& ctx, Float scale, Float phase, Float intensity) {
  return static_cast<Float>(0.5) *
         (static_cast<Float>(1) +
          std::sin(scale * ctx.p.xyz.y + intensity * noise.turb(scale * ctx.p) + phase));
}

Float TextureAmount(
  const std::variant<Float, FloatTexturePtr>& amount,
  const TextureEvalContext& ctx,
  const UniversalTextureEvaluator& evaluator
) {
  if (std::holds_alternative<Float>(amount)) {
    return std::get<Float>(amount);
  }
  const FloatTexturePtr& texture = std::get<FloatTexturePtr>(amount);
  if (!texture) {
    throw std::runtime_error("Mix texture amount is missing");
  }
  return evaluator(*texture, ctx);
}

base::SampledSpectrum ReconstructLinearRGB(
  const base::RGB& rgb,
  TextureSemanticRole role,
  const base::RGBColorSpace& colorSpace,
  const base::SampledWavelengths& lambda
) {
  if (!rgb.IsFinite() || rgb.r < 0 || rgb.g < 0 || rgb.b < 0) {
    throw std::invalid_argument("RGB texture reconstruction requires finite non-negative linear RGB values");
  }

  if (NearlyNeutral(rgb) && !colorSpace.HasSpectralReconstruction()) {
    Float value = rgb.r;
    if (role == TextureSemanticRole::Albedo) {
      if (value > 1) {
        throw std::invalid_argument("Albedo RGB texture reconstruction requires values in [0, 1]");
      }
      return base::SampledSpectrum(value);
    }
    if (role == TextureSemanticRole::Illuminant) {
      if (!colorSpace.illuminant) {
        throw std::invalid_argument("Illuminant RGB texture reconstruction requires a color-space illuminant");
      }
      return colorSpace.illuminant->Sample(lambda) * value;
    }
    if (role == TextureSemanticRole::Unbounded) {
      return base::SampledSpectrum(value);
    }
  }

  if (role == TextureSemanticRole::Albedo) {
    return base::RGBAlbedoSpectrum(colorSpace, rgb).Sample(lambda);
  }
  if (role == TextureSemanticRole::Illuminant) {
    return base::RGBIlluminantSpectrum(colorSpace, rgb).Sample(lambda);
  }
  if (role == TextureSemanticRole::Unbounded) {
    return base::RGBUnboundedSpectrum(colorSpace, rgb).Sample(lambda);
  }
  throw std::invalid_argument("Scalar image data cannot be reconstructed as a spectrum");
}

Float ValidateFiniteScale(Float value, const char* context) {
  if (!std::isfinite(value)) {
    throw std::invalid_argument(std::string(context) + " must be finite");
  }
  return value;
}

FloatTexturePtr MakeFloatTexturePtr(FloatTexture texture) {
  return std::make_shared<const FloatTexture>(std::move(texture));
}

SpectrumTexturePtr MakeSpectrumTexturePtr(SpectrumTexture texture) {
  return std::make_shared<const SpectrumTexture>(std::move(texture));
}

} // namespace

const char* TextureValueTypeName(TextureValueType valueType) {
  switch (valueType) {
  case TextureValueType::Float:
    return "float";
  case TextureValueType::Spectrum:
    return "spectrum";
  }
  return "unknown";
}

const char* TextureKindName(TextureKind kind) {
  switch (kind) {
  case TextureKind::Constant:
    return "constant";
  case TextureKind::Scale:
    return "scale";
  case TextureKind::Mix:
    return "mix";
  case TextureKind::Checker:
    return "checker";
  case TextureKind::Noise:
    return "noise";
  case TextureKind::Image:
    return "image";
  }
  return "unknown";
}

const char* TextureImageEncodingName(TextureImageEncoding encoding) {
  switch (encoding) {
  case TextureImageEncoding::Auto:
    return "auto";
  case TextureImageEncoding::Linear:
    return "linear";
  case TextureImageEncoding::SRGB:
    return "srgb";
  }
  return "unknown";
}

const char* TextureSemanticRoleName(TextureSemanticRole role) {
  switch (role) {
  case TextureSemanticRole::Scalar:
    return "scalar";
  case TextureSemanticRole::Albedo:
    return "albedo";
  case TextureSemanticRole::Illuminant:
    return "illuminant";
  case TextureSemanticRole::Unbounded:
    return "unbounded";
  }
  return "unknown";
}

const char* TextureWrapModeName(TextureWrapMode wrap) {
  switch (wrap) {
  case TextureWrapMode::Repeat:
    return "repeat";
  case TextureWrapMode::Clamp:
    return "clamp";
  case TextureWrapMode::Black:
    return "black";
  }
  return "unknown";
}

const char* TextureFilterModeName(TextureFilterMode filter) {
  switch (filter) {
  case TextureFilterMode::Nearest:
    return "nearest";
  case TextureFilterMode::Bilinear:
    return "bilinear";
  case TextureFilterMode::EWA:
    return "ewa";
  }
  return "unknown";
}

const char* TextureScalarChannelName(TextureScalarChannel channel) {
  switch (channel) {
  case TextureScalarChannel::Luminance:
    return "luminance";
  case TextureScalarChannel::Red:
    return "red";
  case TextureScalarChannel::Green:
    return "green";
  case TextureScalarChannel::Blue:
    return "blue";
  case TextureScalarChannel::Alpha:
    return "alpha";
  }
  return "unknown";
}

TextureFilterFootprint TextureFilterFootprintFromContext(
  const TextureEvalContext& ctx,
  int imageWidth,
  int imageHeight
) {
  TextureFilterFootprint footprint;
  footprint.dudx = ctx.dudx * imageWidth;
  footprint.dvdx = ctx.dvdx * imageHeight;
  footprint.dudy = ctx.dudy * imageWidth;
  footprint.dvdy = ctx.dvdy * imageHeight;
  Float xWidth = std::sqrt(footprint.dudx * footprint.dudx + footprint.dvdx * footprint.dvdx);
  Float yWidth = std::sqrt(footprint.dudy * footprint.dudy + footprint.dvdy * footprint.dvdy);
  footprint.width = std::max(xWidth, yWidth);
  return footprint;
}

ImageCacheKey::ImageCacheKey(
  std::string filename,
  TextureValueType valueType,
  TextureImageEncoding encoding,
  std::string colorSpace,
  TextureSemanticRole role,
  TextureWrapMode wrap,
  TextureFilterMode filter,
  TextureScalarChannel scalarChannel,
  Float maxAnisotropy
)
  : filename(std::move(filename)),
    valueType(valueType),
    encoding(encoding),
    colorSpace(std::move(colorSpace)),
    role(role),
    wrap(wrap),
    filter(filter),
    scalarChannel(scalarChannel),
    maxAnisotropy(maxAnisotropy) {
  if (this->filename.empty()) {
    throw std::invalid_argument("Image cache key filename must not be empty");
  }
  if (!std::isfinite(maxAnisotropy) || maxAnisotropy <= 0) {
    throw std::invalid_argument("Image cache key maxAnisotropy must be finite and positive");
  }
  if (valueType == TextureValueType::Float) {
    this->role = TextureSemanticRole::Scalar;
    this->colorSpace.clear();
  } else if (role == TextureSemanticRole::Scalar) {
    throw std::invalid_argument("Spectrum image cache keys require a non-scalar semantic role");
  }
}

ImageCacheKey ImageCacheKey::Color(
  std::string filename,
  TextureSemanticRole role,
  TextureImageEncoding encoding,
  std::string colorSpace,
  TextureWrapMode wrap,
  TextureFilterMode filter,
  Float maxAnisotropy
) {
  return ImageCacheKey(
    std::move(filename),
    TextureValueType::Spectrum,
    encoding,
    std::move(colorSpace),
    role,
    wrap,
    filter,
    TextureScalarChannel::Luminance,
    maxAnisotropy
  );
}

ImageCacheKey ImageCacheKey::Scalar(
  std::string filename,
  TextureScalarChannel channel,
  TextureImageEncoding encoding,
  TextureWrapMode wrap,
  TextureFilterMode filter,
  Float maxAnisotropy
) {
  return ImageCacheKey(
    std::move(filename),
    TextureValueType::Float,
    encoding,
    "",
    TextureSemanticRole::Scalar,
    wrap,
    filter,
    channel,
    maxAnisotropy
  );
}

bool ImageCacheKey::operator==(const ImageCacheKey& other) const {
  return filename == other.filename && valueType == other.valueType && encoding == other.encoding &&
         colorSpace == other.colorSpace && role == other.role && wrap == other.wrap &&
         filter == other.filter && scalarChannel == other.scalarChannel &&
         maxAnisotropy == other.maxAnisotropy;
}

std::string ImageCacheKey::ToString() const {
  std::ostringstream out;
  out << filename
      << "|type=" << TextureValueTypeName(valueType)
      << "|encoding=" << TextureImageEncodingName(encoding)
      << "|color_space=" << colorSpace
      << "|role=" << TextureSemanticRoleName(role)
      << "|wrap=" << TextureWrapModeName(wrap)
      << "|filter=" << TextureFilterModeName(filter)
      << "|channel=" << TextureScalarChannelName(scalarChannel)
      << "|max_anisotropy=" << maxAnisotropy;
  return out.str();
}

std::size_t ImageCacheKeyHash::operator()(const ImageCacheKey& key) const {
  std::size_t seed = 0;
  HashCombine(seed, key.filename);
  HashCombine(seed, static_cast<int>(key.valueType));
  HashCombine(seed, static_cast<int>(key.encoding));
  HashCombine(seed, key.colorSpace);
  HashCombine(seed, static_cast<int>(key.role));
  HashCombine(seed, static_cast<int>(key.wrap));
  HashCombine(seed, static_cast<int>(key.filter));
  HashCombine(seed, static_cast<int>(key.scalarChannel));
  HashCombine(seed, key.maxAnisotropy);
  return seed;
}

std::shared_ptr<const ImageTextureData> ImageTextureData::FromInterleaved(
  ImageCacheKey key,
  int width,
  int height,
  int channels,
  std::vector<Float> values,
  const base::RGBColorSpace&
) {
  if (width <= 0 || height <= 0 || channels <= 0) {
    throw std::invalid_argument("ImageTextureData dimensions must be positive");
  }
  const std::size_t expected =
    static_cast<std::size_t>(width) * static_cast<std::size_t>(height) * static_cast<std::size_t>(channels);
  if (values.size() != expected) {
    throw std::invalid_argument("ImageTextureData value count does not match dimensions");
  }

  TextureImageEncoding resolvedEncoding = ResolveEncoding(key);
  bool decodeColor = key.valueType == TextureValueType::Spectrum &&
                     resolvedEncoding == TextureImageEncoding::SRGB;
  if (decodeColor) {
    for (int y = 0; y < height; ++y) {
      for (int x = 0; x < width; ++x) {
        std::size_t index =
          (static_cast<std::size_t>(y) * static_cast<std::size_t>(width) + static_cast<std::size_t>(x)) *
          static_cast<std::size_t>(channels);
        Float r = values[index];
        Float g = channels > 1 ? values[index + 1] : r;
        Float b = channels > 2 ? values[index + 2] : r;
        base::RGB decoded = DecodeRGB(base::RGB(r, g, b), TextureImageEncoding::SRGB);
        values[index] = decoded.r;
        if (channels > 1) {
          values[index + 1] = decoded.g;
        }
        if (channels > 2) {
          values[index + 2] = decoded.b;
        }
      }
    }
  }

  for (Float value : values) {
    if (!std::isfinite(value)) {
      throw std::invalid_argument("ImageTextureData values must be finite");
    }
  }

  return std::shared_ptr<const ImageTextureData>(
    new ImageTextureData(std::move(key), width, height, channels, std::move(values), decodeColor)
  );
}

ImageTextureData::ImageTextureData(
  ImageCacheKey key,
  int width,
  int height,
  int channels,
  std::vector<Float> values,
  bool decodedColorToLinear
)
  : key_(std::move(key)),
    width_(width),
    height_(height),
    channels_(channels),
    texels_(std::move(values)),
    decodedColorToLinear_(decodedColorToLinear) {}

base::RGB ImageTextureData::SampleRGB(const TextureEvalContext& ctx) const {
  if (key_.valueType != TextureValueType::Spectrum) {
    throw std::invalid_argument("Scalar image data cannot be sampled as RGB");
  }
  if (key_.filter == TextureFilterMode::Nearest) {
    return SampleRGBNearest(ctx);
  }
  return SampleRGBBilinear(ctx);
}

Float ImageTextureData::SampleScalar(const TextureEvalContext& ctx) const {
  if (key_.valueType != TextureValueType::Float) {
    throw std::invalid_argument("Color image data cannot be sampled as scalar");
  }
  if (key_.filter == TextureFilterMode::Nearest) {
    return SampleScalarNearest(ctx);
  }
  return SampleScalarBilinear(ctx);
}

TextureFilterFootprint ImageTextureData::Footprint(const TextureEvalContext& ctx) const {
  return TextureFilterFootprintFromContext(ctx, width_, height_);
}

Float ImageTextureData::TexelChannel(int x, int y, int channel) const {
  std::size_t index =
    (static_cast<std::size_t>(y) * static_cast<std::size_t>(width_) + static_cast<std::size_t>(x)) *
    static_cast<std::size_t>(channels_);
  if (channel < channels_) {
    return texels_[index + static_cast<std::size_t>(channel)];
  }
  if (channel == 3) {
    return static_cast<Float>(1);
  }
  if (channels_ == 1 && channel < 3) {
    return texels_[index];
  }
  return static_cast<Float>(0);
}

int ImageTextureData::WrapIndex(int index, int size, bool& outside) const {
  if (index >= 0 && index < size) {
    return index;
  }
  if (key_.wrap == TextureWrapMode::Black) {
    outside = true;
    return 0;
  }
  if (key_.wrap == TextureWrapMode::Clamp) {
    return std::clamp(index, 0, size - 1);
  }
  int wrapped = index % size;
  return wrapped < 0 ? wrapped + size : wrapped;
}

Float ImageTextureData::TexelChannelWrapped(int x, int y, int channel) const {
  bool outside = false;
  int wrappedX = WrapIndex(x, width_, outside);
  int wrappedY = WrapIndex(y, height_, outside);
  if (outside) {
    return static_cast<Float>(0);
  }
  return TexelChannel(wrappedX, wrappedY, channel);
}

base::RGB ImageTextureData::TexelRGBWrapped(int x, int y) const {
  return base::RGB(
    TexelChannelWrapped(x, y, 0),
    TexelChannelWrapped(x, y, 1),
    TexelChannelWrapped(x, y, 2)
  );
}

Float ImageTextureData::TexelScalarWrapped(int x, int y) const {
  switch (key_.scalarChannel) {
  case TextureScalarChannel::Luminance: {
    base::RGB rgb = TexelRGBWrapped(x, y);
    return static_cast<Float>(0.2126) * rgb.r +
           static_cast<Float>(0.7152) * rgb.g +
           static_cast<Float>(0.0722) * rgb.b;
  }
  case TextureScalarChannel::Red:
    return TexelChannelWrapped(x, y, 0);
  case TextureScalarChannel::Green:
    return TexelChannelWrapped(x, y, 1);
  case TextureScalarChannel::Blue:
    return TexelChannelWrapped(x, y, 2);
  case TextureScalarChannel::Alpha:
    return TexelChannelWrapped(x, y, 3);
  }
  return 0;
}

base::RGB ImageTextureData::SampleRGBNearest(const TextureEvalContext& ctx) const {
  int x = static_cast<int>(std::floor(ctx.uv[0] * width_));
  int y = static_cast<int>(std::floor((static_cast<Float>(1) - ctx.uv[1]) * height_));
  return TexelRGBWrapped(x, y);
}

base::RGB ImageTextureData::SampleRGBBilinear(const TextureEvalContext& ctx) const {
  Float s = ctx.uv[0] * width_ - static_cast<Float>(0.5);
  Float t = (static_cast<Float>(1) - ctx.uv[1]) * height_ - static_cast<Float>(0.5);
  int x0 = static_cast<int>(std::floor(s));
  int y0 = static_cast<int>(std::floor(t));
  Float dx = s - x0;
  Float dy = t - y0;

  base::RGB c00 = TexelRGBWrapped(x0, y0);
  base::RGB c10 = TexelRGBWrapped(x0 + 1, y0);
  base::RGB c01 = TexelRGBWrapped(x0, y0 + 1);
  base::RGB c11 = TexelRGBWrapped(x0 + 1, y0 + 1);
  base::RGB c0 = (static_cast<Float>(1) - dx) * c00 + dx * c10;
  base::RGB c1 = (static_cast<Float>(1) - dx) * c01 + dx * c11;
  return (static_cast<Float>(1) - dy) * c0 + dy * c1;
}

Float ImageTextureData::SampleScalarNearest(const TextureEvalContext& ctx) const {
  int x = static_cast<int>(std::floor(ctx.uv[0] * width_));
  int y = static_cast<int>(std::floor((static_cast<Float>(1) - ctx.uv[1]) * height_));
  return TexelScalarWrapped(x, y);
}

Float ImageTextureData::SampleScalarBilinear(const TextureEvalContext& ctx) const {
  Float s = ctx.uv[0] * width_ - static_cast<Float>(0.5);
  Float t = (static_cast<Float>(1) - ctx.uv[1]) * height_ - static_cast<Float>(0.5);
  int x0 = static_cast<int>(std::floor(s));
  int y0 = static_cast<int>(std::floor(t));
  Float dx = s - x0;
  Float dy = t - y0;

  Float c00 = TexelScalarWrapped(x0, y0);
  Float c10 = TexelScalarWrapped(x0 + 1, y0);
  Float c01 = TexelScalarWrapped(x0, y0 + 1);
  Float c11 = TexelScalarWrapped(x0 + 1, y0 + 1);
  Float c0 = (static_cast<Float>(1) - dx) * c00 + dx * c10;
  Float c1 = (static_cast<Float>(1) - dx) * c01 + dx * c11;
  return (static_cast<Float>(1) - dy) * c0 + dy * c1;
}

std::shared_ptr<const ImageTextureData> SpectralImageCache::LookupOrCreate(
  const ImageCacheKey& key,
  const Loader& loader
) {
  {
    std::lock_guard<std::mutex> lock(mutex_);
    auto found = images_.find(key);
    if (found != images_.end()) {
      return found->second;
    }
  }

  std::shared_ptr<const ImageTextureData> image = loader(key);
  if (!image) {
    throw std::runtime_error("Image cache loader returned null data");
  }
  if (image->Key() != key) {
    throw std::runtime_error("Image cache loader returned data for a different key");
  }

  std::lock_guard<std::mutex> lock(mutex_);
  auto inserted = images_.emplace(key, image);
  return inserted.first->second;
}

void SpectralImageCache::Store(std::shared_ptr<const ImageTextureData> image) {
  if (!image) {
    throw std::invalid_argument("Cannot store null image data in spectral image cache");
  }
  std::lock_guard<std::mutex> lock(mutex_);
  images_[image->Key()] = std::move(image);
}

std::size_t SpectralImageCache::Size() const {
  std::lock_guard<std::mutex> lock(mutex_);
  return images_.size();
}

FloatTexture::FloatTexture() = default;

FloatTexture::FloatTexture(Variant data) : data_(std::move(data)) {}

FloatTexture FloatTexture::Constant(Float value) {
  return FloatTexture(ConstantFloatTexture{ValidateFiniteScale(value, "Float texture constant")});
}

FloatTexture FloatTexture::Scale(FloatTexture texture, Float scale) {
  return FloatTexture(ScaleFloatTexture{
    MakeFloatTexturePtr(std::move(texture)),
    ValidateFiniteScale(scale, "Float texture scale")
  });
}

FloatTexture FloatTexture::Mix(FloatTexture tex1, FloatTexture tex2, Float amount) {
  return FloatTexture(MixFloatTexture{
    MakeFloatTexturePtr(std::move(tex1)),
    MakeFloatTexturePtr(std::move(tex2)),
    ValidateFiniteScale(amount, "Float texture mix amount")
  });
}

FloatTexture FloatTexture::Mix(FloatTexture tex1, FloatTexture tex2, FloatTexture amount) {
  return FloatTexture(MixFloatTexture{
    MakeFloatTexturePtr(std::move(tex1)),
    MakeFloatTexturePtr(std::move(tex2)),
    MakeFloatTexturePtr(std::move(amount))
  });
}

FloatTexture FloatTexture::Checker(
  FloatTexture tex1,
  FloatTexture tex2,
  int dimension,
  Float uscale,
  Float vscale
) {
  if (dimension != 2 && dimension != 3) {
    throw std::invalid_argument("Checker texture dimension must be 2 or 3");
  }
  if (!(uscale > 0) || !(vscale > 0) || !std::isfinite(uscale) || !std::isfinite(vscale)) {
    throw std::invalid_argument("Checker texture scales must be finite and positive");
  }
  return FloatTexture(CheckerFloatTexture{
    MakeFloatTexturePtr(std::move(tex1)),
    MakeFloatTexturePtr(std::move(tex2)),
    dimension,
    uscale,
    vscale
  });
}

FloatTexture FloatTexture::Noise(Float value1, Float value2, Float scale, Float phase, Float intensity) {
  return FloatTexture(NoiseFloatTexture{
    ValidateFiniteScale(value1, "Noise texture value1"),
    ValidateFiniteScale(value2, "Noise texture value2"),
    ValidateFiniteScale(scale, "Noise texture scale"),
    ValidateFiniteScale(phase, "Noise texture phase"),
    ValidateFiniteScale(intensity, "Noise texture intensity"),
    perlin()
  });
}

FloatTexture FloatTexture::Image(std::shared_ptr<const ImageTextureData> image, Float scale, bool invert) {
  if (!image) {
    throw std::invalid_argument("Float image texture requires image data");
  }
  if (image->ValueType() != TextureValueType::Float) {
    throw std::invalid_argument("Float image texture requires scalar image data");
  }
  return FloatTexture(ImageFloatTexture{
    std::move(image),
    ValidateFiniteScale(scale, "Float image texture scale"),
    invert
  });
}

Float FloatTexture::Evaluate(const TextureEvalContext& ctx) const {
  return UniversalTextureEvaluator()(*this, ctx);
}

TextureKind FloatTexture::Kind() const {
  return std::visit(
    [](const auto& texture) {
      using T = std::decay_t<decltype(texture)>;
      if constexpr (std::is_same<T, ConstantFloatTexture>::value) {
        return TextureKind::Constant;
      } else if constexpr (std::is_same<T, ScaleFloatTexture>::value) {
        return TextureKind::Scale;
      } else if constexpr (std::is_same<T, MixFloatTexture>::value) {
        return TextureKind::Mix;
      } else if constexpr (std::is_same<T, CheckerFloatTexture>::value) {
        return TextureKind::Checker;
      } else if constexpr (std::is_same<T, NoiseFloatTexture>::value) {
        return TextureKind::Noise;
      } else {
        return TextureKind::Image;
      }
    },
    data_
  );
}

SpectrumTexture::SpectrumTexture()
  : data_(ConstantSpectrumTexture{base::Spectrum(base::ConstantSpectrum(0))}) {}

SpectrumTexture::SpectrumTexture(Variant data) : data_(std::move(data)) {}

SpectrumTexture SpectrumTexture::Constant(base::Spectrum spectrum) {
  if (!spectrum.IsValid()) {
    throw std::invalid_argument("Spectrum texture constant requires a valid spectrum");
  }
  return SpectrumTexture(ConstantSpectrumTexture{std::move(spectrum)});
}

SpectrumTexture SpectrumTexture::Constant(Float value) {
  return Constant(base::Spectrum(base::ConstantSpectrum(value)));
}

SpectrumTexture SpectrumTexture::Scale(SpectrumTexture texture, Float scale) {
  return SpectrumTexture(ScaleSpectrumTexture{
    MakeSpectrumTexturePtr(std::move(texture)),
    ValidateFiniteScale(scale, "Spectrum texture scale")
  });
}

SpectrumTexture SpectrumTexture::Mix(SpectrumTexture tex1, SpectrumTexture tex2, Float amount) {
  return SpectrumTexture(MixSpectrumTexture{
    MakeSpectrumTexturePtr(std::move(tex1)),
    MakeSpectrumTexturePtr(std::move(tex2)),
    ValidateFiniteScale(amount, "Spectrum texture mix amount")
  });
}

SpectrumTexture SpectrumTexture::Mix(SpectrumTexture tex1, SpectrumTexture tex2, FloatTexture amount) {
  return SpectrumTexture(MixSpectrumTexture{
    MakeSpectrumTexturePtr(std::move(tex1)),
    MakeSpectrumTexturePtr(std::move(tex2)),
    MakeFloatTexturePtr(std::move(amount))
  });
}

SpectrumTexture SpectrumTexture::Checker(
  SpectrumTexture tex1,
  SpectrumTexture tex2,
  int dimension,
  Float uscale,
  Float vscale
) {
  if (dimension != 2 && dimension != 3) {
    throw std::invalid_argument("Checker texture dimension must be 2 or 3");
  }
  if (!(uscale > 0) || !(vscale > 0) || !std::isfinite(uscale) || !std::isfinite(vscale)) {
    throw std::invalid_argument("Checker texture scales must be finite and positive");
  }
  return SpectrumTexture(CheckerSpectrumTexture{
    MakeSpectrumTexturePtr(std::move(tex1)),
    MakeSpectrumTexturePtr(std::move(tex2)),
    dimension,
    uscale,
    vscale
  });
}

SpectrumTexture SpectrumTexture::Noise(
  SpectrumTexture tex1,
  SpectrumTexture tex2,
  Float scale,
  Float phase,
  Float intensity
) {
  return SpectrumTexture(NoiseSpectrumTexture{
    MakeSpectrumTexturePtr(std::move(tex1)),
    MakeSpectrumTexturePtr(std::move(tex2)),
    ValidateFiniteScale(scale, "Spectrum noise texture scale"),
    ValidateFiniteScale(phase, "Spectrum noise texture phase"),
    ValidateFiniteScale(intensity, "Spectrum noise texture intensity"),
    perlin()
  });
}

SpectrumTexture SpectrumTexture::Image(
  std::shared_ptr<const ImageTextureData> image,
  TextureSemanticRole role,
  base::RGBColorSpace colorSpace,
  Float scale
) {
  if (!image) {
    throw std::invalid_argument("Spectrum image texture requires image data");
  }
  if (image->ValueType() != TextureValueType::Spectrum) {
    throw std::invalid_argument("Spectrum image texture requires color image data");
  }
  if (image->Key().role != role) {
    throw std::invalid_argument("Spectrum image texture role must match the image cache key role");
  }
  if (role == TextureSemanticRole::Scalar) {
    throw std::invalid_argument("Spectrum image texture requires a color semantic role");
  }
  return SpectrumTexture(ImageSpectrumTexture{
    std::move(image),
    role,
    std::move(colorSpace),
    ValidateFiniteScale(scale, "Spectrum image texture scale")
  });
}

base::SampledSpectrum SpectrumTexture::Evaluate(
  const TextureEvalContext& ctx,
  const base::SampledWavelengths& lambda
) const {
  return UniversalTextureEvaluator()(*this, ctx, lambda);
}

TextureKind SpectrumTexture::Kind() const {
  return std::visit(
    [](const auto& texture) {
      using T = std::decay_t<decltype(texture)>;
      if constexpr (std::is_same<T, ConstantSpectrumTexture>::value) {
        return TextureKind::Constant;
      } else if constexpr (std::is_same<T, ScaleSpectrumTexture>::value) {
        return TextureKind::Scale;
      } else if constexpr (std::is_same<T, MixSpectrumTexture>::value) {
        return TextureKind::Mix;
      } else if constexpr (std::is_same<T, CheckerSpectrumTexture>::value) {
        return TextureKind::Checker;
      } else if constexpr (std::is_same<T, NoiseSpectrumTexture>::value) {
        return TextureKind::Noise;
      } else {
        return TextureKind::Image;
      }
    },
    data_
  );
}

Float UniversalTextureEvaluator::operator()(const FloatTexture& texture, const TextureEvalContext& ctx) const {
  return std::visit(
    [&](const auto& node) -> Float {
      using T = std::decay_t<decltype(node)>;
      if constexpr (std::is_same<T, ConstantFloatTexture>::value) {
        return node.value;
      } else if constexpr (std::is_same<T, ScaleFloatTexture>::value) {
        if (!node.texture) {
          throw std::runtime_error("Scale float texture child is missing");
        }
        return (*this)(*node.texture, ctx) * node.scale;
      } else if constexpr (std::is_same<T, MixFloatTexture>::value) {
        if (!node.tex1 || !node.tex2) {
          throw std::runtime_error("Mix float texture child is missing");
        }
        Float amount = TextureAmount(node.amount, ctx, *this);
        return (static_cast<Float>(1) - amount) * (*this)(*node.tex1, ctx) +
               amount * (*this)(*node.tex2, ctx);
      } else if constexpr (std::is_same<T, CheckerFloatTexture>::value) {
        if (!node.tex1 || !node.tex2) {
          throw std::runtime_error("Checker float texture child is missing");
        }
        return UseSecondCheckerTexel(ctx, node.dimension, node.uscale, node.vscale) ?
                 (*this)(*node.tex2, ctx) :
                 (*this)(*node.tex1, ctx);
      } else if constexpr (std::is_same<T, NoiseFloatTexture>::value) {
        Float amount = NoiseWeight(node.noise, ctx, node.scale, node.phase, node.intensity);
        return (static_cast<Float>(1) - amount) * node.value1 + amount * node.value2;
      } else {
        Float value = node.image->SampleScalar(ctx);
        value = node.invert ? static_cast<Float>(1) - value : value;
        return value * node.scale;
      }
    },
    texture.data_
  );
}

base::SampledSpectrum UniversalTextureEvaluator::operator()(
  const SpectrumTexture& texture,
  const TextureEvalContext& ctx,
  const base::SampledWavelengths& lambda
) const {
  return std::visit(
    [&](const auto& node) -> base::SampledSpectrum {
      using T = std::decay_t<decltype(node)>;
      if constexpr (std::is_same<T, ConstantSpectrumTexture>::value) {
        return node.spectrum.Sample(lambda);
      } else if constexpr (std::is_same<T, ScaleSpectrumTexture>::value) {
        if (!node.texture) {
          throw std::runtime_error("Scale spectrum texture child is missing");
        }
        return (*this)(*node.texture, ctx, lambda) * node.scale;
      } else if constexpr (std::is_same<T, MixSpectrumTexture>::value) {
        if (!node.tex1 || !node.tex2) {
          throw std::runtime_error("Mix spectrum texture child is missing");
        }
        Float amount = TextureAmount(node.amount, ctx, *this);
        return (static_cast<Float>(1) - amount) * (*this)(*node.tex1, ctx, lambda) +
               amount * (*this)(*node.tex2, ctx, lambda);
      } else if constexpr (std::is_same<T, CheckerSpectrumTexture>::value) {
        if (!node.tex1 || !node.tex2) {
          throw std::runtime_error("Checker spectrum texture child is missing");
        }
        return UseSecondCheckerTexel(ctx, node.dimension, node.uscale, node.vscale) ?
                 (*this)(*node.tex2, ctx, lambda) :
                 (*this)(*node.tex1, ctx, lambda);
      } else if constexpr (std::is_same<T, NoiseSpectrumTexture>::value) {
        if (!node.tex1 || !node.tex2) {
          throw std::runtime_error("Noise spectrum texture child is missing");
        }
        Float amount = NoiseWeight(node.noise, ctx, node.scale, node.phase, node.intensity);
        return (static_cast<Float>(1) - amount) * (*this)(*node.tex1, ctx, lambda) +
               amount * (*this)(*node.tex2, ctx, lambda);
      } else {
        base::RGB linearRGB = node.image->SampleRGB(ctx) * node.scale;
        return ReconstructLinearRGB(linearRGB, node.role, node.colorSpace, lambda);
      }
    },
    texture.data_
  );
}

base::TextureHandle TextureHandleTable::Add(FloatTexture texture) {
  Entry entry;
  entry.valueType = TextureValueType::Float;
  entry.generation = 1;
  entry.texture = std::move(texture);
  entries_.push_back(std::move(entry));
  return base::TextureHandle::FromIndex(static_cast<base::TextureHandle::IndexType>(entries_.size() - 1), 1);
}

base::TextureHandle TextureHandleTable::Add(SpectrumTexture texture) {
  Entry entry;
  entry.valueType = TextureValueType::Spectrum;
  entry.generation = 1;
  entry.texture = std::move(texture);
  entries_.push_back(std::move(entry));
  return base::TextureHandle::FromIndex(static_cast<base::TextureHandle::IndexType>(entries_.size() - 1), 1);
}

TextureValueType TextureHandleTable::ValueType(base::TextureHandle handle) const {
  return GetEntry(handle).valueType;
}

const FloatTexture& TextureHandleTable::GetFloat(base::TextureHandle handle) const {
  const Entry& entry = GetEntry(handle);
  if (entry.valueType != TextureValueType::Float) {
    throw std::invalid_argument("Texture handle does not reference a FloatTexture");
  }
  return std::get<FloatTexture>(entry.texture);
}

const SpectrumTexture& TextureHandleTable::GetSpectrum(base::TextureHandle handle) const {
  const Entry& entry = GetEntry(handle);
  if (entry.valueType != TextureValueType::Spectrum) {
    throw std::invalid_argument("Texture handle does not reference a SpectrumTexture");
  }
  return std::get<SpectrumTexture>(entry.texture);
}

const TextureHandleTable::Entry& TextureHandleTable::GetEntry(base::TextureHandle handle) const {
  if (!handle.IsValid() || handle.Index() >= entries_.size()) {
    throw std::out_of_range("Texture handle is invalid");
  }
  const Entry& entry = entries_[handle.Index()];
  if (entry.generation != handle.Generation()) {
    throw std::out_of_range("Texture handle generation is stale");
  }
  return entry;
}

TextureDebugEvaluation DebugEvaluateFloatTexture(
  const FloatTexture& texture,
  const TextureEvalContext& ctx
) {
  TextureDebugEvaluation evaluation;
  evaluation.valueType = TextureValueType::Float;
  evaluation.kind = texture.Kind();
  evaluation.scalar = texture.Evaluate(ctx);
  return evaluation;
}

TextureDebugEvaluation DebugEvaluateSpectrumTexture(
  const SpectrumTexture& texture,
  const TextureEvalContext& ctx,
  const base::SampledWavelengths& lambda
) {
  TextureDebugEvaluation evaluation;
  evaluation.valueType = TextureValueType::Spectrum;
  evaluation.kind = texture.Kind();
  evaluation.spectrum = texture.Evaluate(ctx, lambda);
  return evaluation;
}

} // namespace materials
} // namespace rayrender
