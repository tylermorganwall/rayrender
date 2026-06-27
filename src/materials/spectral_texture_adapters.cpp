#include "../materials/spectral_texture.h"

#include "../hitables/hitable.h"
#include "../materials/texturecache.h"

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <utility>
#include <vector>

namespace rayrender {
namespace materials {
namespace {

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

} // namespace

TextureEvalContext TextureEvalContextFromHitRecord(const hit_record& hit) {
  return TextureEvalContextFromHitRecord(hit, TextureRayDifferentials{});
}

TextureEvalContext TextureEvalContextFromHitRecord(
  const hit_record& hit,
  const TextureRayDifferentials& differentials
) {
  TextureEvalContext ctx;
  ctx.p = hit.p;
  ctx.uv = point2f(hit.u, hit.v);
  ctx.dpdu = hit.dpdu;
  ctx.dpdv = hit.dpdv;
  ctx.dpdx = differentials.dpdx;
  ctx.dpdy = differentials.dpdy;
  ctx.dudx = differentials.dudx;
  ctx.dvdx = differentials.dvdx;
  ctx.dudy = differentials.dudy;
  ctx.dvdy = differentials.dvdy;
  return ctx;
}

std::shared_ptr<const ImageTextureData> ImageTextureData::LoadFromFile(
  const ImageCacheKey& key,
  const base::RGBColorSpace& inputColorSpace
) {
  int width = 0;
  int height = 0;
  int channels = 0;
  constexpr int desiredChannels = 4;
  std::vector<Float> values;

  if (IsHighDynamicRangeFile(key.filename)) {
    TextureCache cache;
    Float* data = cache.LookupFloat(key.filename, width, height, channels, desiredChannels);
    values.assign(data, data + static_cast<std::size_t>(width) * height * channels);
  } else {
    unsigned char* data = stbi_load(key.filename.c_str(), &width, &height, &channels, desiredChannels);
    if (!data) {
      throw std::runtime_error(
        "Loading of '" + key.filename + "' failed due to error: " + stbi_failure_reason()
      );
    }
    channels = desiredChannels;
    values.resize(static_cast<std::size_t>(width) * height * channels);
    constexpr Float inv255 = static_cast<Float>(1) / static_cast<Float>(255);
    for (std::size_t i = 0; i < values.size(); ++i) {
      values[i] = static_cast<Float>(data[i]) * inv255;
    }
    stbi_image_free(data);
  }

  return FromInterleaved(key, width, height, channels, std::move(values), inputColorSpace);
}

std::shared_ptr<const ImageTextureData> SpectralImageCache::LookupOrLoad(
  const ImageCacheKey& key,
  const base::RGBColorSpace& inputColorSpace
) {
  return LookupOrCreate(key, [&](const ImageCacheKey& loadKey) {
    return ImageTextureData::LoadFromFile(loadKey, inputColorSpace);
  });
}

} // namespace materials
} // namespace rayrender
