#include "../materials/texturecache.h"
#include "../materials/texture.h"
#include <vector>
#include <memory>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <cctype>
#include <cstdlib>
#ifndef STBIMAGEH
#define STBIMAGEH
#include "../ext/stb/stb_image.h"
#endif

#include <ImathBox.h>
#include <ImfArray.h>
#include <ImfRgbaFile.h>

using namespace OPENEXR_IMF_NAMESPACE;
using namespace IMATH_NAMESPACE;

#include <filesystem>
namespace fs = std::filesystem;

TextureCache::~TextureCache() {
  for (size_t i = 0; i < rawDataFloat.size(); ++i) {
    if (loadedBySTB[i]) {
      stbi_image_free(rawDataFloat[i]);
    } else {
      std::free(rawDataFloat[i]);
    }
  }
  for (size_t i = 0; i < rawDataChar.size(); ++i) {
    stbi_image_free(rawDataChar[i]);
  }
}

Float* TextureCache::LookupFloat(const std::string& filename,
                                 int& nx, int& ny, int& nn, int desired_channels) {
  std::string standardizedFilename = StandardizeFilename(filename);
  auto it = hashTableFloat.find(standardizedFilename);
  if (it != hashTableFloat.end()) {
    auto itDim = hashTableDims.find(standardizedFilename);
    nx = std::get<0>(itDim->second);
    ny = std::get<1>(itDim->second);
    nn = std::get<2>(itDim->second);
    return it->second;
  }
  
  Float* data = LoadImageFloat(filename, nx, ny, nn, desired_channels);
  if (!data) {
    throw std::runtime_error("Failed to load image: " + filename);
  }
  
  hashTableFloat[standardizedFilename] = data;
  hashTableDims[standardizedFilename] = std::make_tuple(nx, ny, nn);
  rawDataFloat.push_back(data);
  return data;
}

unsigned char * TextureCache::LookupChar(const std::string& filename,
                                         int& nx, int& ny, int& nn, int desired_channels) {
  std::string standardizedFilename = StandardizeFilename(filename);
  auto it = hashTableChar.find(standardizedFilename);
  if (it != hashTableChar.end()) {
    auto itDim = hashTableDims.find(standardizedFilename);
    nx = std::get<0>(itDim->second);
    ny = std::get<1>(itDim->second);
    nn = std::get<2>(itDim->second);
    return it->second;
  }
  
  unsigned char * data = LoadImageChar(filename, nx, ny, nn, desired_channels);
  if (!data) {
    throw std::runtime_error("Failed to load image: " + filename);
  }
  
  hashTableChar[standardizedFilename] = data;
  hashTableDims[standardizedFilename] = std::make_tuple(nx, ny, nn);
  rawDataChar.push_back(data);
  return data;
}

std::string TextureCache::StandardizeFilename(const std::string& filename) {
  std::string result = filename;
  std::transform(result.begin(), result.end(), result.begin(), [](unsigned char c) {
    return std::tolower(c);
  });
  return result;
}

float* TextureCache::LoadImageFloat(const std::string& filename, int& width, int& height,
                                    int& channels, int desired_channels) {
  std::string standardizedFilename = StandardizeFilename(filename);
  fs::path filepath(standardizedFilename);
  bool is_exr = filepath.extension() == ".exr";

  float* data = nullptr;
  if (is_exr) {
    try {
      RgbaInputFile file(filename.c_str());
      Box2i dw = file.dataWindow();
      width = dw.max.x - dw.min.x + 1;
      height = dw.max.y - dw.min.y + 1;

      OPENEXR_IMF_NAMESPACE::Array2D<Rgba> px;
      px.resizeErase(height, width);

      file.setFrameBuffer(&px[0][0] - dw.min.x - dw.min.y * width, 1, width);
      file.readPixels(dw.min.y, dw.max.y);

      constexpr int exr_channels = 4;
      const size_t pixel_count = static_cast<size_t>(width) * static_cast<size_t>(height);
      data = static_cast<float*>(std::malloc(sizeof(float) * pixel_count * exr_channels));
      if (!data) {
        throw std::runtime_error("Failed to allocate memory for EXR image: " + filename);
      }

      for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
          const Rgba& p = px[y][x];
          const size_t base = (static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x)) * exr_channels;
          data[base + 0] = p.r;
          data[base + 1] = p.g;
          data[base + 2] = p.b;
          data[base + 3] = p.a;
        }
      }

      channels = exr_channels;
      loadedBySTB.push_back(false);
    } catch (const std::exception& e) {
      throw std::runtime_error("Failed to load EXR image '" + filename + "': " + e.what());
    }
  } else {
    data = stbi_loadf(filename.c_str(), &width, &height, &channels, desired_channels);
    if (!data) {
      throw std::runtime_error("Loading of '" + filename  +
                               "' (float) failed due to error: " + stbi_failure_reason() +
                               "-- nx/ny/channels :"  + std::to_string(width)  +  "/"  +  std::to_string(height)  +  "/"  +  std::to_string(channels));
    }
    if (desired_channels != 0) {
      channels = desired_channels;
    }
    loadedBySTB.push_back(true);
  }

  if (width == 0 || height == 0 || channels == 0) {
    throw std::runtime_error("Could not find " + filename);
  }

  return data;
}

unsigned char * TextureCache::LoadImageChar(const std::string& filename, int& width, int& height, 
                                            int& channels, int desired_channels) {
  unsigned char * data = nullptr;
  data = stbi_load(filename.c_str(), &width, &height, &channels, desired_channels);
  if (!data) {
    throw std::runtime_error("Loading of '" + filename  +
                             "' (char) failed due to error: " + stbi_failure_reason() +
                             "-- nx/ny/channels :"  + std::to_string(width)  +  "/"  +  std::to_string(height)  +  "/"  +  std::to_string(channels));
  }
  if (desired_channels != 0) {
    channels = desired_channels;
  }
  
  if (width == 0 || height == 0 || channels == 0) {
    throw std::runtime_error("Could not find " + filename);
  }
  
  return data;
}
