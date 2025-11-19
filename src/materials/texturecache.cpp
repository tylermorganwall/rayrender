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
#include <ImfHeader.h>
#include <ImfInputFile.h>
#include <ImfTiledInputFile.h>
#include <ImfChannelList.h>
#include <ImfFrameBuffer.h>

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
  const bool is_exr = filepath.extension() == ".exr";

  float* data = nullptr;
  if (is_exr) {
    try {
      // Open header (works for both scanline and tiled via their specific files below)
      Imf::InputFile* scan = nullptr;
      Imf::TiledInputFile* tiled = nullptr;

      // Peek header once to get the window and channel list
      Imf::Header hdr;
      {
        // Prefer scanline path; if it fails with tiled exception, fall back
        try {
          scan = new Imf::InputFile(filename.c_str());
          hdr = scan->header();
        } catch (...) {
          delete scan; scan = nullptr;
          tiled = new Imf::TiledInputFile(filename.c_str());
          hdr = tiled->header();
        }
      }

      const Imath::Box2i dw = hdr.dataWindow();
      width  = dw.max.x - dw.min.x + 1;
      height = dw.max.y - dw.min.y + 1;

      // Decide channel count (3 or 4). We’ll try R,G,B,A; if A missing, we’ll fill 1.0f.
      const bool hasR = hdr.channels().findChannel("R") != nullptr;
      const bool hasG = hdr.channels().findChannel("G") != nullptr;
      const bool hasB = hdr.channels().findChannel("B") != nullptr;
      const bool hasA = hdr.channels().findChannel("A") != nullptr;

      // If only Y exists, you could special-case it here; for now assume RGB present.
      if (!(hasR && hasG && hasB)) {
        throw std::runtime_error("EXR missing R/G/B channels (unsupported minimal loader).");
      }

      const bool wantAlpha = (desired_channels == 4) || (desired_channels == 0 && hasA);
      channels = wantAlpha ? 4 : 3;

      const size_t pixCount = (size_t)width * (size_t)height;
      data = (float*)std::malloc(sizeof(float) * pixCount * channels);
      if (!data) throw std::runtime_error("malloc failed for EXR buffer");

      // Zero/init (so missing A can be set to 1.0f after read)
      std::fill(data, data + pixCount * channels, 0.0f);
      if (!wantAlpha) {
        // nothing
      } else if (!hasA) {
        // Fill A=1 if we’re outputting 4 channels but file lacks A
        for (size_t i = 0; i < pixCount; ++i) data[i*4 + 3] = 1.0f;
      }

      // Build framebuffer slices into interleaved [y*w + x]*C + c layout
      Imf::FrameBuffer fb;
      const size_t xStride = sizeof(float) * (size_t)channels;
      const size_t yStride = xStride * (size_t)width;

      auto insertSlice = [&](const char* name, int cIndex, bool present) {
        if (!present) return;
        char* base = reinterpret_cast<char*>(data)
                   + cIndex * sizeof(float)
                   - (dw.min.x * (ptrdiff_t)xStride + dw.min.y * (ptrdiff_t)yStride);
        fb.insert(name, Imf::Slice(Imf::FLOAT, base,
                                   (size_t)xStride, (size_t)yStride,
                                   1, 1, /*fill*/ 0.0f));
      };

      insertSlice("R", 0, hasR);
      insertSlice("G", 1, hasG);
      insertSlice("B", 2, hasB);
      if (channels == 4) insertSlice("A", 3, hasA); // if A missing, stays as init (1.0f set above)

      // Read pixels
      if (scan) {
        scan->setFrameBuffer(fb);
        scan->readPixels(dw.min.y, dw.max.y);
        delete scan;
      } else {
        tiled->setFrameBuffer(fb);
        // Read all tiles at level 0
        const int tx0 = 0, ty0 = 0;
        const int tx1 = tiled->numXTiles(0) - 1;
        const int ty1 = tiled->numYTiles(0) - 1;
        tiled->readTiles(tx0, tx1, ty0, ty1, /*lx*/0, /*ly*/0);
        delete tiled;
      }

      // Sanitize NaN/Inf (optional but wise)
      for (size_t i = 0, n = pixCount * (size_t)channels; i < n; ++i) {
        if (!std::isfinite(data[i])) data[i] = 0.0f;
      }

      loadedBySTB.push_back(false);
    } catch (const std::exception& e) {
      if (data) std::free(data);
      throw std::runtime_error("Failed to load EXR (float) '" + filename + "': " + e.what());
    }
  } else {
    // Non-EXR: stb float path
    data = stbi_loadf(filename.c_str(), &width, &height, &channels, desired_channels);
    if (!data) {
      throw std::runtime_error(
        "Loading of '" + filename + "' (float) failed: " + stbi_failure_reason()
        + " -- nx/ny/channels: " + std::to_string(width) + "/" + std::to_string(height)
        + "/" + std::to_string(channels));
    }
    if (desired_channels != 0) channels = desired_channels;
    loadedBySTB.push_back(true);
  }

  if (width == 0 || height == 0 || channels == 0) {
    if (data) std::free(data);
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
