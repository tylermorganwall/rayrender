#include "texturecache.h"
#include "texture.h"
#include <vector>
#include <memory>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <cctype>
#ifndef STBIMAGEH
#define STBIMAGEH
#include "stb/stb_image.h"
#endif
#include "tinyobj/tinyexr.h"
#include <filesystem>
namespace fs = std::filesystem;

TextureCache::~TextureCache() {
  for (size_t i = 0; i < rawDataFloat.size(); ++i) {
    if (loadedBySTB[i]) {
      stbi_image_free(rawDataFloat[i]);
    } else {
      free(rawDataFloat[i]);
    }
  }
  for (size_t i = 0; i < rawDataChar.size(); ++i) {
    stbi_image_free(rawDataChar[i]);
  }
}

Float* TextureCache::LookupFloat(const std::string& filename,
                                 int& nx, int& ny, int& nn, int desired_channels) {
  std::string standardizedFilename = StandardizeFilename(filename);;
  auto it = hashTableFloat.find(standardizedFilename);
  if (it != hashTableFloat.end()) {
    auto itDim = hashTableDims.find(standardizedFilename);
    nx = std::get<0>(itDim->second);
    ny = std::get<1>(itDim->second);
    nn = std::get<2>(itDim->second);
    return it->second;
  }
  
  Float* data = LoadImageFloat(standardizedFilename, nx, ny, nn, desired_channels);
  if (!data) {
    throw std::runtime_error("Failed to load image: " + filename);
  }
  
  hashTableFloat[standardizedFilename] = data;
  hashTableDims[standardizedFilename] = std::make_tuple(nx,ny,desired_channels);
  rawDataFloat.push_back(data);
  return data;
}

unsigned char * TextureCache::LookupChar(const std::string& filename,
                                         int& nx, int& ny, int& nn, int desired_channels) {
  std::string standardizedFilename = StandardizeFilename(filename);;
  auto it = hashTableChar.find(standardizedFilename);
  if (it != hashTableChar.end()) {
    auto itDim = hashTableDims.find(standardizedFilename);
    nx = std::get<0>(itDim->second);
    ny = std::get<1>(itDim->second);
    nn = std::get<2>(itDim->second);
    return it->second;
  }
  
  unsigned char * data = LoadImageChar(standardizedFilename, nx, ny, nn, desired_channels);
  if (!data) {
    throw std::runtime_error("Failed to load image: " + filename);
  }
  
  hashTableChar[standardizedFilename] = data;
  hashTableDims[standardizedFilename] = std::make_tuple(nx,ny,desired_channels);
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
  const char* input = standardizedFilename.c_str();
  fs::path filepath(input);
  bool is_exr = filepath.extension() == ".exr";
  
  float* data = nullptr;
  const char* err = nullptr;
  if (is_exr) {
    if (IsEXR(input) != TINYEXR_SUCCESS) {
      throw std::runtime_error("Not an EXR file.");
    }
    EXRVersion exr_version;
    ParseEXRVersionFromFile(&exr_version, input);
    
    EXRHeader header;
    InitEXRHeader(&header);
    
    int header_ret = ParseEXRHeaderFromFile(&header, &exr_version, input, &err);
    if (header_ret != TINYEXR_SUCCESS) {
      if (err) {
      //   Rcpp::Rcout<< "Loading of '" + standardizedFilename << 
      //     "' failed due to: " << err << std::endl;
        FreeEXRErrorMessage(err);
        FreeEXRHeader(&header);
        throw std::runtime_error("Error loading EXR header");
      }
    }
    int ret = LoadEXR(&data, &width, &height, input, &err);
    // channels = header.num_channels;
    channels = 4;
    
    FreeEXRHeader(&header);
    if (err) {
      Rcpp::Rcout<< "Loading of '" + standardizedFilename << "' failed due to: " << err << 
        " -- nx/ny/channels :"  + std::to_string(width)  +  "/"  +  std::to_string(height)  +  "/"  +  std::to_string(channels) << 
        std::endl;
      FreeEXRErrorMessage(err);
    }
    
    if (ret != TINYEXR_SUCCESS) {
      throw std::runtime_error("Failed to load EXR image: " + filename);
    }
    
    loadedBySTB.push_back(false);
  } else {
    data = stbi_loadf(standardizedFilename.c_str(), &width, &height, &channels, desired_channels);
    if (!data) {
      throw std::runtime_error("Loading of '" + standardizedFilename  +
                               "' (float) failed due to error: " + stbi_failure_reason() +
                               "-- nx/ny/channels :"  + std::to_string(width)  +  "/"  +  std::to_string(height)  +  "/"  +  std::to_string(channels));
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
  std::string standardizedFilename = StandardizeFilename(filename);
  const char* input = standardizedFilename.c_str();
  fs::path filepath(input);

  unsigned char * data = nullptr;
  data = stbi_load(standardizedFilename.c_str(), &width, &height, &channels, desired_channels);
  if (!data) {
    throw std::runtime_error("Loading of '" + standardizedFilename  +
                             "' (char) failed due to error: " + stbi_failure_reason() +
                             "-- nx/ny/channels :"  + std::to_string(width)  +  "/"  +  std::to_string(height)  +  "/"  +  std::to_string(channels));
  }
  
  if (width == 0 || height == 0 || channels == 0) {
    throw std::runtime_error("Could not find " + filename);
  }
  
  return data;
}