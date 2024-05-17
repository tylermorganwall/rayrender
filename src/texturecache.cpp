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
  for (size_t i = 0; i < rawData.size(); ++i) {
    if (loadedBySTB[i]) {
      stbi_image_free(rawData[i]);
    } else {
      free(rawData[i]);
    }
  }
}

std::shared_ptr<image_texture_float> TextureCache::Lookup(const std::string& filename,
                                                          int& nx, int& ny,
                                                          Float repeatu = 1.0f, 
                                                          Float repeatv = 1.0f, 
                                                          Float intensity = 1.0f) {
  std::string standardizedFilename = StandardizeFilename(filename);
  std::string standardizedFilenameInfo = standardizedFilename + 
    std::to_string(repeatu) + 
    std::to_string(repeatv) + 
    std::to_string(intensity);
  auto it = hashTable.find(standardizedFilenameInfo);
  if (it != hashTable.end()) {
    auto itDim = hashTableDims.find(standardizedFilenameInfo);
    nx = itDim->second.first;
    ny = itDim->second.second;
    return it->second;
  }
  
  int channels;
  Float* data = LoadImage(standardizedFilename, nx, ny, channels);
  if (!data) {
    throw std::runtime_error("Failed to load image: " + filename);
  }
  
  auto tex = std::make_shared<image_texture_float>(data, nx, ny, channels,
                                                   repeatu, repeatv, intensity);
  hashTable[standardizedFilenameInfo] = tex;
  hashTableDims[standardizedFilenameInfo] = std::pair<int, int>(nx,ny);
  rawData.push_back(data);
  return tex;
}

std::string TextureCache::StandardizeFilename(const std::string& filename) {
  std::string result = filename;
  std::transform(result.begin(), result.end(), result.begin(), [](unsigned char c) {
    return std::tolower(c);
  });
  return result;
}

float* TextureCache::LoadImage(const std::string& filename, int& width, int& height, int& channels) {
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
        std::cerr << "Error loading EXR header: " << err << std::endl;
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
      std::cerr << err << std::endl;
      FreeEXRErrorMessage(err);
    }
    
    if (ret != TINYEXR_SUCCESS) {
      throw std::runtime_error("Failed to load EXR image: " + filename);
    }
    
    loadedBySTB.push_back(false);
  } else {
    data = stbi_loadf(standardizedFilename.c_str(), &width, &height, &channels, 3);
    if (!data) {
      std::cerr << "Load failed: " << stbi_failure_reason() << std::endl;
      throw std::runtime_error("Loading failed: " + standardizedFilename);
    }
    loadedBySTB.push_back(true);
  }
  
  if (width == 0 || height == 0 || channels == 0) {
    throw std::runtime_error("Could not find " + filename);
  }
  
  return data;
}