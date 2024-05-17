#ifndef TEXTURECACHE
#define TEXTURECACHE

//This should include three containers for data: an std::vector<Float* > that holds the raw pointers from
//loading the data via std_image/tinyexr, a std::vector<bool> that tracks whether it was stb_image 
//or tinyexr that loaded it (to properly free the memory at the end),
//and a std::vector<std::shared_ptr<texture> > hashTable to hold the rayrender texture classes that will be looked up in the hash table,
//This should have a lookup function that returns a raw ptr (using ,get()) for a given image filename.
//If the filename isn't in the hash table, then it will load the file and put it in the hash table, and then return the pointer.
//It should ensure the filename is standardized (so as to avoid different hashes from capitalization
//and so on) so that when it encounters the same filename twice, it won't load the image again.
//It should track whether the files were loaded via stb_image or 

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

class TextureCache {
public:
  TextureCache() = default;
  ~TextureCache();
  
  std::shared_ptr<image_texture_float> Lookup(const std::string& filename,
                                              int& nx, int& ny,
                                              Float repeatu, Float repeatv, Float intensity);
  
private:
  std::vector<float*> rawData;
  std::vector<bool> loadedBySTB; // true for stb_image, false for tinyexr
  std::unordered_map<std::string, std::shared_ptr<image_texture_float>> hashTable;
  std::unordered_map<std::string, std::pair<int, int>> hashTableDims;
  
  static std::string StandardizeFilename(const std::string& filename);
  float* LoadImage(const std::string& filename, int& width, int& height, int& channels);
};

#endif // TEXTURECACHE
