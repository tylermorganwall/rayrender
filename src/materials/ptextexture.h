#ifndef RAYRENDER_PTEX_TEXTURE_H
#define RAYRENDER_PTEX_TEXTURE_H

#include "texture.h"
#include <ptex/ptex_api.h>
#include <atomic>
#include <mutex>
#include <string>

// The provider owns its implementation. Only opaque handles and the versioned
// C function table cross the package boundary (no Ptex C++ ABI or static link).
class PtexRuntime {
public:
  PtexRuntime();
  ~PtexRuntime();
  PtexRuntime(const PtexRuntime&) = delete;
  PtexRuntime& operator=(const PtexRuntime&) = delete;
  void RecordFailure(const std::string& filename, const std::string& message) const;
  void Report(bool verbose); // Main thread, after rendering workers have joined.
  const ptex_api_v1* api = nullptr;
  ptex_cache* cache = nullptr;
private:
  Rcpp::RObject api_handle; // Roots the provider namespace until all handles die.
  mutable std::atomic<uint64_t> failures{0};
  mutable std::mutex failure_mutex;
  mutable std::string first_failure;
};

class PtexTextureResource {
public:
  PtexTextureResource(std::shared_ptr<PtexRuntime> runtime, std::string filename);
  ~PtexTextureResource();
  PtexTextureResource(const PtexTextureResource&) = delete;
  PtexTextureResource& operator=(const PtexTextureResource&) = delete;
  point3f Sample(const TextureEvalContext& context, int filter) const;
private:
  std::shared_ptr<PtexRuntime> runtime;
  std::string filename;
  ptex_texture* handle = nullptr;
  ptex_info_v1 info{};
};

#endif
