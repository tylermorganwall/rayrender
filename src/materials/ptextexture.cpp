#include "ptextexture.h"
#include "texturecache.h"
#include <ptex/ptex_r.h>
#include <algorithm>
#include <cmath>

PtexRuntime::PtexRuntime() {
  Rcpp::Environment base = Rcpp::Environment::base_env();
  Rcpp::Environment provider = Rcpp::Environment::namespace_env("ptex");
  Rcpp::Function get_api = provider["ptex_api"];
  api_handle = get_api();
  api = ptex_api_from_R_v1(api_handle);

  Rcpp::Function get_option = base["getOption"];
  double files = Rcpp::as<double>(get_option("rayrender.ptex_cache_files", 100));
  double bytes = Rcpp::as<double>(get_option("rayrender.ptex_cache_memory", 256.0 * 1024 * 1024));
  if (!std::isfinite(files) || files < 1 || files > INT32_MAX || files != std::floor(files) ||
      !std::isfinite(bytes) || bytes < 1 || bytes > 9007199254740991.0 || bytes != std::floor(bytes))
    Rcpp::stop("Ptex cache options must be positive integers (files <= INT32_MAX, bytes <= 2^53-1).");
  ptex_error error{};
  if (api->cache_create(int32_t(files), uint64_t(bytes), 1, &cache, &error) != PTEX_OK)
    Rcpp::stop("Cannot create Ptex cache: %s", error.message);
}

PtexRuntime::~PtexRuntime() {
  if (cache) api->cache_destroy(cache);
}

void PtexRuntime::RecordFailure(const std::string& filename, const std::string& message) const {
  // No R API calls on workers. Keep a bounded diagnostic and count every failed
  // lookup; a corrupt block or invalid face must not throw across the pool.
  if (failures.fetch_add(1, std::memory_order_relaxed) == 0) {
    std::lock_guard<std::mutex> lock(failure_mutex);
    first_failure = filename + ": " + message;
  }
}

void PtexRuntime::Report(bool verbose) {
  const uint64_t count = failures.exchange(0, std::memory_order_relaxed);
  if (count) {
    std::string message;
    {
      std::lock_guard<std::mutex> lock(failure_mutex);
      message.swap(first_failure);
    }
    Rcpp::warning("Ptex: %.0f texture lookups failed and returned zero. First failure: %s",
                  double(count), message.c_str());
  }
  if (verbose) {
    ptex_stats_v1 stats{};
    ptex_error error{};
    if (api->cache_stats(cache, &stats, &error) == PTEX_OK)
      Rcpp::Rcout << "Ptex cache: " << stats.files_accessed << " files accessed, "
                  << stats.block_reads << " block reads, " << stats.file_reopens
                  << " reopens, peak " << double(stats.peak_memory_used) / (1024 * 1024)
                  << " MiB, " << stats.peak_files_open << " open files\n";
  }
}

PtexTextureResource::PtexTextureResource(std::shared_ptr<PtexRuntime> runtime, std::string filename)
  : runtime(std::move(runtime)), filename(std::move(filename)) {
  ptex_error error{};
  if (this->runtime->api->texture_open(this->runtime->cache, this->filename.c_str(), &handle, &error) != PTEX_OK)
    Rcpp::stop("Cannot open Ptex texture '%s': %s", this->filename.c_str(), error.message);
  if (this->runtime->api->texture_info(handle, &info, &error) != PTEX_OK ||
      (info.channels != 1 && info.channels != 3)) {
    this->runtime->api->texture_destroy(handle);
    handle = nullptr;
    Rcpp::stop("Ptex texture '%s' must have one or three channels.", this->filename.c_str());
  }
}

PtexTextureResource::~PtexTextureResource() {
  if (handle) runtime->api->texture_destroy(handle);
}

point3f PtexTextureResource::Sample(const TextureEvalContext& context, int filter) const {
  if (context.face_index < 0) {
    runtime->RecordFailure(filename, "Ptex requires a mesh face and face-local UV coordinates");
    return point3f(0);
  }
  float derivatives[4] = {};
  if (context.footprint.valid) {
    derivatives[0] = context.footprint.dudx;
    derivatives[1] = context.footprint.dvdx;
    derivatives[2] = context.footprint.dudy;
    derivatives[3] = context.footprint.dvdy;
    // The runtime ABI bounds derivatives to one face. Scale the complete
    // footprint uniformly so its orientation/aspect ratio survive saturation.
    float extent = 1;
    for (float x : derivatives) extent = std::max(extent, std::abs(x));
    for (float& x : derivatives) x /= extent;
  }
  float u = std::clamp(float(context.u), 0.f, 1.f);
  float v = std::clamp(float(context.v), 0.f, 1.f);
  // Finite-difference bump probes may just cross the face edge. Bound their
  // center to the face; Ptex still filters across its stored adjacency.
  if (info.mesh_type == 0 && u + v > 1) { const float sum = u + v; u /= sum; v /= sum; }
  ptex_filter_v1 options{filter, 0, 0, 0}; // PBRT: B-spline, no mip-level lerp.
  float result[3] = {};
  ptex_error error{};
  if (runtime->api->sample(handle, &options, context.face_index, u, v,
                           derivatives[0], derivatives[1], derivatives[2], derivatives[3],
                           0, info.channels, result, &error) != PTEX_OK) {
    runtime->RecordFailure(filename, error.message);
    return point3f(0);
  }
  for (int c = 0; c < info.channels; ++c) {
    if (!std::isfinite(result[c])) {
      runtime->RecordFailure(filename, "nonfinite texel value");
      return point3f(0);
    }
  }
  return info.channels == 1 ? point3f(result[0]) : point3f(result[0], result[1], result[2]);
}

std::shared_ptr<const PtexTextureResource> TextureCache::LookupPtex(const std::string& filename) {
  const std::string key = fs::weakly_canonical(filename).string();
  auto found = ptexTextures.find(key);
  if (found != ptexTextures.end()) return found->second;
  if (!ptexRuntime) ptexRuntime = std::make_shared<PtexRuntime>();
  auto resource = std::make_shared<PtexTextureResource>(ptexRuntime, key);
  ptexTextures.emplace(key, resource);
  return resource;
}

void TextureCache::ReportPtex(bool verbose) {
  if (ptexRuntime) ptexRuntime->Report(verbose);
}
