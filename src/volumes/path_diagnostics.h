#ifndef RAYRENDER_PATH_DIAGNOSTICS_H
#define RAYRENDER_PATH_DIAGNOSTICS_H
#include "../core/ray.h"
#include <Rcpp.h>
#include <array>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

// Only explicitly classified transport failures are recoverable. Allocation,
// cancellation and scene construction errors must retain their normal behavior.
enum class PathFailureKind { RepeatedEntry, InvalidExit, Density, Majorant,
                             Roulette, PositionPrecision, Radiance, DiffusionInterior, Count };
class PathFailure : public std::runtime_error {
public:
  PathFailure(PathFailureKind kind, const std::string &message)
      : std::runtime_error(message), kind(kind) {}
  const PathFailureKind kind;
};

class PathDiagnostics {
public:
  // Called only on failed paths. Successful paths do not lock or increment counters.
  void Record(const PathFailure &, const Ray &camera, const Ray &current,
              size_t depth, uint64_t internal_events, uint64_t active_boundary,
              const char *stage);
  // Called on the R thread after workers join; drains per-frame diagnostics.
  Rcpp::List Take();
  // Called on the R thread between sample passes. Consumes only the notice,
  // retaining all counts/examples for the final diagnostics and optional log.
  bool TakeDiffusionInteriorNotice();
private:
  static constexpr size_t kinds = size_t(PathFailureKind::Count);
  static constexpr size_t examples_per_kind = 8;
  std::mutex mutex;
  std::array<uint64_t, kinds> counts{};
  std::array<std::vector<std::string>, kinds> examples;
  bool diffusion_interior_reported = false;
};
#endif
