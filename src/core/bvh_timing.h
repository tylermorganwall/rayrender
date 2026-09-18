#ifndef RAYRENDER_BVH_TIMING_H
#define RAYRENDER_BVH_TIMING_H

#include <chrono>

// Construction runs on the calling thread. Capture constructor-body wall time,
// including bounds, tree construction, flattening and BVH4 conversion. Nested
// constructors are counted but their intervals are not counted twice.
struct BVHBuildTiming {
  double seconds = 0;
  unsigned int count = 0;
  unsigned int depth = 0;
  inline static thread_local BVHBuildTiming* active = nullptr;
  BVHBuildTiming* previous;

  explicit BVHBuildTiming(bool enabled = true) : previous(active) {
    if (enabled) active = this;
  }
  ~BVHBuildTiming() { active = previous; }
  BVHBuildTiming(const BVHBuildTiming&) = delete;
  BVHBuildTiming& operator=(const BVHBuildTiming&) = delete;
};

struct ScopedBVHBuildTiming {
  using Clock = std::chrono::steady_clock;
  BVHBuildTiming* timing = BVHBuildTiming::active;
  Clock::time_point start;
  ScopedBVHBuildTiming() {
    if (timing && timing->depth++ == 0) start = Clock::now();
  }
  ~ScopedBVHBuildTiming() {
    if (!timing) return;
    ++timing->count;
    if (--timing->depth == 0)
      timing->seconds += std::chrono::duration<double>(Clock::now() - start).count();
  }
};

#endif
