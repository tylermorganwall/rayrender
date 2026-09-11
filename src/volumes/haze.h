#ifndef RAYRENDER_VOLUME_HAZE_H
#define RAYRENDER_VOLUME_HAZE_H
#include "medium.h"
#include <atomic>

struct MediumHazeSegment {
  double t_min, t_max;
  bool integrate;
};

// Traverse the actual trilinear density threshold, independently of stochastic
// scattering candidates. Adjacent intervals with the same selection are merged.
// The input ray is in medium coordinates; its parameter remains world distance.
class MediumHazeIterator {
public:
  MediumHazeIterator(const Medium *medium, const Ray &local, double end,
                     bool enabled, const std::atomic<bool> *cancel = nullptr);
  // Tracking callers disable merging to avoid scanning beyond their next event.
  std::optional<MediumHazeSegment> Next(bool merge = true);

private:
  const Medium *medium;
  const std::atomic<bool> *cancel;
  double end, cursor = 0, threshold = 0;
  bool uniform = true, enabled;
  DensityIndexRay index_ray;
  std::optional<RayMajorantIterator> majorants;
  std::optional<RayMajorantSegment> majorant;
  std::array<MediumHazeSegment, 4> pieces;
  size_t piece = 0, piece_count = 0;
  std::optional<MediumHazeSegment> lookahead;
  std::optional<MediumHazeSegment> Piece();
  void SplitCell(double a, double b);
};
#endif
