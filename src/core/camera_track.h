#ifndef RAYRENDER_CAMERA_TRACK_H
#define RAYRENDER_CAMERA_TRACK_H

#include "camera.h"
#include <array>
#include <vector>

namespace animation {
// Columns match ray_camera_motion_columns(). Times are always scene seconds.
struct CameraKey {
  double time = 0;
  std::array<double, 14> values{};
};

class CameraTrack {
public:
  explicit CameraTrack(std::vector<CameraKey> keys);
  CameraKey At(double time) const;
  void Apply(const CameraKey &key, RayCamera &camera) const;
  int Projection() const;

private:
  std::vector<CameraKey> keys;
};
} // namespace animation
#endif
