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
  std::unique_ptr<RayCamera> MakeCamera(const CameraKey &key, Float aspect, Float iso,
                                      Float open, Float close, bool flip,
                                      const RealisticCamera *physical = nullptr) const;
  Ray SampleRay(double scene_time, Float ray_time, Float s, Float t,
                point3f lens, Float aspect, Float iso, bool flip) const;
  const std::vector<CameraKey> &Keys() const { return keys; }

private:
  std::vector<CameraKey> keys;
};
} // namespace animation
#endif
