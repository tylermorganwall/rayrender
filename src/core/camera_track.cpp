#include "camera_track.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace animation {
namespace {
int ProjectionOf(double fov) { return fov < 0 ? -1 : fov == 0 ? 0 : fov == 360 ? 360 : 1; }
}

CameraTrack::CameraTrack(std::vector<CameraKey> input) : keys(std::move(input)) {
  if (keys.empty()) throw std::runtime_error("Camera track requires at least one key.");
  if (Projection() < 0)
    throw std::runtime_error("Realistic-camera tracks in preview are not implemented yet; use offline explicit times.");
  for (size_t i = 0; i < keys.size(); ++i) {
    const auto &key = keys[i];
    if (!std::isfinite(key.time) || (i && key.time <= keys[i-1].time))
      throw std::runtime_error("Camera key times must be finite and strictly increasing.");
    for (double value : key.values)
      if (!std::isfinite(value) || !std::isfinite(Float(value)))
        throw std::runtime_error("Camera key values must be finite in renderer precision.");
    if (key.values[6] < 0 || (Projection() == 1 &&
        (key.values[7] >= 180 || key.values[8] <= 0)) ||
        (Projection() == 0 && (key.values[9] <= 0 || key.values[10] <= 0)))
      throw std::runtime_error("Invalid camera aperture, field of view, focus, or orthographic size.");
    if (ProjectionOf(key.values[7]) != Projection())
      throw std::runtime_error("Interactive camera tracks cannot change projection yet; use offline explicit times.");
  }
}

int CameraTrack::Projection() const { return ProjectionOf(keys.front().values[7]); }

CameraKey CameraTrack::At(double time) const {
  if (!std::isfinite(time)) throw std::runtime_error("Camera time must be finite.");
  auto upper = std::upper_bound(keys.begin(), keys.end(), time,
    [](double value, const CameraKey &key) { return value < key.time; });
  if (upper == keys.begin()) return keys.front();
  if (upper == keys.end()) return keys.back();
  const auto &left = *(upper-1);
  double fraction = (time-left.time)/(upper->time-left.time);
  CameraKey result{time, {}};
  for (size_t i = 0; i < result.values.size(); ++i)
    result.values[i] = left.values[i]*(1-fraction) + upper->values[i]*fraction;
  return result;
}

void CameraTrack::Apply(const CameraKey &key, RayCamera &cam) const {
  const auto &v = key.values;
  const point3f origin(v[0],v[1],v[2]), target(v[3],v[4],v[5]);
  const vec3f up(v[11],v[12],v[13]);
  if ((target-origin).squared_length() < 1e-16 || up.squared_length() < 1e-16 ||
      cross(target-origin,up).squared_length() < 1e-16)
    throw std::runtime_error("Camera track evaluated a degenerate view direction or up vector.");
  if (ProjectionOf(cam.get_fov()) != Projection())
    throw std::runtime_error("Camera track projection does not match the preview camera.");
  // Validate before touching the live view. The renderer calls this only at a
  // completed-batch handoff, never while workers are reading the camera.
  auto replace = [&](auto &live) {
    auto next = live;
    next.update_position_absolute(origin);
    next.update_lookat(target);
    next.update_up(up);
    next.update_focal_absolute(v[8]);
    if (Projection() == 1) {
      next.update_fov_absolute(v[7]);
      next.update_aperture_absolute(v[6]);
    }
    if (Projection() == 0) next.update_ortho_absolute(vec2f(v[9],v[10]));
    for (int axis = 0; axis < 3; ++axis)
      if (!std::isfinite(next.get_origin()[axis]) || !std::isfinite(next.get_w()[axis]) ||
          !std::isfinite(next.get_u()[axis]) || !std::isfinite(next.get_v()[axis]))
        throw std::runtime_error("Camera track exceeds renderer geometry precision.");
    live = std::move(next);
  };
  if (auto *value = dynamic_cast<camera *>(&cam)) replace(*value);
  else if (auto *value = dynamic_cast<ortho_camera *>(&cam)) replace(*value);
  else if (auto *value = dynamic_cast<environment_camera *>(&cam)) replace(*value);
  else throw std::runtime_error("Unsupported timeline camera class.");
}
} // namespace animation
