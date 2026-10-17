#include "camera_track.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace animation {
namespace {
int ProjectionOf(double fov) { return fov < 0 ? -1 : fov == 0 ? 0 : fov == 360 ? 360 : 1; }

bool SeparatedViewNormal(const std::array<vec3<double>,3> &control, double margin, int depth = 0) {
  // cross(target-origin, up) is quadratic for linear camera rows. Its Bezier
  // control hull must exclude a small box about zero throughout the interval.
  // Subdivision certifies difficult curves; uncertainty is rejected, not sampled.
  for (int axis = 0; axis < 3; ++axis) {
    double lo = std::min({control[0][axis],control[1][axis],control[2][axis]});
    double hi = std::max({control[0][axis],control[1][axis],control[2][axis]});
    if (lo > margin || hi < -margin) return true;
  }
  if (depth == 20) return false;
  auto a = (control[0]+control[1])*.5, b = (control[1]+control[2])*.5;
  auto mid = (a+b)*.5;
  return SeparatedViewNormal({control[0],a,mid},margin,depth+1) &&
         SeparatedViewNormal({mid,b,control[2]},margin,depth+1);
}
}

CameraTrack::CameraTrack(std::vector<CameraKey> input) : keys(std::move(input)) {
  if (keys.empty()) throw std::runtime_error("Camera track requires at least one key.");
  for (size_t i = 0; i < keys.size(); ++i) {
    const auto &key = keys[i];
    const int projection = ProjectionOf(key.values[7]);
    if (projection < 0 || ProjectionOf(keys.front().values[7]) < 0) {
      if (projection != ProjectionOf(keys.front().values[7]))
        throw std::runtime_error("Physical-camera tracks cannot switch lens projection.");
      if (key.values[6] != keys.front().values[6] || key.values[8] != keys.front().values[8])
        throw std::runtime_error("Physical-camera playback requires constant aperture and focus; use offline zero-width snapshots for varying optics.");
      if (key.values[8] <= 0)
        throw std::runtime_error("Physical-camera focus must be positive.");
    }
    if (!std::isfinite(key.time) || (i && key.time <= keys[i-1].time))
      throw std::runtime_error("Camera key times must be finite and strictly increasing.");
    for (double value : key.values)
      if (!std::isfinite(value) || !std::isfinite(Float(value)))
        throw std::runtime_error("Camera key values must be finite in renderer precision.");
    if (key.values[6] < 0 || (projection == 1 &&
        (key.values[7] >= 180 || key.values[8] <= 0)) ||
        (projection == 0 && (key.values[9] <= 0 || key.values[10] <= 0)))
      throw std::runtime_error("Invalid camera aperture, field of view, focus, or orthographic size.");
    const auto &v = key.values;
    vec3<double> direction(v[3]-v[0],v[4]-v[1],v[5]-v[2]), up(v[11],v[12],v[13]);
    vec3<double> before = direction, before_up = up;
    if (i && ProjectionOf(keys[i-1].values[7]) == projection) {
      const auto &p = keys[i-1].values;
      before = vec3<double>(p[3]-p[0],p[4]-p[1],p[5]-p[2]);
      before_up = vec3<double>(p[11],p[12],p[13]);
    }
    std::array<vec3<double>,3> normal = {cross(before,before_up),
      (cross(before,up)+cross(direction,before_up))*.5, cross(direction,up)};
    double scale = std::max({normal[0].length(),normal[1].length(),normal[2].length()});
    if (!std::isfinite(scale) || scale == 0 || !SeparatedViewNormal(normal,scale*1e-6))
      throw std::runtime_error("Camera track view direction and up become degenerate or nearly parallel.");
  }
}

std::unique_ptr<RayCamera> CameraTrack::MakeCamera(const CameraKey &key, Float aspect, Float iso,
                                                 Float open, Float close, bool flip,
                                                 const RealisticCamera *physical) const {
  const auto &v = key.values;
  point3f origin(v[0],v[1],v[2]), target(v[3],v[4],v[5]);
  vec3f up(v[11],v[12],v[13]);
  std::unique_ptr<RayCamera> result;
  const int projection = ProjectionOf(v[7]);
  if (projection < 0) {
    if (!physical) throw std::runtime_error("Physical-camera playback requires initialized lens optics.");
    result = std::make_unique<RealisticCamera>(*physical);
    static_cast<RealisticCamera *>(result.get())->SetShutterInterval(open,close);
    result->scene_camera_track.reset();
  } else if (projection == 0)
    result = std::make_unique<ortho_camera>(origin,target,up,v[9],v[10],open,close,iso);
  else if (projection == 360)
    result = std::make_unique<environment_camera>(origin,target,up,open,close,iso);
  else
    result = std::make_unique<camera>(origin,target,up,v[7],aspect,v[6],v[8],open,close,iso);
  result->set_camera_flip_x(flip);
  result->set_shutter_speed(1);
  result->scene_camera_aspect = aspect;
  Apply(key,*result);
  return result;
}

Ray CameraTrack::SampleRay(double scene_time, Float ray_time, Float s, Float t,
                          point3f lens, Float aspect, Float iso, bool flip) const {
  const auto key = At(scene_time);
  const auto &v = key.values;
  const int projection = ProjectionOf(v[7]);
  point3f origin(v[0], v[1], v[2]), target(v[3], v[4], v[5]);
  vec3f up(v[11], v[12], v[13]);
  if (projection < 0)
    throw std::runtime_error("Physical rays require the initialized lens camera's GenerateRay method.");
  // Stack-local cameras reuse the existing projection/lens conventions. They
  // have no attached track, no mutable shared pose, and no R or heap allocation.
  if (projection == 0) {
    ortho_camera sample(origin, target, up, v[9], v[10], ray_time, ray_time, iso);
    sample.set_camera_flip_x(flip);
    return sample.get_ray(s, t, lens, 0);
  }
  if (projection == 360) {
    environment_camera sample(origin, target, up, ray_time, ray_time, iso);
    sample.set_camera_flip_x(flip);
    return sample.get_ray(s, t, lens, 0);
  }
  camera sample(origin, target, up, v[7], aspect, v[6], v[8], ray_time, ray_time, iso);
  sample.set_camera_flip_x(flip);
  return sample.get_ray(s, t, lens, 0);
}

CameraKey CameraTrack::At(double time) const {
  if (!std::isfinite(time)) throw std::runtime_error("Camera time must be finite.");
  auto upper = std::upper_bound(keys.begin(), keys.end(), time,
    [](double value, const CameraKey &key) { return value < key.time; });
  if (upper == keys.begin()) return keys.front();
  if (upper == keys.end()) return keys.back();
  const auto &left = *(upper-1);
  // A projection switch is a cut: hold the entire left camera until the key.
  if (ProjectionOf(left.values[7]) != ProjectionOf(upper->values[7])) return left;
  double fraction = (time-left.time)/(upper->time-left.time);
  CameraKey result{time, {}};
  for (size_t i = 0; i < result.values.size(); ++i)
    result.values[i] = left.values[i]*(1-fraction) + upper->values[i]*fraction;
  return result;
}

void CameraTrack::Apply(const CameraKey &key, RayCamera &cam) const {
  const auto &v = key.values;
  const int projection = ProjectionOf(v[7]);
  const point3f origin(v[0],v[1],v[2]), target(v[3],v[4],v[5]);
  const vec3f up(v[11],v[12],v[13]);
  if ((target-origin).squared_length() < 1e-16 || up.squared_length() < 1e-16 ||
      cross(target-origin,up).squared_length() < 1e-16)
    throw std::runtime_error("Camera track evaluated a degenerate view direction or up vector.");
  if (ProjectionOf(cam.get_fov()) != projection)
    throw std::runtime_error("Camera track projection does not match the preview camera.");
  // Validate before touching the live view. The renderer calls this only at a
  // completed-batch handoff, never while workers are reading the camera.
  auto replace = [&](auto &live) {
    auto next = live;
    next.update_position_absolute(origin);
    next.update_lookat(target);
    next.update_up(up);
    next.update_focal_absolute(v[8]);
    if (projection == 1) {
      next.update_fov_absolute(v[7]);
      next.update_aperture_absolute(v[6]);
    }
    if (projection == 0) next.update_ortho_absolute(vec2f(v[9],v[10]));
    for (int axis = 0; axis < 3; ++axis)
      if (!std::isfinite(next.get_origin()[axis]) || !std::isfinite(next.get_w()[axis]) ||
          !std::isfinite(next.get_u()[axis]) || !std::isfinite(next.get_v()[axis]))
        throw std::runtime_error("Camera track exceeds renderer geometry precision.");
    live = std::move(next);
  };
  if (auto *value = dynamic_cast<RealisticCamera *>(&cam)) {
    if (Float(v[8]) != value->get_focal_distance() || Float(v[6]) != value->get_aperture())
      throw std::runtime_error("Physical-camera lens aperture or focus does not match the authored track.");
    value->SetScenePose(origin,target,up);
  } else if (auto *value = dynamic_cast<camera *>(&cam)) replace(*value);
  else if (auto *value = dynamic_cast<ortho_camera *>(&cam)) replace(*value);
  else if (auto *value = dynamic_cast<environment_camera *>(&cam)) replace(*value);
  else throw std::runtime_error("Unsupported timeline camera class.");
}
} // namespace animation
