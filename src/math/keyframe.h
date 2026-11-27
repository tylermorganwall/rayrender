#ifndef RAYRENDER_KEYFRAME_H
#define RAYRENDER_KEYFRAME_H

#include "transform.h"
#include <array>
#include <memory>
#include <vector>

namespace animation {
// Schema 1 uses seconds, column-vector transforms, and xyzw quaternions.
constexpr int schema_version = 1;
using Value = std::array<double, 4>;
enum class Property { Position, Rotation, Scale, Scalar };
enum class Interpolation { Step, Linear, Cubic, Smooth };
enum class Loop { Hold = 1, Repeat = 2, PingPong = 3 };

// Maps scene seconds to clip seconds without maintaining a playback cursor.
// The start is an anchor, not a visibility gate; repeat is half-open [a,b).
struct ClipClock {
  ClipClock(double begin, double end, double start, double offset, double speed, Loop loop);
  double Phase(double scene_time) const;
  double At(double scene_time) const;
  std::vector<std::pair<double, double>> Intervals(double begin, double end) const;
  const double begin, end, start, offset, speed;
  const Loop loop;
};

// Construct on the owning thread, then publish as const. Evaluation performs no
// allocations, R access, or cursor updates, including after backward seeks.
class Track {
public:
  Track(Property property, Interpolation interpolation, std::vector<double> times,
        std::vector<Value> values, std::vector<Value> incoming = {},
        std::vector<Value> outgoing = {});
  Value At(double time) const;
  std::array<Value, 2> Range(double begin, double end) const;
  void ValidateFloatPacking() const;
  void ValidateYaw(double relative_margin = 1e-8) const;
  int Dimensions() const { return property == Property::Rotation ? 4 : property == Property::Scalar ? 1 : 3; }
  Property property;
  Interpolation interpolation;
  std::vector<double> times;
  std::vector<Value> values, incoming, outgoing;

private:
  std::array<double, 4> Polynomial(size_t segment, int component) const;
};

// Remove only explicitly selected locomotion components, relative to the clip's
// first pose. Yaw is the quaternion's twist about asset +Y; swing is retained.
struct RootMotion {
  std::array<bool, 3> translation{};
  Value position{}, reference_yaw{0, 0, 0, 1};
  bool yaw = false;
};
struct TransformTrack {
  std::shared_ptr<const Track> position, rotation, scale;
  Value pivot{};
  std::shared_ptr<const ClipClock> clock;
  std::shared_ptr<const RootMotion> root_motion;
  Transform At(double time) const;
  aabb Bounds(const aabb &box, double begin, double end) const;
private:
  Transform AtLocal(double time) const;
  aabb LocalBounds(const aabb &box, double begin, double end) const;
};
} // namespace animation
#endif
