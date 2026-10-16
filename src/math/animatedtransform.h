#ifndef ANIMATEDTRANSFORMH
#define ANIMATEDTRANSFORMH

#include "../math/mathinline.h"
#include "../math/transform.h"
#include "../math/interval.h"
#include "../math/quaternion.h"
#include "keyframe.h"

// Published only between completed render batches. Rays retain their existing
// float times; subtracting the exact float frame anchor before adding the
// double epoch preserves snapshots at discontinuous clip/STEP boundaries.
struct SceneTimeContext {
  double epoch = 0;
  Float ray_anchor = 0;
};

// AnimatedTransform Declarations
class AnimatedTransform {
public:
  // AnimatedTransform Public Methods
  AnimatedTransform(Transform* startTransform, Float startTime,
                    Transform*  endTransform, Float endTime,
                    std::shared_ptr<const animation::TransformTrack> track = nullptr,
                    double time_origin = 0,
                    std::shared_ptr<SceneTimeContext> scene_time = nullptr);
  static void Decompose(const Matrix4x4 &m, vec3f *T, Quaternion *R,
                        Matrix4x4 *S);
  void Interpolate(Float time, Transform *t) const;
  Float StartTime() const { return startTime; }
  Float EndTime() const { return endTime; }
  const Transform &StartTransform() const { return *startTransform; }
  const Transform &EndTransform() const { return *endTransform; }
  bool SameMotion(const AnimatedTransform &other) const {
    return track == other.track && scene_time == other.scene_time && time_origin == other.time_origin &&
           startTime == other.startTime && endTime == other.endTime &&
           *startTransform == *other.startTransform && *endTransform == *other.endTransform;
  }
  Ray operator()(const Ray &r) const;
  // rayDifferential operator()(const rayDifferential &r) const;
  point3f operator()(Float time, const point3f &p) const;
  vec3f operator()(Float time, const vec3f &v) const;
  bool HasScale() const {
    if (track) return bool(track->scale);
    return startTransform->HasScale() || endTransform->HasScale();
  }
  aabb MotionBounds(const aabb &b) const;
  aabb MotionBounds(const aabb &b, double begin, double end) const;
  bool UsesSceneTime() const { return bool(track) || bool(scene_time) || time_origin != 0; }
  const std::shared_ptr<const animation::TransformTrack> &Track() const { return track; }
  double TimeOrigin() const { return time_origin; }
  double FrameEpoch() const { return scene_time ? scene_time->epoch : time_origin; }
  Float FrameAnchor() const { return scene_time ? scene_time->ray_anchor : 0; }
  double SceneTime(Float time) const { return FrameEpoch() + (double(time) - double(FrameAnchor())); }
  aabb BoundPointMotion(const point3f &p) const;
  Transform* GetStartTransform() {return(startTransform);};
  Transform* GetEndTransform() {return(endTransform);};
  
  vec3f w();
  vec3f u();
  vec3f v();
  
  
private:
  std::shared_ptr<const animation::TransformTrack> track;
  double time_origin = 0;
  std::shared_ptr<SceneTimeContext> scene_time;
  // AnimatedTransform Private Data
  Transform *startTransform, *endTransform;
  const Float startTime, endTime;
  const bool actuallyAnimated;
  vec3f T[2];
  Quaternion R[2];
  Matrix4x4 S[2];
  bool hasRotation;
  struct DerivativeTerm {
    DerivativeTerm() {}
    DerivativeTerm(Float c, Float x, Float y, Float z)
      : kc(c), kx(x), ky(y), kz(z) {}
    Float kc, kx, ky, kz;
    Float Eval(const point3f &p) const {
      return kc + kx * p.xyz.x + ky * p.xyz.y + kz * p.xyz.z;
    }
  };
  DerivativeTerm c1[3], c2[3], c3[3], c4[3], c5[3];
};

#endif
