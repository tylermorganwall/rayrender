#ifndef ANIMATEDTRANSFORMH
#define ANIMATEDTRANSFORMH

#include "mathinline.h"
#include "transform.h"
#include "interval.h"
#include "quaternion.h"

// AnimatedTransform Declarations
class AnimatedTransform {
public:
  // AnimatedTransform Public Methods
  AnimatedTransform(Transform* startTransform, Float startTime,
                    Transform*  endTransform, Float endTime);
  static void Decompose(const Matrix4x4 &m, vec3f *T, Quaternion *R,
                        Matrix4x4 *S);
  void Interpolate(Float time, Transform *t) const;
  ray operator()(const ray &r) const;
  // rayDifferential operator()(const rayDifferential &r) const;
  point3f operator()(Float time, const point3f &p) const;
  vec3f operator()(Float time, const vec3f &v) const;
  bool HasScale() const {
    return startTransform->HasScale() || endTransform->HasScale();
  }
  aabb MotionBounds(const aabb &b) const;
  aabb BoundPointMotion(const point3f &p) const;
  Transform* GetStartTransform() {return(startTransform);};
  Transform* GetEndTransform() {return(endTransform);};
  
  vec3f w();
  vec3f u();
  vec3f v();
  
  
private:
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
      return kc + kx * p.x() + ky * p.y() + kz * p.z();
    }
  };
  DerivativeTerm c1[3], c2[3], c3[3], c4[3], c5[3];
};

#endif