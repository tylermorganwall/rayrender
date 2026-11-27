#ifndef RAYRENDER_PREVIEW_NAVIGATION_H
#define RAYRENDER_PREVIEW_NAVIGATION_H

#include "camera.h"
#include "../math/aabb.h"
#include <functional>

class hitable;
class VolumeScene;

// Navigation scale is independent of lens focus and the camera's unit forward
// vector. Window backends supply elapsed seconds and dimensionless input axes.
class PreviewNavigation {
public:
  void Initialize(RayCamera& camera, hitable* world, const VolumeScene* volumes,
                  const aabb& finite_bounds, const std::function<bool()>& cancel = {});
  void SetDistance(Float value);
  Float Distance() const { return distance; }
  Float ViewHeight(RayCamera& camera) const;
  void Reset(RayCamera& camera);
  bool Move(RayCamera& camera, Float horizontal, Float vertical, Float forward,
            bool shifted, bool orbit, Float multiplier, double seconds);
private:
  Float distance = 1;
  Float initial_distance = 1;
};
#endif
