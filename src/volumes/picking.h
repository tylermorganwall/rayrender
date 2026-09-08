#ifndef RAYRENDER_VOLUME_PICKING_H
#define RAYRENDER_VOLUME_PICKING_H
#include "boundary.h"

struct RayPick {
  point3f p;
  bool volume;
  bool background = false;
};

// Select the first point with 1 - mean(exp(-optical_depth)) >= opacity,
// or the first ordinary surface. Environment hits are marked as background so
// camera controls can select their direction without using the proxy sphere's
// distance. Misses and cancelled queries return no point.
// No transport samples or renderer RNG state are consumed.
std::optional<RayPick> PickRay(const Ray &ray, hitable *world, const VolumeScene *scene,
                              double opacity = 0.15,
                              const std::function<bool()> &cancel = {});
#endif
