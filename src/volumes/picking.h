#ifndef RAYRENDER_VOLUME_PICKING_H
#define RAYRENDER_VOLUME_PICKING_H
#include "boundary.h"

struct RayPick {
  point3f p;
  bool volume;
};

// Select the first point with 1 - mean(exp(-optical_depth)) >= opacity,
// or the first ordinary surface. Misses and cancelled queries return no point.
// No transport samples or renderer RNG state are consumed.
std::optional<RayPick> PickRay(const Ray &ray, hitable *world, const VolumeScene *scene,
                              double opacity = 0.15,
                              const std::function<bool()> &cancel = {});
#endif
