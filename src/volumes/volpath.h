#ifndef RAYRENDER_VOLPATH_H
#define RAYRENDER_VOLPATH_H
#include "../core/color.h"
class VolumeScene;

// Independent straight-camera-ray transmittance, shared by CPU rendering and
// the Metal hybrid atmosphere service. Visible finite surfaces have zero
// transparency; invisible emitters and boundary-only media preserve the ray.
Float VolumePrimaryTransparency(const Ray &, const VolumeScene &, hitable *, random_gen &,
                                const std::atomic<bool> *cancel);

void color_volume(const Ray &ray, hitable *world, hitable_list *lights, size_t max_depth,
                  size_t roulette_depth, random_gen &rng, Sampler *sampler, Float &transparency,
                  point3f &radiance, normal3f &normal, point3f &albedo,
                  const std::atomic<bool> *cancel);
#endif
