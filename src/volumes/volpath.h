#ifndef RAYRENDER_VOLPATH_H
#define RAYRENDER_VOLPATH_H
#include "../core/color.h"

void color_volume(const Ray &ray, hitable *world, hitable_list *lights, size_t max_depth,
                  size_t roulette_depth, random_gen &rng, Sampler *sampler, Float &transparency,
                  point3f &radiance, normal3f &normal, point3f &albedo,
                  const std::atomic<bool> *cancel);
#endif
