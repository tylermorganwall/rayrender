#ifndef COLORH
#define COLORH

#include <cfloat>

#include "ray.h"
#include "hitablelist.h"
#include "material.h"

enum class IntegratorType {
    ShadowRays = 1,
    BasicPathGuiding = 2,
    Basic = 3
};

point3f color(const ray& r, hitable *world, hitable_list *hlist,
              size_t max_depth, size_t roulette_activate, random_gen& rng, Sampler* sampler,
              bool& alpha, IntegratorType type);

#endif
