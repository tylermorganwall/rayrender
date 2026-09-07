#ifndef COLORH
#define COLORH

#include <atomic>
#include <cfloat>

#include "../core/ray.h"
#include "../hitables/hitablelist.h"
#include "../materials/material.h"

enum class IntegratorType {
    ShadowRays = 1,
    BasicPathGuiding = 2,
    Basic = 3
};

void color(const Ray& r, hitable *world, hitable_list *hlist,
           size_t max_depth, size_t roulette_activate, random_gen& rng, Sampler* sampler,
           Float& transparency, IntegratorType type,
           point3f& color, normal3f& normal, point3f& albedo,
           const std::atomic<bool>* cancel = nullptr);

#endif
