#include "../math/bounds.h"
#include "../math/rng.h"
#include "../math/sampler.h"
#include "../core/ray.h"

template <typename T>
bool Bounds3<T>::HitP(const Ray &r, Float tMax, Float *hitt0, Float *hitt1) const {
    Float t0 = 0, t1 = tMax;
    for (int i = 0; i < 3; ++i) {
        // Update interval for _i_th bounding box slab
        Float invRayDir = r.inv_dir_pad[i];
        Float tNear = (pMin[i] - r.o[i]) * invRayDir;
        Float tFar = (pMax[i] - r.o[i]) * invRayDir;
        // Update parametric interval from slab intersection $t$ values
        if (tNear > tFar)
            std::swap(tNear, tFar);
        // Update _tFar_ to ensure robust ray--bounds intersection
        tFar *= 1 + 2 * gamma(3);

        t0 = tNear > t0 ? tNear : t0;
        t1 = tFar < t1 ? tFar : t1;
        if (t0 > t1)
            return false;
    }
    if (hitt0)
        *hitt0 = t0;
    if (hitt1)
        *hitt1 = t1;
    return true;
}
