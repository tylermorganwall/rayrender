#ifndef HITABLELISTH
#define HITABLELISTH

#include "hitable.h"
#include "sampler.h"

class hitable_list: public hitable {
  public:
    hitable_list() {}
    hitable_list(hitable **l, int n) {list = l; list_size = n;}
    virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng);
    virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
    virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng);
    virtual vec3 random(const vec3& o, random_gen& rng);
    virtual vec3 random(const vec3& o, Sampler* sampler);
    hitable **list;
    int list_size;
};

#endif
