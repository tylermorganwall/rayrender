#ifndef HITABLELISTH
#define HITABLELISTH

#include "hitable.h"
#include "sampler.h"
#include <memory>

class hitable_list: public hitable {
  public:
    hitable_list() {}
    hitable_list(std::shared_ptr<hitable> object) {add(object);}
    virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng);
    virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, Sampler* sampler);
    
    virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
    virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng, Float time = 0);
    virtual Float pdf_value(const vec3& o, const vec3& v, Sampler* sampler, Float time = 0);
    virtual vec3 random(const vec3& o, random_gen& rng, Float time = 0);
    virtual vec3 random(const vec3& o, Sampler* sampler, Float time = 0);
    void add(std::shared_ptr<hitable> object) { objects.push_back(object); }
    int size() {return(objects.size());}
    std::vector<std::shared_ptr<hitable>> objects;
};

#endif
