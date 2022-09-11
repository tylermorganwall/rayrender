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
    virtual Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0);
    virtual Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0);
    virtual vec3f random(const point3f& o, random_gen& rng, Float time = 0);
    virtual vec3f random(const point3f& o, Sampler* sampler, Float time = 0);
    void add(std::shared_ptr<hitable> object) { objects.push_back(object); }
    std::shared_ptr<hitable> back() {return(objects.back());}
    
    int size() {return(objects.size());}
    std::vector<std::shared_ptr<hitable>> objects;
    std::string GetName() const;
    size_t GetSize();
};

#endif
