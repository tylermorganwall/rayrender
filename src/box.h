#ifndef BOXH
#define BOXH

#include "hitablelist.h"
#include "rectangle.h"
#include <memory>

class box : public hitable {
public:
  box() {}
  box(const vec3& p0, const vec3& p1, std::shared_ptr<material> ptr, 
      std::shared_ptr<alpha_texture> alpha_mask, std::shared_ptr<bump_texture> bump_tex);
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler);
  
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const {
    box = aabb(pmin, pmax);
    return(true);
  }
  virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng, Float time = 0) {
    return(list.pdf_value(o,v, rng, time));
  }
  virtual Float pdf_value(const vec3& o, const vec3& v, Sampler* sampler, Float time = 0) {
    return(list.pdf_value(o,v, sampler, time));
  }
  virtual vec3 random(const vec3& o, random_gen& rng, Float time = 0) {
    return(list.random(o, rng, time));
  }
  virtual vec3 random(const vec3& o, Sampler* sampler, Float time = 0) {
    return(list.random(o, sampler, time));
  }
  vec3 pmin, pmax;
  hitable_list list;
};


#endif
