#ifndef BOXH
#define BOXH

#include "hitablelist.h"
#include "rectangle.h"
#include <memory>

class box : public hitable {
public:
  box() {}
  box(const vec3& p0, const vec3& p1, material *ptr, 
      std::shared_ptr<alpha_texture> alpha_mask, std::shared_ptr<bump_texture> bump_tex);
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const {
    box = aabb(pmin, pmax);
    return(true);
  }
  virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
    return(list.pdf_value(o,v, rng));
  }
  virtual vec3 random(const vec3& o, random_gen& rng) {
    return(list.random(o, rng));
  }
  virtual vec3 random(const vec3& o, Sampler* sampler) {
    return(list.random(o, sampler));
  }
  vec3 pmin, pmax;
  hitable_list list;
};


#endif
