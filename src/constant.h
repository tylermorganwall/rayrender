#ifndef CONSTANTH
#define CONSTANTH

#include "hitable.h"
#include "material.h"
#include <float.h>
#include <memory>

class material;

class constant_medium : public hitable {
public:
  constant_medium(std::shared_ptr<hitable> b, Float d, std::shared_ptr<texture> a ) : boundary(b), density(d) {
    phase_function = std::make_shared<isotropic>(a);
  }
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler);
  
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  Float pdf_value(const vec3& o, const vec3& v, random_gen& rng, Float time = 0);
  Float pdf_value(const vec3& o, const vec3& v, Sampler* sampler, Float time = 0);
  vec3 random(const vec3& o, random_gen& rng, Float time = 0);
  vec3 random(const vec3& o, Sampler* sampler, Float time = 0);
  
  std::shared_ptr<hitable> boundary;
  Float density;
  std::shared_ptr<material> phase_function;
};

#endif
