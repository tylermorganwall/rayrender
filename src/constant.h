#ifndef CONSTANTH
#define CONSTANTH

#include "hitable.h"
#include "material.h"
#include <float.h>

class material;

class constant_medium : public hitable {
public:
  constant_medium(hitable *b, Float d, texture *a ) : boundary(b), density(d) {
    phase_function = new isotropic(a);
  }
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const {
    return(boundary->bounding_box(t0,t1,box));
  }
  Float pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
    return(boundary->pdf_value(o,v, rng));
  }
  vec3 random(const vec3& o, random_gen& rng) {
    return(boundary->random(o, rng));
  }
  hitable *boundary;
  Float density;
  material *phase_function;
};

bool constant_medium::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  hit_record rec1, rec2;
  if(boundary->hit(r, -FLT_MAX, FLT_MAX,rec1, rng)) {
    if(boundary->hit(r, rec1.t + 0.0001, FLT_MAX, rec2, rng)) {
      if(rec1.t < t_min) {
        rec1.t = t_min;
      }
      if(rec2.t > t_max) {
        rec2.t = t_max;
      }
      if(rec1.t >= rec2.t) {
        return(false);
      }
      if(rec1.t < 0) {
        rec1.t = 0;
      }
      Float distance_inside_boundary = (rec2.t - rec1.t) * r.direction().length();
      Float hit_distance = -(1/density) * log(rng.unif_rand());
      if(hit_distance < distance_inside_boundary) {
        rec.t = rec1.t + hit_distance / r.direction().length();
        rec.p = r.point_at_parameter(rec.t);
        rec.normal = vec3(1,0,0);
        rec.mat_ptr = phase_function;
        return(true);
      }
    }
  }
  return(false);
}

#endif
