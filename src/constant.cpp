#include "constant.h"


bool constant_medium::bounding_box(Float t0, Float t1, aabb& box) const {
  return(boundary->bounding_box(t0,t1,box));
}

Float constant_medium::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) {
  return(boundary->pdf_value(o,v, rng, time));
}

Float constant_medium::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) {
  return(boundary->pdf_value(o,v, sampler, time));
}

vec3f constant_medium::random(const point3f& o, random_gen& rng, Float time) {
  return(boundary->random(o, rng, time));
}

vec3f constant_medium::random(const point3f& o, Sampler* sampler, Float time) {
  return(boundary->random(o, sampler, time));
}


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
        rec.normal = vec3f(1,0,0);
        rec.mat_ptr = phase_function.get();
        return(true);
      }
    }
  }
  return(false);
}


bool constant_medium::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
  hit_record rec1, rec2;
  if(boundary->hit(r, -FLT_MAX, FLT_MAX,rec1, sampler)) {
    if(boundary->hit(r, rec1.t + 0.0001, FLT_MAX, rec2, sampler)) {
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
      Float hit_distance = -(1/density) * log(sampler->Get1D());
      if(hit_distance < distance_inside_boundary) {
        rec.t = rec1.t + hit_distance / r.direction().length();
        rec.p = r.point_at_parameter(rec.t);
        rec.normal = vec3f(1,0,0);
        rec.mat_ptr = phase_function.get();
        return(true);
      }
    }
  }
  return(false);
}
