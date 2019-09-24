#ifndef DISKH
#define DISKH

#include "hitable.h"
#include "onbh.h"

class disk : public hitable {
public:
  disk() {}
  disk(vec3 cen, Float r, Float i_r, material *mat) : center(cen), radius(r), inner_radius(i_r), mat_ptr(mat) {};
  virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng);
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng);
  virtual vec3 random(const vec3& o, random_gen& rng);
  vec3 center;
  Float radius;
  Float inner_radius;
  material *mat_ptr;
};

bool disk::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  vec3 n(0.0, 1.0, 0.0);
  // First we intersect with the plane containing the disk
  Float t = -r.origin().y() / r.direction().y();
  if(t < t_min || t > t_max) {
    return(false);
  }
  Float x = r.origin().x() + t*r.direction().x();
  Float z = r.origin().z() + t*r.direction().z();
  if(x*x + z*z >= radius * radius || x*x + z*z <= inner_radius * inner_radius) {
    return(false);
  }
  vec3 p = r.point_at_parameter(t);
  p.e[1] = 0;
  rec.p = p;
  rec.normal = n;
  rec.t = t;
  rec.mat_ptr = mat_ptr;
  rec.u = p.x() / (2.0 * radius) + 0.5;
  rec.v = p.z() / (2.0 * radius) + 0.5;
  return(true);
}

Float disk::pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, rng)) {
    Float area =  M_PI * (radius * radius - inner_radius * inner_radius);
    Float distance_squared = rec.t * rec.t * v.squared_length();
    Float cosine = fabs(dot(v,rec.normal)/v.length());
    return(distance_squared / (cosine * area));
  } else {
    return(0);
  }
}

vec3 disk::random(const vec3& o, random_gen& rng) {
  Float r1 = rng.unif_rand();
  Float r2 = sqrt(rng.unif_rand());
  Float phi = 2 * M_PI * r1;
  Float x = ((radius - inner_radius) * r2 + inner_radius) * cos(phi);
  Float z = ((radius - inner_radius) * r2 + inner_radius) * sin(phi);
  return(vec3(x,0,z)+center-o);
}

bool disk::bounding_box(Float t0, Float t1, aabb& box) const {
  box = aabb(-vec3(radius,0.001,radius), vec3(radius,0.001,radius));
  return(true);
}


#endif
