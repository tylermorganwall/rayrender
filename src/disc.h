#ifndef DISCH
#define DISCH

#include "hitable.h"
#include "onbh.h"

class disc : public hitable {
public:
  disc() {}
  disc(vec3 cen, float r, float i_r, material *mat) : center(cen), radius(r), inner_radius(i_r), mat_ptr(mat) {};
  virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec, random_gen& rng);
  virtual bool bounding_box(float t0, float t1, aabb& box) const;
  virtual float pdf_value(const vec3& o, const vec3& v, random_gen& rng);
  virtual vec3 random(const vec3& o, random_gen& rng);
  vec3 center;
  float radius;
  float inner_radius;
  material *mat_ptr;
};

bool disc::hit(const ray& r, float t_min, float t_max, hit_record& rec, random_gen& rng) {
  vec3 n(0.0, 1.0, 0.0);
  // First we intersect with the plane containing the disc
  float t = -r.origin().y() / r.direction().y();
  if(t < t_min || t > t_max) {
    return(false);
  }
  float x = r.origin().x() + t*r.direction().x();
  float z = r.origin().z() + t*r.direction().z();
  if(x*x + z*z >= radius * radius || x*x + z*z <= inner_radius * inner_radius) {
    return(false);
  }
  vec3 p = r.point_at_parameter(t);
  rec.p = p;
  rec.normal = (dot(r.origin(), n) < 0) ? -n : n;
  rec.t = t;
  rec.mat_ptr = mat_ptr;
  rec.u = p.x() / (2.0 * radius) + 0.5;
  rec.v = p.z() / (2.0 * radius) + 0.5;
  return(true);
}

float disc::pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, rng)) {
    float cos_theta_max = sqrt(1 - radius * radius/(center - o).squared_length());
    float solid_angle = 2 * M_PI * (1-cos_theta_max) * dot(rec.normal,v) ;
    return(1/solid_angle);
  } else {
    return(0);
  }
}

vec3 disc::random(const vec3& o, random_gen& rng) {
  vec3 direction = center - o;
  float distance_squared = direction.squared_length();
  onb uvw;
  uvw.build_from_w(direction);
  return(uvw.local(rng.random_to_sphere(radius,distance_squared)));
}

bool disc::bounding_box(float t0, float t1, aabb& box) const {
  box = aabb(-vec3(radius,0.001,radius), vec3(radius,0.001,radius));
  return(true);
}


#endif