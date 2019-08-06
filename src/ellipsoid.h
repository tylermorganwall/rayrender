#ifndef ELLIPSOIDH
#define ELLIPSOIDH

#include "sphere.h"
#include "material.h"

class ellipsoid: public hitable {
  public:
    ellipsoid() {}
    ellipsoid(vec3 cen, float r, vec3 axes, material *mat) : 
      center(cen), radius(r), axes(axes),  mat_ptr(mat) {
      inv_axes = vec3(1.0f/axes.x(), 1.0f/axes.y(), 1.0f/axes.z());
    };
    virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec, random_gen& rng);
    virtual bool bounding_box(float t0, float t1, aabb& box) const;
    virtual float pdf_value(const vec3& o, const vec3& v, random_gen& rng);
    virtual vec3 random(const vec3& o, random_gen& rng);
    vec3 center;
    float radius;
    vec3 axes;
    vec3 inv_axes;
    material *mat_ptr;
};


bool ellipsoid::hit(const ray& r, float t_min, float t_max, hit_record& rec, random_gen& rng) {
  vec3 oc = r.origin() * inv_axes - center;
  float a = dot(r.direction() * inv_axes, r.direction() * inv_axes);
  float b = dot(oc, r.direction() * inv_axes); 
  float c = dot(oc,oc) - radius * radius;
  float discriminant = b * b -  a * c; 
  if(discriminant > 0) {
    float temp = (-b - sqrt(discriminant))/a;
    if(temp < t_max && temp > t_min) {
      rec.t = temp;
      rec.p = r.point_at_parameter(rec.t);
      rec.normal = (rec.p - center) / radius;
      get_sphere_uv(rec.normal, rec.u, rec.v);
      rec.mat_ptr = mat_ptr;
      rec.normal *= inv_axes;
      rec.normal.make_unit_vector();
      return(true);
    }
    temp = (-b + sqrt(discriminant))/a;
    if(temp < t_max && temp > t_min) {
      rec.t = temp;
      rec.p = r.point_at_parameter(rec.t);
      rec.normal = (rec.p - center) / radius * inv_axes;
      get_sphere_uv(rec.normal, rec.u, rec.v);
      rec.mat_ptr = mat_ptr;
      rec.normal *= inv_axes;
      rec.normal.make_unit_vector();
      return(true);
    }
  }
  return(false);
}

float ellipsoid::pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, rng)) {
    float cos_theta_max = sqrt(1 - radius * radius/(center - o).squared_length());
    float solid_angle = 2 * M_PI * (1-cos_theta_max) * axes.x() * axes.y() * axes.z();
    return(1/solid_angle);
  } else {
    return(0);
  }
}

vec3 ellipsoid::random(const vec3& o, random_gen& rng) {
  vec3 direction = center - o;
  float distance_squared = direction.squared_length();
  onb uvw;
  uvw.build_from_w(direction);
  return(uvw.local(rng.random_to_sphere(radius,distance_squared) * inv_axes));
}

bool ellipsoid::bounding_box(float t0, float t1, aabb& box) const {
  box = aabb(center - vec3(radius,radius,radius) * axes, center + vec3(radius,radius,radius) * axes);
  return(true);
}

#endif
