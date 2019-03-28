#ifndef SPHEREH
#define SPHEREH

#include "hitable.h"
#include "material.h"

inline vec3 random_to_sphere(float radius, float distance_squared) {
  float r1 = drand48();
  float r2 = drand48();
  float z = 1 + r2 * (sqrt(1-radius * radius / distance_squared) - 1);
  float phi = 2 * M_PI * r1;
  float x = cos(phi) * sqrt(1-z*z);
  float y = sin(phi) * sqrt(1-z*z);
  return(vec3(x,y,z));
}

class sphere: public hitable {
  public:
    sphere() {}
    sphere(vec3 cen, float r, material *mat) : center(cen), radius(r), mat_ptr(mat) {};
    virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec, random_gen& rng);
    virtual bool bounding_box(float t0, float t1, aabb& box) const;
    virtual float pdf_value(const vec3& o, const vec3& v, random_gen& rng);
    virtual vec3 random(const vec3& o, random_gen& rng);
    vec3 center;
    float radius;
    material *mat_ptr;
};

bool sphere::hit(const ray& r, float t_min, float t_max, hit_record& rec, random_gen& rng) {
  vec3 oc = r.origin() - center;
  float a = dot(r.direction(), r.direction());
  float b = dot(oc, r.direction()); 
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
      return(true);
    }
    temp = (-b + sqrt(discriminant))/a;
    if(temp < t_max && temp > t_min) {
      rec.t = temp;
      rec.p = r.point_at_parameter(rec.t);
      rec.normal = (rec.p - center) / radius;
      get_sphere_uv(rec.normal, rec.u, rec.v);
      rec.mat_ptr = mat_ptr;
      return(true);
    }
  }
  return(false);
}

float sphere::pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, rng)) {
    float cos_theta_max = sqrt(1 - radius * radius/(center - o).squared_length());
    float solid_angle = 2 * M_PI * (1-cos_theta_max);
    return(1/solid_angle);
  } else {
    return(0);
  }
}

vec3 sphere::random(const vec3& o, random_gen& rng) {
  vec3 direction = center - o;
  float distance_squared = direction.squared_length();
  onb uvw;
  uvw.build_from_w(direction);
  return(uvw.local(random_to_sphere(radius,distance_squared)));
}

bool sphere::bounding_box(float t0, float t1, aabb& box) const {
  box = aabb(center - vec3(radius,radius,radius), center + vec3(radius,radius,radius));
  return(true);
}

class moving_sphere: public hitable {
  public:
    moving_sphere() {}
    moving_sphere(vec3 cen0, vec3 cen1, float t0, float t1, float r, material *mat) : 
    center0(cen0), center1(cen1), time0(t0), time1(t1), radius(r), mat_ptr(mat) {};
    virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec, random_gen& rng);
    virtual bool bounding_box(float t0, float t1, aabb& box) const;
    vec3 center(float time) const;
    vec3 center0, center1;
    float time0, time1;
    float radius;
    material *mat_ptr;
};

vec3 moving_sphere::center(float time) const {
  return(center0 + ((time - time0) / (time1 - time0)) * (center1 - center0));
}

bool moving_sphere::bounding_box(float t0, float t1, aabb& box) const {
  aabb box0(center(t0) - vec3(radius,radius,radius), center(t0) + vec3(radius,radius,radius));
  aabb box1(center(t1) - vec3(radius,radius,radius), center(t1) + vec3(radius,radius,radius));
  box = surrounding_box(box0, box1);
  return(true);
}

bool moving_sphere::hit(const ray& r, float t_min, float t_max, hit_record& rec, random_gen& rng) {
  vec3 oc = r.origin() - center(r.time());
  float a = dot(r.direction(), r.direction());
  float b = dot(oc, r.direction()); 
  float c = dot(oc,oc) - radius * radius;
  float discriminant = b * b -  a * c; 
  if(discriminant > 0) {
    float temp = (-b - sqrt(discriminant))/a;
    if(temp < t_max && temp > t_min) {
      rec.t = temp;
      rec.p = r.point_at_parameter(rec.t);
      rec.normal = (rec.p - center(r.time())) / radius;
      get_sphere_uv(rec.normal, rec.u, rec.v);
      rec.mat_ptr = mat_ptr;
      return(true);
    }
    temp = (-b + sqrt(discriminant))/a;
    if(temp < t_max && temp > t_min) {
      rec.t = temp;
      rec.p = r.point_at_parameter(rec.t);
      rec.normal = (rec.p - center(r.time())) / radius;
      get_sphere_uv(rec.normal, rec.u, rec.v);
      rec.mat_ptr = mat_ptr;
      return(true);
    }
  }
  return(false);
}

#endif
