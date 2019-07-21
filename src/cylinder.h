#ifndef CYLINDERH
#define CYLINDERH

#include "hitable.h"
#include "material.h"

class cylinder: public hitable {
public:
  cylinder() {}
  cylinder(float r, float len, material *mat) : radius(r), length(len), mat_ptr(mat) {};
  virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec, random_gen& rng);
  virtual bool bounding_box(float t0, float t1, aabb& box) const;
  virtual float pdf_value(const vec3& o, const vec3& v, random_gen& rng);
  virtual vec3 random(const vec3& o, random_gen& rng);
  void get_cylinder_uv(const vec3& p, float& u, float& v) {
    float phi = atan2(p.z(),p.x());
    // if (phi < 0) phi += 2 * M_PI;
    u = 1 - (phi + M_PI) / (2*M_PI);
    v = (p.y() + length/2)/length;
  };
  float radius;
  float length;
  material *mat_ptr;
};

bool cylinder::hit(const ray& r, float t_min, float t_max, hit_record& rec, random_gen& rng) {
  vec3 oc = r.origin();
  vec3 dir = r.direction();
  dir.e[1] = 0;
  oc.e[1] = 0;
  float a = dot(dir, dir);
  float b = dot(oc, dir); 
  float c = dot(oc, oc) - radius * radius;
  float discriminant = b * b -  a * c; 
  if(discriminant > 0) {
    float temp = (-b - sqrt(discriminant))/a;
    vec3 temppoint = r.point_at_parameter(temp);
    if(temp < t_max && temp > t_min && temppoint.y() > -length/2 && temppoint.y() < length/2){
      rec.t = temp;
      rec.p = temppoint;
      temppoint.e[1] = 0;
      rec.normal = temppoint / radius;
      get_cylinder_uv(rec.p, rec.u, rec.v);
      rec.mat_ptr = mat_ptr;
      return(true);
    }
    temp = (-b + sqrt(discriminant))/a;
    temppoint = r.point_at_parameter(temp);
    if(temp < t_max && temp > t_min && temppoint.y() > -length/2 && temppoint.y() < length/2) {
      rec.t = temp;
      rec.p = temppoint;
      temppoint.e[1] = 0;
      rec.normal = temppoint / radius;
      get_cylinder_uv(rec.p, rec.u, rec.v);
      rec.mat_ptr = mat_ptr;
      return(true);
    }
  }
  return(false);
}

float cylinder::pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, rng)) {
    float cos_theta_max = sqrt(1 - radius * radius/o.squared_length());
    float solid_angle = 2 * M_PI * (1-cos_theta_max);
    return(1/solid_angle);
  } else {
    return(0);
  }
}

vec3 cylinder::random(const vec3& o, random_gen& rng) {
  vec3 direction = -o;
  float distance_squared = direction.squared_length();
  onb uvw;
  uvw.build_from_w(direction);
  return(uvw.local(rng.random_to_sphere(radius,distance_squared)));
}

bool cylinder::bounding_box(float t0, float t1, aabb& box) const {
  box = aabb(-vec3(radius,length/2,radius), vec3(radius,length/2,radius));
  return(true);
}

#endif
