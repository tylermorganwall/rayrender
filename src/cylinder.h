#ifndef CYLINDERH
#define CYLINDERH

#include "hitable.h"
#include "material.h"

class cylinder: public hitable {
public:
  cylinder() {}
  cylinder(Float r, Float len, Float phi_min, Float phi_max, 
           material *mat, alpha_texture *alpha_mask) : 
  radius(r), length(len), phi_min(phi_min), phi_max(phi_max), mat_ptr(mat), alpha_mask(alpha_mask) {};
  ~cylinder() {
    delete mat_ptr;
    delete alpha_mask;
  }
  virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng);
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng);
  virtual vec3 random(const vec3& o, random_gen& rng);
  void get_cylinder_uv(const vec3& p, Float& u, Float& v) {
    Float phi = atan2(p.z(),p.x());
    // if (phi < 0) phi += 2 * M_PI;
    u = 1 - (phi + M_PI) / (2*M_PI);
    v = (p.y() + length/2)/length;
  };
  Float radius;
  Float length;
  Float phi_min;
  Float phi_max;
  material *mat_ptr;
  alpha_texture *alpha_mask;
};

bool cylinder::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  vec3 oc = r.origin();
  vec3 dir = r.direction();
  dir.e[1] = 0;
  oc.e[1] = 0;
  Float a = dot(dir, dir);
  Float b = dot(oc, dir); 
  Float c = dot(oc, oc) - radius * radius;
  Float discriminant = b * b -  a * c; 
  if(discriminant > 0) {
    bool is_hit = true;
    bool second_is_hit = true;
    if(alpha_mask) {
      Float temp = (-b - sqrt(discriminant))/a;
      vec3 temppoint = r.point_at_parameter(temp);
      Float phi = atan2(temppoint.z(),temppoint.x());
      phi = phi < 0 ? phi + 2 * M_PI : phi;
      Float u;
      Float v;
      if(temp < t_max && temp > t_min && 
         temppoint.y() > -length/2 && temppoint.y() < length/2 && phi <= phi_max && phi >= phi_min) {
        Float hitRad = std::sqrt(temppoint.x() * temppoint.x() + temppoint.z() * temppoint.z());
        temppoint.e[0] *= radius / hitRad;
        temppoint.e[2] *= radius / hitRad;
        get_cylinder_uv(temppoint, u, v);
        if(alpha_mask->value(u, v, rec.p).x() < rng.unif_rand()) {
          is_hit = false;
        }
      }
      temp = (-b + sqrt(discriminant))/a;
      temppoint = r.point_at_parameter(temp);
      phi = atan2(temppoint.z(),temppoint.x());
      phi = phi < 0 ? phi + 2 * M_PI : phi;
      if(temp < t_max && temp > t_min && 
         temppoint.y() > -length/2 && temppoint.y() < length/2 && phi <= phi_max && phi >= phi_min) {
        Float hitRad = std::sqrt(temppoint.x() * temppoint.x() + temppoint.z() * temppoint.z());
        temppoint.e[0] *= radius / hitRad;
        temppoint.e[2] *= radius / hitRad;
        get_cylinder_uv(temppoint, u, v);
        if(alpha_mask->value(u, v, rec.p).x() < rng.unif_rand()) {
          if(!is_hit) {
            return(false);
          }
          second_is_hit = false;
        } 
      }
    }
    Float temp = (-b - sqrt(discriminant))/a;
    vec3 temppoint = r.point_at_parameter(temp);
    Float phi = atan2(temppoint.z(),temppoint.x());
    phi = phi < 0 ? phi + 2 * M_PI : phi;
    if(is_hit && temp < t_max && temp > t_min && 
       temppoint.y() > -length/2 && temppoint.y() < length/2 && phi <= phi_max && phi >= phi_min) {
      Float hitRad = std::sqrt(temppoint.x() * temppoint.x() + temppoint.z() * temppoint.z());
      temppoint.e[0] *= radius / hitRad;
      temppoint.e[2] *= radius / hitRad;
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
    phi = atan2(temppoint.z(),temppoint.x());
    phi = phi < 0 ? phi + 2 * M_PI : phi;
    if(second_is_hit && temp < t_max && temp > t_min && 
       temppoint.y() > -length/2 && temppoint.y() < length/2 && phi <= phi_max && phi >= phi_min) {
      Float hitRad = std::sqrt(temppoint.x() * temppoint.x() + temppoint.z() * temppoint.z());
      temppoint.e[0] *= radius / hitRad;
      temppoint.e[2] *= radius / hitRad;
      rec.t = temp;
      rec.p = temppoint;
      temppoint.e[1] = 0;
      rec.normal = -temppoint / radius;
      get_cylinder_uv(rec.p, rec.u, rec.v);
      rec.mat_ptr = mat_ptr;
      return(true);
    }
  }
  return(false);
}

Float cylinder::pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, rng)) {
    Float area = length * radius * (phi_max - phi_min);
    Float distance_squared = rec.t * rec.t * v.squared_length();
    Float cosine = fabs(dot(v,rec.normal)/v.length());
    return(distance_squared / (cosine * area));
  } else {
    return(0);
  }
}

vec3 cylinder::random(const vec3& o, random_gen& rng) {
  Float r1 = rng.unif_rand();
  Float y1 = length*(rng.unif_rand()-0.5);
  Float phi = (phi_max - phi_min) * r1 + phi_min;
  Float x = radius * cos(phi);
  Float z = radius * sin(phi);
  return(vec3(x,y1,z)-o);
}

bool cylinder::bounding_box(Float t0, Float t1, aabb& box) const {
  box = aabb(-vec3(radius,length/2,radius), vec3(radius,length/2,radius));
  return(true);
}

#endif
