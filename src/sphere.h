#ifndef SPHEREH
#define SPHEREH

#include "hitable.h"
#include "material.h"
#include "mathinline.h"

class sphere: public hitable {
  public:
    sphere() {}
    sphere(vec3 cen, Float r, material *mat) : center(cen), radius(r), mat_ptr(mat) {};
    virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng);
    virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
    virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng);
    virtual vec3 random(const vec3& o, random_gen& rng);
    vec3 center;
    Float radius;
    material *mat_ptr;
};

bool sphere::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  vec3 oc = r.origin() - center;
  Float a = dot(r.direction(), r.direction());
  Float b = 2 * dot(oc, r.direction()); 
  Float c = dot(oc,oc) - radius * radius;
  Float temp1, temp2;
  if (!quadratic(a, b, c, &temp1, &temp2)) {
    return(false);
  }
  if(temp1 < t_max && temp1 > t_min) {
    rec.t = temp1;
    rec.p = r.point_at_parameter(rec.t);
    rec.p *= radius / rec.p.length(); 
    rec.normal = (rec.p - center) / radius;
    get_sphere_uv(rec.normal, rec.u, rec.v);
    rec.mat_ptr = mat_ptr;
    return(true);
  }
  if(temp2 < t_max && temp2 > t_min) {
    rec.t = temp2;
    rec.p = r.point_at_parameter(rec.t);
    rec.p *= radius / rec.p.length(); 
    rec.normal = (rec.p - center) / radius;
    get_sphere_uv(rec.normal, rec.u, rec.v);
    rec.mat_ptr = mat_ptr;
    return(true);
  }
  return(false);
}

Float sphere::pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, rng)) {
    Float cos_theta_max = sqrt(1 - radius * radius/(center - o).squared_length());
    Float solid_angle = 2 * M_PI * (1-cos_theta_max);
    return(1/solid_angle);
  } else {
    return(0);
  }
}

vec3 sphere::random(const vec3& o, random_gen& rng) {
  vec3 direction = center - o;
  Float distance_squared = direction.squared_length();
  onb uvw;
  uvw.build_from_w(direction);
  return(uvw.local(rng.random_to_sphere(radius,distance_squared)));
}

bool sphere::bounding_box(Float t0, Float t1, aabb& box) const {
  box = aabb(center - vec3(radius,radius,radius), center + vec3(radius,radius,radius));
  return(true);
}

class moving_sphere: public hitable {
  public:
    moving_sphere() {}
    moving_sphere(vec3 cen0, vec3 cen1, Float t0, Float t1, Float r, material *mat) : 
    center0(cen0), center1(cen1), time0(t0), time1(t1), radius(r), mat_ptr(mat) {};
    virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng);
    virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
    vec3 center(Float time) const;
    vec3 center0, center1;
    Float time0, time1;
    Float radius;
    material *mat_ptr;
};

vec3 moving_sphere::center(Float time) const {
  return(center0 + ((time - time0) / (time1 - time0)) * (center1 - center0));
}

bool moving_sphere::bounding_box(Float t0, Float t1, aabb& box) const {
  aabb box0(center(t0) - vec3(radius,radius,radius), center(t0) + vec3(radius,radius,radius));
  aabb box1(center(t1) - vec3(radius,radius,radius), center(t1) + vec3(radius,radius,radius));
  box = surrounding_box(box0, box1);
  return(true);
}

bool moving_sphere::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  vec3 oc = r.origin() - center(r.time());
  Float a = dot(r.direction(), r.direction());
  Float b = 2 * dot(oc, r.direction()); 
  Float c = dot(oc,oc) - radius * radius;
  Float temp1, temp2;
  if (!quadratic(a, b, c, &temp1, &temp2)) {
    return(false);
  }
  if(temp1 < t_max && temp1 > t_min) {
    rec.t = temp1;
    rec.p = r.point_at_parameter(rec.t);
    rec.p *= radius / rec.p.length(); 
    rec.normal = (rec.p - center(r.time())) / radius;
    get_sphere_uv(rec.normal, rec.u, rec.v);
    rec.mat_ptr = mat_ptr;
    return(true);
  }
  if(temp2 < t_max && temp2 > t_min) {
    rec.t = temp2;
    rec.p = r.point_at_parameter(rec.t);
    rec.p *= radius / rec.p.length(); 
    rec.normal = (rec.p - center(r.time())) / radius;
    get_sphere_uv(rec.normal, rec.u, rec.v);
    rec.mat_ptr = mat_ptr;
    return(true);
  }
  return(false);
}

#endif
