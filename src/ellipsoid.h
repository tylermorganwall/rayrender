#ifndef ELLIPSOIDH
#define ELLIPSOIDH

#include "sphere.h"
#include "material.h"
#include "mathinline.h"

class ellipsoid: public hitable {
  public:
    ellipsoid() {}
    ellipsoid(vec3 cen, Float r, vec3 axes, material *mat, alpha_texture *alpha_mask) : 
      center(cen), radius(r), axes(axes),  mat_ptr(mat), alpha_mask(alpha_mask) {
      inv_axes = vec3(1.0f/axes.x(), 1.0f/axes.y(), 1.0f/axes.z());
      largest_proj_axis = axes.x() * axes.y() * axes.z() / ffmin(axes.x(), ffmin(axes.y(), axes.z()));
    };
    ~ellipsoid() {
      delete mat_ptr;
      delete alpha_mask;
    }
    virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng);
    virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
    virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng);
    virtual vec3 random(const vec3& o, random_gen& rng);
    vec3 center;
    Float radius;
    vec3 axes;
    vec3 inv_axes;
    Float largest_proj_axis;
    material *mat_ptr;
    alpha_texture *alpha_mask;
};


bool ellipsoid::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  ray scaled_ray(r.origin() * inv_axes - center, r.direction() * inv_axes);
  Float a = dot(scaled_ray.direction(), scaled_ray.direction());
  Float b = 2 * dot(scaled_ray.origin(), scaled_ray.direction()); 
  Float c = dot(scaled_ray.origin(),scaled_ray.origin()) - 1;
  Float temp1, temp2;
  if (!quadratic(a, b, c, &temp1, &temp2)) {
    return(false);
  }
  bool is_hit = true;
  bool second_is_hit = true;
  if(alpha_mask) {
    Float u;
    Float v;
    if(temp1 < t_max && temp1 > t_min) {
      vec3 p1 = scaled_ray.point_at_parameter(temp1) ;
      p1 *= 1/p1.length() * axes;; 
      vec3 normal = (p1 - center) * inv_axes;
      normal.make_unit_vector();
      get_sphere_uv(normal, u, v);
      if(alpha_mask->value(u, v, rec.p).x() < rng.unif_rand()) {
        is_hit = false;
      }
    }
    if(temp2 < t_max && temp2 > t_min) {
      vec3 p2 = scaled_ray.point_at_parameter(temp2) ;
      p2 *= 1/p2.length() * axes;; 
      vec3 normal = (p2 - center) * inv_axes;
      normal.make_unit_vector();
      get_sphere_uv(normal, u, v);
      if(alpha_mask->value(u, v, rec.p).x() < rng.unif_rand()) {
        if(!is_hit) {
          return(false);
        }
        second_is_hit = false;
      } 
    }
  }
  if(temp1 < t_max && temp1 > t_min && is_hit) {
    rec.t = temp1;
    rec.p = scaled_ray.point_at_parameter(rec.t) ;
    rec.normal = (rec.p - center);
    rec.p *= 1/rec.p.length() * axes;
    rec.mat_ptr = mat_ptr;
    rec.normal *= inv_axes ;
    rec.normal.make_unit_vector();
    get_sphere_uv(rec.normal, rec.u, rec.v);
    return(true);
  }
  if(temp2 < t_max && temp2 > t_min && second_is_hit) {
    rec.t = temp2;
    rec.p = scaled_ray.point_at_parameter(rec.t) ;
    rec.normal = (rec.p - center);
    rec.p *= 1/rec.p.length() * axes;
    rec.mat_ptr = mat_ptr;
    rec.normal *= inv_axes ;
    if(alpha_mask) {
      rec.normal = -rec.normal;
    }
    rec.normal.make_unit_vector();
    get_sphere_uv(rec.normal, rec.u, rec.v);
    return(true);
  }
  return(false);
}

Float ellipsoid::pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, rng)) {
    Float cos_theta_max = sqrt(1 - 1/(center - o).squared_length());
    Float solid_angle = 2 * M_PI * (1-cos_theta_max) * largest_proj_axis ;
    return(1/solid_angle);
  } else {
    return(0);
  }
}

vec3 ellipsoid::random(const vec3& o, random_gen& rng) {
  vec3 direction = center - o;
  Float distance_squared = direction.squared_length();
  onb uvw;
  uvw.build_from_w(direction);
  return(uvw.local(rng.random_to_sphere(radius,distance_squared) * inv_axes));
}

bool ellipsoid::bounding_box(Float t0, Float t1, aabb& box) const {
  box = aabb(center - vec3(radius,radius,radius) * axes, center + vec3(radius,radius,radius) * axes);
  return(true);
}

#endif
