#ifndef RECTH
#define RECTH

#include "hitable.h"
#include "material.h"

class xy_rect : public hitable {
public:
  xy_rect() {}
  xy_rect(Float _x0, Float _x1, Float _y0, Float _y1, Float _k, 
          std::shared_ptr<material> mat, 
          std::shared_ptr<alpha_texture> alpha_mask, 
          std::shared_ptr<bump_texture> bump_tex, bool flipped) :
    x0(_x0), x1(_x1), y0(_y0), y1(_y1), k(_k), mp(mat), alpha_mask(alpha_mask),
    bump_tex(bump_tex), flipped(flipped) {};
  ~xy_rect() {
    // delete bump_tex;
    // delete alpha_mask;
    // delete mp;
  }
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler);
  
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng, Float time = 0);
  virtual Float pdf_value(const vec3& o, const vec3& v, Sampler* sampler, Float time = 0);
  
  virtual vec3 random(const vec3& o, random_gen& rng, Float time = 0) {
    vec3 random_point = vec3(x0 + rng.unif_rand() * (x1 - x0), y0 + rng.unif_rand() * (y1-y0),k);
    return(random_point - o);
  }
  virtual vec3 random(const vec3& o, Sampler* sampler, Float time = 0) {
    vec2 u = sampler->Get2D();
    vec3 random_point = vec3(x0 + u.x() * (x1 - x0), y0 + u.y() * (y1-y0),k);
    return(random_point - o);
  }
  Float x0, x1, y0, y1, k;
  std::shared_ptr<material> mp;
  std::shared_ptr<alpha_texture> alpha_mask;
  std::shared_ptr<bump_texture> bump_tex;
  bool flipped;
};

class xz_rect : public hitable {
public:
  xz_rect() {}
  xz_rect(Float _x0, Float _x1, Float _z0, Float _z1, Float _k, 
          std::shared_ptr<material> mat, 
          std::shared_ptr<alpha_texture> alpha_mask, 
          std::shared_ptr<bump_texture> bump_tex, bool flipped) :
  x0(_x0), x1(_x1), z0(_z0), z1(_z1), k(_k), mp(mat), alpha_mask(alpha_mask), 
  bump_tex(bump_tex), flipped(flipped) {};
  ~xz_rect() {}
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler);
  
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const {
    box = aabb(vec3(x0,k-0.001,z0), vec3(x1,k+0.001,z1));
    return(true);
  }
  virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng, Float time = 0) {
    hit_record rec;
    if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, rng)) {
      Float area = (x1-x0)*(z1-z0);
      Float distance_squared = rec.t * rec.t * v.squared_length();
      Float cosine = fabs(dot(v, rec.normal)/v.length());
      return(distance_squared / (cosine * area));
    } else {
      return(0);
    }
  }
  virtual Float pdf_value(const vec3& o, const vec3& v, Sampler* sampler, Float time = 0) {
    hit_record rec;
    if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, sampler)) {
      Float area = (x1-x0)*(z1-z0);
      Float distance_squared = rec.t * rec.t * v.squared_length();
      Float cosine = fabs(dot(v, rec.normal)/v.length());
      return(distance_squared / (cosine * area));
    } else {
      return(0);
    }
  }
  virtual vec3 random(const vec3& o, random_gen& rng, Float time = 0) {
    vec3 random_point = vec3(x0 + rng.unif_rand() * (x1 - x0), k, z0 + rng.unif_rand() * (z1-z0));
    return(random_point - o);
  }
  virtual vec3 random(const vec3& o, Sampler* sampler, Float time = 0) {
    vec2 u = sampler->Get2D();
    vec3 random_point = vec3(x0 + u.x() * (x1 - x0), k, z0 + u.y()  * (z1-z0));
    return(random_point - o);
  }
  Float x0, x1, z0, z1, k;
  std::shared_ptr<material> mp;
  std::shared_ptr<alpha_texture> alpha_mask;
  std::shared_ptr<bump_texture> bump_tex;
  bool flipped;
};

class yz_rect : public hitable {
public:
  yz_rect() {}
  yz_rect(Float _y0, Float _y1, Float _z0, Float _z1, Float _k, 
          std::shared_ptr<material> mat, 
          std::shared_ptr<alpha_texture> alpha_mask, std::shared_ptr<bump_texture> bump_tex, bool flipped) :
  y0(_y0), y1(_y1), z0(_z0), z1(_z1), k(_k), mp(mat), alpha_mask(alpha_mask), 
  bump_tex(bump_tex), flipped(flipped) {};
  ~yz_rect() {
    // delete bump_tex;
    // delete alpha_mask;
    // delete mp;
  }
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler);
  
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const {
    box = aabb(vec3(k-0.001,y0,z0), vec3(k+0.001,y1,z1));
    return(true);
  }
  virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng, Float time = 0) {
    hit_record rec;
    if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, rng)) {
      Float area = (y1-y0)*(z1-z0);
      Float distance_squared = rec.t * rec.t * v.squared_length();
      Float cosine = fabs(dot(v,rec.normal)/v.length());
      return(distance_squared / (cosine * area));
    } else {
      return(0);
    }
  }
  virtual Float pdf_value(const vec3& o, const vec3& v, Sampler* sampler, Float time = 0) {
    hit_record rec;
    if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, sampler)) {
      Float area = (y1-y0)*(z1-z0);
      Float distance_squared = rec.t * rec.t * v.squared_length();
      Float cosine = fabs(dot(v,rec.normal)/v.length());
      return(distance_squared / (cosine * area));
    } else {
      return(0);
    }
  }
  virtual vec3 random(const vec3& o, random_gen& rng, Float time = 0) {
    vec3 random_point = vec3(k, y0 + rng.unif_rand() * (y1 - y0), z0 + rng.unif_rand() * (z1-z0));
    return(random_point-o);
  }
  virtual vec3 random(const vec3& o, Sampler* sampler, Float time = 0) {
    vec2 u = sampler->Get2D();
    vec3 random_point = vec3(k, y0 + u.x() * (y1 - y0), z0 + u.y() * (z1-z0));
    return(random_point - o);
  }
  Float y0, y1, z0, z1, k;
  std::shared_ptr<material> mp;
  std::shared_ptr<alpha_texture> alpha_mask;
  std::shared_ptr<bump_texture> bump_tex;
  bool flipped;
};

#endif
