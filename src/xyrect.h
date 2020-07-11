#ifndef RECTH
#define RECTH

#include "hitable.h"
#include "material.h"

class xy_rect : public hitable {
public:
  xy_rect() {}
  xy_rect(Float _x0, Float _x1, Float _y0, Float _y1, Float _k, 
          material *mat, alpha_texture* alpha_mask, bump_texture* bump_tex, bool flipped) :
    x0(_x0), x1(_x1), y0(_y0), y1(_y1), k(_k), mp(mat), alpha_mask(alpha_mask),
    bump_tex(bump_tex), flipped(flipped) {};
  ~xy_rect() {
    delete bump_tex;
    delete alpha_mask;
    delete mp;
  }
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const {
    box = aabb(vec3(x0,y0,k-0.001), vec3(x1,y1,k+0.001));
    return(true);
  }
  virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
    hit_record rec;
    if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, rng)) {
      Float area = (x1-x0)*(y1-y0);
      Float distance_squared = rec.t * rec.t * v.squared_length();
      Float cosine = fabs(dot(v,rec.normal)/v.length());
      return(distance_squared / (cosine * area));
    } else {
      return(0);
    }
  }
  virtual vec3 random(const vec3& o, random_gen& rng) {
    vec3 random_point = vec3(x0 + rng.unif_rand() * (x1 - x0), y0 + rng.unif_rand() * (y1-y0),k);
    return(random_point - o);
  }
  virtual vec3 random(const vec3& o, Sampler* sampler) {
    vec2 u = sampler->Get2D();
    vec3 random_point = vec3(x0 + u.x() * (x1 - x0), y0 + u.y() * (y1-y0),k);
    return(random_point - o);
  }
  Float x0, x1, y0, y1, k;
  material *mp;
  alpha_texture *alpha_mask;
  bump_texture *bump_tex;
  bool flipped;
};

class xz_rect : public hitable {
public:
  xz_rect() {}
  xz_rect(Float _x0, Float _x1, Float _z0, Float _z1, Float _k, 
          material *mat, alpha_texture* alpha_mask, bump_texture* bump_tex, bool flipped) :
  x0(_x0), x1(_x1), z0(_z0), z1(_z1), k(_k), mp(mat), alpha_mask(alpha_mask), 
  bump_tex(bump_tex), flipped(flipped) {};
  ~xz_rect() {
    delete bump_tex;
    delete alpha_mask;
    delete mp;
  }
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const {
    box = aabb(vec3(x0,k-0.001,z0), vec3(x1,k+0.001,z1));
    return(true);
  }
  virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
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
  virtual vec3 random(const vec3& o, random_gen& rng) {
    vec3 random_point = vec3(x0 + rng.unif_rand() * (x1 - x0), k, z0 + rng.unif_rand() * (z1-z0));
    return(random_point - o);
  }
  virtual vec3 random(const vec3& o, Sampler* sampler) {
    vec2 u = sampler->Get2D();
    vec3 random_point = vec3(x0 + u.x() * (x1 - x0), k, z0 + u.y()  * (z1-z0));
    return(random_point - o);
  }
  Float x0, x1, z0, z1, k;
  material *mp;
  alpha_texture *alpha_mask;
  bump_texture *bump_tex;
  bool flipped;
};

class yz_rect : public hitable {
public:
  yz_rect() {}
  yz_rect(Float _y0, Float _y1, Float _z0, Float _z1, Float _k, 
          material *mat, alpha_texture* alpha_mask, bump_texture* bump_tex, bool flipped) :
  y0(_y0), y1(_y1), z0(_z0), z1(_z1), k(_k), mp(mat), alpha_mask(alpha_mask), 
  bump_tex(bump_tex), flipped(flipped) {};
  ~yz_rect() {
    delete bump_tex;
    delete alpha_mask;
    delete mp;
  }
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const {
    box = aabb(vec3(k-0.001,y0,z0), vec3(k+0.001,y1,z1));
    return(true);
  }
  virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
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
  virtual vec3 random(const vec3& o, random_gen& rng) {
    vec3 random_point = vec3(k, y0 + rng.unif_rand() * (y1 - y0), z0 + rng.unif_rand() * (z1-z0));
    return(random_point-o);
  }
  virtual vec3 random(const vec3& o, Sampler* sampler) {
    vec2 u = sampler->Get2D();
    vec3 random_point = vec3(k, y0 + u.x() * (y1 - y0), z0 + u.y() * (z1-z0));
    return(random_point - o);
  }
  Float y0, y1, z0, z1, k;
  material *mp;
  alpha_texture *alpha_mask;
  bump_texture *bump_tex;
  bool flipped;
};

bool xy_rect::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  Float t = (k-r.origin().z()) / r.direction().z();
  if(t < t_min || t > t_max) {
    return(false);
  }
  Float x = r.origin().x() + t*r.direction().x();
  Float y = r.origin().y() + t*r.direction().y();
  if(x < x0 || x > x1 || y < y0 || y > y1) {
    return(false);
  }
  Float u = (x-x0)/(x1-x0);
  Float v = (y-y0)/(y1-y0);
  if(flipped) {
    u = 1 - u;
  }
  if(alpha_mask) {
    if(alpha_mask->value(u, v, rec.p).x() < rng.unif_rand()) {
      return(false);
    }
    rec.normal = dot(r.direction(),vec3(0,0,1)) < 0 ? vec3(0,0,1) : vec3(0,0,-1);
  } else {
    rec.normal = flipped ? vec3(0,0,-1) : vec3(0,0,1);
  }
  rec.u = u;
  rec.v = v;
  rec.t = t;
  
  //Interaction information
  rec.dpdu = vec3(-1, 0, 0);
  rec.dpdv = vec3(0, 1, 0);
  rec.has_bump = bump_tex ? true : false;
  
  if(bump_tex) {
    vec3 bvbu = bump_tex->value(u,v, rec.p);
    bvbu.e[0] *= flipped ? -1 : 1;
    rec.bump_normal = rec.normal + bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv; 
    rec.bump_normal.make_unit_vector();
  }
  
  rec.mat_ptr = mp;
  rec.p = r.point_at_parameter(t);
  rec.p.e[2] = k;

  return(true);
}

bool xz_rect::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  Float t = (k-r.origin().y()) / r.direction().y();
  if(t < t_min || t > t_max) {
    return(false);
  }
  Float x = r.origin().x() + t*r.direction().x();
  Float z = r.origin().z() + t*r.direction().z();
  if(x < x0 || x > x1 || z < z0 || z > z1) {
    return(false);
  }
  Float u = 1-(x-x0)/(x1-x0);
  Float v = (z-z0)/(z1-z0);
  if(flipped) {
    u = 1 - u;
  }
  if(alpha_mask) {
    if(alpha_mask->value(u, v, rec.p).x() < rng.unif_rand()) {
      return(false);
    }
    rec.normal =  dot(r.direction(),vec3(0,1,0)) < 0 ? vec3(0,1,0) : vec3(0,-1,0);
  } else {
    rec.normal =  flipped ? vec3(0,-1,0) : vec3(0,1,0);
  }
  rec.u = u;
  rec.v = v;
  rec.t = t;
  
  //Interaction information
  rec.dpdu = vec3(1, 0, 0);
  rec.dpdv = vec3(0, 0, 1);
  rec.has_bump = bump_tex ? true : false;
  
  if(bump_tex) {
    vec3 bvbu = bump_tex->value(u,v, rec.p);
    bvbu.e[0] *= flipped ? -1 : 1;
    rec.bump_normal = rec.normal + bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv; 
    rec.bump_normal.make_unit_vector();
  }
  
  rec.mat_ptr = mp;
  rec.p = r.point_at_parameter(t);
  rec.p.e[1] = k;

  return(true);
}

bool yz_rect::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  Float t = (k-r.origin().x()) / r.direction().x();
  if(t < t_min || t > t_max) {
    return(false);
  }
  Float z = r.origin().z() + t*r.direction().z();
  Float y = r.origin().y() + t*r.direction().y();
  if(z < z0 || z > z1 || y < y0 || y > y1) {
    return(false);
  }
  Float u = 1-(z-z0)/(z1-z0);
  Float v = (y-y0)/(y1-y0);
  if(flipped) {
    u = 1 - u;
  }
  if(alpha_mask) {
    if(alpha_mask->value(u, v, rec.p).x() < rng.unif_rand()) {
      return(false);
    }
    rec.normal =  dot(r.direction(),vec3(1,0,0)) < 0 ? vec3(1,0,0) : vec3(-1,0,0);
  } else {
    rec.normal =  flipped ? vec3(-1,0,0) : vec3(1,0,0);
  }
  rec.u = u;
  rec.v = v;
  rec.t = t;
  
  //Interaction information
  rec.dpdu = vec3(0, 0, 1);
  rec.dpdv = vec3(0, 1, 0);
  rec.has_bump = bump_tex ? true : false;
  if(bump_tex) {
    vec3 bvbu = bump_tex->value(u,v, rec.p);
    bvbu.e[0] *= flipped ? -1 : 1;
    rec.bump_normal = rec.normal + bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv; 
    rec.bump_normal.make_unit_vector();
  }
  
  rec.mat_ptr = mp;
  rec.p = r.point_at_parameter(t);
  rec.p.e[0] = k;

  return(true);
}

#endif
