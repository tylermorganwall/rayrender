#ifndef DISKH
#define DISKH

#include "hitable.h"
#include "onbh.h"
#include "material.h"

class disk : public hitable {
public:
  disk() {}
  disk(vec3 cen, Float r, Float i_r, material *mat, alpha_texture *alpha_mask, 
       bump_texture* bump_tex) : center(cen), radius(r), 
       inner_radius(i_r), mat_ptr(mat), alpha_mask(alpha_mask), bump_tex(bump_tex) {};
  ~disk() {
    delete mat_ptr;
    delete alpha_mask;
    delete bump_tex;
  }
  virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng);
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng);
  virtual vec3 random(const vec3& o, random_gen& rng);
  virtual vec3 random(const vec3& o, Sampler* sampler);
  
  vec3 center;
  Float radius;
  Float inner_radius;
  material *mat_ptr;
  alpha_texture *alpha_mask;
  bump_texture *bump_tex;
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
  Float radHit2 = x*x + z*z;
  if(radHit2 >= radius * radius || radHit2 <= inner_radius * inner_radius) {
    return(false);
  }
  
  
  vec3 p = r.point_at_parameter(t);
  p.e[1] = 0;
  
  Float u = p.x() / (2.0 * radius) + 0.5;
  Float v = p.z() / (2.0 * radius) + 0.5;
  u = 1 - u;
  if(alpha_mask) {
    if(alpha_mask->value(u, v, rec.p).x() < 1) {
      return(false);
    }
  }
  rec.p = p;
  rec.normal = n;
  rec.t = t;
  rec.mat_ptr = mat_ptr;
  rec.u = u;
  rec.v = v;
  
  //Interaction information
  rec.dpdu = vec3(1, 0, 0);
  rec.dpdv = vec3(0, 0, 1);
  rec.has_bump = bump_tex ? true : false;
  
  if(bump_tex) {
    vec3 bvbu = bump_tex->value(rec.u,rec.v, rec.p);
    rec.bump_normal = rec.normal + bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv; 
    rec.bump_normal.make_unit_vector();
  }
  
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

vec3 disk::random(const vec3& o, Sampler* sampler) {
  vec2 u = sampler->Get2D();
  Float r1 = u.x();
  Float r2 = sqrt(u.y());
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
