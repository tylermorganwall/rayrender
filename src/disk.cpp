#include "disk.h"
#include "raylog.h"

const bool disk::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("Disk");
  
  ray r2 = (*WorldToObject)(r);
  // First we intersect with the plane containing the disk
  Float t = -r2.origin().y() / r2.direction().y();
  bool alpha_miss = false;
  
  if(t < t_min || t > t_max) {
    return(false);
  }
  Float x = r2.origin().x() + t*r2.direction().x();
  Float z = r2.origin().z() + t*r2.direction().z();
  Float radHit2 = x*x + z*z;
  if(radHit2 >= radius * radius || radHit2 <= inner_radius * inner_radius) {
    return(false);
  }
  
  
  point3f p = r2.point_at_parameter(t);
  p.e[1] = 0;

  Float u = p.x() / (2.0 * radius) + 0.5;
  Float v = p.z() / (2.0 * radius) + 0.5;
  u = 1 - u;
  if(alpha_mask) {
    if(alpha_mask->value(u, v, rec.p) < rng.unif_rand()) {
      alpha_miss = true;
    }
    rec.normal =  dot(r2.direction(),normal3f(0,1,0)) < 0 ? normal3f(0,1,0) : normal3f(0,-1,0);
  } else {
    rec.normal = normal3f(0,1,0);  
  }
  rec.p = p;
  
  rec.t = t;
  rec.mat_ptr = mat_ptr.get();
  rec.u = u;
  rec.v = v;

  //Interaction information
  rec.dpdu = vec3f(1, 0, 0);
  rec.dpdv = vec3f(0, 0, 1);
  rec.has_bump = bump_tex ? true : false;
  
  if(bump_tex) {
    point3f bvbu = bump_tex->value(rec.u,rec.v, rec.p);
    rec.bump_normal = rec.normal + convert_to_normal3(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv); 
    rec.bump_normal.make_unit_vector();

  }
  rec.pError = vec3f(0,0,0);
  rec = (*ObjectToWorld)(rec);
  if(!alpha_mask) {
    rec.normal *= reverseOrientation  ? -1 : 1;
    rec.bump_normal *= reverseOrientation  ? -1 : 1;
  }
  rec.shape = this;
  rec.alpha_miss = alpha_miss;
  
  return(true);
}


const bool disk::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("Disk");
  
  ray r2 = (*WorldToObject)(r);
  // First we intersect with the plane containing the disk
  Float t = -r2.origin().y() / r2.direction().y();
  bool alpha_miss = false;
  
  if(t < t_min || t > t_max) {
    return(false);
  }
  Float x = r2.origin().x() + t*r2.direction().x();
  Float z = r2.origin().z() + t*r2.direction().z();
  Float radHit2 = x*x + z*z;
  if(radHit2 >= radius * radius || radHit2 <= inner_radius * inner_radius) {
    return(false);
  }
  
  
  point3f p = r2.point_at_parameter(t);
  p.e[1] = 0;

  Float u = p.x() / (2.0 * radius) + 0.5;
  Float v = p.z() / (2.0 * radius) + 0.5;
  u = 1 - u;
  if(alpha_mask) {
    if(alpha_mask->value(u, v, rec.p) < sampler->Get1D()) {
      alpha_miss = true;
    }
    rec.normal =  dot(r2.direction(),normal3f(0,1,0)) < 0 ? normal3f(0,1,0) : normal3f(0,-1,0);
  } else {
    rec.normal = normal3f(0,1,0);  
  }
  rec.p = p;

  rec.t = t;
  rec.mat_ptr = mat_ptr.get();
  rec.u = u;
  rec.v = v;
  
  //Interaction information
  rec.dpdu = vec3f(1, 0, 0);
  rec.dpdv = vec3f(0, 0, 1);
  rec.has_bump = bump_tex ? true : false;
  
  if(bump_tex) {
    point3f bvbu = bump_tex->value(rec.u,rec.v, rec.p);
    rec.bump_normal = rec.normal + convert_to_normal3(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv); 
    rec.bump_normal.make_unit_vector();
  }
  rec.pError = vec3f(0,0,0);
  rec = (*ObjectToWorld)(rec);
  if(!alpha_mask) {
    rec.normal *= reverseOrientation  ? -1 : 1;
    rec.bump_normal *= reverseOrientation  ? -1 : 1;
  }
  rec.shape = this;
  rec.alpha_miss = alpha_miss;
  
  return(true);
}


bool disk::HitP(const ray& r, Float t_min, Float t_max, random_gen& rng) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("Disk");
  
  ray r2 = (*WorldToObject)(r);
  // First we intersect with the plane containing the disk
  Float t = -r2.origin().y() / r2.direction().y();
  
  if(t < t_min || t > t_max) {
    return(false);
  }
  Float x = r2.origin().x() + t*r2.direction().x();
  Float z = r2.origin().z() + t*r2.direction().z();
  Float radHit2 = x*x + z*z;
  if(radHit2 >= radius * radius || radHit2 <= inner_radius * inner_radius) {
    return(false);
  }
  return(true);
}


bool disk::HitP(const ray& r, Float t_min, Float t_max, Sampler* sampler) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("Disk");
  
  ray r2 = (*WorldToObject)(r);
  // First we intersect with the plane containing the disk
  Float t = -r2.origin().y() / r2.direction().y();
  
  if(t < t_min || t > t_max) {
    return(false);
  }
  Float x = r2.origin().x() + t*r2.direction().x();
  Float z = r2.origin().z() + t*r2.direction().z();
  Float radHit2 = x*x + z*z;
  if(radHit2 >= radius * radius || radHit2 <= inner_radius * inner_radius) {
    return(false);
  }
  return(true);
}



Float disk::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, rng)) {
    Float area =  static_cast<Float>(M_PI) * (radius * radius - inner_radius * inner_radius);
    Float distance_squared = rec.t * rec.t * v.squared_length();
    Float cosine = fabs(dot(v,rec.normal)/v.length());
    return(distance_squared / (cosine * area));
  } else {
    return(0);
  }
}

Float disk::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, sampler)) {
    Float area =  static_cast<Float>(M_PI) * (radius * radius - inner_radius * inner_radius);
    Float distance_squared = rec.t * rec.t * v.squared_length();
    Float cosine = fabs(dot(v,rec.normal)/v.length());
    return(distance_squared / (cosine * area));
  } else {
    return(0);
  }
}

vec3f disk::random(const point3f& o, random_gen& rng, Float time) {
  Float r1 = rng.unif_rand();
  Float r2 = sqrt(rng.unif_rand());
  Float phi = 2 * static_cast<Float>(M_PI) * r1;
  Float x = ((radius - inner_radius) * r2 + inner_radius) * cos(phi);
  Float z = ((radius - inner_radius) * r2 + inner_radius) * sin(phi);
  return((*ObjectToWorld)(point3f(x,0,z))+center-o);
}

vec3f disk::random(const point3f& o, Sampler* sampler, Float time) {
  vec2f u = sampler->Get2D();
  Float r1 = u.x();
  Float r2 = sqrt(u.y());
  Float phi = 2 * static_cast<Float>(M_PI) * r1;
  Float x = ((radius - inner_radius) * r2 + inner_radius) * cos(phi);
  Float z = ((radius - inner_radius) * r2 + inner_radius) * sin(phi);
  return((*ObjectToWorld)(point3f(x,0,z))+center-o);
}

bool disk::bounding_box(Float t0, Float t1, aabb& box) const {
  box = (*ObjectToWorld)(aabb(-point3f(radius,0.001,radius), point3f(radius,0.001,radius)));
  return(true);
}
