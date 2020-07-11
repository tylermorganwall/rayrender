#ifndef HITABLEH
#define HITABLEH

#include "aabb.h"
#include "vec3.h"
#include "texture.h"
#include "rng.h"
#include "sampler.h"
#include <Rcpp.h>

class material;

void get_sphere_uv(const vec3& p, Float& u, Float& v) {
  Float phi = atan2(p.z(),p.x());
  Float theta = asin(p.y());
  u = 1 - (phi + M_PI) / (2*M_PI);
  v = (theta + M_PI/2) / M_PI;
}

struct hit_record {
  Float t;
  Float u;
  Float v;
#ifdef DEBUGBVH
  Float bvh_nodes;
#endif
  vec3 p;
  vec3 normal;
  vec3 dpdu, dpdv;
  vec3 bump_normal;
  bool has_bump;
  material *mat_ptr;
};

class hitable {
  public:
    virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) = 0;
    virtual bool bounding_box(Float t0, Float t1, aabb& box) const = 0;
    virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
      return(0.0);
    }
    virtual vec3 random(const vec3& o, random_gen& rng) {
      return(vec3(0,1,0));
    }
    virtual vec3 random(const vec3& o, Sampler* sampler) {
      return(vec3(0,1,0));
    }
    virtual ~hitable() {}
};

class flip_normals : public hitable {
public:
  flip_normals(hitable *p) : ptr(p) {}
  ~flip_normals() {
    //Rcpp::Rcout << "normal delete " << typeid(*ptr).name() << "\n";
    if(ptr) delete ptr;
  }
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
    if(ptr->hit(r, t_min, t_max, rec, rng)) {
      rec.normal = -rec.normal;
      if(rec.has_bump) {
        rec.bump_normal = -rec.bump_normal;
      }
      return(true);
    } else {
      return(false);
    }
  }
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const {
    return(ptr->bounding_box(t0,t1,box));
  }
  Float pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
    return(ptr->pdf_value(o,v, rng));
  }
  vec3 random(const vec3& o, random_gen& rng) {
    return(ptr->random(o, rng));
  }
  vec3 random(const vec3& o, Sampler* sampler) {
    return(ptr->random(o, sampler));
  }
  
  hitable *ptr;
};



class translate : public hitable {
public: 
  translate(hitable *p, const vec3& displacement) : ptr(p), offset(displacement) {}
  ~translate() {
    //Rcpp::Rcout << "translate delete " << typeid(*ptr).name() << "\n";
    if(ptr) delete ptr;
  }
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  Float pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
    return(ptr->pdf_value(o - offset, v, rng));
  }
  vec3 random(const vec3& o, random_gen& rng) {
    return(ptr->random(o - offset, rng));
  }
  vec3 random(const vec3& o, Sampler* sampler) {
    return(ptr->random(o - offset, sampler));
  }
  hitable *ptr;
  vec3 offset;
};

bool translate::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  ray moved_r(r.origin()-offset, r.direction(), r.time());
  if(ptr->hit(moved_r, t_min, t_max, rec, rng)) {
    rec.p += offset;
    return(true);
  } else {
    return(false);
  }
}

bool translate::bounding_box(Float t0, Float t1, aabb& box) const {
  if(ptr->bounding_box(t0,t1,box)) {
    box = aabb(box.min() + offset, box.max() + offset);
    return(true);
  } else {
    return(false);
  }
}

class scale : public hitable {
public: 
  scale(hitable *p, const vec3& scale_factor) : ptr(p), scale_factor(scale_factor) {
    inv_scale.e[0] = 1 / scale_factor.x();
    inv_scale.e[1] = 1 / scale_factor.y();
    inv_scale.e[2] = 1 / scale_factor.z();
  }
  ~scale() {
    //Rcpp::Rcout << "delete scale" << "\n";
    if(ptr) delete ptr;
  }
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  Float pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
    return(ptr->pdf_value(o * inv_scale, v, rng));
  }
  vec3 random(const vec3& o, random_gen& rng) {
    return(ptr->random(o * inv_scale, rng));
  }
  vec3 random(const vec3& o, Sampler* sampler) {
    return(ptr->random(o * inv_scale, sampler));
  }
  hitable *ptr;
  vec3 scale_factor;
  vec3 inv_scale;
};

bool scale::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  ray scaled_r(r.origin() * inv_scale, r.direction() * inv_scale, r.time());
  if(ptr->hit(scaled_r, t_min, t_max, rec, rng)) {
    rec.p *= scale_factor;
    rec.normal *= scale_factor;
    rec.normal.make_unit_vector();
    if(rec.has_bump) {
      rec.bump_normal *= scale_factor;
      rec.bump_normal.make_unit_vector();
    }
    return(true);
  } else {
    return(false);
  }
}

bool scale::bounding_box(Float t0, Float t1, aabb& box) const {
  if(ptr->bounding_box(t0,t1,box)) {
    box = aabb(box.min() * scale_factor, box.max() * scale_factor);
    return(true);
  } else {
    return(false);
  }
}

class rotate_y : public hitable {
public:
  rotate_y(hitable *p, Float angle);
  ~rotate_y() {
    //Rcpp::Rcout << "rotate_y delete " << typeid(*ptr).name() << "\n";
    if(ptr) delete ptr;
  }
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const {
    box = bbox; 
    return(hasbox);
  }
  Float pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
    vec3 v2 = v;
    v2.e[0] = cos_theta*v.x() - sin_theta*v.z();
    v2.e[2] = sin_theta*v.x() + cos_theta*v.z();
    vec3 o2 = o;
    o2.e[0] = cos_theta*o.x() - sin_theta*o.z();
    o2.e[2] = sin_theta*o.x() + cos_theta*o.z();
    return(ptr->pdf_value(o2,v2, rng));
  }
  vec3 random(const vec3& o, random_gen& rng) {
    vec3 o2 = o;
    o2.e[0] = cos_theta*o.x() - sin_theta*o.z();
    o2.e[2] = sin_theta*o.x() + cos_theta*o.z();
    vec3 temp = ptr->random(o2, rng);
    vec3 temp2 = temp;
    temp2.e[0] = cos_theta*temp.x() + sin_theta*temp.z();
    temp2.e[2] = -sin_theta*temp.x() + cos_theta*temp.z(); 
    return(temp2);
  }
  vec3 random(const vec3& o, Sampler* sampler) {
    vec3 o2 = o;
    o2.e[0] = cos_theta*o.x() - sin_theta*o.z();
    o2.e[2] = sin_theta*o.x() + cos_theta*o.z();
    vec3 temp = ptr->random(o2, sampler);
    vec3 temp2 = temp;
    temp2.e[0] = cos_theta*temp.x() + sin_theta*temp.z();
    temp2.e[2] = -sin_theta*temp.x() + cos_theta*temp.z(); 
    return(temp2);
  }
  hitable *ptr;
  Float sin_theta;
  Float cos_theta;
  bool hasbox;
  aabb bbox;
};

rotate_y::rotate_y(hitable *p, Float angle) : ptr(p) {
  Float radians = (M_PI / 180.0) * angle;
  sin_theta = sin(radians);
  cos_theta = cos(radians);
  hasbox = ptr->bounding_box(0,1,bbox);
  vec3 min(FLT_MAX, FLT_MAX, FLT_MAX);
  vec3 max(-FLT_MAX, -FLT_MAX, -FLT_MAX);
  for(int i = 0; i < 2; i++) {
    for(int j = 0; j < 2; j++) {
      for(int k = 0; k < 2; k++) {
        Float x = i*bbox.max().x() + (1-i)*bbox.min().x();
        Float y = j*bbox.max().y() + (1-j)*bbox.min().y();
        Float z = k*bbox.max().z() + (1-k)*bbox.min().z();
        Float newx = cos_theta*x + sin_theta*z;
        Float newz = -sin_theta*x + cos_theta*z;
        vec3 tester(newx,y,newz);
        for(int c = 0; c < 3; c++) {
          if(tester[c] > max[c]) {
            max.e[c] = tester[c];
          }
          if(tester[c] < min[c]) {
            min.e[c] = tester[c];
          }
        }
      } 
    }
  }
  bbox = aabb(min, max);
}

bool rotate_y::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  vec3 origin = r.origin();
  vec3 direction = r.direction();
  origin.e[0] = cos_theta*r.origin()[0] - sin_theta*r.origin()[2];
  origin.e[2] = sin_theta*r.origin()[0] + cos_theta*r.origin()[2];
  direction.e[0] = cos_theta*r.direction()[0] - sin_theta*r.direction()[2];
  direction.e[2] = sin_theta*r.direction()[0] + cos_theta*r.direction()[2];
  ray rotated_r(origin, direction, r.time());
  if(ptr->hit(rotated_r, t_min, t_max, rec, rng)) {
    vec3 p = rec.p;
    vec3 normal = rec.normal;
    p.e[0] = cos_theta*rec.p.e[0] + sin_theta*rec.p.e[2];
    p.e[2] = -sin_theta*rec.p.e[0] + cos_theta*rec.p.e[2]; 
    normal.e[0] = cos_theta*rec.normal.e[0] + sin_theta*rec.normal.e[2];
    normal.e[2] = -sin_theta*rec.normal.e[0] + cos_theta*rec.normal.e[2]; 
    rec.p = p;
    rec.normal = normal;
    if(rec.has_bump) {
      normal = rec.bump_normal;
      normal.e[0] = cos_theta*rec.bump_normal.e[0] + sin_theta*rec.bump_normal.e[2];
      normal.e[2] = -sin_theta*rec.bump_normal.e[0] + cos_theta*rec.bump_normal.e[2]; 
      rec.bump_normal = normal;
    }
    return(true);
  } else {
    return(false);
  }
}

class rotate_x : public hitable {
public:
  rotate_x(hitable *p, Float angle);
  ~rotate_x() {
    //Rcpp::Rcout << "rotate_x delete " << typeid(*ptr).name() << "\n";
    if(ptr) delete ptr;
  }
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const {
    box = bbox; 
    return(hasbox);
  }
  Float pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
    vec3 v2 = v;
    v2.e[1] = cos_theta*v.y() - sin_theta*v.z();
    v2.e[2] = sin_theta*v.y() + cos_theta*v.z();
    vec3 o2 = o;
    o2.e[1] = cos_theta*o.y() - sin_theta*o.z();
    o2.e[2] = sin_theta*o.y() + cos_theta*o.z();
    return(ptr->pdf_value(o2,v2, rng));
  }
  vec3 random(const vec3& o, random_gen& rng) {
    vec3 o2 = o;
    o2.e[1] = cos_theta*o.y() - sin_theta*o.z();
    o2.e[2] = sin_theta*o.y() + cos_theta*o.z();
    vec3 temp = ptr->random(o2, rng);
    vec3 temp2 = temp;
    temp2.e[1] = cos_theta*temp.y() + sin_theta*temp.z();
    temp2.e[2] = -sin_theta*temp.y() + cos_theta*temp.z(); 
    return(temp2);
  }
  vec3 random(const vec3& o, Sampler* sampler) {
    vec3 o2 = o;
    o2.e[1] = cos_theta*o.y() - sin_theta*o.z();
    o2.e[2] = sin_theta*o.y() + cos_theta*o.z();
    vec3 temp = ptr->random(o2, sampler);
    vec3 temp2 = temp;
    temp2.e[1] = cos_theta*temp.y() + sin_theta*temp.z();
    temp2.e[2] = -sin_theta*temp.y() + cos_theta*temp.z(); 
    return(temp2);
  }
  hitable *ptr;
  Float sin_theta;
  Float cos_theta;
  bool hasbox;
  aabb bbox;
};

rotate_x::rotate_x(hitable *p, Float angle) : ptr(p) {
  Float radians = (M_PI / 180.0) * angle;
  sin_theta = sin(radians);
  cos_theta = cos(radians);
  hasbox = ptr->bounding_box(0,1,bbox);
  vec3 min(FLT_MAX, FLT_MAX, FLT_MAX);
  vec3 max(-FLT_MAX, -FLT_MAX, -FLT_MAX);
  for(int i = 0; i < 2; i++) {
    for(int j = 0; j < 2; j++) {
      for(int k = 0; k < 2; k++) {
        Float x = i*bbox.max().x() + (1-i)*bbox.min().x();
        Float y = j*bbox.max().y() + (1-j)*bbox.min().y();
        Float z = k*bbox.max().z() + (1-k)*bbox.min().z();
        Float newy = cos_theta*y + sin_theta*z;
        Float newz = -sin_theta*y + cos_theta*z;
        vec3 tester(x,newy,newz);
        for(int c = 0; c < 3; c++) {
          if(tester[c] > max[c]) {
            max.e[c] = tester[c];
          }
          if(tester[c] < min[c]) {
            min.e[c] = tester[c];
          }
        }
      } 
    }
  }
  bbox = aabb(min, max);
}

bool rotate_x::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  vec3 origin = r.origin();
  vec3 direction = r.direction();
  origin.e[1] = cos_theta*r.origin()[1] - sin_theta*r.origin()[2];
  origin.e[2] = sin_theta*r.origin()[1] + cos_theta*r.origin()[2];
  direction.e[1] = cos_theta*r.direction()[1] - sin_theta*r.direction()[2];
  direction.e[2] = sin_theta*r.direction()[1] + cos_theta*r.direction()[2];
  ray rotated_r(origin, direction, r.time());
  if(ptr->hit(rotated_r, t_min, t_max, rec, rng)) {
    vec3 p = rec.p;
    vec3 normal = rec.normal;
    p.e[1] = cos_theta*rec.p.e[1] + sin_theta*rec.p.e[2];
    p.e[2] = -sin_theta*rec.p.e[1] + cos_theta*rec.p.e[2]; 
    normal.e[1] = cos_theta*rec.normal.e[1] + sin_theta*rec.normal.e[2];
    normal.e[2] = -sin_theta*rec.normal.e[1] + cos_theta*rec.normal.e[2]; 
    rec.p = p;
    rec.normal = normal;
    if(rec.has_bump) {
      normal = rec.bump_normal;
      normal.e[1] = cos_theta*rec.bump_normal.e[1] + sin_theta*rec.bump_normal.e[2];
      normal.e[2] = -sin_theta*rec.bump_normal.e[1] + cos_theta*rec.bump_normal.e[2]; 
      rec.bump_normal = normal;
    }
    return(true);
  } else {
    return(false);
  }
}
class rotate_z : public hitable {
public:
  rotate_z(hitable *p, Float angle);
  ~rotate_z() {
    //Rcpp::Rcout << "rotate_z delete " << typeid(*ptr).name() << "\n";
    if(ptr) delete ptr;
  }
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const {
    box = bbox; 
    return(hasbox);
  }
  Float pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
    vec3 v2 = v;
    v2.e[0] = cos_theta*v.x() - sin_theta*v.y();
    v2.e[1] = sin_theta*v.x() + cos_theta*v.y();
    vec3 o2 = o;
    o2.e[0] = cos_theta*o.x() - sin_theta*o.y();
    o2.e[1] = sin_theta*o.x() + cos_theta*o.y();
    return(ptr->pdf_value(o2,v2, rng));
  }
  vec3 random(const vec3& o, random_gen& rng) {
    vec3 o2 = o;
    o2.e[0] = cos_theta*o.x() - sin_theta*o.y();
    o2.e[1] = sin_theta*o.x() + cos_theta*o.y();
    vec3 temp = ptr->random(o2, rng);
    vec3 temp2 = temp;
    temp2.e[0] = cos_theta*temp.x() + sin_theta*temp.y();
    temp2.e[1] = -sin_theta*temp.x() + cos_theta*temp.y(); 
    return(temp2);
  }
  vec3 random(const vec3& o, Sampler* sampler) {
    vec3 o2 = o;
    o2.e[0] = cos_theta*o.x() - sin_theta*o.y();
    o2.e[1] = sin_theta*o.x() + cos_theta*o.y();
    vec3 temp = ptr->random(o2, sampler);
    vec3 temp2 = temp;
    temp2.e[0] = cos_theta*temp.x() + sin_theta*temp.y();
    temp2.e[1] = -sin_theta*temp.x() + cos_theta*temp.y(); 
    return(temp2);
  }
  hitable *ptr;
  Float sin_theta;
  Float cos_theta;
  bool hasbox;
  aabb bbox;
};

rotate_z::rotate_z(hitable *p, Float angle) : ptr(p) {
  Float radians = (M_PI / 180.0) * angle;
  sin_theta = sin(radians);
  cos_theta = cos(radians);
  hasbox = ptr->bounding_box(0,1,bbox);
  vec3 min(FLT_MAX, FLT_MAX, FLT_MAX);
  vec3 max(-FLT_MAX, -FLT_MAX, -FLT_MAX);
  for(int i = 0; i < 2; i++) {
    for(int j = 0; j < 2; j++) {
      for(int k = 0; k < 2; k++) {
        Float x = i*bbox.max().x() + (1-i)*bbox.min().x();
        Float y = j*bbox.max().y() + (1-j)*bbox.min().y();
        Float z = k*bbox.max().z() + (1-k)*bbox.min().z();
        Float newx = cos_theta*x + sin_theta*y;
        Float newy = -sin_theta*x + cos_theta*y;
        vec3 tester(newx,newy,z);
        for(int c = 0; c < 3; c++) {
          if(tester[c] > max[c]) {
            max.e[c] = tester[c];
          }
          if(tester[c] < min[c]) {
            min.e[c] = tester[c];
          }
        }
      } 
    }
  }
  bbox = aabb(min, max);
}

bool rotate_z::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  vec3 origin = r.origin();
  vec3 direction = r.direction();
  origin.e[0] = cos_theta*r.origin()[0] - sin_theta*r.origin()[1];
  origin.e[1] = sin_theta*r.origin()[0] + cos_theta*r.origin()[1];
  direction.e[0] = cos_theta*r.direction()[0] - sin_theta*r.direction()[1];
  direction.e[1] = sin_theta*r.direction()[0] + cos_theta*r.direction()[1];
  ray rotated_r(origin, direction, r.time());
  if(ptr->hit(rotated_r, t_min, t_max, rec, rng)) {
    vec3 p = rec.p;
    vec3 normal = rec.normal;
    p.e[0] = cos_theta*rec.p.e[0] + sin_theta*rec.p.e[1];
    p.e[1] = -sin_theta*rec.p.e[0] + cos_theta*rec.p.e[1]; 
    normal.e[0] = cos_theta*rec.normal.e[0] + sin_theta*rec.normal.e[1];
    normal.e[1] = -sin_theta*rec.normal.e[0] + cos_theta*rec.normal.e[1]; 
    rec.p = p;
    rec.normal = normal;
    if(rec.has_bump) {
      normal = rec.bump_normal;
      normal.e[0] = cos_theta*rec.bump_normal.e[0] + sin_theta*rec.bump_normal.e[1];
      normal.e[1] = -sin_theta*rec.bump_normal.e[0] + cos_theta*rec.bump_normal.e[1]; 
      rec.bump_normal = normal;
    }
    return(true);
  } else {
    return(false);
  }
}
#endif
