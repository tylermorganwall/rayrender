#ifndef HITABLEH
#define HITABLEH

#include "aabb.h"
#include "vec3.h"
#include "texture.h"
#include "rng.h"
#include "sampler.h"
#include <Rcpp.h>

class material;

void get_sphere_uv(const vec3& p, Float& u, Float& v);

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

#endif
