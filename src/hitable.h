#ifndef HITABLEH
#define HITABLEH

#include "aabb.h"
#include "vec3.h"
#include "texture.h"
#include "rng.h"
#include "sampler.h"
#include <Rcpp.h>
#include <memory>

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
  material* mat_ptr;
};

class hitable {
  public:
    virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) = 0;
    virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, Sampler* sampler) = 0;
    virtual bool bounding_box(Float t0, Float t1, aabb& box) const = 0;
    virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng, Float time = 0) {
      return(0.0);
    }
    virtual Float pdf_value(const vec3& o, const vec3& v, Sampler* sampler, Float time = 0) {
      return(0.0);
    }
    virtual vec3 random(const vec3& o, random_gen& rng, Float time = 0) {
      return(vec3(0,1,0));
    }
    virtual vec3 random(const vec3& o, Sampler* sampler, Float time = 0) {
      return(vec3(0,1,0));
    }
    virtual ~hitable() {}
};

class flip_normals : public hitable {
  public:
    flip_normals(std::shared_ptr<hitable> p) : ptr(p) {}
    ~flip_normals() {}
    virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
      if(ptr->hit(r, t_min, t_max, rec, rng)) {
        return(true);
      } else {
        return(false);
      }
    }
    virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
      if(ptr->hit(r, t_min, t_max, rec, sampler)) {
        return(true);
      } else {
        return(false);
      }
    }
    virtual bool bounding_box(Float t0, Float t1, aabb& box) const {
      return(ptr->bounding_box(t0,t1,box));
    }
    Float pdf_value(const vec3& o, const vec3& v, random_gen& rng, Float time = 0) {
      return(ptr->pdf_value(o,v, rng, time));
    }
    Float pdf_value(const vec3& o, const vec3& v, Sampler* sampler, Float time = 0) {
      return(ptr->pdf_value(o,v, sampler, time));
    }
    vec3 random(const vec3& o, random_gen& rng, Float time = 0) {
      return(ptr->random(o, rng, time));
    }
    vec3 random(const vec3& o, Sampler* sampler, Float time = 0) {
      return(ptr->random(o, sampler, time));
    }
    
    std::shared_ptr<hitable> ptr;
};



class translate : public hitable {
public: 
  translate(std::shared_ptr<hitable> p, const vec3& displacement) : ptr(p), offset(displacement) {}
  ~translate() {}
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler);
  
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  Float pdf_value(const vec3& o, const vec3& v, random_gen& rng, Float time = 0) {
    return(ptr->pdf_value(o - offset, v, rng, time));
  }
  Float pdf_value(const vec3& o, const vec3& v, Sampler* sampler, Float time = 0) {
    return(ptr->pdf_value(o - offset, v, sampler, time));
  }
  vec3 random(const vec3& o, random_gen& rng, Float time = 0) {
    return(ptr->random(o - offset, rng, time));
  }
  vec3 random(const vec3& o, Sampler* sampler, Float time = 0) {
    return(ptr->random(o - offset, sampler, time));
  }
  std::shared_ptr<hitable> ptr;
  vec3 offset;
};


class scale : public hitable {
public: 
  scale(std::shared_ptr<hitable> p, const vec3& scale_factor) : ptr(p), scale_factor(scale_factor) {
    inv_scale.e[0] = 1 / scale_factor.x();
    inv_scale.e[1] = 1 / scale_factor.y();
    inv_scale.e[2] = 1 / scale_factor.z();
  }
  ~scale() {}
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler);
  
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  Float pdf_value(const vec3& o, const vec3& v, random_gen& rng, Float time = 0) {
    return(ptr->pdf_value(o * inv_scale, v, rng, time));
  }
  Float pdf_value(const vec3& o, const vec3& v, Sampler* sampler, Float time = 0) {
    return(ptr->pdf_value(o * inv_scale, v, sampler, time));
  }
  vec3 random(const vec3& o, random_gen& rng, Float time = 0) {
    return(ptr->random(o * inv_scale, rng, time));
  }
  vec3 random(const vec3& o, Sampler* sampler, Float time = 0) {
    return(ptr->random(o * inv_scale, sampler, time));
  }
  std::shared_ptr<hitable> ptr;
  vec3 scale_factor;
  vec3 inv_scale;
};

class rotate_y : public hitable {
public:
  rotate_y(std::shared_ptr<hitable> p, Float angle);
  ~rotate_y() {}
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler);
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const {
    box = bbox; 
    return(hasbox);
  }
  Float pdf_value(const vec3& o, const vec3& v, random_gen& rng, Float time = 0) {
    vec3 v2 = v;
    v2.e[0] = cos_theta*v.x() - sin_theta*v.z();
    v2.e[2] = sin_theta*v.x() + cos_theta*v.z();
    vec3 o2 = o;
    o2.e[0] = cos_theta*o.x() - sin_theta*o.z();
    o2.e[2] = sin_theta*o.x() + cos_theta*o.z();
    return(ptr->pdf_value(o2,v2, rng, time));
  }
  Float pdf_value(const vec3& o, const vec3& v, Sampler* sampler, Float time = 0) {
    vec3 v2 = v;
    v2.e[0] = cos_theta*v.x() - sin_theta*v.z();
    v2.e[2] = sin_theta*v.x() + cos_theta*v.z();
    vec3 o2 = o;
    o2.e[0] = cos_theta*o.x() - sin_theta*o.z();
    o2.e[2] = sin_theta*o.x() + cos_theta*o.z();
    return(ptr->pdf_value(o2,v2, sampler, time));
  }
  vec3 random(const vec3& o, random_gen& rng, Float time = 0) {
    vec3 o2 = o;
    o2.e[0] = cos_theta*o.x() - sin_theta*o.z();
    o2.e[2] = sin_theta*o.x() + cos_theta*o.z();
    vec3 temp = ptr->random(o2, rng, time);
    vec3 temp2 = temp;
    temp2.e[0] = cos_theta*temp.x() + sin_theta*temp.z();
    temp2.e[2] = -sin_theta*temp.x() + cos_theta*temp.z(); 
    return(temp2);
  }
  vec3 random(const vec3& o, Sampler* sampler, Float time = 0) {
    vec3 o2 = o;
    o2.e[0] = cos_theta*o.x() - sin_theta*o.z();
    o2.e[2] = sin_theta*o.x() + cos_theta*o.z();
    vec3 temp = ptr->random(o2, sampler, time);
    vec3 temp2 = temp;
    temp2.e[0] = cos_theta*temp.x() + sin_theta*temp.z();
    temp2.e[2] = -sin_theta*temp.x() + cos_theta*temp.z(); 
    return(temp2);
  }
  std::shared_ptr<hitable> ptr;
  Float sin_theta;
  Float cos_theta;
  bool hasbox;
  aabb bbox;
};

class rotate_x : public hitable {
public:
  rotate_x(std::shared_ptr<hitable> p, Float angle);
  ~rotate_x() {}
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler);
  
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const {
    box = bbox; 
    return(hasbox);
  }
  Float pdf_value(const vec3& o, const vec3& v, random_gen& rng, Float time = 0) {
    vec3 v2 = v;
    v2.e[1] = cos_theta*v.y() - sin_theta*v.z();
    v2.e[2] = sin_theta*v.y() + cos_theta*v.z();
    vec3 o2 = o;
    o2.e[1] = cos_theta*o.y() - sin_theta*o.z();
    o2.e[2] = sin_theta*o.y() + cos_theta*o.z();
    return(ptr->pdf_value(o2,v2, rng, time));
  }
  Float pdf_value(const vec3& o, const vec3& v, Sampler* sampler, Float time = 0) {
    vec3 v2 = v;
    v2.e[1] = cos_theta*v.y() - sin_theta*v.z();
    v2.e[2] = sin_theta*v.y() + cos_theta*v.z();
    vec3 o2 = o;
    o2.e[1] = cos_theta*o.y() - sin_theta*o.z();
    o2.e[2] = sin_theta*o.y() + cos_theta*o.z();
    return(ptr->pdf_value(o2,v2, sampler, time));
  }
  vec3 random(const vec3& o, random_gen& rng, Float time = 0) {
    vec3 o2 = o;
    o2.e[1] = cos_theta*o.y() - sin_theta*o.z();
    o2.e[2] = sin_theta*o.y() + cos_theta*o.z();
    vec3 temp = ptr->random(o2, rng, time);
    vec3 temp2 = temp;
    temp2.e[1] = cos_theta*temp.y() + sin_theta*temp.z();
    temp2.e[2] = -sin_theta*temp.y() + cos_theta*temp.z(); 
    return(temp2);
  }
  vec3 random(const vec3& o, Sampler* sampler, Float time = 0) {
    vec3 o2 = o;
    o2.e[1] = cos_theta*o.y() - sin_theta*o.z();
    o2.e[2] = sin_theta*o.y() + cos_theta*o.z();
    vec3 temp = ptr->random(o2, sampler, time);
    vec3 temp2 = temp;
    temp2.e[1] = cos_theta*temp.y() + sin_theta*temp.z();
    temp2.e[2] = -sin_theta*temp.y() + cos_theta*temp.z(); 
    return(temp2);
  }
  std::shared_ptr<hitable> ptr;
  Float sin_theta;
  Float cos_theta;
  bool hasbox;
  aabb bbox;
};

class rotate_z : public hitable {
public:
  rotate_z(std::shared_ptr<hitable> p, Float angle);
  ~rotate_z() {}
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler);
  
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const {
    box = bbox; 
    return(hasbox);
  }
  Float pdf_value(const vec3& o, const vec3& v, random_gen& rng, Float time = 0) {
    vec3 v2 = v;
    v2.e[0] = cos_theta*v.x() - sin_theta*v.y();
    v2.e[1] = sin_theta*v.x() + cos_theta*v.y();
    vec3 o2 = o;
    o2.e[0] = cos_theta*o.x() - sin_theta*o.y();
    o2.e[1] = sin_theta*o.x() + cos_theta*o.y();
    return(ptr->pdf_value(o2,v2, rng, time));
  }
  Float pdf_value(const vec3& o, const vec3& v, Sampler* sampler, Float time = 0) {
    vec3 v2 = v;
    v2.e[0] = cos_theta*v.x() - sin_theta*v.y();
    v2.e[1] = sin_theta*v.x() + cos_theta*v.y();
    vec3 o2 = o;
    o2.e[0] = cos_theta*o.x() - sin_theta*o.y();
    o2.e[1] = sin_theta*o.x() + cos_theta*o.y();
    return(ptr->pdf_value(o2,v2, sampler, time));
  }
  vec3 random(const vec3& o, random_gen& rng, Float time = 0) {
    vec3 o2 = o;
    o2.e[0] = cos_theta*o.x() - sin_theta*o.y();
    o2.e[1] = sin_theta*o.x() + cos_theta*o.y();
    vec3 temp = ptr->random(o2, rng, time);
    vec3 temp2 = temp;
    temp2.e[0] = cos_theta*temp.x() + sin_theta*temp.y();
    temp2.e[1] = -sin_theta*temp.x() + cos_theta*temp.y(); 
    return(temp2);
  }
  vec3 random(const vec3& o, Sampler* sampler, Float time = 0) {
    vec3 o2 = o;
    o2.e[0] = cos_theta*o.x() - sin_theta*o.y();
    o2.e[1] = sin_theta*o.x() + cos_theta*o.y();
    vec3 temp = ptr->random(o2, sampler, time);
    vec3 temp2 = temp;
    temp2.e[0] = cos_theta*temp.x() + sin_theta*temp.y();
    temp2.e[1] = -sin_theta*temp.x() + cos_theta*temp.y(); 
    return(temp2);
  }
  std::shared_ptr<hitable> ptr;
  Float sin_theta;
  Float cos_theta;
  bool hasbox;
  aabb bbox;
};

#endif
