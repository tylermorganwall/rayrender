#ifndef CSGH
#define CSGH

#include "hitable.h"
#include "material.h"
#include "mathinline.h"


class ImplicitShape { 
  public: 
    virtual float getDistance(const vec3& from) const = 0; 
    virtual ~ImplicitShape() {} 
}; 

class csg_sphere : public ImplicitShape {
  public: 
    csg_sphere(const vec3& c, const float& r) : 
      center(c), radius(r) {} 
    float getDistance(const vec3& from) const { 
      return((from - center).length() - radius); 
    } 
    vec3 center;
    float radius; 
}; 

class Plane : public ImplicitShape { 
  public: 
    Plane(const vec3& nn = vec3(0, 1, 0), const vec3& pp = vec3(0)) : 
      n(nn), pointOnPlane(pp) {} 
    float getDistance(const vec3& from) const { 
      return dot(n, from - pointOnPlane); 
    } 
    vec3 n, pointOnPlane; 
}; 


class csg: public hitable {
  public:
    csg() {}
    ~csg() {
      delete mat_ptr;
    }
    csg(vec3 cen, Float r, material *mat) : 
      center(cen), radius(r), mat_ptr(mat) {};
    virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng);
    virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
    virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
      return(1);
    }
    virtual vec3 random(const vec3& o, random_gen& rng) {
      return(vec3(0,1,0));
    }
    virtual vec3 random(const vec3& o, Sampler* sampler) {
      return(vec3(0,1,0));
    }
    vec3 center;
    Float radius;
    material *mat_ptr;
};

#endif