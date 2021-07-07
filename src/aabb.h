#ifndef AABBH
#define AABBH

#include "ray.h"
#include "rng.h"
#include "mathinline.h"
#include "sampler.h"

class aabb {
  public: 
    aabb() {
      Float minNum = std::numeric_limits<Float>::lowest();
      Float maxNum = std::numeric_limits<Float>::max();
      bounds[0] = vec3f(maxNum, maxNum, maxNum);
      bounds[1] = vec3f(minNum, minNum, minNum);
      centroid = vec3f(0,0,0);
      diag = vec3f(0);
    }
    aabb(vec3f a) {
      bounds[0] = a;
      bounds[1] = a;
      centroid = a;
      diag = vec3f(0);
    }
    aabb(const vec3f& a, const vec3f& b) { 
      bounds[0] = vec3f(fmin(a.x(), b.x()), fmin(a.y(), b.y()),fmin(a.z(), b.z()));
      bounds[1] = vec3f(fmax(a.x(), b.x()), fmax(a.y(), b.y()),fmax(a.z(), b.z()));
      centroid = (a + b)/2;
      diag = b - a;
    }
    aabb(const aabb &box) {
      bounds[0] = box.bounds[0]; 
      bounds[1] = box.bounds[1];
      centroid = box.centroid;
      diag = box.diag;
    } 
    
    vec3f min() const {return(bounds[0]);}
    vec3f max() const {return(bounds[1]);}
    
    bool hit(const ray& r, Float tmin, Float tmax, random_gen& rng);
    bool hit(const ray& r, Float tmin, Float tmax, Sampler* sampler);
    
    const vec3f offset(const vec3f o);
    Float surface_area();
    Float Volume();
    
    void Expand(Float delta);
    vec3f bounds[2];
    vec3f centroid;
    vec3f diag;
};

inline aabb surrounding_box(aabb box0, aabb box1) {
  vec3f small(fmin(box0.min().x(), box1.min().x()),
             fmin(box0.min().y(), box1.min().y()),
             fmin(box0.min().z(), box1.min().z()));
  vec3f big(fmax(box0.max().x(), box1.max().x()),
           fmax(box0.max().y(), box1.max().y()),
           fmax(box0.max().z(), box1.max().z()));
  return(aabb(small,big));
}

inline aabb Expand(aabb box, Float delta) {
  return(aabb(box.min() - vec3f(delta, delta, delta),
              box.max() + vec3f(delta, delta, delta)));
}

inline aabb Expand(aabb box, vec3f delta) {
  return(aabb(box.min() - delta,
              box.max() + delta));
}

#endif
