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
      bounds[0] = vec3(maxNum, maxNum, maxNum);
      bounds[1] = vec3(minNum, minNum, minNum);
      centroid = vec3(0,0,0);
      diag = vec3(0);
    }
    aabb(vec3 a) {
      bounds[0] = a;
      bounds[1] = a;
      centroid = a;
      diag = vec3(0);
    }
    aabb(const vec3& a, const vec3& b) { 
      bounds[0] = vec3(fmin(a.x(), b.x()), fmin(a.y(), b.y()),fmin(a.z(), b.z()));
      bounds[1] = vec3(fmax(a.x(), b.x()), fmax(a.y(), b.y()),fmax(a.z(), b.z()));
      centroid = (a + b)/2;
      diag = b - a;
    }
    aabb(const aabb &box) {
      bounds[0] = box.bounds[0]; 
      bounds[1] = box.bounds[1];
      centroid = box.centroid;
      diag = box.diag;
    } 
    
    vec3 min() const {return(bounds[0]);}
    vec3 max() const {return(bounds[1]);}
    
    bool hit(const ray& r, Float tmin, Float tmax, random_gen& rng);
    bool hit(const ray& r, Float tmin, Float tmax, Sampler* sampler);
    
    const vec3 offset(const vec3 o);
    Float surface_area();
    Float Volume();
    
    void Expand(Float delta);
    vec3 bounds[2];
    vec3 centroid;
    vec3 diag;
};

inline aabb surrounding_box(aabb box0, aabb box1) {
  vec3 small(fmin(box0.min().x(), box1.min().x()),
             fmin(box0.min().y(), box1.min().y()),
             fmin(box0.min().z(), box1.min().z()));
  vec3 big(fmax(box0.max().x(), box1.max().x()),
           fmax(box0.max().y(), box1.max().y()),
           fmax(box0.max().z(), box1.max().z()));
  return(aabb(small,big));
}

inline aabb Expand(aabb box, Float delta) {
  return(aabb(box.min() - vec3(delta, delta, delta),
              box.max() + vec3(delta, delta, delta)));
}

inline aabb Expand(aabb box, vec3 delta) {
  return(aabb(box.min() - delta,
              box.max() + delta));
}

#endif
