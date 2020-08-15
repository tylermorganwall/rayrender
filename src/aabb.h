#ifndef AABBH
#define AABBH

#include "ray.h"
#include "rng.h"
#include "mathinline.h"

class aabb {
  public: 
    aabb() {}
    aabb(const vec3& a, const vec3& b) { 
      bounds[0] = vec3(ffmin(a.x(), b.x()), ffmin(a.y(), b.y()),ffmin(a.z(), b.z()));
      bounds[1] = vec3(ffmax(a.x(), b.x()), ffmax(a.y(), b.y()),ffmax(a.z(), b.z()));
      centroid = (a + b)/2;
      diagonal = b - a;
    }
    
    vec3 min() const {return(bounds[0]);}
    vec3 max() const {return(bounds[1]);}
    
    bool hit(const ray& r, Float tmin, Float tmax, random_gen& rng);
    const vec3 offset(const vec3 o);
    Float surface_area();
    void Expand(Float delta);
    vec3 bounds[2];
    vec3 centroid;
    vec3 diagonal;
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

#endif
