#ifndef AABBH
#define AABBH
#include "ray.h"
#include "rng.h"

inline float ffmin(float a, float b) { return(a < b ? a : b);}
inline float ffmax(float a, float b) { return(a > b ? a : b);}


class aabb {
  public: 
    aabb() {}
    aabb(const vec3& a, const vec3& b) { 
      _min = a; 
      _max = b;
      centroid = (a + b)/2;
      diagonal = b - a;
    }
    
    vec3 min() const {return(_min);}
    vec3 max() const {return(_max);}
    
    bool hit(const ray& r, float tmin, float tmax, random_gen& rng);
    const vec3 offset(const vec3 o);
    float surface_area();
    vec3 _min;
    vec3 _max;
    vec3 centroid;
    vec3 diagonal;
};

float aabb::surface_area() {
  return(2*(diagonal.x() * diagonal.y() + diagonal.x() * diagonal.z() + diagonal.y()*diagonal.z()));
}

bool aabb::hit(const ray& r, float tmin, float tmax, random_gen& rng) {
  for(int a = 0; a < 3; a++) {
    float invD = 1.0f / r.direction()[a];
    float t0 = (min()[a] - r.origin()[a]) * invD;
    float t1 = (max()[a] - r.origin()[a]) * invD;
    if(invD < 0.0f) {
      std::swap(t0,t1);
    }
    tmin = t0 > tmin ? t0 : tmin;
    tmax = t1 < tmax ? t1 : tmax;
    if(tmax <= tmin) {
      return(false);
    }
  }
  return(true);
}

aabb surrounding_box(aabb box0, aabb box1) {
  vec3 small(fmin(box0.min().x(), box1.min().x()),
             fmin(box0.min().y(), box1.min().y()),
             fmin(box0.min().z(), box1.min().z()));
  vec3 big(fmax(box0.max().x(), box1.max().x()),
           fmax(box0.max().y(), box1.max().y()),
           fmax(box0.max().z(), box1.max().z()));
  return(aabb(small,big));
}

const vec3 aabb::offset(const vec3 p) {
  vec3 o = p - _min;
  if (_max.x() > _min.x()) o.e[0] /= _max.x() - _min.x();
  if (_max.y() > _min.y()) o.e[1] /= _max.y() - _min.y();
  if (_max.z() > _min.z()) o.e[2] /= _max.z() - _min.z();
  return o;
}

#endif
