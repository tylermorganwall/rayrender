#ifndef AABBH
#define AABBH

#include "ray.h"
#include "rng.h"
#include "mathinline.h"

class aabb {
  public: 
    aabb() {}
    aabb(const vec3& a, const vec3& b) { 
      bounds[0] = a;
      bounds[1] = b;
      centroid = (a + b)/2;
      diagonal = b - a;
    }
    
    vec3 min() const {return(bounds[0]);}
    vec3 max() const {return(bounds[1]);}
    
    bool hit(const ray& r, Float tmin, Float tmax, random_gen& rng);
    const vec3 offset(const vec3 o);
    Float surface_area();
    vec3 bounds[2];
    vec3 centroid;
    vec3 diagonal;
};

Float aabb::surface_area() {
  return(2*(diagonal.x() * diagonal.y() + diagonal.x() * diagonal.z() + diagonal.y()*diagonal.z()));
}

bool aabb::hit(const ray &r, Float tmin, Float tmax, random_gen& rng) {
  Float txmin, txmax, tymin, tymax, tzmin, tzmax;
  txmin = (bounds[  r.sign[0]].x()-r.origin().x()) * r.inv_dir.x();
  txmax = (bounds[1-r.sign[0]].x()-r.origin().x()) * r.inv_dir_pad.x();
  tymin = (bounds[  r.sign[1]].y()-r.origin().y()) * r.inv_dir.y();
  tymax = (bounds[1-r.sign[1]].y()-r.origin().y()) * r.inv_dir_pad.y();
  tzmin = (bounds[  r.sign[2]].z()-r.origin().z()) * r.inv_dir.z();
  tzmax = (bounds[1-r.sign[2]].z()-r.origin().z()) * r.inv_dir_pad.z();
  tmin = ffmax(tzmin, ffmax(tymin, ffmax(txmin, tmin)));
  tmax = ffmin(tzmax, ffmin(tymax, ffmin(txmax, tmax)));
  return(tmin <= tmax);
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
  vec3 o = p - min();
  if (max().x() > min().x()) o.e[0] /= max().x() - min().x();
  if (max().y() > min().y()) o.e[1] /= max().y() - min().y();
  if (max().z() > min().z()) o.e[2] /= max().z() - min().z();
  return o;
}

#endif
