#ifndef AABBH
#define AABBH

#include "ray.h"
#include "rng.h"
#include "point3.h"
#include "mathinline.h"
#include "sampler.h"

class aabb {
  public: 
    aabb() {
      Float minNum = std::numeric_limits<Float>::lowest();
      Float maxNum = std::numeric_limits<Float>::max();
      bounds[0] = point3f(maxNum, maxNum, maxNum);
      bounds[1] = point3f(minNum, minNum, minNum);
      centroid = vec3f(0,0,0);
      diag = vec3f(0);
    }
    aabb(vec3f a) {
      bounds[0] = point3f(a.x(),a.y(),a.z());
      bounds[1] = point3f(a.x(),a.y(),a.z());
      centroid  = a;
      diag = vec3f(0);
    }
    aabb(point3f a) {
      bounds[0] = a;
      bounds[1] = a;
      centroid = vec3f(a.x(),a.y(),a.z());
      diag = vec3f(0);
    }
    aabb(const vec3f& a, const vec3f& b) { 
      bounds[0] = point3f(fmin(a.x(), b.x()), fmin(a.y(), b.y()),fmin(a.z(), b.z()));
      bounds[1] = point3f(fmax(a.x(), b.x()), fmax(a.y(), b.y()),fmax(a.z(), b.z()));
      centroid = (a + b)/2;
      diag = b - a;
    }
    aabb(const point3f& a, const point3f& b) { 
      bounds[0] = point3f(fmin(a.x(), b.x()), fmin(a.y(), b.y()),fmin(a.z(), b.z()));
      bounds[1] = point3f(fmax(a.x(), b.x()), fmax(a.y(), b.y()),fmax(a.z(), b.z()));
      point3f temp = (a + b)/2;
      centroid = vec3f(temp.x(),temp.y(),temp.z());
      diag = b - a;
    }
    aabb(const aabb &box) {
      bounds[0] = box.bounds[0]; 
      bounds[1] = box.bounds[1];
      centroid = box.centroid;
      diag = box.diag;
    } 
    
    point3f min() const {return(bounds[0]);}
    point3f max() const {return(bounds[1]);}
    
    
    const Float surface_area() const;
    const Float Volume() const;
    
    bool hit(const ray& r, Float tmin, Float tmax, random_gen& rng);
    bool hit(const ray& r, Float tmin, Float tmax, Sampler* sampler);
    
    const point3f offset(const point3f o);
    const point3f offset(const vec3f o);
    
    const point3f Corner(int corner) const;
    const point3f Lerp(const point3f &t) const;
    const point2f Lerp(const point2f &t) const;
    
    void Expand(Float delta);
    point3f bounds[2];
    vec3f centroid;
    vec3f diag;
};

inline aabb surrounding_box(aabb box0, aabb box1) {
  point3f small(fmin(box0.min().x(), box1.min().x()),
             fmin(box0.min().y(), box1.min().y()),
             fmin(box0.min().z(), box1.min().z()));
  point3f big(fmax(box0.max().x(), box1.max().x()),
           fmax(box0.max().y(), box1.max().y()),
           fmax(box0.max().z(), box1.max().z()));
  return(aabb(small,big));
}

inline aabb surrounding_box(aabb box0, point3f point1) {
  point3f small(fmin(box0.min().x(), point1.x()),
                fmin(box0.min().y(), point1.y()),
                fmin(box0.min().z(), point1.z()));
  point3f big(fmax(box0.max().x(), point1.x()),
              fmax(box0.max().y(), point1.y()),
              fmax(box0.max().z(), point1.z()));
  return(aabb(small,big));
}

inline aabb surrounding_box(aabb box0, vec3f point1) {
  point3f small(fmin(box0.min().x(), point1.x()),
                fmin(box0.min().y(), point1.y()),
                fmin(box0.min().z(), point1.z()));
  point3f big(fmax(box0.max().x(), point1.x()),
              fmax(box0.max().y(), point1.y()),
              fmax(box0.max().z(), point1.z()));
  return(aabb(small,big));
}

inline aabb Expand(aabb box, Float delta) {
  return(aabb(box.min() + -point3f(delta, delta, delta),
              box.max() + point3f(delta, delta, delta)));
}

inline aabb Expand(aabb box, vec3f delta) {
  return(aabb(box.min() - delta,
              box.max() + delta));
}


inline bool Inside(const point3f &p, const aabb &b) {
  return (p.x() >= b.min().x() && p.x() <= b.max().x() && 
          p.y() >= b.min().y() && p.y() <= b.max().y() && 
          p.z() >= b.min().z() && p.z() <= b.max().z());
}

inline bool InsideExclusive(const point3f &p, const aabb &b) {
  return (p.x() >= b.min().x() && p.x()  < b.max().x() && p.y() >= b.min().y() &&
          p.y() <  b.max().y() && p.z() >= b.min().z() && p.z()  < b.max().z());
}

inline std::ostream& operator<<(std::ostream &os, const aabb &t) {
  os << "Low: " << t.bounds[0] <<  " High: " << t.bounds[1];
  return os;
}


#endif
