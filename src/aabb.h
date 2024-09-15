#ifndef AABBH
#define AABBH

#include "ray.h"
#include "rng.h"
#include "point3.h"
#include "mathinline.h"
#include "sampler.h"
#include "simd.h"
#include <algorithm>

class aabb {
  public: 
    aabb() {
      Float minNum = std::numeric_limits<Float>::lowest();
      Float maxNum = std::numeric_limits<Float>::max();
      bounds[0] = point3f(maxNum, maxNum, maxNum);
      bounds[1] = point3f(minNum, minNum, minNum);
    }
    aabb(vec3f a) {
      bounds[0] = point3f(a.x(),a.y(),a.z());
      bounds[1] = point3f(a.x(),a.y(),a.z());
    }
    aabb(point3f a) {
      bounds[0] = a;
      bounds[1] = a;
    }
    aabb(const vec3f& a, const vec3f& b) { 
      bounds[0] = point3f(fmin(a.x(), b.x()), fmin(a.y(), b.y()),fmin(a.z(), b.z()));
      bounds[1] = point3f(fmax(a.x(), b.x()), fmax(a.y(), b.y()),fmax(a.z(), b.z()));
    }
    aabb(const point3f& a, const point3f& b) { 
      bounds[0] = point3f(fmin(a.x(), b.x()), fmin(a.y(), b.y()),fmin(a.z(), b.z()));
      bounds[1] = point3f(fmax(a.x(), b.x()), fmax(a.y(), b.y()),fmax(a.z(), b.z()));
    }
    aabb(const aabb &box) {
      bounds[0] = box.bounds[0]; 
      bounds[1] = box.bounds[1];
    } 
    
    point3f min() const {return(bounds[0]);}
    point3f max() const {return(bounds[1]);}
    int MaxDimension() const;
    
    const Float surface_area() const;
    const Float Volume() const;
    const point3f Centroid() const;
    const point3f Diag() const;
    
    const bool hit(const ray& r, Float tmin, Float tmax, random_gen& rng) const;
    const bool hit(const ray& r, Float tmin, Float tmax, Sampler* sampler) const;
    
    const point3f offset(const point3f o) const;
    const point3f offset(const vec3f o) const;
    
    const point3f Corner(int corner) const;
    const point3f Lerp(const point3f &t) const;
    const point2f Lerp(const point2f &t) const;
    
    void Expand(Float delta);
    point3f bounds[2];
};

struct BBox4 {
    union {
        FVec4 corners[6];             // order: minX, minY, minZ, maxX, maxY, maxZ
        float cornersFloat[2][3][4];  // indexed as corner[minOrMax][XYZ][bboxNumber]
        float cornersFloatAlt[6][4];
    };

#if defined(__x86_64__)
    inline const __m128* minCornerSSE() const { return &corners[0].v; }
    inline const __m128* maxCornerSSE() const { return &corners[3].v; }
#elif defined(__aarch64__)
    inline const float32x4_t* minCornerNeon() const { return &corners[0].v; }
    inline const float32x4_t* maxCornerNeon() const { return &corners[3].v; }
#endif

    inline void setBBox(int boxNum, const FVec4& minCorner, const FVec4& maxCorner) {
        cornersFloat[0][0][boxNum] = fmin(minCorner[0], maxCorner[0]);
        cornersFloat[0][1][boxNum] = fmin(minCorner[1], maxCorner[1]);
        cornersFloat[0][2][boxNum] = fmin(minCorner[2], maxCorner[2]);
        cornersFloat[1][0][boxNum] = fmax(minCorner[0], maxCorner[0]);
        cornersFloat[1][1][boxNum] = fmax(minCorner[1], maxCorner[1]);
        cornersFloat[1][2][boxNum] = fmax(minCorner[2], maxCorner[2]);
    }
    BBox4() {}
    
    BBox4(const aabb& a, const aabb& b, const aabb& c, const aabb& d) {
        setBBox(0, a.min(), a.max());
        setBBox(1, b.min(), b.max());
        setBBox(2, c.min(), c.max());
        setBBox(3, d.min(), d.max());
    }
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

void rayBBoxIntersect4(const ray& ray,
                       const BBox4& bbox4,
                       Float tMin,
                       Float tMax,
                       IVec4& hits,
                       FVec4& tMins,
                       FVec4& tMaxs);

#endif
