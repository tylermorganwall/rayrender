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
      bounds[0] = point3f(ffmin(a.x(), b.x()), ffmin(a.y(), b.y()),ffmin(a.z(), b.z()));
      bounds[1] = point3f(ffmax(a.x(), b.x()), ffmax(a.y(), b.y()),ffmax(a.z(), b.z()));
    }
    aabb(const point3f& a, const point3f& b) { 
      bounds[0] = point3f(ffmin(a.x(), b.x()), ffmin(a.y(), b.y()),ffmin(a.z(), b.z()));
      bounds[1] = point3f(ffmax(a.x(), b.x()), ffmax(a.y(), b.y()),ffmax(a.z(), b.z()));
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

struct alignas(16) BBox4 {
    union {
        FVec4 corners[6];             // order: minX, maxX, minY, maxY, minZ, maxZ
        float cornersFloat[2][3][4];  // indexed as corner[minOrMax][XYZ][bboxNumber]
    };

    inline FVec4 getMinX() const { return corners[0]; }
    inline FVec4 getMinY() const { return corners[2]; }
    inline FVec4 getMinZ() const { return corners[4]; }
    inline FVec4 getMaxX() const { return corners[1]; }
    inline FVec4 getMaxY() const { return corners[3]; }
    inline FVec4 getMaxZ() const { return corners[5]; }

    inline void setBBox(int boxNum, const point3f& minCorner, const point3f& maxCorner) {
      // Set minX and maxX
      corners[0].xyzw[boxNum] = ffmin(minCorner.x(), maxCorner.x()); // minX
      corners[1].xyzw[boxNum] = ffmax(minCorner.x(), maxCorner.x()); // maxX

      // Set minY and maxY
      corners[2].xyzw[boxNum] = ffmin(minCorner.y(), maxCorner.y()); // minY
      corners[3].xyzw[boxNum] = ffmax(minCorner.y(), maxCorner.y()); // maxY

      // Set minZ and maxZ
      corners[4].xyzw[boxNum] = ffmin(minCorner.z(), maxCorner.z()); // minZ
      corners[5].xyzw[boxNum] = ffmax(minCorner.z(), maxCorner.z()); // maxZ
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
  point3f small(ffmin(box0.min().x(), box1.min().x()),
             ffmin(box0.min().y(), box1.min().y()),
             ffmin(box0.min().z(), box1.min().z()));
  point3f big(ffmax(box0.max().x(), box1.max().x()),
           ffmax(box0.max().y(), box1.max().y()),
           ffmax(box0.max().z(), box1.max().z()));
  return(aabb(small,big));
}

inline aabb surrounding_box(aabb box0, point3f point1) {
  point3f small(ffmin(box0.min().x(), point1.x()),
                ffmin(box0.min().y(), point1.y()),
                ffmin(box0.min().z(), point1.z()));
  point3f big(ffmax(box0.max().x(), point1.x()),
              ffmax(box0.max().y(), point1.y()),
              ffmax(box0.max().z(), point1.z()));
  return(aabb(small,big));
}

inline aabb surrounding_box(aabb box0, vec3f point1) {
  point3f small(ffmin(box0.min().x(), point1.x()),
                ffmin(box0.min().y(), point1.y()),
                ffmin(box0.min().z(), point1.z()));
  point3f big(ffmax(box0.max().x(), point1.x()),
              ffmax(box0.max().y(), point1.y()),
              ffmax(box0.max().z(), point1.z()));
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
                       FVec4& tEnters
                       );


void rayBBoxIntersect4Serial(const ray& ray,
                       const BBox4& bbox4,
                       Float tMin,
                       Float tMax,
                       IVec4& hits,
                       FVec4& tEnters);
#endif
