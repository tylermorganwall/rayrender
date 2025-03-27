#ifndef BOUNDSH
#define BOUNDSH

#include "../math/point2.h"
#include "../math/point3.h"
#include "../math/mathinline.h"

class Ray;
class random_gen;
class Sampler;

// Bounds Declarations
template <typename T>
class Bounds2 {
public:
  // Bounds2 Public Methods
  Bounds2() {
    T minNum = std::numeric_limits<T>::lowest();
    T maxNum = std::numeric_limits<T>::max();
    pMin = point2<T>(maxNum, maxNum);
    pMax = point2<T>(minNum, minNum);
  }
  explicit Bounds2(const point2<T> &p) : pMin(p), pMax(p) {}
  Bounds2(const point2<T> &p1, const point2<T> &p2) {
    pMin = point2<T>(ffmin(p1.xy.x, p2.xy.x), ffmin(p1.xy.y, p2.xy.y));
    pMax = point2<T>(ffmax(p1.xy.x, p2.xy.x), ffmax(p1.xy.y, p2.xy.y));
  }
  template <typename U>
  explicit operator Bounds2<U>() const {
    return Bounds2<U>((point2<U>)pMin, (point2<U>)pMax);
  }
  
  vec2<T> Diagonal() const { return pMax - pMin; }
  T Area() const {
    vec2<T> d = pMax - pMin;
    return (d.xy.x * d.xy.y);
  }
  int MaximumExtent() const {
    vec2<T> diag = Diagonal();
    if (diag.xy.x > diag.xy.y)
      return 0;
    else
      return 1;
  }
  inline const point2<T> &operator[](int i) const {
    // DCHECK(i == 0 || i == 1);
    return (i == 0) ? pMin : pMax;
  }
  inline point2<T> &operator[](int i) {
    // DCHECK(i == 0 || i == 1);
    return (i == 0) ? pMin : pMax;
  }
  bool operator==(const Bounds2<T> &b) const {
    return b.pMin == pMin && b.pMax == pMax;
  }
  bool operator!=(const Bounds2<T> &b) const {
    return b.pMin != pMin || b.pMax != pMax;
  }
  point2<T> Lerp(const point2f &t) const {
    return point2<T>(lerp(t.xy.x, pMin.xy.x, pMax.xy.x),
                     lerp(t.xy.y, pMin.xy.y, pMax.xy.y));
  }
  vec2<T> Offset(const point2<T> &p) const {
    vec2<T> o = p - pMin;
    if (pMax.xy.x > pMin.xy.x) o.xy.x /= pMax.xy.x - pMin.xy.x;
    if (pMax.xy.y > pMin.xy.y) o.xy.y /= pMax.xy.y - pMin.xy.y;
    return o;
  }
  void BoundingSphere(point2<T> *c, Float *rad) const {
    *c = (pMin + pMax) / 2;
    *rad = Inside(*c, *this) ? Distance(*c, pMax) : 0;
  }
  friend std::ostream &operator<<(std::ostream &os, const Bounds2<T> &b) {
    os << "[ " << b.pMin << " - " << b.pMax << " ]";
    return os;
  }
  
  // Bounds2 Public Data
  point2<T> pMin, pMax;
};

template <typename T>
class Bounds3 {
public:
  // Bounds3 Public Methods
  Bounds3() {
    T minNum = std::numeric_limits<T>::lowest();
    T maxNum = std::numeric_limits<T>::max();
    pMin = point3<T>(maxNum, maxNum, maxNum);
    pMax = point3<T>(minNum, minNum, minNum);
  }
  explicit Bounds3(const point3<T> &p) : pMin(p), pMax(p) {}
  Bounds3(const point3<T> &p1, const point3<T> &p2)
    : pMin(ffmin(p1.xyz.x, p2.xyz.x), ffmin(p1.xyz.y, p2.xyz.y), ffmin(p1.xyz.z, p2.xyz.z)),
      pMax(ffmax(p1.xyz.x, p2.xyz.x), ffmax(p1.xyz.y, p2.xyz.y), ffmax(p1.xyz.z, p2.xyz.z)) {}
  const point3<T> &operator[](int i) const;
  point3<T> &operator[](int i);
  bool operator==(const Bounds3<T> &b) const {
    return b.pMin == pMin && b.pMax == pMax;
  }
  bool operator!=(const Bounds3<T> &b) const {
    return b.pMin != pMin || b.pMax != pMax;
  }
  point3<T> Corner(int corner) const {
    // DCHECK(corner >= 0 && corner < 8);
    return point3<T>((*this)[(corner & 1)].xyz.x,
                     (*this)[(corner & 2) ? 1 : 0].xyz.y,
                     (*this)[(corner & 4) ? 1 : 0].xyz.z);
  }
  vec3<T> Diagonal() const { return pMax - pMin; }
  T SurfaceArea() const {
    vec3<T> d = Diagonal();
    return 2 * (d.xyz.x * d.xyz.y + d.xyz.x * d.xyz.z + d.xyz.y * d.xyz.z);
  }
  T Volume() const {
    vec3<T> d = Diagonal();
    return d.xyz.x * d.xyz.y * d.xyz.z;
  }
  int MaximumExtent() const {
    vec3<T> d = Diagonal();
    if (d.xyz.x > d.xyz.y && d.xyz.x > d.xyz.z)
      return 0;
    else if (d.xyz.y > d.xyz.z)
      return 1;
    else
      return 2;
  }
  point3<T> Lerp(const point3f &t) const {
    return point3<T>(lerp(t.xyz.x, pMin.xyz.x, pMax.xyz.x),
                     lerp(t.xyz.y, pMin.xyz.y, pMax.xyz.y),
                     lerp(t.xyz.z, pMin.xyz.z, pMax.xyz.z));
  }
  vec3<T> Offset(const point3<T> &p) const {
    vec3<T> o = p - pMin;
    if (pMax.xyz.x > pMin.xyz.x) o.xyz.x /= pMax.xyz.x - pMin.xyz.x;
    if (pMax.xyz.y > pMin.xyz.y) o.xyz.y /= pMax.xyz.y - pMin.xyz.y;
    if (pMax.xyz.z > pMin.xyz.z) o.xyz.z /= pMax.xyz.z - pMin.xyz.z;
    return o;
  }
  void BoundingSphere(point3<T> *center, Float *radius) const {
    *center = (pMin + pMax) / 2;
    *radius = Inside(*center, *this) ? Distance(*center, pMax) : 0;
  }
  template <typename U>
  explicit operator Bounds3<U>() const {
    return Bounds3<U>((point3<U>)pMin, (point3<U>)pMax);
  }
  // bool IntersectP(const Ray &ray, Float *hitt0 = nullptr,
  //                 Float *hitt1 = nullptr) const;
  // inline bool IntersectP(const Ray &ray, const vec3f &invDir,
  //                        const int dirIsNeg[3]) const;
  friend std::ostream &operator<<(std::ostream &os, const Bounds3<T> &b) {
    os << "[ " << b.pMin << " - " << b.pMax << " ]";
    return os;
  }
                    
  bool HitP(const Ray &r, Float tMax = Infinity, Float *hitt0 = nullptr,
            Float *hitt1 = nullptr) const;
                  
  // Bounds3 Public Data
  point3<T> pMin, pMax;
};

typedef Bounds2<Float> Bounds2f;
typedef Bounds2<int> Bounds2i;
typedef Bounds3<Float> Bounds3f;
typedef Bounds3<int> Bounds3i;

template <typename T>
point2<T> Min(const point2<T> &pa, const point2<T> &pb) {
  return point2<T>(ffmin(pa.xy.x, pb.xy.x), ffmin(pa.xy.y, pb.xy.y));
}

template <typename T>
point2<T> Max(const point2<T> &pa, const point2<T> &pb) {
  return point2<T>(ffmax(pa.xy.x, pb.xy.x), ffmax(pa.xy.y, pb.xy.y));
}


template <typename T>
Bounds2<T> UnionB(const Bounds2<T> &b, const point2<T> &p) {
  Bounds2<T> ret;
  ret.pMin = Min(b.pMin, p);
  ret.pMax = Max(b.pMax, p);
  return ret;
}

template <typename T>
Bounds2<T> UnionB(const Bounds2<T> &b, const Bounds2<T> &b2) {
  Bounds2<T> ret;
  ret.pMin = Min(b.pMin, b2.pMin);
  ret.pMax = Max(b.pMax, b2.pMax);
  return ret;
}

template <typename T>
Bounds2<T> IntersectB(const Bounds2<T> &b1, const Bounds2<T> &b2) {
  // Important: assign to pMin/pMax directly and don't run the Bounds2()
  // constructor, since it takes min/max of the points passed to it.  In
  // turn, that breaks returning an invalid bound for the case where we
  // intersect non-overlapping bounds (as we'd like to happen).
  Bounds2<T> ret;
  ret.pMin = Max(b1.pMin, b2.pMin);
  ret.pMax = Min(b1.pMax, b2.pMax);
  return ret;
}

template <typename T>
bool Overlaps(const Bounds2<T> &ba, const Bounds2<T> &bb) {
  bool x = (ba.pMax.xyz.x >= bb.pMin.xyz.x) && (ba.pMin.xyz.x <= bb.pMax.xyz.x);
  bool y = (ba.pMax.xyz.y >= bb.pMin.xyz.y) && (ba.pMin.xyz.y <= bb.pMax.xyz.y);
  return (x && y);
}

template <typename T>
bool Inside(const point2<T> &pt, const Bounds2<T> &b) {
  return (pt.xy.x >= b.pMin.xy.x && pt.xy.x <= b.pMax.xy.x && pt.xy.y >= b.pMin.xy.y &&
          pt.xy.y <= b.pMax.xy.y);
}

template <typename T>
bool InsideExclusive(const point2<T> &pt, const Bounds2<T> &b) {
  return (pt.xy.x >= b.pMin.xy.x && pt.xy.x < b.pMax.xy.x && pt.xy.y >= b.pMin.xy.y &&
          pt.xy.y < b.pMax.xy.y);
}

template <typename T>
bool InsideExclusive(const point3<T> &pt, const Bounds3<T> &b) {
  return (pt.xyz.x >= b.pMin.xyz.x && pt.xyz.x < b.pMax.xyz.x && 
          pt.xyz.y >= b.pMin.xyz.y && pt.xyz.y < b.pMax.xyz.y &&
          pt.xyz.z >= b.pMin.xyz.z && pt.xyz.z < b.pMax.xyz.z);
}

template <typename T, typename U>
Bounds2<T> Expand(const Bounds2<T> &b, U delta) {
  return Bounds2<T>(b.pMin - vec2<T>(delta, delta),
                    b.pMax + vec2<T>(delta, delta));
}

#endif