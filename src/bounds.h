#ifndef BOUNDSH
#define BOUNDSH

#include "point2.h"
#include "point3.h"
#include "mathinline.h"

class ray;
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
    pMin = point2<T>(ffmin(p1.x, p2.x), ffmin(p1.y, p2.y));
    pMax = point2<T>(ffmax(p1.x, p2.x), ffmax(p1.y, p2.y));
  }
  template <typename U>
  explicit operator Bounds2<U>() const {
    return Bounds2<U>((point2<U>)pMin, (point2<U>)pMax);
  }
  
  vec2<T> Diagonal() const { return pMax - pMin; }
  T Area() const {
    vec2<T> d = pMax - pMin;
    return (d.x * d.y);
  }
  int MaximumExtent() const {
    vec2<T> diag = Diagonal();
    if (diag.x > diag.y)
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
    return point2<T>(lerp(t.x, pMin.x, pMax.x),
                     lerp(t.y, pMin.y, pMax.y));
  }
  vec2<T> Offset(const point2<T> &p) const {
    vec2<T> o = p - pMin;
    if (pMax.x > pMin.x) o.x /= pMax.x - pMin.x;
    if (pMax.y > pMin.y) o.y /= pMax.y - pMin.y;
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
    : pMin(ffmin(p1.x, p2.x), ffmin(p1.y, p2.y), ffmin(p1.z, p2.z)),
      pMax(ffmax(p1.x, p2.x), ffmax(p1.y, p2.y), ffmax(p1.z, p2.z)) {}
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
    return point3<T>((*this)[(corner & 1)].x,
                     (*this)[(corner & 2) ? 1 : 0].y,
                     (*this)[(corner & 4) ? 1 : 0].z);
  }
  vec3<T> Diagonal() const { return pMax - pMin; }
  T SurfaceArea() const {
    vec3<T> d = Diagonal();
    return 2 * (d.x * d.y + d.x * d.z + d.y * d.z);
  }
  T Volume() const {
    vec3<T> d = Diagonal();
    return d.x * d.y * d.z;
  }
  int MaximumExtent() const {
    vec3<T> d = Diagonal();
    if (d.x > d.y && d.x > d.z)
      return 0;
    else if (d.y > d.z)
      return 1;
    else
      return 2;
  }
  point3<T> Lerp(const point3f &t) const {
    return point3<T>(lerp(t.x, pMin.x, pMax.x),
                     lerp(t.y, pMin.y, pMax.y),
                     lerp(t.z, pMin.z, pMax.z));
  }
  vec3<T> Offset(const point3<T> &p) const {
    vec3<T> o = p - pMin;
    if (pMax.x > pMin.x) o.x /= pMax.x - pMin.x;
    if (pMax.y > pMin.y) o.y /= pMax.y - pMin.y;
    if (pMax.z > pMin.z) o.z /= pMax.z - pMin.z;
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
                    
  bool HitP(const ray &r, Float tMax = Infinity, Float *hitt0 = nullptr,
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
  return point2<T>(ffmin(pa.x, pb.x), ffmin(pa.y, pb.y));
}

template <typename T>
point2<T> Max(const point2<T> &pa, const point2<T> &pb) {
  return point2<T>(ffmax(pa.x, pb.x), ffmax(pa.y, pb.y));
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
  bool x = (ba.pMax.x >= bb.pMin.x) && (ba.pMin.x <= bb.pMax.x);
  bool y = (ba.pMax.y >= bb.pMin.y) && (ba.pMin.y <= bb.pMax.y);
  return (x && y);
}

template <typename T>
bool Inside(const point2<T> &pt, const Bounds2<T> &b) {
  return (pt.x >= b.pMin.x && pt.x <= b.pMax.x && pt.y >= b.pMin.y &&
          pt.y <= b.pMax.y);
}

template <typename T>
bool InsideExclusive(const point2<T> &pt, const Bounds2<T> &b) {
  return (pt.x >= b.pMin.x && pt.x < b.pMax.x && pt.y >= b.pMin.y &&
          pt.y < b.pMax.y);
}

template <typename T>
bool InsideExclusive(const point3<T> &pt, const Bounds3<T> &b) {
  return (pt.x >= b.pMin.x && pt.x < b.pMax.x && 
          pt.y >= b.pMin.y && pt.y < b.pMax.y &&
          pt.z >= b.pMin.z && pt.z < b.pMax.z);
}

template <typename T, typename U>
Bounds2<T> Expand(const Bounds2<T> &b, U delta) {
  return Bounds2<T>(b.pMin - vec2<T>(delta, delta),
                    b.pMax + vec2<T>(delta, delta));
}

#endif