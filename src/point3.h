#ifndef POINT3H
#define POINT3H

#include <iostream>
#include <cmath>

#include "float.h"
#include "vec3.h" 
#include "normal.h"

template <typename T> class point3 {
public:
  point3() {}
  point3(T e0, T e1, T e2) {e[0] = e0; e[1] = e1; e[2] = e2;}
  point3(vec3<T> e1) {e[0] = e1.x(); e[1] = e1.y(); e[2] = e1.z();}
  point3(T e0) {e[0] = e0; e[1] = e0; e[2] = e0;}
  template <typename U> operator vec3<U>() const {
    return vec3<U>(x(), y(), z());
  }
  point3<T> operator+(const vec3<T> &v) const {
    return point3<T>(e[0] + v.x(), e[1] + v.y(), e[2] + v.z());
  }
  point3<T> &operator+=(const vec3<T> &v) {
    e[0] += v.x(); e[1] += v.y(); e[2] += v.z();
    return *this;
  }
  vec3<T> operator-(const point3<T> &p) const {
    return vec3<T>(x() - p.x(), y() - p.y(), z() - p.z());
  }
  point3<T> operator-(const vec3<T> &v) const {
    return point3<T>(x() - v.x(), y() - v.y(), z() - v.z());
  }
  point3<T> &operator-=(const vec3<T> &v) {
    e[0] -= v.x(); e[1] -= v.y(); e[2] -= v.z();
    return *this;
  }
  inline T x() const { return e[0]; }
  inline T y() const { return e[1]; }
  inline T z() const { return e[2]; }
  inline T r() const { return e[0]; }
  inline T g() const { return e[1]; }
  inline T b() const { return e[2]; }
  
  inline const point3<T>& operator+() const { return *this; }
  inline point3<T> operator-() const { return point3<T>(-e[0], -e[1], -e[2]); }
  inline T operator[](int i) const { return e[i]; }
  inline T& operator[](int i) { return e[i]; }
  
  
  inline point3<T>& operator+=(const point3<T> &v2);
  inline point3<T>& operator-=(const point3<T> &v2);
  inline point3<T>& operator*=(const point3<T> &v2);
  inline point3<T>& operator/=(const point3<T> &v2);
  
  inline point3<T>& operator+=(const Float& t);
  inline point3<T>& operator-=(const Float& t);
  inline point3<T>& operator*=(const Float& t);
  inline point3<T>& operator/=(const Float& t);
  
  bool HasNaNs() const {
    return(std::isnan(e[0]) || std::isnan(e[1]) || std::isnan(e[2]));
  }
  
  inline Float length() const { return std::sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]); }
  inline Float squared_length() const { return e[0]*e[0] + e[1]*e[1] + e[2]*e[2]; }
  inline point3<T> pow(Float exponent) const {
    return(point3<T>(std::pow(e[0],exponent),std::pow(e[1],exponent),std::pow(e[2],exponent)));
  }
  inline void make_unit_vector();
  
  T e[3];
};

template<typename T> 
inline void point3<T>::make_unit_vector() {
  Float k = 1.0 / std::sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
  e[0] *= k; e[1] *= k; e[2] *= k; 
}

template<typename T> 
inline std::istream& operator>>(std::istream &is, point3<T> &t) {
  is >> t.e[0] >> t.e[1] >> t.e[2];
  return is;
}

template<typename T> 
inline std::ostream& operator<<(std::ostream &os, const point3<T> &t) {
  os << t.e[0] << ", " << t.e[1] << ", " << t.e[2];
  return os;
}

template<typename T> 
inline point3<T> operator+(const point3<T> &v1, const point3<T> &v2) {
  return point3<T>(v1.e[0] + v2.e[0],v1.e[1] + v2.e[1],v1.e[2] + v2.e[2]);
}

template<typename T> 
inline point3<T> operator-(const point3<T> &v1, const point3<T> &v2) {
  return point3<T>(v1.e[0] - v2.e[0],v1.e[1] - v2.e[1],v1.e[2] - v2.e[2]);
}

template<typename T> 
inline point3<T> operator*(const point3<T> &v1, const point3<T> &v2) {
  return point3<T>(v1.e[0] * v2.e[0],v1.e[1] * v2.e[1],v1.e[2] * v2.e[2]);
}

template<typename T> 
inline point3<T> operator/(const point3<T> &v1, const point3<T> &v2) {
  return point3<T>(v1.e[0] / v2.e[0],v1.e[1] / v2.e[1],v1.e[2] / v2.e[2]);
}

template<typename T> 
inline point3<T> operator+(const point3<T> &v, Float t) {
  return point3<T>(v.e[0] + t,v.e[1] + t,v.e[2] + t);
}

template<typename T> 
inline point3<T> operator-(const point3<T> &v, Float t) {
  return point3<T>(v.e[0] - t,v.e[1] - t,v.e[2] - t);
}

template<typename T> 
inline point3<T> operator*(Float t, const point3<T> &v) {
  return point3<T>(t*v.e[0], t*v.e[1], t*v.e[2]);
}

template<typename T> 
inline point3<T> operator*(const point3<T> &v, Float t) {
  return point3<T>(t*v.e[0], t*v.e[1], t*v.e[2]);
}

template<typename T> 
inline point3<T> operator/(const point3<T> &v, Float t) {
  return point3<T>(v.e[0]/t, v.e[1]/t, v.e[2]/t);
}

template<typename T> 
inline point3<T> operator/(Float t,const point3<T> &v) {
  return point3<T>(t/v.e[0], t/v.e[1], t/v.e[2]);
}

template<typename T> 
inline Float dot(const point3<T> &v1, const point3<T> &v2) {
  return (v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2]);
}

template<typename T> 
inline point3<T> cross(const point3<T> &v1, const point3<T> &v2) {
  return(point3<T>(DifferenceOfProducts(v1.y(), v2.z(), v1.z(), v2.y()),
                 DifferenceOfProducts(v1.z(), v2.x(), v1.x(), v2.z()),
                 DifferenceOfProducts(v1.x(), v2.y(), v1.y(), v2.x())));
}

template<typename T> 
inline point3<T>& point3<T>::operator+=(const point3<T> &v) {
  e[0] += v.e[0];
  e[1] += v.e[1];
  e[2] += v.e[2];
  return(*this);
}

template<typename T> 
inline point3<T>& point3<T>::operator-=(const point3<T> &v) {
  e[0] -= v.e[0];
  e[1] -= v.e[1];
  e[2] -= v.e[2];
  return(*this);
}

template<typename T> 
inline point3<T>& point3<T>::operator*=(const point3<T> &v) {
  e[0] *= v.e[0];
  e[1] *= v.e[1];
  e[2] *= v.e[2];
  return(*this);
}

template<typename T> 
inline point3<T>& point3<T>::operator/=(const point3<T> &v) {
  e[0] /= v.e[0];
  e[1] /= v.e[1];
  e[2] /= v.e[2];
  return(*this);
}

template<typename T> 
inline point3<T>& point3<T>::operator-=(const Float& t) {
  e[0] -= t;
  e[1] -= t;
  e[2] -= t;
  return(*this);
}

template<typename T> 
inline point3<T>& point3<T>::operator+=(const Float& t) {
  e[0] += t;
  e[1] += t;
  e[2] += t;
  return(*this);
}

template<typename T> 
inline point3<T>& point3<T>::operator*=(const Float& t) {
  e[0] *= t;
  e[1] *= t;
  e[2] *= t;
  return(*this);
}

template<typename T> 
inline point3<T>& point3<T>::operator/=(const Float& t) {
  Float k = 1.0/t;
  
  e[0] *= k;
  e[1] *= k;
  e[2] *= k;
  return(*this);
}

template<typename T> 
inline point3<T> unit_vector(point3<T> v) {
  return(v/v.length());
}

template<typename T> 
inline Float MinComponent(const point3<T> &v) {
  return(std::min(v.x(), std::min(v.y(), v.z())));
}

template<typename T> 
inline Float MaxComponent(const point3<T> &v) {
  return(std::fmax(v.x(), std::fmax(v.y(), v.z())));
}

template<typename T> 
inline int MaxDimension(const point3<T> &v) {
  return((v.x() > v.y()) ? ((v.x() > v.z()) ? 0 : 2) : ((v.y() > v.z()) ? 1 : 2));
}


template<typename T> 
inline point3<T> Min(const point3<T> &p1, const point3<T> &p2) {
  return(point3<T>(fmin(p1.x(), p2.x()), fmin(p1.y(), p2.y()),fmin(p1.z(), p2.z())));
}

template<typename T> 
inline point3<T> Max(const point3<T> &p1, const point3<T> &p2) {
  return(point3<T>(fmax(p1.x(), p2.x()), fmax(p1.y(), p2.y()),fmax(p1.z(), p2.z())));
}

template<typename T> 
inline point3<T> Permute(const point3<T> &v, int x, int y, int z) {
  return(point3<T>(v.e[x], v.e[y], v.e[z]));
}

template<typename T> 
inline point3<T> Abs(const point3<T> &v) {
  return(point3<T>(fabs(v.x()), fabs(v.y()), fabs(v.z())));
}

template <typename T> inline Float
Distance(const point3<T> &p1, const point3<T> &p2) {
  return (p1 - p2).length();
}

template <typename T> inline Float
DistanceSquared(const point3<T> &p1, const point3<T> &p2) {
  return (p1 - p2).squared_length();
}

template <typename T> point3<T>
Lerp(Float t, const point3<T> &p0, const point3<T> &p1) {
  return (1 - t) * p0 + t * p1;
}

template <typename T> point3<T> Floor(const point3<T> &p) {
  return point3<T>(std::floor(p.x()), std::floor(p.y()), std::floor(p.z()));
}

template <typename T> point3<T> Ceil(const point3<T> &p) {
  return point3<T>(std::ceil(p.x()), std::ceil(p.y()), std::ceil(p.z()));
}


template<typename T> 
inline Float dot(const vec3<T> &v1, const point3<T> &v2) {
  return (v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2]);
}

template<typename T> 
inline Float dot(const point3<T> &v1, const vec3<T> &v2) {
  return (v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2]);
}



typedef point3<Float> point3f;
typedef point3<int>   point3i;


#endif
