#ifndef POINT2H
#define POINT2H

#include <iostream>
#include <cmath>

#include "../math/float.h"
#include "../math/vec3.h"
#include "../math/point3.h"
#include "../math/vec2.h"

template <typename T> class point2 {
public:
  point2() {}
  point2(T e0, T e1) {e[0] = e0; e[1] = e1;}
  explicit point2(const vec3<T> &p) { e[0] = p.xy.x; e[1] = p.xy.y; }
  explicit point2(const point3<T> &p) { e[0] = p.xy.x; e[1] = p.xy.y; }
  
  inline const point2<T>& operator+() const { return *this; }
  inline point2<T> operator-() const { return point2<T>(-e[0], -e[1]); }

  inline T operator[](int i) const { return e[i]; }
  
  inline point2<T>& operator+=(const point2<T> &v2);
  inline point2<T>& operator-=(const point2<T> &v2);
  inline point2<T>& operator*=(const point2<T> &v2);
  inline point2<T>& operator/=(const point2<T> &v2);
  inline point2<T>& operator*=(const Float t);
  inline point2<T>& operator/=(const Float t);
  inline bool operator==(const point2<T> &v2) const;

  bool HasNaNs() const {
    return(std::isnan(e[0]) || std::isnan(e[1]));
  }
  inline Float length() const { return sqrt(e[0]*e[0] + e[1]*e[1]); }
  inline Float squared_length() const { return e[0]*e[0] + e[1]*e[1]; }
  inline void make_unit_vector();
  
  union {
    Float e[2];
    XYstruct<T> xy;
  };
};


template<typename T> 
inline std::istream& operator>>(std::istream &is, point2<T> &t) {
  is >> t.e[0] >> t.e[1];
  return is;
}

template<typename T> 
inline std::ostream& operator<<(std::ostream &os, const point2<T> &t) {
  os << t.e[0] << ", " << t.e[1];
  return os;
}

template<typename T> 
inline bool point2<T>::operator==(const point2<T> &t) const {
  return xy.x == t.xy.x && xy.y == t.xy.y;
}

template<typename T> 
inline void point2<T>::make_unit_vector() {
  Float k = 1.0 / sqrt(e[0]*e[0] + e[1]*e[1]);
  e[0] *= k; e[1] *= k; 
}

template<typename T> 
inline point2<T> operator+(const point2<T> &v1, const point2<T> &v2) {
  return point2<T>(v1.e[0] + v2.e[0],v1.e[1] + v2.e[1]);
}

template<typename T> 
inline point2<T> operator+(const point2<T> &v1, const vec2<T> &v2) {
  return point2<T>(v1.e[0] + v2.e[0],v1.e[1] + v2.e[1]);
}

template<typename T> 
inline vec2<T> operator-(const point2<T> &v1, const point2<T> &v2) {
  return vec2<T>(v1.e[0] - v2.e[0],v1.e[1] - v2.e[1]);
}

template<typename T> 
inline point2<T> operator-(const point2<T> &v1, const vec2<T> &v2) {
  return point2<T>(v1.e[0] - v2.e[0],v1.e[1] - v2.e[1]);
}

template<typename T> 
inline point2<T> operator/(const point2<T> &v1, const point2<T> &v2) {
  return point2<T>(v1.e[0] / v2.e[0],v1.e[1] / v2.e[1]);
}



template<typename T> 
inline point2<T> operator*(const point2<T> &v, const Float t) {
  return point2<T>(t*v.e[0], t*v.e[1]);
}


template<typename T> 
inline point2<T> operator/(const point2<T> &v, Float t) {
  return point2<T>(v.e[0]/t, v.e[1]/t);
}

template<typename T> 
inline point2<T> operator*(const Float t,const point2<T> &v) {
  return point2<T>(t*v.e[0], t*v.e[1]);
}


template<typename T> 
inline point2<T> operator/(Float t,const point2<T> &v) {
  return point2<T>(t/v.e[0], t/v.e[1]);
}

template<typename T> 
inline Float dot(const point2<T> &v1, const point2<T> &v2) {
  return (v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1]);
}

template<typename T> 
inline point2<T>& point2<T>::operator+=(const point2<T> &v) {
  e[0] += v.e[0];
  e[1] += v.e[1];
  return(*this);
}

template<typename T> 
inline point2<T>& point2<T>::operator*=(const point2<T> &v) {
  e[0] *= v.e[0];
  e[1] *= v.e[1];
  return(*this);
}

template<typename T> 
inline point2<T>& point2<T>::operator/=(const point2<T> &v) {
  e[0] /= v.e[0];
  e[1] /= v.e[1];
  return(*this);
}

template<typename T> 
inline point2<T>& point2<T>::operator-=(const point2<T> &v) {
  e[0] -= v.e[0];
  e[1] -= v.e[1];
  return(*this);
}

template<typename T> 
inline point2<T>& point2<T>::operator*=(const Float t) {
  e[0] *= t;
  e[1] *= t;
  return(*this);
}

template<typename T> 
inline point2<T>& point2<T>::operator/=(const Float t) {
  Float k = 1.0/t;
  
  e[0] *= k;
  e[1] *= k;
  return(*this);
}

template<typename T> 
inline point2<T> unit_vector(point2<T> v) {
  return(v/v.length());
}


typedef point2<Float> point2f;
typedef point2<int>   point2i;

#endif
