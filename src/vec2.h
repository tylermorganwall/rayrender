#ifndef VEC2
#define VEC2

#include <iostream>
#include <cmath>
#include "float.h"
#include "vec3.h"
#include "point3.h"

template <typename T> class vec2 {
public:
  vec2() {}
  vec2(T e0, T e1) {e[0] = e0; e[1] = e1;}
  explicit vec2<T>(const vec3<T> &p) { e[0] = p.x(); e[1] = p.y(); }
  explicit vec2<T>(const point3<T> &p) { e[0] = p.x(); e[1] = p.y(); }
  
  inline T x() const { return e[0]; }
  inline T y() const { return e[1]; }
  inline T u() const { return e[0]; }
  inline T v() const { return e[1]; }

  inline const vec2<T>& operator+() const { return *this; }
  inline vec2<T> operator-() const { return vec2<T>(-e[0], -e[1]); }
  inline T operator[](int i) const { return e[i]; }
  
  inline vec2<T>& operator+=(const vec2<T> &v2);
  inline vec2<T>& operator-=(const vec2<T> &v2);
  inline vec2<T>& operator*=(const vec2<T> &v2);
  inline vec2<T>& operator/=(const vec2<T> &v2);
  inline vec2<T>& operator*=(const Float t);
  inline vec2<T>& operator/=(const Float t);
  bool HasNaNs() const {
    return(std::isnan(e[0]) || std::isnan(e[1]));
  }
  inline Float length() const { return sqrt(e[0]*e[0] + e[1]*e[1]); }
  inline Float squared_length() const { return e[0]*e[0] + e[1]*e[1]; }
  inline void make_unit_vector();
  
  Float e[2];
};


template<typename T> 
inline std::istream& operator>>(std::istream &is, vec2<T> &t) {
  is >> t.e[0] >> t.e[1];
  return is;
}

template<typename T> 
inline std::ostream& operator<<(std::ostream &os, const vec2<T> &t) {
  os << t.e[0] << ", " << t.e[1];
  return os;
}

template<typename T> 
inline void vec2<T>::make_unit_vector() {
  Float k = 1.0 / sqrt(e[0]*e[0] + e[1]*e[1]);
  e[0] *= k; e[1] *= k; 
}

template<typename T> 
inline vec2<T> operator+(const vec2<T> &v1, const vec2<T> &v2) {
  return vec2<T>(v1.e[0] + v2.e[0],v1.e[1] + v2.e[1]);
}

template<typename T> 
inline vec2<T> operator-(const vec2<T> &v1, const vec2<T> &v2) {
  return vec2<T>(v1.e[0] - v2.e[0],v1.e[1] - v2.e[1]);
}

template<typename T> 
inline vec2<T> operator*(const vec2<T> &v1, const vec2<T> &v2) {
  return vec2<T>(v1.e[0] * v2.e[0],v1.e[1] * v2.e[1]);
}

template<typename T> 
inline vec2<T> operator/(const vec2<T> &v1, const vec2<T> &v2) {
  return vec2<T>(v1.e[0] / v2.e[0],v1.e[1] / v2.e[1]);
}

template<typename T> 
inline vec2<T> operator*(Float t, const vec2<T> &v) {
  return vec2<T>(t*v.e[0], t*v.e[1]);
}

template<typename T> 
inline vec2<T> operator*(const vec2<T> &v, Float t) {
  return vec2<T>(t*v.e[0], t*v.e[1]);
}

template<typename T> 
inline vec2<T> operator/(const vec2<T> &v, Float t) {
  return vec2<T>(v.e[0]/t, v.e[1]/t);
}

template<typename T> 
inline Float dot(const vec2<T> &v1, const vec2<T> &v2) {
  return (v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1]);
}

template<typename T> 
inline vec2<T>& vec2<T>::operator+=(const vec2<T> &v) {
  e[0] += v.e[0];
  e[1] += v.e[1];
  return(*this);
}

template<typename T> 
inline vec2<T>& vec2<T>::operator*=(const vec2<T> &v) {
  e[0] *= v.e[0];
  e[1] *= v.e[1];
  return(*this);
}

template<typename T> 
inline vec2<T>& vec2<T>::operator/=(const vec2<T> &v) {
  e[0] /= v.e[0];
  e[1] /= v.e[1];
  return(*this);
}

template<typename T> 
inline vec2<T>& vec2<T>::operator-=(const vec2<T> &v) {
  e[0] -= v.e[0];
  e[1] -= v.e[1];
  return(*this);
}

template<typename T> 
inline vec2<T>& vec2<T>::operator*=(const Float t) {
  e[0] *= t;
  e[1] *= t;
  return(*this);
}

template<typename T> 
inline vec2<T>& vec2<T>::operator/=(const Float t) {
  Float k = 1.0/t;
  
  e[0] *= k;
  e[1] *= k;
  return(*this);
}

template<typename T> 
inline vec2<T> unit_vector(vec2<T> v) {
  return(v/v.length());
}


typedef vec2<Float> vec2f;
typedef vec2<int>   vec2i;

#endif
