#ifndef VEC2
#define VEC2

#include <iostream>
#include <cmath>
#include "../math/float.h"
#include "../math/vec3.h"
#include "../math/point3.h"

template <typename T> class vec2 {
public:
  vec2() {}
  vec2(T e0, T e1) {xy.x = e0; xy.y = e1;}
  explicit vec2(const vec3<T> &p) { xy.x = p.xyz.x; xy.y = p.xyz.y; }
  explicit vec2(const point3<T> &p) { xy.x = p.xyz.x; xy.y = p.xyz.y; }


  inline const vec2<T>& operator+() const { return *this; }
  inline vec2<T> operator-() const { return vec2<T>(-xy.x, -xy.y); }
  inline T operator[](int i) const { return e[i]; }
  
  inline vec2<T>& operator+=(const vec2<T> &v2);
  inline vec2<T>& operator-=(const vec2<T> &v2);
  inline vec2<T>& operator*=(const vec2<T> &v2);
  inline vec2<T>& operator/=(const vec2<T> &v2);
  inline vec2<T>& operator*=(const Float t);
  inline vec2<T>& operator/=(const Float t);
  bool HasNaNs() const {
    return(std::isnan(xy.x) || std::isnan(xy.y));
  }
  inline Float length() const { return sqrt(xy.x*xy.x + xy.y*xy.y); }
  inline Float squared_length() const { return xy.x*xy.x + xy.y*xy.y; }
  inline void make_unit_vector();
  
  union {
    Float e[2];
    XYstruct<T> xy;
  };
};


template<typename T> 
inline std::istream& operator>>(std::istream &is, vec2<T> &t) {
  is >> t.xy.x >> t.xy.y;
  return is;
}

template<typename T> 
inline std::ostream& operator<<(std::ostream &os, const vec2<T> &t) {
  os << t.xy.x << ", " << t.xy.y;
  return os;
}

template<typename T> 
inline void vec2<T>::make_unit_vector() {
  Float k = 1.0 / sqrt(xy.x*xy.x + xy.y*xy.y);
  xy.x *= k; xy.y *= k; 
}

template<typename T> 
inline vec2<T> operator+(const vec2<T> &v1, const vec2<T> &v2) {
  return vec2<T>(v1.xy.x + v2.xy.x,v1.xy.y + v2.xy.y);
}

template<typename T> 
inline vec2<T> operator-(const vec2<T> &v1, const vec2<T> &v2) {
  return vec2<T>(v1.xy.x - v2.xy.x,v1.xy.y - v2.xy.y);
}

template<typename T> 
inline vec2<T> operator*(const vec2<T> &v1, const vec2<T> &v2) {
  return vec2<T>(v1.xy.x * v2.xy.x,v1.xy.y * v2.xy.y);
}

template<typename T> 
inline vec2<T> operator/(const vec2<T> &v1, const vec2<T> &v2) {
  return vec2<T>(v1.xy.x / v2.xy.x,v1.xy.y / v2.xy.y);
}

template<typename T> 
inline vec2<T> operator*(Float t, const vec2<T> &v) {
  return vec2<T>(t*v.xy.x, t*v.xy.y);
}

template<typename T> 
inline vec2<T> operator*(const vec2<T> &v, Float t) {
  return vec2<T>(t*v.xy.x, t*v.xy.y);
}

template<typename T> 
inline vec2<T> operator/(const vec2<T> &v, Float t) {
  return vec2<T>(v.xy.x/t, v.xy.y/t);
}

template<typename T> 
inline Float dot(const vec2<T> &v1, const vec2<T> &v2) {
  return (v1.xy.x * v2.xy.x + v1.xy.y * v2.xy.y);
}

template<typename T> 
inline vec2<T>& vec2<T>::operator+=(const vec2<T> &v) {
  xy.x += v.xy.x;
  xy.y += v.xy.y;
  return(*this);
}

template<typename T> 
inline vec2<T>& vec2<T>::operator*=(const vec2<T> &v) {
  xy.x *= v.xy.x;
  xy.y *= v.xy.y;
  return(*this);
}

template<typename T> 
inline vec2<T>& vec2<T>::operator/=(const vec2<T> &v) {
  xy.x /= v.xy.x;
  xy.y /= v.xy.y;
  return(*this);
}

template<typename T> 
inline vec2<T>& vec2<T>::operator-=(const vec2<T> &v) {
  xy.x -= v.xy.x;
  xy.y -= v.xy.y;
  return(*this);
}

template<typename T> 
inline vec2<T>& vec2<T>::operator*=(const Float t) {
  xy.x *= t;
  xy.y *= t;
  return(*this);
}

template<typename T> 
inline vec2<T>& vec2<T>::operator/=(const Float t) {
  Float k = 1.0/t;
  xy.x *= k;
  xy.y *= k;
  return(*this);
}

template<typename T> 
inline vec2<T> unit_vector(vec2<T> v) {
  return(v/v.length());
}


typedef vec2<Float> vec2f;
typedef vec2<int>   vec2i;

#endif
