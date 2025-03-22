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
  vec2(T e0, T e1) {x = e0; y = e1;}
  explicit vec2(const vec3<T> &p) { x = p.x; y = p.y; }
  explicit vec2(const point3<T> &p) { x = p.x; y = p.y; }


  inline const vec2<T>& operator+() const { return *this; }
  inline vec2<T> operator-() const { return vec2<T>(-x, -y); }
  inline T operator[](int i) const { return e[i]; }
  
  inline vec2<T>& operator+=(const vec2<T> &v2);
  inline vec2<T>& operator-=(const vec2<T> &v2);
  inline vec2<T>& operator*=(const vec2<T> &v2);
  inline vec2<T>& operator/=(const vec2<T> &v2);
  inline vec2<T>& operator*=(const Float t);
  inline vec2<T>& operator/=(const Float t);
  bool HasNaNs() const {
    return(std::isnan(x) || std::isnan(y));
  }
  inline Float length() const { return sqrt(x*x + y*y); }
  inline Float squared_length() const { return x*x + y*y; }
  inline void make_unit_vector();
  
  union {
    Float e[2];
    struct {
        Float x;
        Float y;
    };
    struct {
        Float u;
        Float v;
    };
  };
};


template<typename T> 
inline std::istream& operator>>(std::istream &is, vec2<T> &t) {
  is >> t.x >> t.y;
  return is;
}

template<typename T> 
inline std::ostream& operator<<(std::ostream &os, const vec2<T> &t) {
  os << t.x << ", " << t.y;
  return os;
}

template<typename T> 
inline void vec2<T>::make_unit_vector() {
  Float k = 1.0 / sqrt(x*x + y*y);
  x *= k; y *= k; 
}

template<typename T> 
inline vec2<T> operator+(const vec2<T> &v1, const vec2<T> &v2) {
  return vec2<T>(v1.x + v2.x,v1.y + v2.y);
}

template<typename T> 
inline vec2<T> operator-(const vec2<T> &v1, const vec2<T> &v2) {
  return vec2<T>(v1.x - v2.x,v1.y - v2.y);
}

template<typename T> 
inline vec2<T> operator*(const vec2<T> &v1, const vec2<T> &v2) {
  return vec2<T>(v1.x * v2.x,v1.y * v2.y);
}

template<typename T> 
inline vec2<T> operator/(const vec2<T> &v1, const vec2<T> &v2) {
  return vec2<T>(v1.x / v2.x,v1.y / v2.y);
}

template<typename T> 
inline vec2<T> operator*(Float t, const vec2<T> &v) {
  return vec2<T>(t*v.x, t*v.y);
}

template<typename T> 
inline vec2<T> operator*(const vec2<T> &v, Float t) {
  return vec2<T>(t*v.x, t*v.y);
}

template<typename T> 
inline vec2<T> operator/(const vec2<T> &v, Float t) {
  return vec2<T>(v.x/t, v.y/t);
}

template<typename T> 
inline Float dot(const vec2<T> &v1, const vec2<T> &v2) {
  return (v1.x * v2.x + v1.y * v2.y);
}

template<typename T> 
inline vec2<T>& vec2<T>::operator+=(const vec2<T> &v) {
  x += v.x;
  y += v.y;
  return(*this);
}

template<typename T> 
inline vec2<T>& vec2<T>::operator*=(const vec2<T> &v) {
  x *= v.x;
  y *= v.y;
  return(*this);
}

template<typename T> 
inline vec2<T>& vec2<T>::operator/=(const vec2<T> &v) {
  x /= v.x;
  y /= v.y;
  return(*this);
}

template<typename T> 
inline vec2<T>& vec2<T>::operator-=(const vec2<T> &v) {
  x -= v.x;
  y -= v.y;
  return(*this);
}

template<typename T> 
inline vec2<T>& vec2<T>::operator*=(const Float t) {
  x *= t;
  y *= t;
  return(*this);
}

template<typename T> 
inline vec2<T>& vec2<T>::operator/=(const Float t) {
  Float k = 1.0/t;
  
  x *= k;
  y *= k;
  return(*this);
}

template<typename T> 
inline vec2<T> unit_vector(vec2<T> v) {
  return(v/v.length());
}


typedef vec2<Float> vec2f;
typedef vec2<int>   vec2i;

#endif
