#ifndef VEC2
#define VEC2

#include <iostream>
#include <cmath>

#ifndef FLOATDEF
#define FLOATDEF
#ifdef RAY_FLOAT_AS_DOUBLE
typedef double Float;
#else
typedef float Float;
#endif 
#endif

#include "vec3.h"

class vec2 {
public:
  vec2() {}
  vec2(Float e0, Float e1) {e[0] = e0; e[1] = e1;}
  explicit vec2(const vec3 &p) { e[0] = p.x(); e[1] = p.y(); }
  inline Float x() const { return e[0]; }
  inline Float y() const { return e[1]; }
  inline Float u() const { return e[0]; }
  inline Float v() const { return e[1]; }

  inline const vec2& operator+() const { return *this; }
  inline vec2 operator-() const { return vec2(-e[0], -e[1]); }
  inline Float operator[](int i) const { return e[i]; }
  
  inline vec2& operator+=(const vec2 &v2);
  inline vec2& operator-=(const vec2 &v2);
  inline vec2& operator*=(const vec2 &v2);
  inline vec2& operator/=(const vec2 &v2);
  inline vec2& operator*=(const Float t);
  inline vec2& operator/=(const Float t);
  bool HasNaNs() const {
    return(std::isnan(e[0]) || std::isnan(e[1]));
  }
  inline Float length() const { return sqrt(e[0]*e[0] + e[1]*e[1]); }
  inline Float squared_length() const { return e[0]*e[0] + e[1]*e[1]; }
  inline void make_unit_vector();
  
  Float e[2];
};


inline std::istream& operator>>(std::istream &is, vec2 &t) {
  is >> t.e[0] >> t.e[1];
  return is;
}

inline std::ostream& operator<<(std::ostream &os, const vec2 &t) {
  os << t.e[0] << ", " << t.e[1];
  return os;
}

inline void vec2::make_unit_vector() {
  Float k = 1.0 / sqrt(e[0]*e[0] + e[1]*e[1]);
  e[0] *= k; e[1] *= k; 
}

inline vec2 operator+(const vec2 &v1, const vec2 &v2) {
  return vec2(v1.e[0] + v2.e[0],v1.e[1] + v2.e[1]);
}

inline vec2 operator-(const vec2 &v1, const vec2 &v2) {
  return vec2(v1.e[0] - v2.e[0],v1.e[1] - v2.e[1]);
}

inline vec2 operator*(const vec2 &v1, const vec2 &v2) {
  return vec2(v1.e[0] * v2.e[0],v1.e[1] * v2.e[1]);
}

inline vec2 operator/(const vec2 &v1, const vec2 &v2) {
  return vec2(v1.e[0] / v2.e[0],v1.e[1] / v2.e[1]);
}

inline vec2 operator*(Float t, const vec2 &v) {
  return vec2(t*v.e[0], t*v.e[1]);
}

inline vec2 operator*(const vec2 &v, Float t) {
  return vec2(t*v.e[0], t*v.e[1]);
}

inline vec2 operator/(const vec2 &v, Float t) {
  return vec2(v.e[0]/t, v.e[1]/t);
}

inline Float dot(const vec2 &v1, const vec2 &v2) {
  return (v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1]);
}

inline vec2& vec2::operator+=(const vec2 &v) {
  e[0] += v.e[0];
  e[1] += v.e[1];
  return(*this);
}

inline vec2& vec2::operator*=(const vec2 &v) {
  e[0] *= v.e[0];
  e[1] *= v.e[1];
  return(*this);
}

inline vec2& vec2::operator/=(const vec2 &v) {
  e[0] /= v.e[0];
  e[1] /= v.e[1];
  return(*this);
}

inline vec2& vec2::operator-=(const vec2 &v) {
  e[0] -= v.e[0];
  e[1] -= v.e[1];
  return(*this);
}

inline vec2& vec2::operator*=(const Float t) {
  e[0] *= t;
  e[1] *= t;
  return(*this);
}

inline vec2& vec2::operator/=(const Float t) {
  Float k = 1.0/t;
  
  e[0] *= k;
  e[1] *= k;
  return(*this);
}

inline vec2 unit_vector(vec2 v) {
  return(v/v.length());
}

#endif