#ifndef VECH
#define VECH

#include <iostream>
#include <cmath>
#include "mathinline.h"

class vec3 {
public:
  vec3() {}
  vec3(Float e0, Float e1, Float e2) {e[0] = e0; e[1] = e1; e[2] = e2;}
  inline Float x() const { return e[0]; }
  inline Float y() const { return e[1]; }
  inline Float z() const { return e[2]; }
  inline Float r() const { return e[0]; }
  inline Float g() const { return e[1]; }
  inline Float b() const { return e[2]; }
  
  inline const vec3& operator+() const { return *this; }
  inline vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); }
  inline Float operator[](int i) const { return e[i]; }
  
  inline vec3& operator+=(const vec3 &v2);
  inline vec3& operator-=(const vec3 &v2);
  inline vec3& operator*=(const vec3 &v2);
  inline vec3& operator/=(const vec3 &v2);
  inline vec3& operator*=(const Float t);
  inline vec3& operator/=(const Float t);
  
  inline Float length() const { return sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]); }
  inline Float squared_length() const { return e[0]*e[0] + e[1]*e[1] + e[2]*e[2]; }
  inline void make_unit_vector();
  
  Float e[3];
};

inline void vec3::make_unit_vector() {
  Float k = 1.0 / sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
  e[0] *= k; e[1] *= k; e[2] *= k; 
}

inline vec3 operator+(const vec3 &v1, const vec3 &v2) {
  return vec3(v1.e[0] + v2.e[0],v1.e[1] + v2.e[1],v1.e[2] + v2.e[2]);
}

inline vec3 operator-(const vec3 &v1, const vec3 &v2) {
  return vec3(v1.e[0] - v2.e[0],v1.e[1] - v2.e[1],v1.e[2] - v2.e[2]);
}

inline vec3 operator*(const vec3 &v1, const vec3 &v2) {
  return vec3(v1.e[0] * v2.e[0],v1.e[1] * v2.e[1],v1.e[2] * v2.e[2]);
}

inline vec3 operator/(const vec3 &v1, const vec3 &v2) {
  return vec3(v1.e[0] / v2.e[0],v1.e[1] / v2.e[1],v1.e[2] / v2.e[2]);
}

inline vec3 operator*(Float t, const vec3 &v) {
  return vec3(t*v.e[0], t*v.e[1], t*v.e[2]);
}

inline vec3 operator*(const vec3 &v, Float t) {
  return vec3(t*v.e[0], t*v.e[1], t*v.e[2]);
}

inline vec3 operator/(const vec3 &v, Float t) {
  return vec3(v.e[0]/t, v.e[1]/t, v.e[2]/t);
}

inline Float dot(const vec3 &v1, const vec3 &v2) {
  return (v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2]);
}

inline vec3 cross(const vec3 &v1, const vec3 &v2) {
  return( vec3( (v1.e[1] * v2.e[2] - v1.e[2] * v2.e[1]),
                (-(v1.e[0] * v2.e[2] - v1.e[2] * v2.e[0])),
                (v1.e[0] * v2.e[1] - v1.e[1] * v2.e[0])));
}

inline vec3& vec3::operator+=(const vec3 &v) {
  e[0] += v.e[0];
  e[1] += v.e[1];
  e[2] += v.e[2];
  return(*this);
}

inline vec3& vec3::operator*=(const vec3 &v) {
  e[0] *= v.e[0];
  e[1] *= v.e[1];
  e[2] *= v.e[2];
  return(*this);
}

inline vec3& vec3::operator/=(const vec3 &v) {
  e[0] /= v.e[0];
  e[1] /= v.e[1];
  e[2] /= v.e[2];
  return(*this);
}

inline vec3& vec3::operator-=(const vec3 &v) {
  e[0] -= v.e[0];
  e[1] -= v.e[1];
  e[2] -= v.e[2];
  return(*this);
}

inline vec3& vec3::operator*=(const Float t) {
  e[0] *= t;
  e[1] *= t;
  e[2] *= t;
  return(*this);
}

inline vec3& vec3::operator/=(const Float t) {
  Float k = 1.0/t;
  
  e[0] *= k;
  e[1] *= k;
  e[2] *= k;
  return(*this);
}

inline vec3 unit_vector(vec3 v) {
  return(v/v.length());
}

#endif
