#ifndef VECH
#define VECH

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

inline float DifferenceOfProducts(float a, float b, float c, float d) {
  float cd = c * d;
  float err = std::fma(-c, d, cd);
  float dop = std::fma(a, b, -cd);
  return(dop + err);
}

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
  inline Float& operator[](int i) { return e[i]; }
  
  inline vec3& operator+=(const vec3 &v2);
  inline vec3& operator-=(const vec3 &v2);
  inline vec3& operator*=(const vec3 &v2);
  inline vec3& operator/=(const vec3 &v2);
  inline vec3& operator*=(const Float t);
  inline vec3& operator/=(const Float t);
  
  bool HasNaNs() const {
    return(std::isnan(e[0]) || std::isnan(e[1]) || std::isnan(e[2]));
  }
  
  inline Float length() const { return sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]); }
  inline Float squared_length() const { return e[0]*e[0] + e[1]*e[1] + e[2]*e[2]; }
  inline vec3 pow(Float exponent) const {
    return(vec3(std::pow(e[0],exponent),std::pow(e[1],exponent),std::pow(e[2],exponent)));
  }
  inline void make_unit_vector();
  
  Float e[3];
};

inline void vec3::make_unit_vector() {
  Float k = 1.0 / sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
  e[0] *= k; e[1] *= k; e[2] *= k; 
}

inline std::istream& operator>>(std::istream &is, vec3 &t) {
  is >> t.e[0] >> t.e[1] >> t.e[2];
  return is;
}

inline std::ostream& operator<<(std::ostream &os, const vec3 &t) {
  os << t.e[0] << " " << t.e[1] << " " << t.e[2];
  return os;
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
  return(vec3(DifferenceOfProducts(v1.y(), v2.z(), v1.z(), v2.y()),
              DifferenceOfProducts(v1.z(), v2.x(), v1.x(), v2.z()),
              DifferenceOfProducts(v1.x(), v2.y(), v1.y(), v2.x())));
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

inline Float MinComponent(const vec3 &v) {
  return(std::min(v.x(), std::min(v.y(), v.z())));
}

inline Float MaxComponent(const vec3 &v) {
  return(std::max(v.x(), std::max(v.y(), v.z())));
}

inline int MaxDimension(const vec3 &v) {
  return((v.x() > v.y()) ? ((v.x() > v.z()) ? 0 : 2) : ((v.y() > v.z()) ? 1 : 2));
}

inline vec3 Min(const vec3 &p1, const vec3 &p2) {
  return(vec3(std::min(p1.x(), p2.x()), std::min(p1.y(), p2.y()), std::min(p1.z(), p2.z())));
}

inline vec3 Max(const vec3 &p1, const vec3 &p2) {
  return(vec3(std::max(p1.x(), p2.x()), std::max(p1.y(), p2.y()),std::max(p1.z(), p2.z())));
}

inline vec3 Permute(const vec3 &v, int x, int y, int z) {
  return(vec3(v.e[x], v.e[y], v.e[z]));
}

inline vec3 Abs(const vec3 &v) {
  return(vec3(v.x(), v.y(), v.z()));
}

#endif
