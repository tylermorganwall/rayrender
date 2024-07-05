#ifndef NORMALH
#define NORMALH

#include <iostream>
#include <cmath>
#include "vec3.h"
#include "point3.h"
#include "float.h"

template <typename T> class normal3 {
public:
  normal3() {}
  normal3(T e0, T e1, T e2) {e[0] = e0; e[1] = e1; e[2] = e2;}
  normal3(T e0) {e[0] = e0; e[1] = e0; e[2] = e0;}
  template <typename U> normal3(const vec3<U> &p) { 
    e[0] = (T)p.x();
    e[1] = (T)p.y();
    e[2] = (T)p.z();
  }
 
  
  
  inline T x() const { return e[0]; }
  inline T y() const { return e[1]; }
  inline T z() const { return e[2]; }
  inline T r() const { return e[0]; }
  inline T g() const { return e[1]; }
  inline T b() const { return e[2]; }
  
  inline const normal3<T>& operator+() const { return *this; }
  
  inline normal3<T> operator-() const { return normal3<T>(-e[0], -e[1], -e[2]); }
  inline T operator[](int i) const { return e[i]; }
  inline T& operator[](int i) { return e[i]; }
  
  inline normal3<T>& operator+=(const normal3<T> &v2);
  inline normal3<T>& operator-=(const normal3<T> &v2);
  inline normal3<T>& operator*=(const normal3<T> &v2);
  inline normal3<T>& operator/=(const normal3<T> &v2);
  inline normal3<T>& operator+=(const Float t);
  inline normal3<T>& operator-=(const Float t);
  inline normal3<T>& operator*=(const Float t);
  inline normal3<T>& operator/=(const Float t);
  
  bool HasNaNs() const {
    return(std::isnan(e[0]) || std::isnan(e[1]) || std::isnan(e[2]));
  }
  
  inline Float length() const { return std::sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]); }
  inline Float squared_length() const { return e[0]*e[0] + e[1]*e[1] + e[2]*e[2]; }
  inline normal3<T> pow(Float exponent) const {
    return(normal3<T>(std::pow(e[0],exponent),std::pow(e[1],exponent),std::pow(e[2],exponent)));
  }
  inline vec3<T> convert_to_vec3() const {
    return(vec3<T>(e[0],e[1],e[2]));
  }
  inline void make_unit_vector();
  
  T e[3];
};

template<typename T> 
inline void normal3<T>::make_unit_vector() {
  Float k = 1.0 / std::sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
  e[0] *= k; e[1] *= k; e[2] *= k; 
}

template<typename T> 
inline std::istream& operator>>(std::istream &is, normal3<T> &t) {
  is >> t.e[0] >> t.e[1] >> t.e[2];
  return is;
}

template<typename T> 
inline std::ostream& operator<<(std::ostream &os, const normal3<T> &t) {
  os << t.e[0] << ", " << t.e[1] << ", " << t.e[2];
  return os;
}

template<typename T> 
inline normal3<T> operator+(const normal3<T> &v1, const normal3<T> &v2) {
  return normal3<T>(v1.e[0] + v2.e[0],v1.e[1] + v2.e[1],v1.e[2] + v2.e[2]);
}

template<typename T> 
inline normal3<T> operator-(const normal3<T> &v1, const normal3<T> &v2) {
  return normal3<T>(v1.e[0] - v2.e[0],v1.e[1] - v2.e[1],v1.e[2] - v2.e[2]);
}

template<typename T> 
inline normal3<T> operator*(const normal3<T> &v1, const normal3<T> &v2) {
  return normal3<T>(v1.e[0] * v2.e[0],v1.e[1] * v2.e[1],v1.e[2] * v2.e[2]);
}

template<typename T> 
inline normal3<T> operator/(const normal3<T> &v1, const normal3<T> &v2) {
  return normal3<T>(v1.e[0] / v2.e[0],v1.e[1] / v2.e[1],v1.e[2] / v2.e[2]);
}

template<typename T> 
inline normal3<T> operator+(const normal3<T> &v, Float t) {
  return normal3<T>(v.e[0] + t,v.e[1] + t,v.e[2] + t);
}

template<typename T> 
inline normal3<T> operator-(const normal3<T> &v, Float t) {
  return normal3<T>(v.e[0] - t,v.e[1] - t,v.e[2] - t);
}

template<typename T> 
inline normal3<T> operator*(Float t, const normal3<T> &v) {
  return normal3<T>(t*v.e[0], t*v.e[1], t*v.e[2]);
}

template<typename T> 
inline normal3<T> operator*(const normal3<T> &v, Float t) {
  return normal3<T>(t*v.e[0], t*v.e[1], t*v.e[2]);
}

template<typename T> 
inline normal3<T> operator/(const normal3<T> &v, Float t) {
  return normal3<T>(v.e[0]/t, v.e[1]/t, v.e[2]/t);
}

template<typename T> 
inline normal3<T> operator/(Float t,const normal3<T> &v) {
  return normal3<T>(t/v.e[0], t/v.e[1], t/v.e[2]);
}

template<typename T> 
inline Float dot(const normal3<T> &v1, const normal3<T> &v2) {
  return (v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2]);
}

template<typename T> 
inline Float dot(const vec3<T> &v1, const normal3<T> &v2) {
  return (v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2]);
}

template<typename T> 
inline Float dot(const normal3<T> &v1, const vec3<T> &v2) {
  return (v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2]);
}

template<typename T> 
inline Float AbsDot(const vec3<T> &v1, const normal3<T> &v2) {
  return std::fabs(v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2]);
}

template<typename T> 
inline Float AbsDot(const normal3<T> &v1, const vec3<T> &v2) {
  return std::fabs(v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2]);
}

template<typename T> 
inline normal3<T> cross(const normal3<T> &v1, const normal3<T> &v2) {
  return(normal3<T>(DifferenceOfProducts(v1.y(), v2.z(), v1.z(), v2.y()),
                 DifferenceOfProducts(v1.z(), v2.x(), v1.x(), v2.z()),
                 DifferenceOfProducts(v1.x(), v2.y(), v1.y(), v2.x())));
}

template<typename T> 
inline normal3<T> cross(const vec3<T> &v1, const normal3<T> &v2) {
  return(normal3<T>(DifferenceOfProducts(v1.y(), v2.z(), v1.z(), v2.y()),
                    DifferenceOfProducts(v1.z(), v2.x(), v1.x(), v2.z()),
                    DifferenceOfProducts(v1.x(), v2.y(), v1.y(), v2.x())));
}

template<typename T> 
inline normal3<T> cross(const normal3<T> &v1, const vec3<T> &v2) {
  return(normal3<T>(DifferenceOfProducts(v1.y(), v2.z(), v1.z(), v2.y()),
                    DifferenceOfProducts(v1.z(), v2.x(), v1.x(), v2.z()),
                    DifferenceOfProducts(v1.x(), v2.y(), v1.y(), v2.x())));
}



template<typename T> 
inline normal3<T>& normal3<T>::operator+=(const normal3<T> &v) {
  e[0] += v.e[0];
  e[1] += v.e[1];
  e[2] += v.e[2];
  return(*this);
}

template<typename T> 
inline normal3<T>& normal3<T>::operator*=(const normal3<T> &v) {
  e[0] *= v.e[0];
  e[1] *= v.e[1];
  e[2] *= v.e[2];
  return(*this);
}

template<typename T> 
inline normal3<T>& normal3<T>::operator/=(const normal3<T> &v) {
  e[0] /= v.e[0];
  e[1] /= v.e[1];
  e[2] /= v.e[2];
  return(*this);
}

template<typename T> 
inline normal3<T>& normal3<T>::operator-=(const normal3<T> &v) {
  e[0] -= v.e[0];
  e[1] -= v.e[1];
  e[2] -= v.e[2];
  return(*this);
}

template<typename T> 
inline normal3<T>& normal3<T>::operator+=(const Float t) {
  e[0] += t;
  e[1] += t;
  e[2] += t;
  return(*this);
}

template<typename T> 
inline normal3<T>& normal3<T>::operator*=(const Float t) {
  e[0] *= t;
  e[1] *= t;
  e[2] *= t;
  return(*this);
}

template<typename T> 
inline normal3<T>& normal3<T>::operator/=(const Float t) {
  Float k = 1.0/t;
  
  e[0] *= k;
  e[1] *= k;
  e[2] *= k;
  return(*this);
}

template<typename T> 
inline normal3<T> unit_vector(normal3<T> v) {
  return(v/v.length());
}

template<typename T> 
inline Float MinComponent(const normal3<T> &v) {
  return(std::min(v.x(), std::min(v.y(), v.z())));
}

template<typename T> 
inline Float MaxComponent(const normal3<T> &v) {
  return(std::fmax(v.x(), std::fmax(v.y(), v.z())));
}

template<typename T> 
inline int MaxDimension(const normal3<T> &v) {
  return((v.x() > v.y()) ? ((v.x() > v.z()) ? 0 : 2) : ((v.y() > v.z()) ? 1 : 2));
}

template<typename T> 
inline normal3<T> Min(const normal3<T> &p1, const normal3<T> &p2) {
  return(normal3<T>(fmin(p1.x(), p2.x()), fmin(p1.y(), p2.y()),fmin(p1.z(), p2.z())));
}

template<typename T> 
inline normal3<T> Max(const normal3<T> &p1, const normal3<T> &p2) {
  return(normal3<T>(fmax(p1.x(), p2.x()), fmax(p1.y(), p2.y()),fmax(p1.z(), p2.z())));
}

template<typename T> 
inline normal3<T> Permute(const normal3<T> &v, int x, int y, int z) {
  return(normal3<T>(v.e[x], v.e[y], v.e[z]));
}

template<typename T> 
inline normal3<T> Abs(const normal3<T> &v) {
  return(normal3<T>(fabs(v.x()), fabs(v.y()), fabs(v.z())));
}

template <typename T> inline normal3<T>
Faceforward(const normal3<T> &n, const vec3<T> &v) {
  return (dot(n, v) < 0.f) ? -n : n;
}


typedef normal3<Float> normal3f;
typedef normal3<int>   normal3i;


#endif
