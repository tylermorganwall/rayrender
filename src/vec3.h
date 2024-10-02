#ifndef VECH
#define VECH

#include <iostream>
#include <cmath>
#include "float.h"

template<class T>
inline T ffmin(T a, T b);

template<class T>
inline T ffmax(T a, T b);


inline Float DifferenceOfProducts(Float a, Float b, Float c, Float d) {
  Float cd = c * d;
#ifndef RAY_FLOAT_AS_DOUBLE
  Float err = std::fmaf(-c, d, cd);
  Float dop = std::fmaf(a, b, -cd);
#else
  Float err = std::fma(-c, d, cd);
  Float dop = std::fma(a, b, -cd);
#endif
  return(dop + err);
}

// inline Float DifferenceOfProducts(Float a, Float b, Float c, Float d) {
//   return(a * b - c * d);
// }

template <typename T> class vec3 {
public:
  vec3() {}
  vec3(T e0, T e1, T e2) {e[0] = e0; e[1] = e1; e[2] = e2;}
  vec3(T e0) {e[0] = e0; e[1] = e0; e[2] = e0;}
  inline T x() const { return e[0]; }
  inline T y() const { return e[1]; }
  inline T z() const { return e[2]; }
  inline T r() const { return e[0]; }
  inline T g() const { return e[1]; }
  inline T b() const { return e[2]; }
  
  inline const vec3<T>& operator+() const { return *this; }
  inline vec3<T> operator-() const { return vec3<T>(-e[0], -e[1], -e[2]); }
  inline T operator[](int i) const { return e[i]; }
  inline T& operator[](int i) { return e[i]; }
  
  inline vec3<T>& operator+=(const vec3<T> &v2);
  inline vec3<T>& operator-=(const vec3<T> &v2);
  inline vec3<T>& operator*=(const vec3<T> &v2);
  inline vec3<T>& operator/=(const vec3<T> &v2);
  
  inline vec3<T>& operator+=(const Float& t);
  inline vec3<T>& operator-=(const Float& t);
  inline vec3<T>& operator*=(const Float& t);
  inline vec3<T>& operator/=(const Float& t);
  
  bool HasNaNs() const {
    return(std::isnan(e[0]) || std::isnan(e[1]) || std::isnan(e[2]));
  }
  
  inline Float length() const { return std::sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]); }
  inline Float squared_length() const { return e[0]*e[0] + e[1]*e[1] + e[2]*e[2]; }
  inline vec3<T> pow(Float exponent) const {
    return(vec3<T>(std::pow(e[0],exponent),std::pow(e[1],exponent),std::pow(e[2],exponent)));
  }
  inline void make_unit_vector();
  
  T e[3];
};

template<typename T> 
inline void vec3<T>::make_unit_vector() {
  Float k = 1.0 / std::sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
  e[0] *= k; e[1] *= k; e[2] *= k; 
}

template<typename T> 
inline std::istream& operator>>(std::istream &is, vec3<T> &t) {
  is >> t.e[0] >> t.e[1] >> t.e[2];
  return is;
}

template<typename T> 
inline std::ostream& operator<<(std::ostream &os, const vec3<T> &t) {
  os << t.e[0] << ", " << t.e[1] << ", " << t.e[2];
  return os;
}

template<typename T> 
inline vec3<T> operator+(const vec3<T> &v1, const vec3<T> &v2) {
  return vec3<T>(v1.e[0] + v2.e[0],v1.e[1] + v2.e[1],v1.e[2] + v2.e[2]);
}

template<typename T> 
inline vec3<T> operator-(const vec3<T> &v1, const vec3<T> &v2) {
  return vec3<T>(v1.e[0] - v2.e[0],v1.e[1] - v2.e[1],v1.e[2] - v2.e[2]);
}

template<typename T> 
inline vec3<T> operator*(const vec3<T> &v1, const vec3<T> &v2) {
  return vec3<T>(v1.e[0] * v2.e[0],v1.e[1] * v2.e[1],v1.e[2] * v2.e[2]);
}

template<typename T> 
inline vec3<T> operator/(const vec3<T> &v1, const vec3<T> &v2) {
  return vec3<T>(v1.e[0] / v2.e[0],v1.e[1] / v2.e[1],v1.e[2] / v2.e[2]);
}

template<typename T> 
inline vec3<T> operator+(const vec3<T> &v, Float t) {
  return vec3<T>(v.e[0] + t,v.e[1] + t,v.e[2] + t);
}

template<typename T> 
inline vec3<T> operator-(const vec3<T> &v, Float t) {
  return vec3<T>(v.e[0] - t,v.e[1] - t,v.e[2] - t);
}

template<typename T> 
inline vec3<T> operator*(Float t, const vec3<T> &v) {
  return vec3<T>(t*v.e[0], t*v.e[1], t*v.e[2]);
}

template<typename T> 
inline vec3<T> operator*(const vec3<T> &v, Float t) {
  return vec3<T>(t*v.e[0], t*v.e[1], t*v.e[2]);
}

template<typename T> 
inline vec3<T> operator/(const vec3<T> &v, Float t) {
  Float inv_t = 1./t;
  return vec3<T>(v.e[0]*inv_t, v.e[1]*inv_t, v.e[2]*inv_t);
}

template<typename T> 
inline vec3<T> operator/(Float t,const vec3<T> &v) {
  return vec3<T>(t/v.e[0], t/v.e[1], t/v.e[2]);
}

template<typename T> 
inline Float dot(const vec3<T> &v1, const vec3<T> &v2) {
  return (v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2]);
}

template<typename T> 
inline vec3<T> cross(const vec3<T> &v1, const vec3<T> &v2) {
  return(vec3<T>(DifferenceOfProducts(v1.y(), v2.z(), v1.z(), v2.y()),
                 DifferenceOfProducts(v1.z(), v2.x(), v1.x(), v2.z()),
                 DifferenceOfProducts(v1.x(), v2.y(), v1.y(), v2.x())));
}

template<typename T> 
inline bool parallelVectors(const vec3<T> &v1, const vec3<T> &v2) {
  Float x = DifferenceOfProducts(v1.y(), v2.z(), v1.z(), v2.y());
  Float y = DifferenceOfProducts(v1.z(), v2.x(), v1.x(), v2.z());
  Float z = DifferenceOfProducts(v1.x(), v2.y(), v1.y(), v2.x());
  if(x == 0 && y == 0 && z == 0) {
    [[unlikely]];
    return(true);
  }
  return(false);
}

template<typename T> 
inline vec3<T>& vec3<T>::operator+=(const vec3<T> &v) {
  e[0] += v.e[0];
  e[1] += v.e[1];
  e[2] += v.e[2];
  return(*this);
}

template<typename T> 
inline vec3<T>& vec3<T>::operator-=(const vec3<T> &v) {
  e[0] -= v.e[0];
  e[1] -= v.e[1];
  e[2] -= v.e[2];
  return(*this);
}

template<typename T> 
inline vec3<T>& vec3<T>::operator*=(const vec3<T> &v) {
  e[0] *= v.e[0];
  e[1] *= v.e[1];
  e[2] *= v.e[2];
  return(*this);
}

template<typename T> 
inline vec3<T>& vec3<T>::operator/=(const vec3<T> &v) {
  e[0] /= v.e[0];
  e[1] /= v.e[1];
  e[2] /= v.e[2];
  return(*this);
}

template<typename T> 
inline vec3<T>& vec3<T>::operator-=(const Float& t) {
  e[0] -= t;
  e[1] -= t;
  e[2] -= t;
  return(*this);
}

template<typename T> 
inline vec3<T>& vec3<T>::operator+=(const Float& t) {
  e[0] += t;
  e[1] += t;
  e[2] += t;
  return(*this);
}

template<typename T> 
inline vec3<T>& vec3<T>::operator*=(const Float& t) {
  e[0] *= t;
  e[1] *= t;
  e[2] *= t;
  return(*this);
}

template<typename T> 
inline vec3<T>& vec3<T>::operator/=(const Float& t) {
  Float k = 1.0/t;
  
  e[0] *= k;
  e[1] *= k;
  e[2] *= k;
  return(*this);
}

template<typename T> 
inline vec3<T> unit_vector(vec3<T> v) {
  return(v/v.length());
}

template<typename T> 
inline Float MinComponent(const vec3<T> &v) {
  return(ffmin(v.x(), ffmin(v.y(), v.z())));
}

template<typename T> 
inline Float MaxComponent(const vec3<T> &v) {
  return(ffmax(v.x(), ffmax(v.y(), v.z())));
}

template<typename T> 
inline int MaxDimension(const vec3<T> &v) {
  return((v.x() > v.y()) ? ((v.x() > v.z()) ? 0 : 2) : ((v.y() > v.z()) ? 1 : 2));
}

template<typename T> 
inline vec3<T> Min(const vec3<T> &p1, const vec3<T> &p2) {
  return(vec3<T>(ffmin(p1.x(), p2.x()), ffmin(p1.y(), p2.y()),ffmin(p1.z(), p2.z())));
}

template<typename T> 
inline vec3<T> Max(const vec3<T> &p1, const vec3<T> &p2) {
  return(vec3<T>(ffmax(p1.x(), p2.x()), ffmax(p1.y(), p2.y()),ffmax(p1.z(), p2.z())));
}

template<typename T> 
inline vec3<T> Permute(const vec3<T> &v, int x, int y, int z) {
  return(vec3<T>(v.e[x], v.e[y], v.e[z]));
}

// In-place Permute function
template<typename T>
inline void PermuteInPlace(vec3<T>& v, int x, int y, int z) {
  T temp[3] = { v.e[x], v.e[y], v.e[z] };
  v.e[0] = temp[0]; v.e[1] = temp[1]; v.e[2] = temp[2];
}

template<typename T> 
inline vec3<T> Abs(const vec3<T> &v) {
  return(vec3<T>(std::fabs(v.x()), std::fabs(v.y()), std::fabs(v.z())));
}




typedef vec3<Float> vec3f;
typedef vec3<int>   vec3i;


#endif
