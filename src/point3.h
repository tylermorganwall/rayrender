#ifndef POINT3H
#define POINT3H

#include <iostream>
#include <cmath>

#include "float.h"
#include "vec3.h"

template <typename T> class alignas(16) point3 {
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
  
  inline point3<T>& operator+=(const T& t);
  inline point3<T>& operator-=(const T& t);
  inline point3<T>& operator*=(const T& t);
  inline point3<T>& operator/=(const T& t);
  
  inline bool HasNaNs() const {
    if constexpr (std::is_floating_point<T>::value) {
        return (std::isnan(e[0]) || std::isnan(e[1]) || std::isnan(e[2]));
    } else {
        return false;
    }
  }
  
  inline T length() const { return std::sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]); }
  inline T squared_length() const { return e[0]*e[0] + e[1]*e[1] + e[2]*e[2]; }
  template<typename U = T>
  inline typename std::enable_if<std::is_floating_point<U>::value, point3<T>>::type
  pow(T exponent) const {
      return point3<T>(std::pow(e[0], exponent), std::pow(e[1], exponent), std::pow(e[2], exponent));
  }
  
  inline void make_unit_vector();
  
  T e[3];
};

template<typename T> 
inline void point3<T>::make_unit_vector() {
  // Float len = length();
  //   if(len - 1 < 1e-8) {
  //     volatile int f = 1;
  //   }
  T k = 1.0 / std::sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
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
inline point3<T> operator+(const point3<T> &v, T t) {
  return point3<T>(v.e[0] + t,v.e[1] + t,v.e[2] + t);
}

template<typename T> 
inline point3<T> operator-(const point3<T> &v, T t) {
  return point3<T>(v.e[0] - t,v.e[1] - t,v.e[2] - t);
}

template<typename T> 
inline point3<T> operator*(T t, const point3<T> &v) {
  return point3<T>(t*v.e[0], t*v.e[1], t*v.e[2]);
}

template<typename T> 
inline point3<T> operator*(const point3<T> &v, T t) {
  return point3<T>(t*v.e[0], t*v.e[1], t*v.e[2]);
}

template<typename T>
inline point3<T> operator*(int t, const point3<T>& v) {
    T t_f = static_cast<T>(t);
    return point3<T>(t_f*v.e[0], t_f*v.e[1], t_f*v.e[2]);
}

template<typename T>
inline point3<T> operator*(const point3<T>& v, int t) {
    T t_f = static_cast<T>(t);
    return point3<T>(t_f*v.e[0], t_f*v.e[1], t_f*v.e[2]);
}

template<typename T> 
inline point3<T> operator/(const point3<T> &v, T t) {
    T inv_t = 1.f / t;
  return point3<T>(v.e[0] * inv_t, v.e[1] * inv_t, v.e[2] * inv_t);
}

template<typename T> 
inline point3<T> operator/(T t,const point3<T> &v) {
  return point3<T>(t/v.e[0], t/v.e[1], t/v.e[2]);
}

template<typename T> 
inline T dot(const point3<T> &v1, const point3<T> &v2) {
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
inline point3<T>& point3<T>::operator-=(const T& t) {
  e[0] -= t;
  e[1] -= t;
  e[2] -= t;
  return(*this);
}

template<typename T> 
inline point3<T>& point3<T>::operator+=(const T& t) {
  e[0] += t;
  e[1] += t;
  e[2] += t;
  return(*this);
}

template<typename T> 
inline point3<T>& point3<T>::operator*=(const T& t) {
  e[0] *= t;
  e[1] *= t;
  e[2] *= t;
  return(*this);
}

template<typename T> 
inline point3<T>& point3<T>::operator/=(const T& t) {
  T k = 1.0/t;
  
  e[0] *= k;
  e[1] *= k;
  e[2] *= k;
  return(*this);
}

template<typename T> 
inline point3<T> unit_vector(point3<T> v) {
    // T len = v.length();
    // if(len - 1 < 1e-8) {
    //   volatile int f = 1;
    // }
  return(v/v.length());
}

template<typename T> 
inline T MinComponent(const point3<T> &v) {
  return(ffmin(v.x(), ffmin(v.y(), v.z())));
}

template<typename T> 
inline T MaxComponent(const point3<T> &v) {
  return(ffmax(v.x(), ffmax(v.y(), v.z())));
}

template<typename T> 
inline int MaxDimension(const point3<T> &v) {
  return((v.x() > v.y()) ? ((v.x() > v.z()) ? 0 : 2) : ((v.y() > v.z()) ? 1 : 2));
}


template<typename T> 
inline point3<T> Min(const point3<T> &p1, const point3<T> &p2) {
  return(point3<T>(ffmin(p1.x(), p2.x()), ffmin(p1.y(), p2.y()),ffmin(p1.z(), p2.z())));
}

template<typename T> 
inline point3<T> Max(const point3<T> &p1, const point3<T> &p2) {
  return(point3<T>(ffmax(p1.x(), p2.x()), ffmax(p1.y(), p2.y()),ffmax(p1.z(), p2.z())));
}

template<typename T> 
inline point3<T> Permute(const point3<T> &v, int x, int y, int z) {
  return(point3<T>(v.e[x], v.e[y], v.e[z]));
}

// In-place Permute function
template<typename T>
inline void PermuteInPlace(point3<T>& v, int x, int y, int z) {
  T temp[3] = { v.e[x], v.e[y], v.e[z] };
  v.e[0] = temp[0]; v.e[1] = temp[1]; v.e[2] = temp[2];
}

template<typename T> 
inline point3<T> Abs(const point3<T> &v) {
  return(point3<T>(std::abs(v.x()), std::abs(v.y()), std::abs(v.z())));
}

template<typename T>
inline point3<T> operator/(const point3<T>& v, int t) {
    T k = 1.0 / static_cast<Float>(t);
    return point3<T>(v.e[0]*k, v.e[1]*k, v.e[2]*k);
}

template <typename T> inline T
Distance(const point3<T> &p1, const point3<T> &p2) {
  return (p1 - p2).length();
}

template <typename T> inline T
DistanceSquared(const point3<T> &p1, const point3<T> &p2) {
  return (p1 - p2).squared_length();
}

template <typename T> point3<T>
Lerp(T t, const point3<T> &p0, const point3<T> &p1) {
  return (1 - t) * p0 + t * p1;
}

template <typename T> point3<T> Floor(const point3<T> &p) {
  return point3<T>(std::floor(p.x()), std::floor(p.y()), std::floor(p.z()));
}

template <typename T> point3<T> Ceil(const point3<T> &p) {
  return point3<T>(std::ceil(p.x()), std::ceil(p.y()), std::ceil(p.z()));
}


template<typename T> 
inline T dot(const vec3<T> &v1, const point3<T> &v2) {
  return (v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2]);
}

template<typename T> 
inline T dot(const point3<T> &v1, const vec3<T> &v2) {
  return (v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2]);
}

#ifdef RAYSIMDVEC

// Specialize for Float
template<>
class alignas(16) point3<Float> {
public:
    FVec4 e;

    // Constructors
    point3() {
        e = simd_set1(0.0f);
    }

    point3(Float e0, Float e1, Float e2) {
        alignas(16) float values[4] = { e0, e1, e2, 0.0f };
        e = simd_load(values);
    }

    point3(Float e0) {
        e = simd_set1(e0);
    }

    // Accessors
    inline Float x() const { return e[0]; }
    inline Float y() const { return e[1]; }
    inline Float z() const { return e[2]; }
    inline Float r() const { return e[0]; }
    inline Float g() const { return e[1]; }
    inline Float b() const { return e[2]; }

    // Operators with vec
    inline point3<Float> operator+(const vec3<Float>& v) const;
    inline point3<Float>& operator+=(const vec3<Float>& v);
    inline point3<Float> operator-(const vec3<Float>& v) const;
    inline point3<Float>& operator-=(const vec3<Float>& v);
    inline vec3<Float> operator-(const point3<Float>& p) const;
    inline point3<Float> operator*(const point3<Float>& p) const;
    inline point3<Float> operator/(const point3<Float>& p) const;

    // Operators 
    inline const point3<Float>& operator+() const { return *this; }
    inline point3<Float> operator-() const {
        point3<Float> result;
        result.e = simd_sub(simd_set1(0.0f), e);
        return result;
    }
    inline Float operator[](int i) const { return e[i]; }
    inline Float& operator[](int i) { return e[i]; }

    inline point3<Float>& operator+=(const point3<Float>& v2) {
        e = simd_add(e, v2.e);
        return *this;
    }

    inline point3<Float>& operator-=(const point3<Float>& v2) {
        e = simd_sub(e, v2.e);
        return *this;
    }

    inline point3<Float>& operator*=(const point3<Float>& v2) {
        e = simd_mul(e, v2.e);
        return *this;
    }

    inline point3<Float>& operator/=(const point3<Float>& v2) {
        e = simd_div(e, v2.e);
        return *this;
    }

    inline point3<Float>& operator+=(const Float& t) {
        e = simd_add(e, simd_set1(t));
        return *this;
    }

    inline point3<Float>& operator-=(const Float& t) {
        e = simd_sub(e, simd_set1(t));
        return *this;
    }

    inline point3<Float>& operator*=(const Float& t) {
        e = simd_mul(e, simd_set1(t));
        return *this;
    }

    inline point3<Float>& operator/=(const Float& t) {
        e = simd_div(e, simd_set1(t));
        return *this;
    }

    inline bool HasNaNs() const {
        return (std::isnan(e[0]) || std::isnan(e[1]) || std::isnan(e[2]));
    }

    inline Float length() const {
        return std::sqrt(squared_length());
    }

    inline Float squared_length() const {
      return(simd_squared_length(e));
        // FVec4 mul = simd_mul(e, e);
        // return mul[0] + mul[1] + mul[2];
    }

    inline point3<Float> pow(Float exponent) const {
        point3<Float> result;
        result.e[0] = std::pow(e[0], exponent);
        result.e[1] = std::pow(e[1], exponent);
        result.e[2] = std::pow(e[2], exponent);
        result.e[3] = 0.0f;
        return result;
    }

    inline point3<Float> sgn() const {
        point3<Float> result;
        result.e = simd_sgn(e);
        return result;
    }

    inline void make_unit_vector() {
          // Float len = length();
          // if(len - 1 < 1e-8) {
          //   volatile int f = 1;
          // }
        Float len = length();
        e = simd_div(e, simd_set1(len));
    }
};

// Operators outside the class
inline point3<Float> operator+(const point3<Float>& v1, const point3<Float>& v2) {
    point3<Float> result;
    result.e = simd_add(v1.e, v2.e);
    return result;
}

inline point3<Float> point3<Float>::operator+(const vec3<Float>& v) const {
    point3<Float> result;
    result.e = simd_add(e, v.e);
    return result;
}

inline point3<Float> operator+(const vec3<Float>& v1, const point3<Float>& v2) {
    point3<Float> result;
    result.e = simd_add(v1.e, v2.e);
    return result;
}

inline point3<Float>& point3<Float>::operator+=(const vec3<Float>& v) {
    e = simd_add(e, v.e);
    return *this;
}

inline point3<Float>& point3<Float>::operator-=(const vec3<Float>& v) {
    e = simd_sub(e, v.e);
    return *this;
}


inline vec3<Float> point3<Float>::operator-(const point3<Float>& p) const {
    vec3<Float> result;
    result.e = simd_sub(e, p.e);
    return result;
}

inline point3<Float> point3<Float>::operator-(const vec3<Float>& v) const {
    point3<Float> result;
    result.e = simd_sub(e, v.e);
    return result;
}

inline point3<Float> point3<Float>::operator*(const point3<Float>& p) const {
    point3<Float> result;
    result.e = simd_mul(e, p.e);
    return result;
}

inline point3<Float> point3<Float>::operator/(const point3<Float>& p) const {
    point3<Float> result;
    result.e = simd_div(e, p.e);
    return result;
}

inline point3<Float> operator+(const point3<Float>& v, Float t) {
    point3<Float> result;
    result.e = simd_add(v.e, simd_set1(t));
    return result;
}

inline point3<Float> operator-(const point3<Float>& v, Float t) {
    point3<Float> result;
    result.e = simd_sub(v.e, simd_set1(t));
    return result;
}


inline point3<Float> operator*(Float t, const point3<Float>& v) {
    point3<Float> result;
    result.e = simd_mul(simd_set1(t), v.e);
    return result;
}

inline point3<Float> operator*(const point3<Float>& v, Float t) {
    point3<Float> result;
    result.e = simd_mul(v.e, simd_set1(t));
    return result;
}

inline point3<Float> operator/(const point3<Float>& v, Float t) {
    point3<Float> result;
    result.e = simd_div(v.e, simd_set1(t));
    return result;
}

inline point3<Float> operator/(Float t, const point3<Float>& v) {
    point3<Float> result;
    result.e = simd_div(simd_set1(t), v.e);
    return result;
}

inline point3<Float> cross(const point3<Float>& v1, const point3<Float>& v2) {
  point3<Float> result;
  result.e = simd_cross(v1.e, v2.e);
  return result;
}

inline Float Distance(const point3<Float>& p1, const point3<Float>& p2) {
    return (p1 - p2).length();
}

inline Float DistanceSquared(const point3<Float>& p1, const point3<Float>& p2) {
    return (p1 - p2).squared_length();
}

inline point3<Float> Lerp(Float t, const point3<Float>& p0, const point3<Float>& p1) {
    point3<Float> result;
    FVec4 temp1 = simd_mul(simd_set1(1 - t), p0.e);
    FVec4 temp2 = simd_mul(simd_set1(t), p1.e);
    result.e = simd_add(temp1, temp2);
    return result;
}

// inline point3<Float> cross(const point3<Float>& v1, const point3<Float>& v2) {
//     point3<Float> result;
// #ifdef HAS_SSE
//     result.e.v = _mm_sub_ps(
//         _mm_mul_ps(_mm_shuffle_ps(v1.e.v, v1.e.v, _MM_SHUFFLE(3, 0, 2, 1)), _mm_shuffle_ps(v2.e.v, v2.e.v, _MM_SHUFFLE(3, 1, 0, 2))),
//         _mm_mul_ps(_mm_shuffle_ps(v1.e.v, v1.e.v, _MM_SHUFFLE(3, 1, 0, 2)), _mm_shuffle_ps(v2.e.v, v2.e.v, _MM_SHUFFLE(3, 0, 2, 1)))
//     );
// #else
//     result.e[0] = DifferenceOfProducts(v1.y(), v2.z(), v1.z(), v2.y());
//     result.e[1] = DifferenceOfProducts(v1.z(), v2.x(), v1.x(), v2.z());
//     result.e[2] = DifferenceOfProducts(v1.x(), v2.y(), v1.y(), v2.x());
//     result.e[3] = 0.0f;
// #endif
//     return result;
// }

inline point3<Float> unit_vector(const point3<Float>& v) {
    // Float len = v.length();
    // if(len - 1 < 1e-8) {
    //   volatile int f = 1;
    // }
    return v / v.length();
}

inline Float MinComponent(const point3<Float>& v) {
    return std::fmin(v.x(), std::fmin(v.y(), v.z()));
}

inline Float MaxComponent(const point3<Float>& v) {
    return std::fmax(v.x(), std::fmax(v.y(), v.z()));
}

inline int MaxDimension(const point3<Float>& v) {
    if (v.x() > v.y()) {
        return (v.x() > v.z()) ? 0 : 2;
    } else {
        return (v.y() > v.z()) ? 1 : 2;
    }
}

inline point3<Float> Min(const point3<Float>& p1, const point3<Float>& p2) {
    point3<Float> result;
    result.e = simd_min(p1.e, p2.e);
    return result;
}

inline point3<Float> Max(const point3<Float>& p1, const point3<Float>& p2) {
    point3<Float> result;
    result.e = simd_max(p1.e, p2.e);
    return result;
}

inline point3<Float> Permute(const point3<Float>& v, int x, int y, int z) {
    point3<Float> result;
    result.e[0] = v.e[x];
    result.e[1] = v.e[y];
    result.e[2] = v.e[z];
    result.e[3] = 0.0f;
    return result;
}

inline void PermuteInPlace(point3<Float>& v, int x, int y, int z) {
    Float temp[3] = { v.e[x], v.e[y], v.e[z] };
    v.e[0] = temp[0];
    v.e[1] = temp[1];
    v.e[2] = temp[2];
    v.e[3] = 0.0f;
}

inline point3<Float> Abs(const point3<Float>& v) {
    point3<Float> result;
#ifdef HAS_SSE
    // SSE implementation using _mm_andnot_ps to clear the sign bit
    result.e.v = _mm_andnot_ps(_mm_set1_ps(-0.0f), v.e.v);
#elif defined(HAS_NEON)
    // NEON implementation using vabsq_f32 to compute absolute values
    result.e.v = vabsq_f32(v.e.v);
#else
    // Fallback to scalar implementation if SIMD is not available
#ifdef RAY_FLOAT_AS_DOUBLE
    result.e[0] = ffabs(v.e[0]);
    result.e[1] = ffabs(v.e[1]);
    result.e[2] = ffabs(v.e[2]);
    result.e[3] = 0.0;
#else 
    result.e[0] = ffabs(v.e[0]);
    result.e[1] = ffabs(v.e[1]);
    result.e[2] = ffabs(v.e[2]);
    result.e[3] = 0.0;
#endif
#endif
    return result;
}

inline Float dot(const point3<Float>& v1, const point3<Float>& v2) {
    return simd_dot(v1.e, v2.e);
}

inline Float dot(const vec3<Float>& v1, const point3<Float>& p2) {
    return simd_dot(v1.e, p2.e);
}

inline Float dot(const point3<Float>& p1, const vec3<Float>& v2) {
    return simd_dot(p1.e, v2.e);
}


// Specialize for int type
template<>
class alignas(16) point3<int> {
public:
    IVec4 e;

    // Constructors
    point3() {
        e = simd_set1(0);
    }

    point3(int e0, int e1, int e2) {
        alignas(16) int values[4] = { e0, e1, e2, 0 };
    #ifdef HAS_SSE
        e.v = _mm_loadu_si128(reinterpret_cast<const __m128i*>(values));
    #elif defined(HAS_NEON)
        e.v = vld1q_s32(values);
    #else
        for (int i = 0; i < 4; ++i) {
            e.xyzw[i] = values[i];
        }
    #endif
    }

    point3(int e0) {
        e = simd_set1(e0);
    }

    // Accessors
    inline int x() const { return e[0]; }
    inline int y() const { return e[1]; }
    inline int z() const { return e[2]; }

    inline int& x() { return e[0]; }
    inline int& y() { return e[1]; }
    inline int& z() { return e[2]; }

    // Operators
    inline const point3<int>& operator+() const { return *this; }
    inline point3<int> operator-() const {
        point3<int> result;
        result.e = simd_sub(simd_set1(0), e);
        return result;
    }
    inline int operator[](int i) const { return e[i]; }
    inline int& operator[](int i) { return e[i]; }

    inline point3<int>& operator+=(const point3<int>& v2) {
        e = simd_add(e, v2.e);
        return *this;
    }

    inline point3<int>& operator-=(const point3<int>& v2) {
        e = simd_sub(e, v2.e);
        return *this;
    }

    inline point3<int>& operator*=(const point3<int>& v2) {
        e = simd_mul(e, v2.e);
        return *this;
    }

    inline point3<int>& operator/=(const point3<int>& v2) {
        e = simd_div(e, v2.e);
        return *this;
    }

    inline point3<int>& operator+=(const int& t) {
        e = simd_add(e, simd_set1(t));
        return *this;
    }

    inline point3<int>& operator-=(const int& t) {
        e = simd_sub(e, simd_set1(t));
        return *this;
    }

    inline point3<int>& operator*=(const int& t) {
        e = simd_mul(e, simd_set1(t));
        return *this;
    }

    inline point3<int>& operator/=(const int& t) {
        e = simd_div(e, simd_set1(t));
        return *this;
    }

    inline bool HasNaNs() const {
        // Integers do not have NaNs, but we can check for some invalid value if needed
        return false;
    }

    inline point3<int> sgn() const {
        point3<int> result;
        result.e = simd_sgn(e);
        return result;
    }
};


// Operators outside the class
inline point3<int> operator+(const point3<int>& v1, const point3<int>& v2) {
    point3<int> result;
    result.e = simd_add(v1.e, v2.e);
    return result;
}

inline point3<int> operator-(const point3<int>& v1, const point3<int>& v2) {
    point3<int> result;
    result.e = simd_sub(v1.e, v2.e);
    return result;
}

inline point3<int> operator*(const point3<int>& v1, const point3<int>& v2) {
    point3<int> result;
    result.e = simd_mul(v1.e, v2.e);
    return result;
}

inline point3<int> operator/(const point3<int>& v1, const point3<int>& v2) {
    point3<int> result;
    result.e = simd_div(v1.e, v2.e);
    return result;
}

inline point3<int> operator+(const point3<int>& v, int t) {
    point3<int> result;
    result.e = simd_add(v.e, simd_set1(t));
    return result;
}

inline point3<int> operator-(const point3<int>& v, int t) {
    point3<int> result;
    result.e = simd_sub(v.e, simd_set1(t));
    return result;
}

inline point3<int> operator*(int t, const point3<int>& v) {
    point3<int> result;
    result.e = simd_mul(simd_set1(t), v.e);
    return result;
}

inline point3<int> operator*(const point3<int>& v, int t) {
    point3<int> result;
    result.e = simd_mul(v.e, simd_set1(t));
    return result;
}

inline point3<int> operator/(const point3<int>& v, int t) {
    point3<int> result;
    result.e = simd_div(v.e, simd_set1(t));
    return result;
}

inline point3<int> operator/(int t, const point3<int>& v) {
    point3<int> result;
    result.e = simd_div(simd_set1(t), v.e);
    return result;
}

inline int dot(const point3<int>& v1, const point3<int>& v2) {
    IVec4 mul = simd_mul(v1.e, v2.e);
    int sum = mul[0] + mul[1] + mul[2];
    return sum;
}

inline point3<int> cross(const point3<int>& v1, const point3<int>& v2) {
    point3<int> result;
    int x = v1.y() * v2.z() - v1.z() * v2.y();
    int y = v1.z() * v2.x() - v1.x() * v2.z();
    int z = v1.x() * v2.y() - v1.y() * v2.x();
    result = point3<int>(x, y, z);
    return result;
}

// inline point3<int>& point3<int>::operator*=(const int& t) {
//     e = simd_mul(e, simd_set1(t));
//     return *this;
// }

// inline point3<int>& point3<int>::operator/=(const int& t) {
//     e = simd_div(e, simd_set1(t));
//     return *this;
// }

inline point3<int> Min(const point3<int>& p1, const point3<int>& p2) {
    point3<int> result;
#if defined(__SSE4_1__)
    result.e.v = _mm_min_epi32(p1.e.v, p2.e.v);
#elif defined(HAS_NEON)
    result.e.v = vminq_s32(p1.e.v, p2.e.v);
#else
    for (int i = 0; i < 3; ++i) {
        result.e[i] = std::min(p1.e[i], p2.e[i]);
    }
#endif
    return result;
}

inline point3<int> Max(const point3<int>& p1, const point3<int>& p2) {
    point3<int> result;
#if defined(__SSE4_1__)
    result.e.v = _mm_max_epi32(p1.e.v, p2.e.v);
#elif defined(HAS_NEON)
    result.e.v = vmaxq_s32(p1.e.v, p2.e.v);
#else
    for (int i = 0; i < 3; ++i) {
        result.e[i] = std::max(p1.e[i], p2.e[i]);
    }
#endif
    return result;
}

inline int MinComponent(const point3<int>& v) {
    int minVal = std::min(v.x(), std::min(v.y(), v.z()));
    return minVal;
}

inline int MaxComponent(const point3<int>& v) {
    int maxVal = std::max(v.x(), std::max(v.y(), v.z()));
    return maxVal;
}

inline int MaxDimension(const point3<int>& v) {
    if (v.x() > v.y()) {
        return (v.x() > v.z()) ? 0 : 2;
    } else {
        return (v.y() > v.z()) ? 1 : 2;
    }
}

inline point3<int> Abs(const point3<int>& v) {
    point3<int> result;
#ifdef __SSE3__
    result.e.v = _mm_abs_epi32(v.e.v);
#elif defined(HAS_NEON)
    result.e.v = vabsq_s32(v.e.v);
#else
    for (int i = 0; i < 3; ++i) {
        result.e[i] = std::abs(v.e[i]);
    }
#endif
    return result;
}

inline point3<int> Permute(const point3<int>& v, int x, int y, int z) {
    point3<int> result;
    result.e[0] = v.e[x];
    result.e[1] = v.e[y];
    result.e[2] = v.e[z];
    result.e[3] = 0;
    return result;
}

inline void PermuteInPlace(point3<int>& v, int x, int y, int z) {
    int temp[3] = { v.e[x], v.e[y], v.e[z] };
    v.e[0] = temp[0];
    v.e[1] = temp[1];
    v.e[2] = temp[2];
    v.e[3] = 0;
}

inline std::istream& operator>>(std::istream& is, point3<int>& t) {
    is >> t.e[0] >> t.e[1] >> t.e[2];
    return is;
}

inline std::ostream& operator<<(std::ostream& os, const point3<int>& t) {
    os << t.e[0] << ", " << t.e[1] << ", " << t.e[2];
    return os;
}

#endif


#ifdef RAY_FLOAT_AS_DOUBLE
typedef point3<double> point3f;
#else 
typedef point3<float> point3f;
#endif
typedef point3<int>   point3i;




#endif
