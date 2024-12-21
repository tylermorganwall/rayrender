#ifndef NORMALH
#define NORMALH

#include <iostream>
#include <cmath>
#include "vec3.h"
#include "point3.h"
#include "float.h"
#include "simd.h" 


#ifndef RAYSIMDVEC

class normal3f {
public:
  normal3f() {}
  normal3f(Float e0, Float e1, Float e2) {e[0] = e0; e[1] = e1; e[2] = e2;}
  normal3f(Float e0) {e[0] = e0; e[1] = e0; e[2] = e0;}
  
  inline Float x() const { return e[0]; }
  inline Float y() const { return e[1]; }
  inline Float z() const { return e[2]; }
  inline Float r() const { return e[0]; }
  inline Float g() const { return e[1]; }
  inline Float b() const { return e[2]; }
  
  inline const normal3f& operator+() const { return *this; }
  
  inline normal3f operator-() const { return normal3f(-e[0], -e[1], -e[2]); }
  inline Float operator[](int i) const { return e[i]; }
  inline Float& operator[](int i) { return e[i]; }
  
  inline normal3f& operator+=(const normal3f &v2);
  inline normal3f& operator-=(const normal3f &v2);
  inline normal3f& operator*=(const normal3f &v2);
  inline normal3f& operator/=(const normal3f &v2);
  inline normal3f& operator+=(const Float t);
  inline normal3f& operator-=(const Float t);
  inline normal3f& operator*=(const Float t);
  inline normal3f& operator/=(const Float t);
  
  bool HasNaNs() const {
    return(std::isnan(e[0]) || std::isnan(e[1]) || std::isnan(e[2]));
  }
  
  inline Float length() const { return std::sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]); }
  inline Float squared_length() const { return e[0]*e[0] + e[1]*e[1] + e[2]*e[2]; }
  inline normal3f pow(Float exponent) const {
    return(normal3f(std::pow(e[0],exponent),std::pow(e[1],exponent),std::pow(e[2],exponent)));
  }
  inline void make_unit_vector();
  
  Float e[3];
};


inline void normal3f::make_unit_vector() {
  // Float len = length();
  // if(len - 1 < 1e-8) {
  //   volatile int f = 1;
  // }
  Float k = 1.0f / std::sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
  e[0] *= k; e[1] *= k; e[2] *= k; 
}


inline std::istream& operator>>(std::istream &is, normal3f &t) {
  is >> t.e[0] >> t.e[1] >> t.e[2];
  return is;
}


inline std::ostream& operator<<(std::ostream &os, const normal3f &t) {
  os << t.e[0] << ", " << t.e[1] << ", " << t.e[2];
  return os;
}


inline normal3f operator+(const normal3f &v1, const normal3f &v2) {
  return normal3f(v1.e[0] + v2.e[0],v1.e[1] + v2.e[1],v1.e[2] + v2.e[2]);
}


inline normal3f operator-(const normal3f &v1, const normal3f &v2) {
  return normal3f(v1.e[0] - v2.e[0],v1.e[1] - v2.e[1],v1.e[2] - v2.e[2]);
}


inline normal3f operator*(const normal3f &v1, const normal3f &v2) {
  return normal3f(v1.e[0] * v2.e[0],v1.e[1] * v2.e[1],v1.e[2] * v2.e[2]);
}


inline normal3f operator/(const normal3f &v1, const normal3f &v2) {
  return normal3f(v1.e[0] / v2.e[0],v1.e[1] / v2.e[1],v1.e[2] / v2.e[2]);
}


inline normal3f operator+(const normal3f &v, Float t) {
  return normal3f(v.e[0] + t,v.e[1] + t,v.e[2] + t);
}


inline normal3f operator-(const normal3f &v, Float t) {
  return normal3f(v.e[0] - t,v.e[1] - t,v.e[2] - t);
}


inline normal3f operator*(Float t, const normal3f &v) {
  return normal3f(t*v.e[0], t*v.e[1], t*v.e[2]);
}


inline normal3f operator*(const normal3f &v, Float t) {
  return normal3f(t*v.e[0], t*v.e[1], t*v.e[2]);
}


inline normal3f operator/(const normal3f &v, Float t) {
  return normal3f(v.e[0]/t, v.e[1]/t, v.e[2]/t);
}


inline normal3f operator/(Float t,const normal3f &v) {
  return normal3f(t/v.e[0], t/v.e[1], t/v.e[2]);
}


inline Float dot(const normal3f &v1, const normal3f &v2) {
  return (v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2]);
}


inline Float dot(const vec3f &v1, const normal3f &v2) {
  return (v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2]);
}


inline Float dot(const normal3f &v1, const vec3f &v2) {
  return (v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2]);
}


inline Float AbsDot(const vec3f &v1, const normal3f &v2) {
  return std::fabs(v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2]);
}


inline Float AbsDot(const normal3f &v1, const vec3f &v2) {
  return std::fabs(v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2]);
}


inline normal3f cross(const normal3f &v1, const normal3f &v2) {
  return(normal3f(DifferenceOfProducts(v1.y(), v2.z(), v1.z(), v2.y()),
                    DifferenceOfProducts(v1.z(), v2.x(), v1.x(), v2.z()),
                    DifferenceOfProducts(v1.x(), v2.y(), v1.y(), v2.x())));
}


inline normal3f cross(const vec3f &v1, const normal3f &v2) {
  return(normal3f(DifferenceOfProducts(v1.y(), v2.z(), v1.z(), v2.y()),
                    DifferenceOfProducts(v1.z(), v2.x(), v1.x(), v2.z()),
                    DifferenceOfProducts(v1.x(), v2.y(), v1.y(), v2.x())));
}


inline normal3f cross(const normal3f &v1, const vec3f &v2) {
  return(normal3f(DifferenceOfProducts(v1.y(), v2.z(), v1.z(), v2.y()),
                    DifferenceOfProducts(v1.z(), v2.x(), v1.x(), v2.z()),
                    DifferenceOfProducts(v1.x(), v2.y(), v1.y(), v2.x())));
}




inline normal3f& normal3f::operator+=(const normal3f &v) {
  e[0] += v.e[0];
  e[1] += v.e[1];
  e[2] += v.e[2];
  return(*this);
}


inline normal3f& normal3f::operator*=(const normal3f &v) {
  e[0] *= v.e[0];
  e[1] *= v.e[1];
  e[2] *= v.e[2];
  return(*this);
}


inline normal3f& normal3f::operator/=(const normal3f &v) {
  e[0] /= v.e[0];
  e[1] /= v.e[1];
  e[2] /= v.e[2];
  return(*this);
}


inline normal3f& normal3f::operator-=(const normal3f &v) {
  e[0] -= v.e[0];
  e[1] -= v.e[1];
  e[2] -= v.e[2];
  return(*this);
}


inline normal3f& normal3f::operator+=(const Float t) {
  e[0] += t;
  e[1] += t;
  e[2] += t;
  return(*this);
}


inline normal3f& normal3f::operator*=(const Float t) {
  e[0] *= t;
  e[1] *= t;
  e[2] *= t;
  return(*this);
}

inline normal3f& normal3f::operator/=(const Float t) {
  Float k = 1.0/t;
  
  e[0] *= k;
  e[1] *= k;
  e[2] *= k;
  return(*this);
}


inline normal3f unit_vector(normal3f v) {
  //   Float len = v.length();
  // if(len - 1 < 1e-8) {
  //   volatile int f = 1;
  // }
  return(v/v.length());
}


inline Float MinComponent(const normal3f &v) {
  return(ffmin(v.x(), ffmin(v.y(), v.z())));
}


inline Float MaxComponent(const normal3f &v) {
  return(ffmax(v.x(), ffmax(v.y(), v.z())));
}


inline int MaxDimension(const normal3f &v) {
  return((v.x() > v.y()) ? ((v.x() > v.z()) ? 0 : 2) : ((v.y() > v.z()) ? 1 : 2));
}


inline normal3f Min(const normal3f &p1, const normal3f &p2) {
  return(normal3f(ffmin(p1.x(), p2.x()), ffmin(p1.y(), p2.y()),ffmin(p1.z(), p2.z())));
}


inline normal3f Max(const normal3f &p1, const normal3f &p2) {
  return(normal3f(fmax(p1.x(), p2.x()), fmax(p1.y(), p2.y()),fmax(p1.z(), p2.z())));
}


inline normal3f Permute(const normal3f &v, int x, int y, int z) {
  return(normal3f(v.e[x], v.e[y], v.e[z]));
}


inline normal3f Abs(const normal3f &v) {
  return(normal3f(fabs(v.x()), fabs(v.y()), fabs(v.z())));
}

inline normal3f Faceforward(const normal3f &n, const vec3f &v) {
  return (dot(n, v) < 0.f) ? -n : n;
}

inline normal3f Faceforward(const normal3f &n, const normal3f &v) {
  return (dot(n, v) < 0.f) ? -n : n;
}

#else

class alignas(16) normal3f {
public:
    FVec4 e;  // SIMD vector for storage

    // Constructors
    normal3f() {
        e = simd_set1(0.0f);
    }

    normal3f(Float e0, Float e1, Float e2) {
        alignas(16) float values[4] = { e0, e1, e2, 0.0f };
        e = simd_load(values);
    }

    normal3f(Float e0) {
        e = simd_set1(e0);
    }

    template <typename U>
    normal3f(const vec3<U>& p) {
        e = simd_set(static_cast<Float>(p.x()), 
                     static_cast<Float>(p.y()), 
                     static_cast<Float>(p.z()), 
                     0.0f);
    }

    // Accessors
    inline Float x() const { return e[0]; }
    inline Float y() const { return e[1]; }
    inline Float z() const { return e[2]; }
    inline Float r() const { return e[0]; }
    inline Float g() const { return e[1]; }
    inline Float b() const { return e[2]; }

    // Unary operators
    inline const normal3f& operator+() const { return *this; }

    inline normal3f operator-() const {
        normal3f result;
        result.e = simd_sub(simd_set1(0.0f), e);
        return result;
    }

    inline Float operator[](int i) const { return e[i]; }
    inline Float& operator[](int i) { return e[i]; }

    // Compound assignment operators
    inline normal3f& operator+=(const normal3f& v);
    inline normal3f& operator-=(const normal3f& v);
    inline normal3f& operator*=(const normal3f& v);
    inline normal3f& operator/=(const normal3f& v);
    inline normal3f& operator+=(const Float t);
    inline normal3f& operator-=(const Float t);
    inline normal3f& operator*=(const Float t);
    inline normal3f& operator/=(const Float t);

    inline normal3f& operator*=(const vec3f& v);
    inline normal3f& operator*=(const point3f& v);

    // Methods
    inline Float length() const;
    inline Float squared_length() const;
    inline normal3f pow(Float exponent) const;
    inline vec3f convert_to_vec3() const;
    inline void make_unit_vector();
    inline bool HasNaNs() const;
};


inline bool normal3f::HasNaNs() const {
    return (std::isnan(e[0]) || std::isnan(e[1]) || std::isnan(e[2]));
}

inline Float normal3f::squared_length() const {
    FVec4 mul = simd_mul(e, e);
    return mul[0] + mul[1] + mul[2];
}

inline Float normal3f::length() const {
    return std::sqrt(squared_length());
}

inline normal3f normal3f::pow(Float exponent) const {
    normal3f result;
    result.e[0] = std::pow(e[0], exponent);
    result.e[1] = std::pow(e[1], exponent);
    result.e[2] = std::pow(e[2], exponent);
    result.e[3] = 0.0f;
    return result;
}

inline vec3f normal3f::convert_to_vec3() const {
    vec3f result;
    result.e = e;  // Copy the SIMD vector
    return result;
}

inline void normal3f::make_unit_vector() {
  //   Float len = length();
  // if(len - 1 < 1e-8) {
  //   volatile int f = 1;
  // }
    Float len = length();
    e = simd_div(e, simd_set1(len));
}

inline normal3f& normal3f::operator+=(const normal3f& v) {
    e = simd_add(e, v.e);
    return *this;
}

inline normal3f& normal3f::operator-=(const normal3f& v) {
    e = simd_sub(e, v.e);
    return *this;
}

inline normal3f& normal3f::operator*=(const normal3f& v) {
    e = simd_mul(e, v.e);
    return *this;
}

inline normal3f& normal3f::operator/=(const normal3f& v) {
    e = simd_div(e, v.e);
    return *this;
}

inline normal3f& normal3f::operator+=(const Float t) {
    e = simd_add(e, simd_set1(t));
    return *this;
}

inline normal3f& normal3f::operator-=(const Float t) {
    e = simd_sub(e, simd_set1(t));
    return *this;
}

inline normal3f& normal3f::operator*=(const Float t) {
    e = simd_mul(e, simd_set1(t));
    return *this;
}

inline normal3f& normal3f::operator/=(const Float t) {
    e = simd_div(e, simd_set1(t));
    return *this;
}

inline normal3f& normal3f::operator*=(const point3f& v) {
    e = simd_mul(e, v.e);
    return *this;
}

inline normal3f& normal3f::operator*=(const vec3f& v) {
    e = simd_mul(e, v.e);
    return *this;
}

inline normal3f Faceforward(const normal3f& n, const vec3f& v) {
    // Compute the dot product (scalar value)
    Float d = simd_dot(n.e, v.e);

    // Use std::copysign to determine the sign
    Float sgn_d = std::copysign(1.0f, d);

    // Multiply the normal vector by the sign
    normal3f result;
    result.e = simd_mul(n.e, simd_set1(sgn_d));

    return result;
}

inline normal3f Faceforward(const normal3f& n, const normal3f& v) {
    // Compute the dot product (scalar value)
    Float d = simd_dot(n.e, v.e);

    // Use std::copysign to determine the sign
    Float sgn_d = std::copysign(1.0f, d);

    // Multiply the normal vector by the sign
    normal3f result;
    result.e = simd_mul(n.e, simd_set1(sgn_d));

    return result;
}

inline normal3f operator+(const normal3f& v1, const normal3f& v2) {
    normal3f result;
    result.e = simd_add(v1.e, v2.e);
    return result;
}

inline normal3f operator-(const normal3f& v1, const normal3f& v2) {
    normal3f result;
    result.e = simd_sub(v1.e, v2.e);
    return result;
}

inline normal3f operator*(const normal3f& v1, const normal3f& v2) {
    normal3f result;
    result.e = simd_mul(v1.e, v2.e);
    return result;
}

inline normal3f operator/(const normal3f& v1, const normal3f& v2) {
    normal3f result;
    result.e = simd_div(v1.e, v2.e);
    return result;
}

inline normal3f operator+(const normal3f& v, Float t) {
    normal3f result;
    result.e = simd_add(v.e, simd_set1(t));
    return result;
}

inline normal3f operator-(const normal3f& v, Float t) {
    normal3f result;
    result.e = simd_sub(v.e, simd_set1(t));
    return result;
}

inline normal3f operator*(const normal3f& v, Float t) {
    normal3f result;
    result.e = simd_mul(v.e, simd_set1(t));
    return result;
}

inline normal3f operator*(Float t, const normal3f& v) {
    normal3f result;
    result.e = simd_mul(simd_set1(t), v.e);
    return result;
}

inline normal3f operator/(const normal3f& v, Float t) {
    normal3f result;
    result.e = simd_div(v.e, simd_set1(t));
    return result;
}

inline normal3f operator/(Float t, const normal3f& v) {
    normal3f result;
    result.e = simd_div(simd_set1(t), v.e);
    return result;
}


inline Float dot(const normal3f& v1, const normal3f& v2) {
    return simd_dot(v1.e, v2.e);
}

inline Float dot(const vec3f& v1, const normal3f& v2) {
    return simd_dot(v1.e, v2.e);
}

inline Float dot(const normal3f& v1, const vec3f& v2) {
    return simd_dot(v1.e, v2.e);
}

inline Float AbsDot(const vec3f& v1, const normal3f& v2) {
    return ffabs(simd_dot(v1.e, v2.e));
}

inline Float AbsDot(const normal3f& v1, const vec3f& v2) {
    return ffabs(simd_dot(v1.e, v2.e));
}

inline normal3f cross(const normal3f& v1, const normal3f& v2) {
    normal3f result;
    result.e = simd_cross(v1.e, v2.e);
    return result;
}

inline normal3f cross(const vec3f& v1, const normal3f& v2) {
    normal3f result;
    result.e = simd_cross(v1.e, v2.e);
    return result;
}

inline normal3f cross(const normal3f& v1, const vec3f& v2) {
    normal3f result;
    result.e = simd_cross(v1.e, v2.e);
    return result;
}

inline normal3f unit_vector(const normal3f& v) {
  //   Float len = v.length();
  // if(len - 1 < 1e-8) {
  //   volatile int f = 1;
  // }
    return v / v.length();
}

inline Float MinComponent(const normal3f& v) {
    return std::fmin(v.x(), std::fmin(v.y(), v.z()));
}

inline Float MaxComponent(const normal3f& v) {
    return std::fmax(v.x(), std::fmax(v.y(), v.z()));
}

inline int MaxDimension(const normal3f& v) {
    if (v.x() > v.y()) {
        return (v.x() > v.z()) ? 0 : 2;
    } else {
        return (v.y() > v.z()) ? 1 : 2;
    }
}

inline normal3f Min(const normal3f& p1, const normal3f& p2) {
    normal3f result;
    result.e = simd_min(p1.e, p2.e);
    return result;
}

inline normal3f Max(const normal3f& p1, const normal3f& p2) {
    normal3f result;
    result.e = simd_max(p1.e, p2.e);
    return result;
}

inline normal3f Permute(const normal3f& v, int x, int y, int z) {
    normal3f result;
    result.e[0] = v.e[x];
    result.e[1] = v.e[y];
    result.e[2] = v.e[z];
    result.e[3] = 0.0f;
    return result;
}

inline std::istream& operator>>(std::istream& is, normal3f& t) {
    is >> t.e[0] >> t.e[1] >> t.e[2];
    return is;
}

inline std::ostream& operator<<(std::ostream& os, const normal3f& t) {
    os << t.e[0] << ", " << t.e[1] << ", " << t.e[2];
    return os;
}

inline normal3f Abs(const normal3f& v) {
    normal3f result;
    result.e = simd_abs(v.e);
    return result;
}

#endif

#endif
