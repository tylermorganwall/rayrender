#ifndef VECTYPESH
#define VECTYPESH

#include "../math/vec3.h"
#include "../math/point3.h"
#include "../math/normal.h"
#include "../math/simd.h"
#include "../math/vec2.h"
#include "../math/point2.h"

#ifdef RAYSIMDVEC

// point3f convert_to_point3(const vec3f& p);
// point3f convert_to_point3(const normal3f& p);
// vec3f convert_to_vec3(const point3f& p);
// vec3f convert_to_vec3(const normal3f& p);
// normal3f convert_to_normal3(const vec3f& p);
// normal3f convert_to_normal3(const point3f& p);


inline point3f convert_to_point3(const vec3f& p) {
    point3f result;
    result.e = p.e;  // Assuming 'e' is publicly accessible or there are public methods to access it
    return result;
}

inline point3f convert_to_point3(const normal3f& p) {
    point3f result;
    result.e = p.e;  // Assuming 'e' is publicly accessible or there are public methods to access it
    return result;
}


inline vec3f convert_to_vec3(const point3f& p) {
    vec3f result;
    result.e = p.e;  // Assuming 'e' is publicly accessible or there are public methods to access it
    return result;
}

inline vec3f convert_to_vec3(const normal3f& p) {
    vec3f result;
    result.e = p.e;  // Assuming 'e' is publicly accessible or there are public methods to access it
    return result;
}

inline normal3f convert_to_normal3(const vec3f& p) {
    normal3f result;
    result.e = p.e;  // Assuming 'e' is publicly accessible or there are public methods to access it
    return result;
}

inline normal3f convert_to_normal3(const point3f& p) {
    normal3f result;
    result.e = p.e;  // Assuming 'e' is publicly accessible or there are public methods to access it
    return result;
}


inline point3<Float> operator*(const point3<Float>& v1, const vec3<Float>& v2) {
    point3<Float> result;
    result.e = simd_mul(v1.e, v2.e);
    return result;
}

inline point3<Float> operator*(const point3<Float>& v1, const normal3f& v2) {
    point3<Float> result;
    result.e = simd_mul(v1.e, v2.e);
    return result;
}

inline normal3f operator*(const normal3f& v1, const vec3<Float>& v2) {
    normal3f result;
    result.e = simd_mul(v1.e, v2.e);
    return result;
}


inline point3f operator*=(const point3f& v1, const point3f& v2) {
    point3f result;
    result.e = simd_mul(v1.e, v2.e);
    return result;
}

inline point3f operator*=(const point3f& v1, const normal3f& v2) {
    point3f result;
    result.e = simd_mul(v1.e, v2.e);
    return result;
}


#else

template <typename T>
inline point3<T> convert_to_point3(const vec3<T>& p) {
    point3<T> result;
    result.e[0] = p.x;result.e[1] = p.y;result.e[2] = p.z;
    return result;
}

inline point3f convert_to_point3(const normal3f& p) {
    point3f result;
    result.e[0] = p.x;result.e[1] = p.y;result.e[2] = p.z;
    return result;
}

template <typename T>
inline vec3<T> convert_to_vec3(const point3<T>& p) {
    vec3<T> result;
    result.e[0] = p.x;result.e[1] = p.y;result.e[2] = p.z;
    return result;
}

inline vec3f convert_to_vec3(const normal3f& p) {
    vec3f result;
    result.e[0] = p.x;result.e[1] = p.y;result.e[2] = p.z;
    return result;
}

template <typename T>
normal3f convert_to_normal3(const vec3<T>& p) {
    normal3f result;
    result.e[0] = p.x;result.e[1] = p.y;result.e[2] = p.z;
    return result;
}

template <typename T>
inline normal3f convert_to_normal3(const point3<T>& p) {
    normal3f result;
    result.e[0] = p.x; result.e[1] = p.y; result.e[2] = p.z;
    return result;
}

template<typename T> 
inline point3<T> operator*(const point3<T>& v1, const vec3<T>& v2) {
    point3<T> result;
    result.e[0] = v1.x * v2.x; result.e[1] = v1.y * v2.y; result.e[2] = v1.z * v2.z; 
    return result;
}

template<typename T> 
inline point3<T> operator*(const point3<T>& v1, const normal3f& v2) {
    point3<T> result;
    result.e[0] = v1.x * v2.x; result.e[1] = v1.y * v2.y; result.e[2] = v1.z * v2.z; 
    return result;
}

template<typename T> 
inline point3<T> operator+(const vec3<T>& v1, const point3<T>& v2) {
    point3<T> result;
    result.e[0] = v1.x + v2.x; result.e[1] = v1.y + v2.y;result.e[2] = v1.z + v2.z;
    return result;
}

template<typename T> 
inline normal3f operator*(const normal3f& v1, const vec3<T>& v2) {
    normal3f result;
    result.e[0] = v1.x * v2.x; result.e[1] = v1.y * v2.y; result.e[2] = v1.z * v2.z; 
    return result;
}


template<typename T> 
inline vec3<T> operator*=(const vec3<T>& v1, const point3<T>& v2) {
    vec3<T> result;
    result.e[0] = v1.x * v2.x; result.e[1] = v1.y * v2.y; result.e[2] = v1.z * v2.z; 
    return result;
}

template<typename T> 
inline vec3<T> operator*=(const vec3<T>& v1, const normal3f& v2) {
    vec3<T> result;
    result.e[0] = v1.x * v2.x; result.e[1] = v1.y * v2.y; result.e[2] = v1.z * v2.z; 
    return result;
}

template<typename T> 
inline point3<T> operator*=(const point3<T>& v1, const point3<T>& v2) {
    point3<T> result;
    result.e[0] = v1.x * v2.x; result.e[1] = v1.y * v2.y; result.e[2] = v1.z * v2.z; 
    return result;
}


template<typename T> 
inline point3<T> operator*=(const point3<T>& v1, const normal3f& v2) {
    point3<T> result;
    result.e[0] = v1.x * v2.x; result.e[1] = v1.y * v2.y; result.e[2] = v1.z * v2.z; 
    return result;
}

template<typename T> 
inline normal3f operator*=(const normal3f& v1, const point3<T>& v2) {
    normal3f result;
    result.e[0] = v1.x * v2.x; result.e[1] = v1.y * v2.y; result.e[2] = v1.z * v2.z; 
    return result;
}


template<typename T> 
inline normal3f operator*=(const normal3f& v1, const vec3<T>& v2) {
    normal3f result;
    result.e[0] = v1.x * v2.x; result.e[1] = v1.y * v2.y; result.e[2] = v1.z * v2.z; 
    return result;
}

#endif

#endif