#ifndef VECTYPESH
#define VECTYPESH

#include "vec3.h"
#include "point3.h"
#include "normal.h"
#include "simd.h"

point3f convert_to_point3f(const vec3f& p);
point3f convert_to_point3f(const normal3f& p);
vec3f convert_to_vec3f(const point3f& p);
vec3f convert_to_vec3f(const normal3f& p);
normal3f convert_to_normal3f(const vec3f& p);
normal3f convert_to_normal3f(const point3f& p);

inline point3<Float> operator*(const point3<Float>& v1, const vec3<Float>& v2) {
    point3<Float> result;
    result.e = simd_mul(v1.e, v2.e);
    return result;
}

inline normal3f operator*(const normal3f& v1, const vec3<Float>& v2) {
    normal3f result;
    result.e = simd_mul(v1.e, v2.e);
    return result;
}

inline vec3f operator*=(const vec3f& v1, const point3f& v) {
    vec3f result;
    result.e = simd_mul(result.e, v.e);
    return result;
}

inline vec3f operator*=(const vec3f& v1, const normal3f& v) {
    vec3f result;
    result.e = simd_mul(result.e, v.e);
    return result;
}

inline point3f operator*=(const point3f& v1, const point3f& v) {
    point3f result;
    result.e = simd_mul(result.e, v.e);
    return result;
}

inline point3f operator*=(const point3f& v1, const normal3f& v) {
    point3f result;
    result.e = simd_mul(result.e, v.e);
    return result;
}

// inline normal3f operator*=(const normal3f& v1, const point3f& v) {
//     normal3f result;
//     result.e = simd_mul(result.e, v.e);
//     return result;
// }

// inline normal3f operator*=(const normal3f& v1, const vec3f& v) {
//     normal3f result;
//     result.e = simd_mul(result.e, v.e);
//     return result;
// }

#endif