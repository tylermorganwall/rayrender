#ifndef MATHINLINEH
#define MATHINLINEH

#include <algorithm>
#include <cmath>
#include <cstring>
#include "vec3.h"

template<class T>
inline T ffmin(T a, T b) { return(a < b ? a : b);}

template<class T>
inline T ffmax(T a, T b) { return(a > b ? a : b);}

template<class T>
inline T lerp(Float t, T v1, T v2) {
  return((1-t) * v1 + t * v2);
}

template<class T>
inline T lerp(T t, T v1, T v2) {
  return((1-t) * v1 + t * v2);
}

inline Float saturate(Float v1) {
  v1 = v1 < 0 ? 0 : v1;
  v1 = v1 > 1 ? 1 : v1;
  return(v1);
}

inline vec3 saturate(vec3 v1) {
  v1.e[0] = saturate(v1.e[0]);
  v1.e[1] = saturate(v1.e[1]);
  v1.e[2] = saturate(v1.e[2]);
  return(v1);
}

inline vec3 clamp(const vec3& c, Float clamplow, Float clamphigh) {
  vec3 temp = c;
  if(c[0] > clamphigh) {
    temp.e[0] = clamphigh;
  } else if(c[0] < clamplow) {
    temp.e[0] = clamplow;
  }
  if(c[1] > clamphigh) {
    temp.e[1] = clamphigh;
  } else if(c[1] < clamplow) {
    temp.e[1] = clamplow;
  }
  if(c[2] > clamphigh) {
    temp.e[2] = clamphigh;
  } else if(c[2] < clamplow) {
    temp.e[2] = clamplow;
  }
  return(temp);
}

inline Float clamp(const Float& c, Float clamplow, Float clamphigh) {
  if(c > clamphigh) {
    return(clamphigh);
  } else if(c < clamplow) {
    return(clamplow);
  }
  return(c);
}

inline Float CosTheta(const vec3 &w) { return w.z(); }
inline Float Cos2Theta(const vec3 &w) { return w.z() * w.z(); }
inline Float AbsCosTheta(const vec3 &w) { return std::abs(w.z()); }

inline Float Sin2Theta(const vec3 &w) {
  return(std::max((Float)0, (Float)1 - Cos2Theta(w)));
}

inline Float SinTheta(const vec3 &w) {
  return(std::sqrt(Sin2Theta(w)));
}

inline Float CosPhi(const vec3 &w) {
  Float sinTheta = SinTheta(w);
  return((sinTheta == 0) ? 1 : clamp(w.x() / sinTheta, -1, 1));
}

inline Float SinPhi(const vec3 &w) {
  Float sinTheta = SinTheta(w);
  return((sinTheta == 0) ? 0 : clamp(w.y() / sinTheta, -1, 1));
}

inline Float Cos2Phi(const vec3 &w) {
  return(CosPhi(w) * CosPhi(w));
}

inline Float Sin2Phi(const vec3 &w) {
  return(SinPhi(w) * SinPhi(w));
}

inline bool quadratic(Float a, Float b, Float c, Float *t0, Float *t1) {
  double discrim = (double)b * (double)b - 4 * (double)a * (double)c;
  if (discrim < 0) {
    return(false);
  }
  double rootDiscrim = std::sqrt(discrim);
  double q = (b < 0) ? -0.5 * (b - rootDiscrim) : -0.5 * (b + rootDiscrim);
  *t0 = q / a;
  *t1 = c / q;
  if (*t0 > *t1) {
    std::swap(*t0, *t1);
  }
  return(true);
}

//interval search, pbrt utilities
template <typename Predicate> int FindInterval(int size, const Predicate &pred) {
  int first = 0, len = size;
  while (len > 0) {
    int half = len >> 1, middle = first + half;
    if (pred(middle)) {
      first = middle + 1;
      len -= half + 1;
    } else {
      len = half;
    }
  }
  return clamp(first - 1, 0, size - 2);
}

inline Float SphericalPhi(const vec3 &v) {
  float p = atan2f(v.x(), v.z());
  return((p < 0.0f) ? p + 2.0f*M_PI : p);
}

inline Float SphericalTheta(const vec3 &v) {
  return(acosf(clamp(v.y(), -1.0f, 1.0f)));
}

template <typename IN_T, typename OUT_T>
inline OUT_T reinterpret_type(const IN_T in) {
  OUT_T out;
  memcpy(&out, &in, sizeof(out));
  return out;
}

// constexpr Float origin() { return 1.0f / 32.0f; }
// constexpr Float float_scale() { return 1.0f / 65536.0f; }
// constexpr Float int_scale() { return 256.0f; }
// 
// vec3 offset_ray(const vec3 p, const vec3 n) {
//   int of_i[3] = {(int)(int_scale() * n.x()), (int)(int_scale() * n.y()), (int)(int_scale() * n.z())};
//   vec3 p_i(reinterpret_type<int,Float>(reinterpret_type<Float,int>(p.x())+((p.x() < 0) ? -of_i[0] : of_i[0])), 
//            reinterpret_type<int,Float>(reinterpret_type<Float,int>(p.y())+((p.y() < 0) ? -of_i[1] : of_i[1])), 
//            reinterpret_type<int,Float>(reinterpret_type<Float,int>(p.z())+((p.z() < 0) ? -of_i[2] : of_i[2])));
//   return(vec3(fabs(p.x()) < origin() ? p.x() + float_scale()*n.x() : p_i.x(), 
//               fabs(p.y()) < origin() ? p.y() + float_scale()*n.y() : p_i.y(),
//               fabs(p.z()) < origin() ? p.z() + float_scale()*n.z() : p_i.z()));
// }

#endif
