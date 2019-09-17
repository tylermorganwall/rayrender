#ifndef MATHINLINEH
#define MATHINLINEH

#include <algorithm>
#include <cmath>

#ifdef RAY_FLOAT_AS_DOUBLE
typedef double Float;
#else
typedef float Float;
#endif 

template<class T>
inline T ffmin(T a, T b) { return(a < b ? a : b);}

template<class T>
inline T ffmax(T a, T b) { return(a > b ? a : b);}

template<class T>
inline T lerp(Float t, T v1, T v2) {
  return((1-t) * v1 + t * v2);
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

#endif
