#ifndef MATHINLINEH
#define MATHINLINEH

template<class T>
inline T lerp(float t, T v1, T v2) {
  return((1-t) * v1 + t * v2);
}

inline bool quadratic(float a, float b, float c, float *t0, float *t1) {
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