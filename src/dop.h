#ifndef DOPH
#define DOPH

#include <cmath>

template<typename T>
inline T DifferenceOfProducts(T a, T b, T c, T d) {
    return a * b - c * d;
}

inline float DifferenceOfProducts(float a, float b, float c, float d) {
    // Use std::fmaf for improved precision
    // Computes (a * b) - (c * d) with reduced rounding error
    float cd = c * d;
    float err = std::fmaf(-c, d, cd); // Error term from c * d
    float dop = std::fmaf(a, b, -cd); // (a * b) - cd
    return dop + err;                 // Corrected result
}

inline double DifferenceOfProducts(double a, double b, double c, double d) {
  double cd = c * d;
  double err = std::fma(-c, d, cd);
  double dop = std::fma(a, b, -cd);
  return(dop + err);
}

#endif