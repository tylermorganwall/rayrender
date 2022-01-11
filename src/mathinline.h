#ifndef MATHINLINEH
#define MATHINLINEH

#include <algorithm>
#include <cmath>
#include <cstring>
#include "float.h"
#include "vec3.h"
#include "vec2.h"
#include "point3.h"
#include "point2.h"

#include "normal.h"
#include <limits>
#include <array>

static const Float mpi_over_180 = M_PI/180;
static const Float SqrtPiOver8 = 0.626657069f;
static const Float ONE_OVER_2_PI = 1 / (2 * M_PI);
static const Float Infinity =  Infinity;

static constexpr Float MachineEpsilon = std::numeric_limits<Float>::epsilon() * 0.5;

inline constexpr Float gamma(int n) {
  return (n * MachineEpsilon) / (1 - n * MachineEpsilon);
}

template<class T>
inline T ffmin(T a, T b) { return(a < b ? a : b);}

template<class T>
inline T ffmax(T a, T b) { return(a > b ? a : b);}

template<class T>
inline T lerp(Float t, T v1, T v2) {
  return((1-t) * v1 + t * v2);
}

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

inline vec3f sgn(vec3f v) {
  return(vec3f(sgn(v.x()),sgn(v.y()),sgn(v.z())));
}

inline Float Radians(Float deg) { 
  return (M_PI / 180) * deg; 
}
inline Float Degrees(Float rad) { 
  return (180 / M_PI) * rad; 
}



inline point3f de_nan(const point3f& c) {
  point3f temp = c;
  if(std::isnan(c[0])) temp.e[0] = 1.0f;
  if(std::isnan(c[1])) temp.e[1] = 0.0f;
  if(std::isnan(c[2])) temp.e[2] = 1.0f;
  return(temp);
}



inline vec3f rand_to_unit(vec2f u) {
  Float a = 2.0*u.x() - 1.0; 
  Float b = 2.0*u.y() - 1.0;
  if (b == 0) {
    b = 1.0;
  }
  Float r, phi;
  if (a*a > b*b) {
    r = a;
    phi = (M_PI/4)*(b/a);
  } else {
    r = b;
    phi = (M_PI_2) - (M_PI_4)*(a/b);
  }
  return(vec3f(r*cos(phi),r*sin(phi),0));
}

inline vec3f rand_cosine_direction(vec2f u) {
  Float r1 = u.x();
  Float r2 = u.y();
  Float z = std::sqrt(1.0-r2);
  Float phi = 2.0 * M_PI * r1;
  Float x = cos(phi) * std::sqrt(r2);
  Float y = sin(phi) * std::sqrt(r2);
  return(vec3f(x, y, z));
}

inline vec3f rand_to_sphere(Float radius, Float distance_squared, vec2f u) {
  Float r1 = u.x();
  Float r2 = u.y();
  Float z = 1.0 + r2 * (std::sqrt(1.0-radius * radius / distance_squared) - 1);
  Float phi = 2.0 * M_PI * r1;
  Float x = std::cos(phi) * std::sqrt(1-z*z);
  Float y = std::sin(phi) * std::sqrt(1-z*z);
  return(vec3f(x,y,z));
}

inline vec3f SphericalDirection(Float sinTheta, 
                                   Float cosTheta, Float phi) {
  return vec3f(sinTheta * std::cos(phi), 
               sinTheta * std::sin(phi),
               cosTheta);
}

inline vec3f SphericalDirection(Float sinTheta, Float cosTheta, 
                                Float phi, const vec3f &x, const vec3f &y,
                                const vec3f &z) {
  return sinTheta * std::cos(phi) * x + sinTheta * std::sin(phi) * y + cosTheta * z;
}

// 
// point2f ConcentricSampleDisk(const point2f &u) {
//   point2f uOffset = 2.f * u - vec2f(1, 1);
//   
//   if (uOffset.x() == 0 && uOffset.y() == 0) {
//     return(point2f(0, 0));
//   }
//   
//   Float theta, r;
//   if (std::fabs(uOffset.x()) > std::fabs(uOffset.y())) {
//     r = uOffset.x();
//     theta = M_PI_4 * (uOffset.y() / uOffset.x());
//   } else {
//     r = uOffset.y();
//     theta = M_PI_2 - M_PI_4 * (uOffset.x() / uOffset.y());
//   }
//   return(r * point2f(std::cos(theta), std::sin(theta)));
//   
// }


// template<class T>
// inline T lerp(T t, T v1, T v2) {
//   return((1-t) * v1 + t * v2);
// }

inline Float saturate(Float v1) {
  v1 = v1 < 0 ? 0 : v1;
  v1 = v1 > 1 ? 1 : v1;
  return(v1);
}

inline vec3f saturate(vec3f v1) {
  v1.e[0] = saturate(v1.e[0]);
  v1.e[1] = saturate(v1.e[1]);
  v1.e[2] = saturate(v1.e[2]);
  return(v1);
}

inline vec3f clamp(const vec3f& c, Float clamplow, Float clamphigh) {
  vec3f temp = c;
  if(c.e[0] > clamphigh) {
    temp.e[0] = clamphigh;
  } else if(c[0] < clamplow) {
    temp.e[0] = clamplow;
  }
  if(c.e[1] > clamphigh) {
    temp.e[1] = clamphigh;
  } else if(c[1] < clamplow) {
    temp.e[1] = clamplow;
  }
  if(c.e[2] > clamphigh) {
    temp.e[2] = clamphigh;
  } else if(c[2] < clamplow) {
    temp.e[2] = clamplow;
  }
  return(temp);
}



inline point3f clamp_point(const point3f& c, Float clamplow, Float clamphigh) {
  point3f temp = c;
  if(c.e[0] > clamphigh) {
    temp.e[0] = clamphigh;
  } else if(c[0] < clamplow) {
    temp.e[0] = clamplow;
  }
  if(c.e[1] > clamphigh) {
    temp.e[1] = clamphigh;
  } else if(c[1] < clamplow) {
    temp.e[1] = clamplow;
  }
  if(c.e[2] > clamphigh) {
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

inline vec3f clamp(vec3f input, vec3f low, vec3f high) {
  vec3f final(clamp(input.x(), low.x(), high.x()),
              clamp(input.y(), low.y(), high.y()),
              clamp(input.z(), low.z(), high.z()));
  return(final);
}

inline point3f clamp(point3f input, point3f low, point3f high) {
  point3f final(clamp(input.x(), low.x(), high.x()),
                clamp(input.y(), low.y(), high.y()),
                clamp(input.z(), low.z(), high.z()));
  return(final);
}



// General Utility Functions
inline Float Sqr(Float v) { return v * v; }
template <int n>
static Float Pow(Float v) {
  static_assert(n > 0, "Power can't be negative");
  Float n2 = Pow<n / 2>(v);
  return n2 * n2 * Pow<n & 1>(v);
}

template <>
inline Float Pow<1>(Float v) {
  return(v);
}

template <>
inline Float Pow<0>(Float v) {
  return(1);
}

inline Float SafeASin(Float x) {
  return(std::asin(clamp(x, -1, 1)));
}

inline Float SafeSqrt(Float x) {
  return(std::sqrt(std::fmax(Float(0), x)));
}

inline vec3f Exp(vec3f a) {
  return(vec3f(std::exp(a.x()),std::exp(a.y()),std::exp(a.z())));
}

inline point3f Exp(point3f a) {
  return(point3f(std::exp(a.x()),std::exp(a.y()),std::exp(a.z())));
}

inline Float Logistic(Float x, Float s) {
  x = std::abs(x);
  return(std::exp(-x / s) / (s * Sqr(1 + std::exp(-x / s))));
}

inline Float LogisticCDF(Float x, Float s) {
  return(1 / (1 + std::exp(-x / s)));
}

inline Float TrimmedLogistic(Float x, Float s, Float a, Float b) {
  return(Logistic(x, s) / (LogisticCDF(b, s) - LogisticCDF(a, s)));
}


inline Float CosTheta(const vec3f &w) { return w.z(); }
inline Float Cos2Theta(const vec3f &w) { return w.z() * w.z(); }
inline Float AbsCosTheta(const vec3f &w) { return std::fabs(w.z()); }

inline Float Sin2Theta(const vec3f &w) {
  return(std::fmax((Float)0, (Float)1 - Cos2Theta(w)));
}

inline Float SinTheta(const vec3f &w) {
  return(std::sqrt(Sin2Theta(w)));
}

inline Float CosPhi(const vec3f &w) {
  Float sinTheta = SinTheta(w);
  return((sinTheta == 0) ? 1 : clamp(w.x() / sinTheta, -1, 1));
}

inline Float SinPhi(const vec3f &w) {
  Float sinTheta = SinTheta(w);
  return((sinTheta == 0) ? 0 : clamp(w.y() / sinTheta, -1, 1));
}

inline Float Cos2Phi(const vec3f &w) {
  return(CosPhi(w) * CosPhi(w));
}

inline Float Sin2Phi(const vec3f &w) {
  return(SinPhi(w) * SinPhi(w));
}

inline Float TanTheta(const vec3f &w) {
  return(SinTheta(w) / CosTheta(w));
}

inline Float Tan2Theta(const vec3f &w) {
  return(Sin2Theta(w) / Cos2Theta(w));
}

inline Float FrDielectric(Float cosThetaI, Float etaI, Float etaT) {
  cosThetaI = clamp(cosThetaI, -1, 1);
  bool entering = cosThetaI > 0.f;
  if (!entering) {
    std::swap(etaI, etaT);
    cosThetaI = std::fabs(cosThetaI);
  }
  Float sinThetaI = std::sqrt(std::fmax((Float)0, 1 - cosThetaI * cosThetaI));
  Float sinThetaT = etaI / etaT * sinThetaI;
  Float cosThetaT = std::sqrt(std::fmax((Float)0, 1 - sinThetaT * sinThetaT));
  Float Rparl = ((etaT * cosThetaI) - (etaI * cosThetaT)) /
    ((etaT * cosThetaI) + (etaI * cosThetaT));
  Float Rperp = ((etaI * cosThetaI) - (etaT * cosThetaT)) /
    ((etaI * cosThetaI) + (etaT * cosThetaT));
  return((Rparl * Rparl + Rperp * Rperp) / 2);
}



inline Float FrDielectric(Float cosThetaI, Float eta) {
  cosThetaI = clamp(cosThetaI, -1, 1);
  // Potentially flip interface orientation for Fresnel equations
  if (cosThetaI < 0) {
    eta = 1 / eta;
    cosThetaI = -cosThetaI;
  }
  
  // Compute $\cos\,\theta_\roman{t}$ for Fresnel equations using Snell's law
  Float sin2Theta_i = 1 - Sqr(cosThetaI);
  Float sin2Theta_t = sin2Theta_i * Sqr(eta);
  if (sin2Theta_t >= 1) return 1.0f;
  Float cosTheta_t = SafeSqrt(1 - sin2Theta_t);
  
  Float r_parl = (cosThetaI - eta * cosTheta_t) / (cosThetaI + eta * cosTheta_t);
  Float r_perp = (eta * cosThetaI - cosTheta_t) / (eta * cosThetaI + cosTheta_t);
  return (r_parl * r_parl + r_perp * r_perp) / 2;
}

inline bool quadratic(Float a, Float b, Float c, Float *t0, Float *t1) {
  Float discrim = DifferenceOfProducts(b, b, 4 * a, c);
  if (discrim < 0) {
    return(false);
  }
  Float rootDiscrim = std::sqrt(discrim);
  Float q = (b < 0) ? -0.5 * (b - rootDiscrim) : -0.5 * (b + rootDiscrim);
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

inline Float SphericalPhi(const vec3f &v) {
  float p = atan2f(v.x(), v.y());
  return((p < 0.0f) ? p + 2.0f*M_PI : p);
}

inline Float SphericalTheta(const vec3f &v) {
  return(acosf(clamp(v.z(), -1.0f, 1.0f)));
}

inline bool SameHemisphere(const vec3f &w, const vec3f &wp) {
  return(w.z() * wp.z() > 0);
}

inline bool SameHemisphere(const normal3f &w, const vec3f &wp) {
  return(w.z() * wp.z() > 0);
}

inline bool SameHemisphere(const vec3f &w, const normal3f &wp) {
  return(w.z() * wp.z() > 0);
}

inline bool refract(const vec3f& v, const vec3f& n, Float ni_over_nt, vec3f& refracted) {
  vec3f uv = unit_vector(v);
  Float dt = dot(uv, n);
  Float discriminant = 1.0 - ni_over_nt * ni_over_nt * (1 - dt * dt);
  if(discriminant > 0) {
    refracted = ni_over_nt * (uv - n * dt) - n * sqrt(discriminant);
    return(true);
  } else {
    return(false);
  }
}

inline bool refract(const vec3f& v, const normal3f& n2, Float ni_over_nt, vec3f& refracted) {
  vec3f uv = unit_vector(v);
  vec3f n = vec3f(n2.x(),n2.y(),n2.z());
  Float dt = dot(uv, n);
  Float discriminant = 1.0 - ni_over_nt * ni_over_nt * (1 - dt * dt);
  if(discriminant > 0) {
    refracted = ni_over_nt * (uv - n * dt) - n * sqrt(discriminant);
    return(true);
  } else {
    return(false);
  }
}

inline vec3f refract(const vec3f& uv, const vec3f& n, Float ni_over_nt) {
  Float cos_theta = dot(-uv, n);
  vec3f r_out_parallel =  ni_over_nt * (uv + cos_theta*n);
  vec3f r_out_perp = -sqrt(1.0 - r_out_parallel.squared_length()) * n;
  return(r_out_parallel + r_out_perp);
}

inline vec3f refract(const vec3f& uv, const normal3f& n, Float ni_over_nt) {
  vec3f n2 = vec3f(n.x(),n.y(),n.z());
  Float cos_theta = dot(-uv, n);
  vec3f r_out_parallel =  ni_over_nt * (uv + cos_theta*n2);
  vec3f r_out_perp = -sqrt(1.0 - r_out_parallel.squared_length()) * n2;
  return(r_out_parallel + r_out_perp);
}


template <typename T>
inline Float AbsDot(const vec3<T> &v1, const vec3<T> &v2) {
  return(std::fabs(dot(v1,v2)));
} 


template <typename IN_T, typename OUT_T>
inline OUT_T reinterpret_type(const IN_T in) {
  OUT_T out;
  memcpy(&out, &in, sizeof(out));
  return out;
}

inline Float int_to_float(const int in) {
  Float out;
  memcpy(&out, &in, sizeof(out));
  return out;
}

inline int float_to_int(const Float in) {
  int out;
  memcpy(&out, &in, sizeof(out));
  return out;
}




inline Float ErfInv(Float x) {
  Float w, p;
  x = clamp(x, -.99999f, .99999f);
  w = -std::log((1 - x) * (1 + x));
  if (w < 5) {
    w = w - 2.5f;
    p = 2.81022636e-08f;
    p = 3.43273939e-07f + p * w;
    p = -3.5233877e-06f + p * w;
    p = -4.39150654e-06f + p * w;
    p = 0.00021858087f + p * w;
    p = -0.00125372503f + p * w;
    p = -0.00417768164f + p * w;
    p = 0.246640727f + p * w;
    p = 1.50140941f + p * w;
  } else {
    w = std::sqrt(w) - 3;
    p = -0.000200214257f;
    p = 0.000100950558f + p * w;
    p = 0.00134934322f + p * w;
    p = -0.00367342844f + p * w;
    p = 0.00573950773f + p * w;
    p = -0.0076224613f + p * w;
    p = 0.00943887047f + p * w;
    p = 1.00167406f + p * w;
    p = 2.83297682f + p * w;
  }
  return p * x;
}

inline Float Erf(Float x) {
  // constants
  Float a1 = 0.254829592f;
  Float a2 = -0.284496736f;
  Float a3 = 1.421413741f;
  Float a4 = -1.453152027f;
  Float a5 = 1.061405429f;
  Float p = 0.3275911f;
  
  // Save the sign of x
  int sign = 1;
  if (x < 0) sign = -1;
  x = std::fabs(x);
  
  // A&S formula 7.1.26
  Float t = 1 / (1 + p * x);
  Float y = 1 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * std::exp(-x * x);
  
  return sign * y;
}

template <typename T>
inline vec3<T> Reflect(const vec3<T> &wo, const normal3<T> &n) {
  normal3<T> n2 = 2 * dot(wo, n) * n;
  return(-wo + vec3<T>(n2.x(),n2.y(),n2.z()));
}

template <typename T>
inline vec3<T> Reflect(const vec3<T> &wo, const vec3<T> &n) {
  return(-wo + 2 * dot(wo, n) * n);
}

inline uint32_t FloatToBits(float f) {
  uint32_t ui;
  memcpy(&ui, &f, sizeof(float));
  return ui;
}

inline float BitsToFloat(uint32_t ui) {
  float f;
  memcpy(&f, &ui, sizeof(uint32_t));
  return f;
}


inline float NextFloatUp(float v) {
  // Handle infinity and negative zero for _NextFloatUp()_
  if (std::isinf(v) && v > 0.) return v;
  if (v == -0.f) v = 0.f;
  
  // Advance _v_ to next higher float
  uint32_t ui = FloatToBits(v);
  if (v >= 0) {
    ++ui;
  } else {
    --ui;
  }
  return BitsToFloat(ui);
}

inline float NextFloatDown(float v) {
  // Handle infinity and positive zero for _NextFloatDown()_
  if (std::isinf(v) && v < 0.) return v;
  if (v == 0.f) v = -0.f;
  uint32_t ui = FloatToBits(v);
  if (v > 0)
    --ui;
  else
    ++ui;
  return BitsToFloat(ui);
}

inline uint64_t FloatToBits(double f) {
  uint64_t ui;
  memcpy(&ui, &f, sizeof(double));
  return ui;
}

inline double BitsToFloat(uint64_t ui) {
  double f;
  memcpy(&f, &ui, sizeof(uint64_t));
  return f;
}

inline double NextFloatUp(double v, int delta = 1) {
  if (std::isinf(v) && v > 0.) return v;
  if (v == -0.f) v = 0.f;
  uint64_t ui = FloatToBits(v);
  if (v >= 0.)
    ui += delta;
  else
    ui -= delta;
  return BitsToFloat(ui);
}

inline double NextFloatDown(double v, int delta = 1) {
  if (std::isinf(v) && v < 0.) return v;
  if (v == 0.f) v = -0.f;
  uint64_t ui = FloatToBits(v);
  if (v > 0.)
    ui -= delta;
  else
    ui += delta;
  return BitsToFloat(ui);
}

static const Float MaxT = NextFloatDown(std::numeric_limits<Float>::infinity());
static const Float MinT = NextFloatUp((Float)0);
constexpr Float origin() { return 1.0f / 32.0f; }
constexpr Float float_scale() { return 1.0f / 65536.0f; }
constexpr Float int_scale() { return 256.0f; }

inline point3f offset_ray(const vec3f p, const normal3f n) {
  int of_i[3] = {(int)(int_scale() * n.x()), (int)(int_scale() * n.y()), (int)(int_scale() * n.z())};
  point3f p_i(int_to_float(float_to_int(p.x())+((p.x() < 0) ? -of_i[0] : of_i[0])),
              int_to_float(float_to_int(p.y())+((p.y() < 0) ? -of_i[1] : of_i[1])),
              int_to_float(float_to_int(p.z())+((p.z() < 0) ? -of_i[2] : of_i[2])));
  return(point3f(std::fabs(p.x()) < origin() ? p.x() + float_scale()*n.x() : p_i.x(),
                 std::fabs(p.y()) < origin() ? p.y() + float_scale()*n.y() : p_i.y(),
                 std::fabs(p.z()) < origin() ? p.z() + float_scale()*n.z() : p_i.z()));
}

inline Float Log2(Float x) {
  const Float invLog2 = 1.442695040888963387004650940071;
  return(std::log(x) * invLog2);
}

inline point3f RGBtoHSV(point3f& rgb) {
  Float max_val = ffmax(ffmax(rgb.r(), rgb.g()), rgb.b());
  Float min_val = ffmin(ffmin(rgb.r(), rgb.g()), rgb.b());
  Float delta_val = max_val - min_val;
  point3f hsv;
  
  if(delta_val > 0) {
    if(max_val == rgb.r()) {
      hsv.e[0] = 60 * (fmod(((rgb.g() - rgb.b()) / delta_val), 6));
    } else if(max_val == rgb.g()) {
      hsv.e[0] = 60 * (((rgb.b() - rgb.r()) / delta_val) + 2);
    } else if(max_val == rgb.b()) {
      hsv.e[0] = 60 * (((rgb.r() - rgb.g()) / delta_val) + 4);
    }
    if(max_val > 0) {
      hsv.e[1] = delta_val / max_val;
    } else {
      hsv.e[1] = 0;
    }
    hsv.e[2] = max_val;
  } else {
    hsv.e[0] = 0;
    hsv.e[1] = 0;
    hsv.e[2] = max_val;
  }
  if(hsv.e[0] < 0) {
    hsv.e[0] = 360 + hsv.e[0];
  }
  return(hsv);
}


inline point3f HSVtoRGB(point3f hsv) {
  Float chroma = hsv.z() * hsv.y(); 
  Float fHPrime = fmod(hsv.x() / 60.0, 6);
  Float x_val = chroma * (1 - std::fabs(fmod(fHPrime, 2) - 1));
  Float m_val = hsv.z() - chroma;

  if(0 <= fHPrime && fHPrime < 1) {
    point3f rgb(chroma,x_val,0);
    rgb += m_val;
    return(rgb);
  } else if(1 <= fHPrime && fHPrime < 2) {
    point3f rgb(x_val,chroma,0);
    rgb += m_val;
    return(rgb);
  } else if(2 <= fHPrime && fHPrime < 3) {
    point3f rgb(0,chroma,x_val);
    rgb += m_val;
    return(rgb);
  } else if(3 <= fHPrime && fHPrime < 4) {
    point3f rgb(0,x_val,chroma);
    rgb += m_val;
    return(rgb);
  } else if(4 <= fHPrime && fHPrime < 5) {
    point3f rgb(x_val,0,chroma);
    rgb += m_val;
    return(rgb);
  } else if(5 <= fHPrime && fHPrime < 6) {
    point3f rgb(chroma,0,x_val);
    rgb += m_val;
    return(rgb);
  } else {
    point3f rgb(0,0,0);
    rgb += m_val;
    return(rgb);
  }
}


//Hair utilities
static const int pMax = 3;


inline uint32_t Compact1By1(uint32_t x) {
  // TODO: as of Haswell, the PEXT instruction could do all this in a
  // single instruction.
  // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
  x &= 0x55555555;
  // x = --fe --dc --ba --98 --76 --54 --32 --10
  x = (x ^ (x >> 1)) & 0x33333333;
  // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
  x = (x ^ (x >> 2)) & 0x0f0f0f0f;
  // x = ---- ---- fedc ba98 ---- ---- 7654 3210
  x = (x ^ (x >> 4)) & 0x00ff00ff;
  // x = ---- ---- ---- ---- fedc ba98 7654 3210
  x = (x ^ (x >> 8)) & 0x0000ffff;
  return(x);
}

inline vec2f DemuxFloat(Float f) {
  uint64_t v = f * (1ull << 32);
  uint32_t bits[2] = {Compact1By1(v), Compact1By1(v >> 1)};
  return(vec2f(bits[0] / Float(1 << 16), bits[1] / Float(1 << 16)));
}

inline Float I0(Float x) {
  Float val = 0;
  Float x2i = 1;
  int64_t ifact = 1;
  int i4 = 1;
  // I0(x) \approx Sum_i x^(2i) / (4^i (i!)^2)
  for (int i = 0; i < 10; ++i) {
    if (i > 1) {
      ifact *= i;
    }
    val += x2i / (i4 * Sqr(ifact));
    x2i *= x * x;
    i4 *= 4;
  }
  return(val);
}

inline Float LogI0(Float x) {
  if (x > 12) {
    return(x + 0.5 * (-std::log(2 * M_PI) + std::log(1 / x) + 1 / (8 * x)));
  } else {
    return(std::log(I0(x)));
  }
}

inline std::array<point3f, pMax + 1> Ap(Float cosThetaO, Float eta, Float h,
                                     const point3f &T) {
  std::array<point3f, pMax + 1> ap;
  // Compute $p=0$ attenuation at initial cylinder intersection
  Float cosGammaO = SafeSqrt(1 - h * h);
  Float cosTheta = cosThetaO * cosGammaO;
  Float f = FrDielectric(cosTheta, 1.f, eta);
  ap[0] = f;
  
  // Compute $p=1$ attenuation term
  ap[1] = Sqr(1 - f) * T;
  
  // Compute attenuation terms up to $p=_pMax_$
  for (int p = 2; p < pMax; ++p) {
    ap[p] = ap[p - 1] * T * f;
  }
  
  // Compute attenuation term accounting for remaining orders of scattering
  ap[pMax] = ap[pMax - 1] * f * T / (point3f(1.0f,1.0f,1.0f) + -T * f);
  return(ap);
}


inline Float Mp(Float cosThetaI, Float cosThetaO, Float sinThetaI,
                Float sinThetaO, Float v) {
  Float a = cosThetaI * cosThetaO / v;
  Float b = sinThetaI * sinThetaO / v;
  Float mp = (v <= .1)
    ? (std::exp(LogI0(a) - b - 1 / v + 0.6931f + std::log(1 / (2 * v))))
    : (std::exp(-b) * I0(a)) / (std::sinh(1 / v) * 2 * v);
  return(mp);
}

inline Float Phi(int p, Float gammaO, Float gammaT) {
  return(2 * p * gammaT - 2 * gammaO + p * M_PI);
}

inline Float Np(Float phi, int p, Float s, Float gammaO, Float gammaT) {
  Float dphi = phi - Phi(p, gammaO, gammaT);
  // Remap _dphi_ to $[-\pi,\pi]$
  while (dphi > M_PI) dphi -= 2 * M_PI;
  while (dphi < -M_PI) dphi += 2 * M_PI;
  return(TrimmedLogistic(dphi, s, -M_PI, M_PI));
}

inline Float SampleTrimmedLogistic(Float u, Float s, Float a, Float b) {
  Float k = LogisticCDF(b, s) - LogisticCDF(a, s);
  Float x = -s * std::log(1 / (u * k + LogisticCDF(a, s)) - 1);
  return(clamp(x, a, b));
}

//Sampling helpers


inline int CountTrailingZeros(uint32_t v) {
// #if defined(PBRT_IS_MSVC)
//   unsigned long index;
//   if (_BitScanForward(&index, v))
//     return index;
//   else
//     return 32;
// #else
  return __builtin_ctz(v);
// #endif
}

// Low Discrepancy Inline Functions
inline uint32_t ReverseBits32(uint32_t n) {
  n = (n << 16) | (n >> 16);
  n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8);
  n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4);
  n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2);
  n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1);
  return n;
}

inline uint64_t ReverseBits64(uint64_t n) {
  uint64_t n0 = ReverseBits32((uint32_t)n);
  uint64_t n1 = ReverseBits32((uint32_t)(n >> 32));
  return (n0 << 32) | n1;
}

template <typename T> inline void
CoordinateSystem(const vec3<T> &v1, vec3<T> *v2, vec3<T> *v3) {
  if (std::fabs(v1.x()) > std::fabs(v1.y())) {
    *v2 = vec3<T>(-v1.z(), 0, v1.x()) /
      std::sqrt(v1.x() * v1.x() + v1.z() * v1.z());
  } else {
    *v2 = vec3<T>(0, v1.z(), -v1.y()) /
      std::sqrt(v1.y() * v1.y() + v1.z() * v1.z());
  }
  *v3 = cross(v1, *v2);
}

inline point3f OffsetRayOrigin(const point3f &p, const vec3f &pError,
                               const normal3f &n, const vec3f &w) {
  Float d = dot(Abs(n), pError);
  vec3f offset = d * vec3f(n.x(),n.y(),n.z());
  if (dot(w, n) < 0) {
    offset = -offset;
  }
  point3f po = p + offset;
  for (int i = 0; i < 3; ++i) {
    if (offset.e[i] > 0) {
      po.e[i] = NextFloatUp(po.e[i]);
    } else if (offset.e[i] < 0)  {
      po.e[i] = NextFloatDown(po.e[i]);
    }
  }
  
  return po;
}


inline Float UniformConePdf(Float cosThetaMax) {
  return 1 / (2 * M_PI * (1 - cosThetaMax));
}

inline vec3f UniformSampleCone(const point2f &u, Float cosThetaMax) {
  Float cosTheta = ((Float)1 - u[0]) + u[0] * cosThetaMax;
  Float sinTheta = std::sqrt((Float)1 - cosTheta * cosTheta);
  Float phi = u[1] * 2 * M_PI;
  return vec3f(std::cos(phi) * sinTheta, std::sin(phi) * sinTheta,
               cosTheta);
}

inline vec3f UniformSampleCone(const point2f &u, Float cosThetaMax,
                        const vec3f &x, const vec3f &y,
                        const vec3f &z) {
  Float cosTheta = lerp(u[0], cosThetaMax, (Float)1.);
  Float sinTheta = std::sqrt((Float)1. - cosTheta * cosTheta);
  Float phi = u[1] * 2 * M_PI;
  return std::cos(phi) * sinTheta * x + std::sin(phi) * sinTheta * y +
    cosTheta * z;
}

inline vec3f UniformSampleHemisphere(const vec2f &u) {
  Float z = u[0];
  Float r = std::sqrt(std::fmax((Float)0, (Float)1. - z * z));
  Float phi = 2 * M_PI * u[1];
  return vec3f(r * std::cos(phi), r * std::sin(phi), z);
}

inline Float UniformHemispherePdf() { return 0.5*M_1_PI; }

inline vec3f UniformSampleSphere(const vec2f &u) {
  Float z = 1 - 2 * u[0];
  Float r = std::sqrt(std::fmax((Float)0, (Float)1 - z * z));
  Float phi = 2 * M_PI * u[1];
  return vec3f(r * std::cos(phi), r * std::sin(phi), z);
}

inline Float UniformSpherePdf() { return 0.25*M_1_PI; }

inline point2f UniformSampleDisk(const vec2f &u) {
  Float r = std::sqrt(u[0]);
  Float theta = 2 * M_PI * u[1];
  return point2f(r * std::cos(theta), r * std::sin(theta));
}

template <typename T>
inline normal3<T> Faceforward(const normal3<T> &n, const normal3<T> &n2) {
  return (dot(n, n2) < 0.f) ? -n : n;
}

template <typename T>
inline vec3<T> Faceforward(const vec3<T> &v, const vec3<T> &v2) {
  return (dot(v, v2) < 0.f) ? -v : v;
}

template <typename T>
inline vec3<T> Faceforward(const vec3<T> &v, const normal3<T> &n2) {
  return (dot(v, n2) < 0.f) ? -v : v;
}

#endif
