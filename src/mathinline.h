#ifndef MATHINLINEH
#define MATHINLINEH

#include <algorithm>
#include <cmath>
#include <cstring>
#include "vec3.h"
#include "vec2.h"

template<class T>
inline T ffmin(T a, T b) { return(a < b ? a : b);}

template<class T>
inline T ffmax(T a, T b) { return(a > b ? a : b);}

template<class T>
inline T lerp(Float t, T v1, T v2) {
  return((1-t) * v1 + t * v2);
}

inline vec3 rand_to_unit(vec2 u) {
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
  return(vec3(r*cos(phi),r*sin(phi),0));
}

inline vec3 rand_cosine_direction(vec2 u) {
  Float r1 = u.x();
  Float r2 = u.y();
  Float z = std::sqrt(1.0-r2);
  Float phi = 2.0 * M_PI * r1;
  Float x = cos(phi) * std::sqrt(r2);
  Float y = sin(phi) * std::sqrt(r2);
  return(vec3(x, y, z));
}

inline vec3 rand_to_sphere(Float radius, Float distance_squared, vec2 u) {
  Float r1 = u.x();
  Float r2 = u.y();
  Float z = 1.0 + r2 * (std::sqrt(1.0-radius * radius / distance_squared) - 1);
  Float phi = 2.0 * M_PI * r1;
  Float x = std::cos(phi) * std::sqrt(1-z*z);
  Float y = std::sin(phi) * std::sqrt(1-z*z);
  return(vec3(x,y,z));
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

inline Float TanTheta(const vec3 &w) {
  return(SinTheta(w) / CosTheta(w));
}

inline Float Tan2Theta(const vec3 &w) {
  return(Sin2Theta(w) / Cos2Theta(w));
}

Float FrDielectric(Float cosThetaI, Float etaI, Float etaT) {
  cosThetaI = clamp(cosThetaI, -1, 1);
  bool entering = cosThetaI > 0.f;
  if (!entering) {
    std::swap(etaI, etaT);
    cosThetaI = std::abs(cosThetaI);
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

inline Float SphericalPhi(const vec3 &v) {
  float p = atan2f(v.x(), v.y());
  return((p < 0.0f) ? p + 2.0f*M_PI : p);
}

inline Float SphericalTheta(const vec3 &v) {
  return(acosf(clamp(v.z(), -1.0f, 1.0f)));
}

inline vec3 SphericalDirection(Float sinTheta, Float cosTheta, Float phi) {
  return(vec3(sinTheta * std::cos(phi),
              sinTheta * std::sin(phi),
              cosTheta));
}

inline bool SameHemisphere(const vec3 &w, const vec3 &wp) {
  return(w.z() * wp.z() > 0);
}


inline Float AbsDot(const vec3 &v1, const vec3 &v2) {
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

inline vec3 Reflect(const vec3 &wo, const vec3 &n) {
  return(-wo + 2 * dot(wo, n) * n);
}

constexpr Float origin() { return 1.0f / 32.0f; }
constexpr Float float_scale() { return 1.0f / 65536.0f; }
constexpr Float int_scale() { return 256.0f; }

vec3 offset_ray(const vec3 p, const vec3 n) {
  int of_i[3] = {(int)(int_scale() * n.x()), (int)(int_scale() * n.y()), (int)(int_scale() * n.z())};
  vec3 p_i(int_to_float(float_to_int(p.x())+((p.x() < 0) ? -of_i[0] : of_i[0])),
           int_to_float(float_to_int(p.y())+((p.y() < 0) ? -of_i[1] : of_i[1])),
           int_to_float(float_to_int(p.z())+((p.z() < 0) ? -of_i[2] : of_i[2])));
  return(vec3(fabs(p.x()) < origin() ? p.x() + float_scale()*n.x() : p_i.x(),
              fabs(p.y()) < origin() ? p.y() + float_scale()*n.y() : p_i.y(),
              fabs(p.z()) < origin() ? p.z() + float_scale()*n.z() : p_i.z()));
}

inline Float Log2(Float x) {
  const Float invLog2 = 1.442695040888963387004650940071;
  return(std::log(x) * invLog2);
}

inline vec3 RGBtoHSV(vec3& rgb) {
  Float max_val = ffmax(ffmax(rgb.r(), rgb.g()), rgb.b());
  Float min_val = ffmin(ffmin(rgb.r(), rgb.g()), rgb.b());
  Float delta_val = max_val - min_val;
  vec3 hsv;
  
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


inline vec3 HSVtoRGB(vec3 hsv) {
  Float chroma = hsv.z() * hsv.y(); 
  Float fHPrime = fmod(hsv.x() / 60.0, 6);
  Float x_val = chroma * (1 - fabs(fmod(fHPrime, 2) - 1));
  Float m_val = hsv.z() - chroma;

  if(0 <= fHPrime && fHPrime < 1) {
    vec3 rgb(chroma,x_val,0);
    rgb += m_val;
    return(rgb);
  } else if(1 <= fHPrime && fHPrime < 2) {
    vec3 rgb(x_val,chroma,0);
    rgb += m_val;
    return(rgb);
  } else if(2 <= fHPrime && fHPrime < 3) {
    vec3 rgb(0,chroma,x_val);
    rgb += m_val;
    return(rgb);
  } else if(3 <= fHPrime && fHPrime < 4) {
    vec3 rgb(0,x_val,chroma);
    rgb += m_val;
    return(rgb);
  } else if(4 <= fHPrime && fHPrime < 5) {
    vec3 rgb(x_val,0,chroma);
    rgb += m_val;
    return(rgb);
  } else if(5 <= fHPrime && fHPrime < 6) {
    vec3 rgb(chroma,0,x_val);
    rgb += m_val;
    return(rgb);
  } else {
    vec3 rgb(0,0,0);
    rgb += m_val;
    return(rgb);
  }
}


#endif
