#ifndef RAYH
#define RAYH

#include "float.h"
#include "vec3.h"
#include "point3.h"
#include <vector>
#include "mathinline.h"
#include "simd.h"

class dielectric;

#ifdef RAY_FLOAT_AS_DOUBLE
inline Float add_ulp_magnitude(Float f, int ulps) {
  if (!std::isfinite(f)) return f;
  const unsigned long long bits = reinterpret_type<Float, unsigned long long>(f);
  return reinterpret_type<unsigned long long, Float>(bits + ulps);
}
#else
inline Float add_ulp_magnitude(Float f, int ulps) {
  if (!std::isfinite(f)) return f;
  const unsigned bits = reinterpret_type<Float, unsigned>(f);
  return reinterpret_type<unsigned, Float>(bits + ulps);
}
#endif

class ray {
  public: 
    ray() {}
    ray(const point3f& a, const vec3f& b, 
        Float ti = 0.0, Float tmax = Infinity) {
      A = a;
      B = b;
      _time = ti;
      tMax = tmax;
      vec3 inv_dir = vec3f(1/b.x(), 1/b.y(), 1/b.z());
      inv_dir_pad.e[0] = add_ulp_magnitude(inv_dir.x(), 2);
      inv_dir_pad.e[1] = add_ulp_magnitude(inv_dir.y(), 2);
      inv_dir_pad.e[2] = add_ulp_magnitude(inv_dir.z(), 2);
      kz = MaxDimension(Abs(B));
      kx = kz + 1;
      if (kx == 3) kx = 0;
      ky = kx + 1;
      if (ky == 3) ky = 0;
      vec3f dPermuted = Permute(B, kx, ky, kz);
      Float Sx = -dPermuted.x() / dPermuted.z();
      Float Sy = -dPermuted.y() / dPermuted.z();
      Float Sz = 1.f / dPermuted.z();
      Svec = vec3f(Sx, Sy, Sz);
#ifdef RAYSIMD
      origin4[0] = simd_set1(a.x());
      origin4[1] = simd_set1(a.y());
      origin4[2] = simd_set1(a.z());
      inv_dir_pad4[0] = simd_set1(inv_dir_pad.x());
      inv_dir_pad4[1] = simd_set1(inv_dir_pad.y());
      inv_dir_pad4[2] = simd_set1(inv_dir_pad.z());
      maskPos4[0] = simd_cmpge(inv_dir_pad4[0], simd_set1(0.0f));
      maskPos4[1] = simd_cmpge(inv_dir_pad4[1], simd_set1(0.0f));
      maskPos4[2] = simd_cmpge(inv_dir_pad4[2], simd_set1(0.0f));
#endif
    }
    ray(const point3f& a, const vec3f& b,  std::vector<dielectric* > *priority2, 
        Float ti = 0.0, Float tmax = Infinity) {
      A = a; 
      B = b; 
      _time = ti;
      tMax = tmax;
      vec3 inv_dir = vec3f(1/b.x(), 1/b.y(), 1/b.z());
      inv_dir_pad.e[0] = add_ulp_magnitude(inv_dir.x(), 2);
      inv_dir_pad.e[1] = add_ulp_magnitude(inv_dir.y(), 2);
      inv_dir_pad.e[2] = add_ulp_magnitude(inv_dir.z(), 2);
      pri_stack = priority2;
      kz = MaxDimension(Abs(B));
      kx = kz + 1;
      if (kx == 3) kx = 0;
      ky = kx + 1;
      if (ky == 3) ky = 0;
      vec3f dPermuted = Permute(B, kx, ky, kz);
      Float Sx = -dPermuted.x() / dPermuted.z();
      Float Sy = -dPermuted.y() / dPermuted.z();
      Float Sz = 1.f / dPermuted.z();
      Svec = vec3f(Sx, Sy, Sz);
#ifdef RAYSIMD
      origin4[0] = simd_set1(a.x());
      origin4[1] = simd_set1(a.y());
      origin4[2] = simd_set1(a.z());
      inv_dir_pad4[0] = simd_set1(inv_dir_pad.x());
      inv_dir_pad4[1] = simd_set1(inv_dir_pad.y());
      inv_dir_pad4[2] = simd_set1(inv_dir_pad.z());
      maskPos4[0] = simd_cmpge(inv_dir_pad4[0], simd_set1(0.0f));
      maskPos4[1] = simd_cmpge(inv_dir_pad4[1], simd_set1(0.0f));
      maskPos4[2] = simd_cmpge(inv_dir_pad4[2], simd_set1(0.0f));
#endif
    }
    point3f operator()(Float t) const { return A + B * t; }
    
    point3f origin() const {return(A);}
    vec3f direction() const {return(B);}
    // vec3f inverse_dir() const {return(inv_dir);}
    Float time() const {return _time;}
    point3f point_at_parameter(Float t) const {return(A + t*B);}
    bool has_priority() const {return(pri_stack ? true : false);}
    size_t get_priority_size() const {return(pri_stack->size());}
    
    point3f A;
    vec3f B;
    // vec3f inv_dir;
    vec3f inv_dir_pad;
  #ifdef RAYSIMD
    FVec4 origin4[3];
    FVec4 inv_dir_pad4[3];
    SimdMask maskPos4[3];
  #endif

    // int sign[3];
    vec3f Svec;
    Float _time;
    // Float Sx, Sy, Sz;
    int kx, ky, kz;
    mutable Float tMax;
    std::vector<dielectric*> *pri_stack;
};

inline std::istream& operator>>(std::istream &is, ray &r) {
  is >> r.A.e[0] >> r.A.e[1] >> r.A.e[2];
  return is;
}

inline std::ostream& operator<<(std::ostream &os, const ray &r) {
  os << "Origin: " << r.A.e[0] << ", " << r.A.e[1] << ", " << r.A.e[2] << " Dir: " << r.B.e[0] << ", " << r.B.e[1] << ", " << r.B.e[2] ;
  return os;
}

#endif
