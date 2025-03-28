#ifndef RAYH
#define RAYH

#include "../math/float.h"
#include "../math/vec3.h"
#include "../math/point3.h"
#include <vector>
#include "../math/mathinline.h"
#include "../math/simd.h"

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

class Ray {
  public: 
    Ray() {}
    Ray(const point3f& a, const vec3f& b, 
        Float ti = 0.0, Float tmax = Infinity) {
      o = a;
      d = b;
      _time = ti;
      tMax = tmax;
      vec3 inv_dir = vec3f(1/b.xyz.x, 1/b.xyz.y, 1/b.xyz.z);
      inv_dir_pad.e[0] = add_ulp_magnitude(inv_dir.xyz.x, 2);
      inv_dir_pad.e[1] = add_ulp_magnitude(inv_dir.xyz.y, 2);
      inv_dir_pad.e[2] = add_ulp_magnitude(inv_dir.xyz.z, 2);
      kz = MaxDimension(Abs(d));
      kx = kz + 1;
      if (kx == 3) kx = 0;
      ky = kx + 1;
      if (ky == 3) ky = 0;
      vec3f dPermuted = Permute(d, kx, ky, kz);
      Float Sx = -dPermuted.xyz.x / dPermuted.xyz.z;
      Float Sy = -dPermuted.xyz.y / dPermuted.xyz.z;
      Float Sz = 1.f / dPermuted.xyz.z;
      Svec = vec3f(Sx, Sy, Sz);
#ifdef RAYSIMD
      origin4[0] = simd_set1(a.xyz.x);
      origin4[1] = simd_set1(a.xyz.y);
      origin4[2] = simd_set1(a.xyz.z);
      inv_dir_pad4[0] = simd_set1(inv_dir_pad.xyz.x);
      inv_dir_pad4[1] = simd_set1(inv_dir_pad.xyz.y);
      inv_dir_pad4[2] = simd_set1(inv_dir_pad.xyz.z);
      maskPos4[0] = simd_cmpge(inv_dir_pad4[0], simd_set1(0.0f));
      maskPos4[1] = simd_cmpge(inv_dir_pad4[1], simd_set1(0.0f));
      maskPos4[2] = simd_cmpge(inv_dir_pad4[2], simd_set1(0.0f));
#endif
    }
    Ray(const point3f& a, const vec3f& b,  std::vector<dielectric* > *priority2, 
        Float ti = 0.0, Float tmax = Infinity) {
      o = a; 
      d = b; 
      _time = ti;
      tMax = tmax;
      vec3 inv_dir = vec3f(1/b.xyz.x, 1/b.xyz.y, 1/b.xyz.z);
      inv_dir_pad.e[0] = add_ulp_magnitude(inv_dir.xyz.x, 2);
      inv_dir_pad.e[1] = add_ulp_magnitude(inv_dir.xyz.y, 2);
      inv_dir_pad.e[2] = add_ulp_magnitude(inv_dir.xyz.z, 2);
      pri_stack = priority2;
      kz = MaxDimension(Abs(d));
      kx = kz + 1;
      if (kx == 3) kx = 0;
      ky = kx + 1;
      if (ky == 3) ky = 0;
      vec3f dPermuted = Permute(d, kx, ky, kz);
      Float Sx = -dPermuted.xyz.x / dPermuted.xyz.z;
      Float Sy = -dPermuted.xyz.y / dPermuted.xyz.z;
      Float Sz = 1.f / dPermuted.xyz.z;
      Svec = vec3f(Sx, Sy, Sz);
#ifdef RAYSIMD
      origin4[0] = simd_set1(a.xyz.x);
      origin4[1] = simd_set1(a.xyz.y);
      origin4[2] = simd_set1(a.xyz.z);
      inv_dir_pad4[0] = simd_set1(inv_dir_pad.xyz.x);
      inv_dir_pad4[1] = simd_set1(inv_dir_pad.xyz.y);
      inv_dir_pad4[2] = simd_set1(inv_dir_pad.xyz.z);
      maskPos4[0] = simd_cmpge(inv_dir_pad4[0], simd_set1(0.0f));
      maskPos4[1] = simd_cmpge(inv_dir_pad4[1], simd_set1(0.0f));
      maskPos4[2] = simd_cmpge(inv_dir_pad4[2], simd_set1(0.0f));
#endif
    }
    point3f operator()(Float t) const { return o + d * t; }
    
    point3f origin() const {return(o);}
    vec3f direction() const {return(d);}
    Float time() const {return _time;}
    bool has_priority() const {return(pri_stack ? true : false);}
    size_t get_priority_size() const {return(pri_stack->size());}
    
    point3f o;
    vec3f d;
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

inline std::istream& operator>>(std::istream &is, Ray &r) {
  is >> r.o.e[0] >> r.o.e[1] >> r.o.e[2];
  return is;
}

inline std::ostream& operator<<(std::ostream &os, const Ray &r) {
  os << "Origin: " << r.o.e[0] << ", " << r.o.e[1] << ", " << r.o.e[2] << " Dir: " << r.d.e[0] << ", " << r.d.e[1] << ", " << r.d.e[2] ;
  return os;
}

#endif
