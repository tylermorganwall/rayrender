#ifndef RAYH
#define RAYH
#include "vec3.h"
#include <cstring>

template <typename IN_T, typename OUT_T>
inline OUT_T reinterpret_type(const IN_T in) {
  OUT_T out;
  memcpy(&out, &in, sizeof(out));
  return out;
}
inline Float add_ulp_magnitude(Float f, int ulps) {
  if (!std::isfinite(f)) return f;
  const unsigned bits = reinterpret_type<Float, unsigned>(f);
  return reinterpret_type<unsigned, Float>(bits + ulps);
}

// static constexpr Float MachineEpsilon =
//   std::numeric_limits<Float>::epsilon() * 0.5;
// 
// inline constexpr Float gamma(int n) {
//   return (n * MachineEpsilon) / (1 - n * MachineEpsilon);
// }

class ray {
  public: 
    ray() {}
    ray(const vec3& a, const vec3& b, Float ti = 0.0) {
      A = a; 
      B = b; 
      _time = ti;
      inv_dir = vec3(1/b.x(), 1/b.y(), 1/b.z());
      inv_dir_pad.e[0] = add_ulp_magnitude(inv_dir.x(), 2);
      inv_dir_pad.e[1] = add_ulp_magnitude(inv_dir.y(), 2);
      inv_dir_pad.e[2] = add_ulp_magnitude(inv_dir.z(), 2);
      sign[0] = (inv_dir.x() < 0);
      sign[1] = (inv_dir.y() < 0);
      sign[2] = (inv_dir.z() < 0);
    }
    
    vec3 origin() const {return(A);}
    vec3 direction() const {return(B);}
    vec3 inverse_dir() const {return(inv_dir);}
    Float time() const {return _time;}
    vec3 point_at_parameter(Float t) const {return(A + t*B);}
    
    vec3 A;
    vec3 B;
    vec3 inv_dir;
    vec3 inv_dir_pad;
    int sign[3];
    Float _time;
};

#endif
