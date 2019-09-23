#ifndef RNGH
#define RNGH

#include <random>
#define extended extended_rng
#include "pcg/pcg_random.hpp"
#undef extended
#include "vec3.h"



class random_gen {
  pcg32 rng;
public:
  random_gen(unsigned int seed) : rng(seed) {}
  random_gen() : rng(pcg_extras::seed_seq_from<std::random_device>{}) { }
  Float unif_rand() {
    return(std::ldexp(rng(),-32)); 
  }
  vec3 random_in_unit_disk() {
    vec3 p;
    do {
      p = 2.0 * vec3(unif_rand(), unif_rand(), 0) - vec3(1, 1, 0);
    } while (dot(p,p) >= 1.0);
    return(p);
  }
  vec3 random_in_unit_sphere() {
    vec3 p;
    do {
      p = 2.0 * vec3(unif_rand(),unif_rand(),unif_rand()) - vec3(1,1,1);
    } while (p.squared_length() >= 1.0);
    return(p);
  }
  vec3 random_cosine_direction() {
    Float r1 = unif_rand();
    Float r2 = unif_rand();
    Float z = std::sqrt(1-r2);
    Float phi = 2*M_PI*r1;
    Float x = cos(phi)*2*std::sqrt(r2);
    Float y = sin(phi)*2*std::sqrt(r2);
    return vec3(x, y, z);
  }
  vec3 random_to_sphere(Float radius, Float distance_squared) {
    Float r1 = unif_rand();
    Float r2 = unif_rand();
    Float z = 1 + r2 * (std::sqrt(1-radius * radius / distance_squared) - 1);
    Float phi = 2 * M_PI * r1;
    Float x = std::cos(phi) * std::sqrt(1-z*z);
    Float y = std::sin(phi) * std::sqrt(1-z*z);
    return(vec3(x,y,z));
  }
};


#endif
