#ifndef RNGH
#define RNGH

#include <random>
#include "pcg_random.hpp"
#include "vec3.h"



class random_gen {
  pcg32 rng;
public:
  random_gen() : rng(pcg_extras::seed_seq_from<std::random_device>{}) { }
  float unif_rand() {
    return float(rand()) * 0x1.0p-32f; 
  }
  vec3 random_in_unit_disk() {
    vec3 p;
    do {
      p = 2.0 * vec3(unif_rand(),unif_rand(),0) - vec3(1,1,0);
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
    float r1 = unif_rand();
    float r2 = unif_rand();
    float z = sqrt(1-r2);
    float phi = 2*M_PI*r1;
    float x = cos(phi)*2*sqrt(r2);
    float y = sin(phi)*2*sqrt(r2);
    return vec3(x, y, z);
  }
};


#endif