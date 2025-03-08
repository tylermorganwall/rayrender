#ifndef RNGH
#define RNGH

#include <random>
#define extended extended_rng
#include "pcg/pcg_random.hpp"
#undef extended
#include "float.h"
#include "vec3.h"
#include "point3.h"

class random_gen {
public:
  random_gen(unsigned int seed) : rng(seed) {}
  random_gen() : rng(pcg_extras::seed_seq_from<std::random_device>{}) { }
  ~random_gen();
  Float unif_rand();
  vec3f random_in_unit_disk();
  vec3f random_in_unit_sphere();
  vec3f random_cosine_direction();
  vec3f random_to_sphere(Float radius, Float distance_squared);

  uint32_t UniformUInt32(uint32_t b);
  void SetSequence(unsigned int seed);
  pcg32 rng;
};


#endif
