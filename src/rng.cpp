#include "rng.h"

Float random_gen::unif_rand() {
  return(std::ldexp(rng(),-32)); 
}

vec3f random_gen::random_in_unit_disk() {
  vec3f p;
  do {
    p = 2.0 * vec3f(unif_rand(), unif_rand(), 0) - vec3f(1, 1, 0);
  } while (dot(p,p) >= 1.0);
  return(p);
}

vec3f random_gen::random_in_unit_sphere() {
  vec3f p;
  do {
    p = 2.0 * vec3f(unif_rand(),unif_rand(),unif_rand()) - vec3f(1,1,1);
  } while (p.squared_length() >= 1.0);
  return(p);
}

vec3f random_gen::random_cosine_direction() {
  Float r1 = unif_rand();
  Float r2 = unif_rand();
  Float z = std::sqrt(1.0-r2);
  Float phi = 2.0 * M_PI * r1;
  Float x = cos(phi) * std::sqrt(r2);
  Float y = sin(phi) * std::sqrt(r2);
  return(vec3f(x, y, z));
}

vec3f random_gen::random_to_sphere(Float radius, Float distance_squared) {
  Float r1 = unif_rand();
  Float r2 = unif_rand();
  Float z = 1.0 + r2 * (std::sqrt(1.0-radius * radius / distance_squared) - 1);
  Float phi = 2.0 * M_PI * r1;
  Float x = std::cos(phi) * std::sqrt(1-z*z);
  Float y = std::sin(phi) * std::sqrt(1-z*z);
  return(vec3f(x,y,z));
}

uint32_t random_gen::UniformUInt32(uint32_t b) {
  uint32_t threshold = (~b + 1u) % b;
  while (true) {
    uint32_t r = rng();
    if (r >= threshold)
      return r % b;
  }
}

void random_gen::SetSequence(unsigned int seed) {
  rng = seed;
}

random_gen::~random_gen() {}
