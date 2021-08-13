#include "perlin.h"


void permute(int *p, int n) {
#ifdef RAY_REPRODUCE_PERLIN
  random_gen rng(1);
#else
  random_gen rng;
#endif
  for(int i = n-1; i > 0; i--) {
    int target = int(rng.unif_rand()*(i+1));
    int tmp = p[i];
    p[i] = p[target];
    p[target] = tmp;
  }
  return;
}
int* perlin_generate_perm() {
  int * p = new int[256];
  for(int i = 0; i < 256; i++) {
    p[i] = i;
  }
  permute(p, 256);
  return(p);
}

point3f* perlin_generate() {
#ifdef RAY_REPRODUCE_PERLIN
  random_gen rng(1);
#else
  random_gen rng;
#endif
  point3f *p  = new point3f[256];
  for(int i = 0; i < 256; i++) {
    p[i] = unit_vector(point3f(-1.0f + 2*rng.unif_rand(),-1.0f + 2*rng.unif_rand(),-1.0f + 2*rng.unif_rand()));
  }
  return(p);
}

point3f *perlin::ranvec = perlin_generate();
int *perlin::perm_x = perlin_generate_perm();
int *perlin::perm_y = perlin_generate_perm();
int *perlin::perm_z = perlin_generate_perm();
