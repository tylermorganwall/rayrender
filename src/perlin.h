#ifndef PERLINH
#define PERLINH

class perlin {
public:
  float noise(const vec3& p) const {
    float u = p.x() - floor(p.x());
    float v = p.y() - floor(p.y());
    float w = p.z() - floor(p.z());
    int i = int(4*p.x()) & 255;
    int j = int(4*p.y()) & 255;
    int k = int(4*p.z()) & 255;
    return(ranfloat[perm_x[i] ^ perm_y[j] ^ perm_z[k]]);
  }
  static float *ranfloat;
  static int *perm_x;
  static int *perm_y;
  static int *perm_z;
};

static float* perlin_generate() {
  float *p  = new float[256];
  for(int i = 0; i < 256; i++) {
    p[i] = drand48();
  }
  return(p);
}

void permute(int *p, int n) {
  for(int i = n-1; i > 0; i--) {
    int target = int(drand48()*(i+1));
    int tmp = p[i];
    p[i] = p[target];
    p[target] = tmp;
  }
  return;
}

static int* perlin_generate_perm() {
  int * p = new int[256];
  for(int i = 0; i < 256; i++) {
    p[i] = i;
  }
  permute(p, 256);
  return(p);
}

float *perlin::ranfloat = perlin_generate();
int *perlin::perm_x = perlin_generate_perm();
int *perlin::perm_y = perlin_generate_perm();
int *perlin::perm_z = perlin_generate_perm();
#endif

