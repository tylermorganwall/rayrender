#ifndef PERLINH
#define PERLINH

#include "rng.h"
#include "point3.h"

inline Float perlin_interp(point3f c[2][2][2], Float u, Float v, Float w) {
  Float uu = u*u*(3-2*u);
  Float vv = v*v*(3-2*v);
  Float ww = w*w*(3-2*w);
  Float accum = 0;
  for(int i = 0; i < 2; i++) {
    for(int j = 0; j < 2; j++) {
      for(int k = 0; k < 2; k++) {
        point3f weight_v(u-i,v-j, w-k);
        accum += (i*uu + (1-i)*(1-uu))*(j*vv + (1-j)*(1-vv))*(k*ww + (1-k)*(1-ww))*dot(c[i][j][k], weight_v);
      } 
    }
  }
  return(accum);
}

class perlin {
public:
  perlin() {};
  Float noise(const point3f& p) const {
    Float u = p.x() - floor(p.x());
    Float v = p.y() - floor(p.y());
    Float w = p.z() - floor(p.z());
    int i = floor(p.x());
    int j = floor(p.y());
    int k = floor(p.z());
    point3f c[2][2][2];
    for(int di = 0; di < 2; di++) {
      for(int dj = 0; dj < 2; dj++) {
        for(int dk = 0; dk < 2; dk++) {
          c[di][dj][dk] = ranvec[perm_x[(i + di) & 255] ^ perm_y[(j + dj) & 255] ^ perm_z[(k + dk) & 255]];
        } 
      }
    }
    return(perlin_interp(c,u,v,w));
  }
  Float turb(const point3f& p, int depth=7) const {
    Float accum = 0;
    point3f temp_p = p;
    Float weight = 0.5;
    for(int i = 0; i < depth; i++) {
      accum += weight * noise(temp_p);
      weight *= 0.5;
      temp_p *= 2;
    }
    return(std::fabs(accum));
  }
  static point3f *ranvec;
  static int *perm_x;
  static int *perm_y;
  static int *perm_z;
  random_gen rng;
};


#endif

