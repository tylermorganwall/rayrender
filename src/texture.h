#ifndef TEXTUREH
#define TEXTUREH

#include "perlin.h"

class texture {
public: 
  virtual vec3 value(float u, float v, const vec3& p) const = 0;
};

class constant_texture : public texture {
public: 
  constant_texture() {}
  constant_texture(vec3 c) : color(c) {}
  virtual vec3 value(float u, float v, const vec3& p) const {
    return(color);
  }
  vec3 color;
};

class checker_texture : public texture {
public:
  checker_texture() {}
  checker_texture(texture *t0, texture *t1, float p) : even(t0), odd(t1), period(p) {}
  virtual vec3 value(float u, float v, const vec3& p) const {
    float invperiod = 1.0/period;
    float sines  = sin(invperiod*p.x()*M_PI) * sin(invperiod*p.y()*M_PI) * sin(invperiod*p.z()*M_PI);
    if(sines < 0) {
      return(odd->value(u,v,p));
    } else {
      return(even->value(u,v,p));
    }
  }
  texture *even;
  texture *odd;
  float period;
};

class noise_texture : public texture {
public:
  noise_texture() {}
  noise_texture(float sc, vec3 c, vec3 c2, float ph, float inten) : scale(sc), color(c), color2(c2), phase(ph), intensity(inten) {}
  virtual vec3 value(float u, float v, const vec3& p) const {
    float weight = 0.5*(1+sin(scale*p.y()  + intensity*noise.turb(scale * p) + phase));
    return(color * (1-weight) + color2 * weight);
  }
  perlin noise;
  float scale;
  vec3 color;
  vec3 color2;
  float phase;
  float intensity;
};

class image_texture : public texture {
public:
  image_texture() {}
  image_texture(unsigned char *pixels, int A, int B) : data(pixels), nx(A), ny(B) {}
  virtual vec3 value(float u, float v, const vec3& p) const;
  unsigned char *data;
  int nx, ny;
};

vec3 image_texture::value(float u, float v, const vec3& p) const {
  int i = u * nx;
  int j = (1-v) * ny - 0.001;
  if (i < 0) i = 0;
  if (j < 0) j = 0;
  if (i > nx-1) i = nx-1;
  if (j > ny-1) j = ny-1;
  float r = int(data[3*i + 3*nx*j])/255.0;
  float g = int(data[3*i + 3*nx*j+1])/255.0;
  float b = int(data[3*i + 3*nx*j+2])/255.0;
  return(vec3(r,g,b));
}

#endif