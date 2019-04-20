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
  noise_texture(float sc, vec3 c, vec3 c2, float ph, float inten, random_gen rng) : 
    scale(sc), color(c), color2(c2), phase(ph), intensity(inten) {
    noise = new perlin(rng);
  }
  ~noise_texture() {
    if(noise != 0) delete noise;
  }
  virtual vec3 value(float u, float v, const vec3& p) const {
    float weight = 0.5*(1+sin(scale*p.y()  + intensity*noise->turb(scale * p) + phase));
    return(color * (1-weight) + color2 * weight);
  }
  perlin *noise;
  float scale;
  vec3 color;
  vec3 color2;
  float phase;
  float intensity;
  random_gen rng;
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

class triangle_texture : public texture {
public:
  triangle_texture() {}
  triangle_texture(vec3 a, vec3 b, vec3 c) : a(a), b(b), c(c) {}
  virtual vec3 value(float u, float v, const vec3& p) const {
    return((1 - u - v) * a + u * b + v * c);
  };
  vec3 a,b,c;
};

class triangle_image_texture : public texture {
public:
  triangle_image_texture() {}
  triangle_image_texture(unsigned char *pixels, int A, int B,
                         float tex_u_a, float tex_v_a,
                         float tex_u_b, float tex_v_b,
                         float tex_u_c, float tex_v_c) : data(pixels),  nx(A), ny(B),
                         a_u(tex_u_a), a_v(tex_v_a), 
                         b_u(tex_u_b), b_v(tex_u_b), 
                         c_u(tex_u_c), c_v(tex_v_c) {}
  virtual vec3 value(float u, float v, const vec3& p) const;
  
  unsigned char *data;
  int nx, ny;
  float a_u, a_v, b_u, b_v,c_u, c_v;
};

vec3 triangle_image_texture::value(float u, float v, const vec3& p) const {
  int i = a_u * nx;
  int j = (1-a_v) * ny - 0.001;
  if (i < 0) i = 0;
  if (j < 0) j = 0;
  if (i > nx-1) i = nx-1;
  if (j > ny-1) j = ny-1;
  vec3 a(data[3*i + 3*nx*j]/255.0, data[3*i + 3*nx*j+1]/255.0, data[3*i + 3*nx*j+2]/255.0);
  i = b_u * nx;
  j = (1-b_v) * ny - 0.001;
  if (i < 0) i = 0;
  if (j < 0) j = 0;
  if (i > nx-1) i = nx-1;
  if (j > ny-1) j = ny-1;
  vec3 b(data[3*i + 3*nx*j]/255.0, data[3*i + 3*nx*j+1]/255.0, data[3*i + 3*nx*j+2]/255.0);
  i = c_u * nx;
  j = (1-c_v) * ny - 0.001;
  if (i < 0) i = 0;
  if (j < 0) j = 0;
  if (i > nx-1) i = nx-1;
  if (j > ny-1) j = ny-1;
  vec3 c(data[3*i + 3*nx*j]/255.0, data[3*i + 3*nx*j+1]/255.0, data[3*i + 3*nx*j+2]/255.0);
  return((1 - u - v) * a + u * b + v * c);
}

#endif