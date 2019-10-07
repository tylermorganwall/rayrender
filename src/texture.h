#ifndef TEXTUREH
#define TEXTUREH

#include "perlin.h"

class texture {
public: 
  virtual vec3 value(Float u, Float v, const vec3& p) const = 0;
};

class constant_texture : public texture {
public: 
  constant_texture() {}
  constant_texture(vec3 c) : color(c) {}
  virtual vec3 value(Float u, Float v, const vec3& p) const {
    return(color);
  }
  vec3 color;
};

class checker_texture : public texture {
public:
  checker_texture() {}
  checker_texture(texture *t0, texture *t1, Float p) : even(t0), odd(t1), period(p) {}
  virtual vec3 value(Float u, Float v, const vec3& p) const {
    Float invperiod = 1.0/period;
    Float sines  = sin(invperiod*p.x()*M_PI) * sin(invperiod*p.y()*M_PI) * sin(invperiod*p.z()*M_PI);
    if(sines < 0) {
      return(odd->value(u,v,p));
    } else {
      return(even->value(u,v,p));
    }
  }
  texture *even;
  texture *odd;
  Float period;
};

class noise_texture : public texture {
public:
  noise_texture() {}
  noise_texture(Float sc, vec3 c, vec3 c2, Float ph, Float inten) : 
    scale(sc), color(c), color2(c2), phase(ph), intensity(inten) {
  }
  ~noise_texture() {
    if(noise != 0) delete noise;
  }
  virtual vec3 value(Float u, Float v, const vec3& p) const {
    Float weight = 0.5*(1+sin(scale*p.y()  + intensity*noise->turb(scale * p) + phase));
    return(color * (1-weight) + color2 * weight);
  }
  perlin *noise;
  Float scale;
  vec3 color;
  vec3 color2;
  Float phase;
  Float intensity;
};

class image_texture : public texture {
public:
  image_texture() {}
  image_texture(unsigned char *pixels, int A, int B, int nn) : data(pixels), nx(A), ny(B), channels(nn) {}
  virtual vec3 value(Float u, Float v, const vec3& p) const;
  unsigned char *data;
  int nx, ny, channels;
};

vec3 image_texture::value(Float u, Float v, const vec3& p) const {
  int i = u * nx;
  int j = (1-v) * ny - 0.001;
  if (i < 0) i = 0;
  if (j < 0) j = 0;
  if (i > nx-1) i = nx-1;
  if (j > ny-1) j = ny-1;
  Float r = int(data[channels*i + channels*nx*j])/255.0;
  Float g = int(data[channels*i + channels*nx*j+1])/255.0;
  Float b = int(data[channels*i + channels*nx*j+2])/255.0;
  return(vec3(r,g,b));
}

class triangle_texture : public texture {
public:
  triangle_texture() {}
  triangle_texture(vec3 a, vec3 b, vec3 c) : a(a), b(b), c(c) {}
  virtual vec3 value(Float u, Float v, const vec3& p) const {
    return((1 - u - v) * a + u * b + v * c);
  };
  vec3 a,b,c;
};

class triangle_image_texture : public texture {
public:
  triangle_image_texture() {}
  triangle_image_texture(unsigned char *pixels, int A, int B, int nn,
                         Float tex_u_a, Float tex_v_a,
                         Float tex_u_b, Float tex_v_b,
                         Float tex_u_c, Float tex_v_c) : data(pixels),  nx(A), ny(B), channels(nn),
                         a_u(tex_u_a), a_v(tex_v_a), 
                         b_u(tex_u_b), b_v(tex_v_b), 
                         c_u(tex_u_c), c_v(tex_v_c) {}
  virtual vec3 value(Float u, Float v, const vec3& p) const;
  
  unsigned char *data;
  int nx, ny, channels;
  Float a_u, a_v, b_u, b_v,c_u, c_v;
};

vec3 triangle_image_texture::value(Float u, Float v, const vec3& p) const {
  Float uu = ((1 - u - v) * a_u + u * b_u + v * c_u);
  Float vv = ((1 - u - v) * a_v + u * b_v + v * c_v);
  int i = uu * nx;
  int j = (1-vv) * ny - 0.001;
  if (i < 0) i = 0;
  if (j < 0) j = 0;
  if (i > nx-1) i = nx-1;
  if (j > ny-1) j = ny-1;
  Float r = int(data[channels*i + channels*nx*j])/255.0;
  Float g = int(data[channels*i + channels*nx*j+1])/255.0;
  Float b = int(data[channels*i + channels*nx*j+2])/255.0;
  return(vec3(r,g,b));
}

#endif
