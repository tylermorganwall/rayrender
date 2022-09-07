#ifndef TEXTUREH
#define TEXTUREH

#include "perlin.h"
#include "Rcpp.h"
#include "point3.h"
#include "mathinline.h"
#include <memory>

class texture {
public: 
  virtual point3f value(Float u, Float v, const point3f& p) const = 0;
  virtual ~texture() {};
};

class constant_texture : public texture {
public: 
  constant_texture() {}
  constant_texture(point3f c) : color(c) {}
  virtual point3f value(Float u, Float v, const point3f& p) const {
    return(color);
  }
  point3f color;
};

class checker_texture : public texture {
public:
  checker_texture() {}
  ~checker_texture() {
    // if(even) delete even;
    // if(odd) delete odd;
  }
  checker_texture(std::shared_ptr<texture> t0, std::shared_ptr<texture> t1, Float p) : even(t0), odd(t1), period(p) {}
  virtual point3f value(Float u, Float v, const point3f& p) const {
    Float invperiod = 1.0/period;
    Float sinx  = sin(invperiod*p.x()*M_PI);
    sinx = sinx == 0 ? 1 : sinx;
    Float siny  = sin(invperiod*p.y()*M_PI);
    siny = siny == 0 ? 1 : siny;
    Float sinz  = sin(invperiod*p.z()*M_PI);
    sinz = sinz == 0 ? 1 : sinz;
    if(sinx * siny * sinz < 0) {
      return(odd->value(u,v,p));
    } else {
      return(even->value(u,v,p));
    }
  }
  std::shared_ptr<texture> even;
  std::shared_ptr<texture> odd;
  Float period;
};

class noise_texture : public texture {
public:
  noise_texture() {}
  noise_texture(Float sc, point3f c, point3f c2, Float ph, Float inten) : 
    scale(sc), color(c), color2(c2), phase(ph), intensity(inten) {
    noise = new perlin();
  }
  ~noise_texture() {
    if(noise) delete noise;
  }
  virtual point3f value(Float u, Float v, const point3f& p) const {
    Float weight = 0.5*(1+sin(scale*p.y()  + intensity*noise->turb(scale * vec3f(p.x(),p.y(),p.z())) + phase));
    return(color * (1-weight) + color2 * weight);
  }
  perlin *noise;
  Float scale;
  point3f color;
  point3f color2;
  Float phase;
  Float intensity;
};

class world_gradient_texture : public texture {
public:
  world_gradient_texture() {}
  world_gradient_texture(point3f p1, point3f p2, point3f c1, point3f c2, bool hsv2) : 
    point1(p1)  {
    gamma_color1 = hsv2 ? RGBtoHSV(c1) : c1;
    gamma_color2 = hsv2 ? RGBtoHSV(c2) : c2;
    dir = p2 - p1;
    inv_trans_length = 1.0/dir.squared_length();
    hsv = hsv2;
  }
  ~world_gradient_texture() {}
  virtual point3f value(Float u, Float v, const point3f& p) const {
    vec3f offsetp = p - point1;
    Float mix = clamp(dot(offsetp, dir)*inv_trans_length, 0, 1);
    point3f color = gamma_color1 * (1-mix) + mix * gamma_color2;
    return(hsv ? HSVtoRGB(color) : color);
  }
  point3f point1;
  point3f gamma_color1, gamma_color2;
  Float inv_trans_length;
  vec3f dir;
  bool hsv;
};


class image_texture_float : public texture {
public:
  image_texture_float() {}
  image_texture_float(Float *pixels, int A, int B, int nn, 
                Float repeatu = 1.0f, Float repeatv = 1.0f, Float intensity = 1.0f) : 
    data(pixels), nx(A), ny(B), channels(nn), repeatu(repeatu), repeatv(repeatv), intensity(intensity) {}
  virtual point3f value(Float u, Float v, const point3f& p) const;
  Float *data;
  int nx, ny, channels;
  Float repeatu, repeatv;
  Float intensity;
};

class image_texture_char : public texture {
public:
  image_texture_char() {}
  image_texture_char(unsigned char * pixels, int A, int B, int nn, 
                Float repeatu = 1.0f, Float repeatv = 1.0f, Float intensity = 1.0f) : 
    data(pixels), nx(A), ny(B), channels(nn), repeatu(repeatu), repeatv(repeatv), intensity(intensity) {}
  virtual point3f value(Float u, Float v, const point3f& p) const;
  unsigned char * data;
  int nx, ny, channels;
  Float repeatu, repeatv;
  Float intensity;
};


class triangle_texture : public texture {
public:
  triangle_texture() {}
  triangle_texture(point3f a, point3f b, point3f c) : a(a), b(b), c(c) {}
  virtual point3f value(Float u, Float v, const point3f& p) const;
  point3f a,b,c;
};


class gradient_texture : public texture {
public: 
  gradient_texture() {}
  gradient_texture(point3f c1, point3f c2, bool v, bool hsv2) : 
    aligned_v(v) {
    gamma_color1 = hsv2 ? RGBtoHSV(c1) : c1;
    gamma_color2 = hsv2 ? RGBtoHSV(c2) : c2;
    hsv = hsv2;
  }
  virtual point3f value(Float u, Float v, const point3f& p) const {
    point3f final_color = aligned_v ? gamma_color1 * (1-u) + u * gamma_color2 : gamma_color1 * (1-v) + v * gamma_color2;
    return(hsv ? HSVtoRGB(final_color) : final_color);
  }
  point3f gamma_color1, gamma_color2;
  bool aligned_v;
  bool hsv;
};

class alpha_texture {
public:
  alpha_texture() {}
  alpha_texture(unsigned char *pixels, int A, int B, int nn) : 
    data(pixels), nx(A), ny(B), channels(nn) {}
  Float value(Float u, Float v, const point3f& p) const;
  unsigned char *data;
  int nx, ny, channels;
};


class bump_texture {
public:
  bump_texture() {}
  bump_texture(unsigned char *pixels, int A, int B, int nn, Float intensity, 
               Float repeatu = 1.f, Float repeatv = 1.f) : 
    data(pixels), nx(A), ny(B), channels(nn), intensity(intensity), 
    repeatu(repeatu), repeatv(repeatv) {}
  point3f value(Float u, Float v, const point3f& p) const;
  Float raw_value(Float u, Float v, const point3f& p) const;
  
  unsigned char *data;
  int nx, ny, channels;
  Float intensity;
  Float repeatu, repeatv;
};

class roughness_texture {
public:
  roughness_texture() {}
  roughness_texture(unsigned char *pixels, int A, int B, int nn) : 
    data(pixels), nx(A), ny(B), channels(nn) {}
  point2f value(Float u, Float v) const;
  static Float RoughnessToAlpha(Float roughness);
  unsigned char *data;
  int nx, ny, channels;
  vec3f u_vec, v_vec;
};


#endif
