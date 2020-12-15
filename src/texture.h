#ifndef TEXTUREH
#define TEXTUREH

#include "perlin.h"
#include "Rcpp.h"
#include "mathinline.h"
#include <memory>

class texture {
public: 
  virtual vec3 value(Float u, Float v, const vec3& p) const = 0;
  virtual ~texture() {};
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
  ~checker_texture() {
    // if(even) delete even;
    // if(odd) delete odd;
  }
  checker_texture(std::shared_ptr<texture> t0, std::shared_ptr<texture> t1, Float p) : even(t0), odd(t1), period(p) {}
  virtual vec3 value(Float u, Float v, const vec3& p) const {
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
  noise_texture(Float sc, vec3 c, vec3 c2, Float ph, Float inten) : 
    scale(sc), color(c), color2(c2), phase(ph), intensity(inten) {
    noise = new perlin();
  }
  ~noise_texture() {
    if(noise) delete noise;
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

class world_gradient_texture : public texture {
public:
  world_gradient_texture() {}
  world_gradient_texture(vec3 p1, vec3 p2, vec3 c1, vec3 c2, bool hsv2) : 
    point1(p1)  {
    gamma_color1 = hsv2 ? RGBtoHSV(c1) : c1;
    gamma_color2 = hsv2 ? RGBtoHSV(c2) : c2;
    dir = p2 - p1;
    inv_trans_length = 1.0/dir.squared_length();
    hsv = hsv2;
  }
  ~world_gradient_texture() {}
  virtual vec3 value(Float u, Float v, const vec3& p) const {
    vec3 offsetp = p - point1;
    Float mix = clamp(dot(offsetp, dir)*inv_trans_length, 0, 1);
    vec3 color = gamma_color1 * (1-mix) + mix * gamma_color2;
    return(hsv ? HSVtoRGB(color) : color);
  }
  vec3 point1;
  vec3 gamma_color1, gamma_color2;
  Float inv_trans_length;
  vec3 dir;
  bool hsv;
};

class image_texture : public texture {
public:
  image_texture() {}
  image_texture(Float *pixels, int A, int B, int nn, Float repeatu, Float repeatv, Float intensity) : 
    data(pixels), nx(A), ny(B), channels(nn), repeatu(repeatu), repeatv(repeatv), intensity(intensity) {}
  virtual vec3 value(Float u, Float v, const vec3& p) const;
  Float *data;
  int nx, ny, channels;
  Float repeatu, repeatv;
  Float intensity;
};


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
  ~triangle_image_texture() {}
  triangle_image_texture(Float *pixels, int A, int B, int nn,
                         Float tex_u_a, Float tex_v_a,
                         Float tex_u_b, Float tex_v_b,
                         Float tex_u_c, Float tex_v_c) : data(pixels),  nx(A), ny(B), channels(nn),
                         a_u(tex_u_a), a_v(tex_v_a), 
                         b_u(tex_u_b), b_v(tex_v_b), 
                         c_u(tex_u_c), c_v(tex_v_c) {}
  virtual vec3 value(Float u, Float v, const vec3& p) const;
  
  Float *data;
  int nx, ny, channels;
  Float a_u, a_v, b_u, b_v, c_u, c_v;
};

class gradient_texture : public texture {
public: 
  gradient_texture() {}
  gradient_texture(vec3 c1, vec3 c2, bool v, bool hsv2) : 
    aligned_v(v) {
    gamma_color1 = hsv2 ? RGBtoHSV(c1) : c1;
    gamma_color2 = hsv2 ? RGBtoHSV(c2) : c2;
    hsv = hsv2;
  }
  virtual vec3 value(Float u, Float v, const vec3& p) const {
    vec3 final_color = aligned_v ? gamma_color1 * (1-u) + u * gamma_color2 : gamma_color1 * (1-v) + v * gamma_color2;
    return(hsv ? HSVtoRGB(final_color) : final_color);
  }
  vec3 gamma_color1, gamma_color2;
  bool aligned_v;
  bool hsv;
};

class alpha_texture {
public:
  alpha_texture() {}
  alpha_texture(Float *pixels, int A, int B, int nn) : data(pixels), nx(A), ny(B), channels(nn) {
    u_vec = vec3(0,1,0);
    v_vec = vec3(0,0,1);
  }
  alpha_texture(Float *pixels, int A, int B, int nn, vec3 u, vec3 v) : 
                data(pixels), nx(A), ny(B), channels(nn), u_vec(u), v_vec(v) {}
  vec3 value(Float u, Float v, const vec3& p) const;
  Float channel_value(Float u, Float v, const vec3& p) const;
  Float *data;
  int nx, ny, channels;
  vec3 u_vec, v_vec;
};


class bump_texture {
public:
  bump_texture() {}
  bump_texture(Float *pixels, int A, int B, int nn, Float intensity) : 
    data(pixels), nx(A), ny(B), channels(nn), intensity(intensity) { 
    u_vec = vec3(0,1,0);
    v_vec = vec3(0,0,1);
  }
  bump_texture(Float *pixels, int A, int B, int nn, vec3 u, vec3 v, Float intensity) : 
    data(pixels), nx(A), ny(B), channels(nn), u_vec(u), v_vec(v), intensity(intensity) {}
  vec3 value(Float u, Float v, const vec3& p) const;
  vec3 mesh_value(Float u, Float v, const vec3& p) const;
  Float *data;
  int nx, ny, channels;
  vec3 u_vec, v_vec;
  Float intensity;
};

#endif
