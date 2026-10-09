#ifndef TEXTUREH
#define TEXTUREH

#include "../math/perlin.h"
#include "Rcpp.h"
#include "../math/point3.h"
#include "../math/mathinline.h"
#include <memory>

struct hit_record;
struct TextureFootprint {
  Float dudx=0, dvdx=0, dudy=0, dvdy=0;
  bool valid=false;
};
struct HeightLevel {
  int width=0, height=0;
  std::vector<Float> pixels;
};
struct HeightImage {
  std::vector<HeightLevel> levels;
  HeightImage(int width, int height, std::vector<Float> pixels);
};
struct TextureEvalContext {
  point3f p{0}, object_p{0};
  normal3f geometric_normal{0, 1, 0}, object_normal{0, 1, 0};
  Float u = 0, v = 0, time = 0;
  vec3f dpdx{0}, dpdy{0};
  TextureFootprint footprint;
  bool has_derivatives = false;
  int face_index = -1; // Ptex source face, not the triangulated primitive index.
  static TextureEvalContext FromHit(const hit_record& hit);
};

class texture {
public:
  virtual point3f value(Float u, Float v, const point3f& p) const = 0;
  virtual point3f value(const TextureEvalContext& context) const {
    return value(context.u, context.v, context.p);
  }
  virtual point3f value(const hit_record& hit) const;
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
    Float sinx  = sin(invperiod*p.xyz.x*M_PI);
    sinx = sinx == 0 ? 1 : sinx;
    Float siny  = sin(invperiod*p.xyz.y*M_PI);
    siny = siny == 0 ? 1 : siny;
    Float sinz  = sin(invperiod*p.xyz.z*M_PI);
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
    Float weight = 0.5*(1+sin(scale*p.xyz.y  + intensity*noise->turb(scale * p) + phase));
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
                Float repeatu = 1.0f, Float repeatv = 1.0f, Float intensity = 1.0f, Float offsetu = 0.f, Float offsetv = 0.f) :
    data(pixels), nx(A), ny(B), channels(nn), repeatu(repeatu), repeatv(repeatv), offsetu(offsetu), offsetv(offsetv), intensity(intensity) {}
  virtual point3f value(Float u, Float v, const point3f& p) const;
  Float *data;
  int nx, ny, channels;
  Float repeatu, repeatv;
  Float offsetu = 0, offsetv = 0;
  Float intensity;
};

// Latitude-longitude images have periodic longitude and clamped latitude.
// Pixel centers, rather than endpoint texels, must partition the full sphere.
class latlong_image_texture final : public image_texture_float {
public:
  latlong_image_texture(Float *pixels, int width, int height, int channels, Float intensity = 1)
      : image_texture_float(pixels, width, height, channels, 1, 1, intensity) {}
  point3f value(Float u, Float v, const point3f& p) const override;
};

class image_texture_char : public texture {
public:
  image_texture_char() {}
  image_texture_char(unsigned char * pixels, int A, int B, int nn,
                Float repeatu = 1.0f, Float repeatv = 1.0f, Float intensity = 1.0f, Float offsetu = 0.f, Float offsetv = 0.f) :
    data(pixels), nx(A), ny(B), channels(nn), repeatu(repeatu), repeatv(repeatv), offsetu(offsetu), offsetv(offsetv), intensity(intensity) {}
  virtual point3f value(Float u, Float v, const point3f& p) const;
  unsigned char * data;
  int nx, ny, channels;
  Float repeatu, repeatv;
  Float offsetu = 0, offsetv = 0;
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
  explicit alpha_texture(Float opacity) : data(nullptr), nx(0), ny(0), channels(0), opacity(opacity) {}
  alpha_texture(unsigned char *pixels, int A, int B, int nn, Float offsetu = 0.f, Float offsetv = 0.f,
                Float repeatu = 1.f, Float repeatv = 1.f) :
    offsetu(offsetu), offsetv(offsetv), repeatu(repeatu), repeatv(repeatv), data(pixels), nx(A), ny(B), channels(nn) {}
  Float offsetu = 0, offsetv = 0;
  Float repeatu = 1, repeatv = 1;
  Float value(Float u, Float v, const point3f& p) const;
  unsigned char *data;
  int nx, ny, channels;
  Float opacity = 1;
};


class bump_texture {
public:
  bump_texture() {}
  bump_texture(std::shared_ptr<const texture> height_texture, Float intensity) :
    height_texture(std::move(height_texture)), nx(0), ny(0), channels(0), intensity(intensity),
    repeatu(1), repeatv(1) {}
  bump_texture(unsigned char *pixels, int A, int B, int nn, Float intensity,
               Float repeatu = 1.f, Float repeatv = 1.f, Float offsetu = 0.f, Float offsetv = 0.f) :
    data(pixels), nx(A), ny(B), channels(nn), intensity(intensity),
    repeatu(repeatu), repeatv(repeatv), offsetu(offsetu), offsetv(offsetv) {}
  bump_texture(std::shared_ptr<const HeightImage> pixels, int A, int B, int nn, Float intensity,
               Float repeatu = 1.f, Float repeatv = 1.f, Float offsetu = 0.f, Float offsetv = 0.f) :
    image(std::move(pixels)), nx(A), ny(B), channels(nn), intensity(intensity),
    repeatu(repeatu), repeatv(repeatv), offsetu(offsetu), offsetv(offsetv) {}
  point3f value(Float u, Float v, const point3f& p, TextureFootprint footprint = {}) const;
  Float raw_value(Float u, Float v, const point3f& p, TextureFootprint footprint = {}) const;
  // PBRT BumpMap with footprint-based steps and filtered heights. Geometry supplies derivatives
  // of its smooth normal; planar surfaces use zero. Tangents are updated too.
  normal3f perturb(Float u, Float v, const point3f& p, normal3f n,
                   vec3f& dpdu, vec3f& dpdv,
                   normal3f dndu = normal3f(0), normal3f dndv = normal3f(0),
                   TextureFootprint footprint = {}, const TextureEvalContext* context = nullptr) const;

  unsigned char *data = nullptr;
  std::shared_ptr<const HeightImage> image;
  std::shared_ptr<const texture> height_texture;
  int nx, ny, channels;
  Float intensity;
  Float repeatu, repeatv;
  Float offsetu = 0, offsetv = 0;
};

class roughness_texture {
public:
  roughness_texture() {}
  roughness_texture(unsigned char *pixels, int A, int B, int nn, Float offsetu = 0.f, Float offsetv = 0.f) :
    offsetu(offsetu), offsetv(offsetv), data(pixels), nx(A), ny(B), channels(nn) {}
  Float offsetu = 0, offsetv = 0;
  point2f value(Float u, Float v) const;
  point2f raw_value(Float u, Float v) const;
  static Float RoughnessToAlpha(Float roughness);
  unsigned char *data;
  int nx, ny, channels;
  vec3f u_vec, v_vec;
};


#endif
