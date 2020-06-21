#ifndef INFINITEAREALIGHTH
#define INFINITEAREALIGHTH

#include "hitable.h"
#include "distributions.h"
#include "mathinline.h"
#include "texture.h"
#include "onbh.h"
#include "vec3.h"

#ifndef FLOATDEF
#define FLOATDEF
#ifdef RAY_FLOAT_AS_DOUBLE
typedef double Float;
#else
typedef float Float;
#endif 
#endif

class InfiniteAreaLight: public hitable {
public:
  InfiniteAreaLight() {}
  ~InfiniteAreaLight() {
    delete distribution;
    delete mat_ptr;
  }
  InfiniteAreaLight(int width, int height, Float r, vec3 center, texture *image,  material *mat);
  virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng);
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng);
  virtual vec3 random(const vec3& o, random_gen& rng);
  virtual vec3 random(const vec3& o, Sampler* sampler);
  
  int width, height;
  Float radius;
  vec3 center;
  material *mat_ptr;
  Distribution2D *distribution;
};

InfiniteAreaLight::InfiniteAreaLight(int width, int height, Float r, vec3 center, 
                                     texture *image, material *mat)
  : width(width), height(height), radius(r), center(center), mat_ptr(mat) {
  //Set up distribution
  std::unique_ptr<Float[]>  img(new Float[width * height]);
  for (int v = 0; v < height; ++v) {
    Float vp = (Float)v / (Float)height;
    Float sinTheta = std::sin(M_PI * Float(v + .5f) / Float(height));
    for (int u = 0; u < width; ++u) {
      Float up = (Float)u / (Float)width;
      vec3 rgb = image->value(up, vp, center);
      img[u + v * width] =  0.212671f*rgb.r() + 0.715160f*rgb.g() + 0.072169f*rgb.b();
      img[u + v * width] *= sinTheta;
    }
  }
  distribution = new Distribution2D(img.get(), width, height);
}

bool InfiniteAreaLight::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  vec3 oc = r.origin() - center;
  Float a = dot(r.direction(), r.direction());
  Float b = 2 * dot(oc, r.direction()); 
  Float c = dot(oc,oc) - radius * radius;
  Float temp1, temp2;
  if (!quadratic(a, b, c, &temp1, &temp2)) {
    return(false);
  }
  if(temp1 < t_max && temp1 > t_min) {
    rec.t = temp1;
    rec.p = r.point_at_parameter(rec.t);
    rec.p *= radius / rec.p.length(); 
    rec.normal = (rec.p - center) / radius;
    vec3 v2(-r.direction().z(),r.direction().y(),r.direction().x());
    get_sphere_uv(unit_vector(v2), rec.u, rec.v);
    rec.u = 1 - rec.u;
    rec.mat_ptr = mat_ptr;
    return(true);
  }
  if(temp2 < t_max && temp2 > t_min) {
    rec.t = temp1;
    rec.p = r.point_at_parameter(rec.t);
    rec.p *= radius / rec.p.length(); 
    rec.normal = (rec.p - center) / radius;
    vec3 v2(-r.direction().z(),r.direction().y(),r.direction().x());
    get_sphere_uv(unit_vector(v2), rec.u, rec.v);
    rec.u = 1 - rec.u;
    rec.mat_ptr = mat_ptr;
    return(true);
  }
  return(false);
}

Float InfiniteAreaLight::pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, rng)) {
    vec3 v2(-v.z(),v.y(),v.x());
    get_sphere_uv(unit_vector(v2), rec.u, rec.v);
    rec.u = 1 - rec.u;
    Float sinTheta = std::sin(rec.v * M_PI);
    if (sinTheta == 0) {
      return(0);
    }
    //u = phi, v = theta
    return(distribution->Pdf(vec2(rec.u, rec.v)) /
      (2 * M_PI * M_PI * sinTheta));
  } else {
    return(0);
  }
}


vec3 InfiniteAreaLight::random(const vec3& o, random_gen& rng) {
  vec2 u(rng.unif_rand(), rng.unif_rand());
  Float mapPdf;
  vec2 uv = distribution->SampleContinuous(u, &mapPdf);
  if (mapPdf == 0) {
    return(vec3(0.f,0.f,0.f));
  }
  //theta vertical, phi horizontal
  Float theta = (1-uv[1]) * M_PI, phi = (1-uv[0]) * 2.0f * M_PI;
  Float cosTheta = std::cos(theta), sinTheta = std::sin(theta);
  Float sinPhi = std::sin(phi), cosPhi = std::cos(phi);
  vec3 d(sinTheta * sinPhi, cosTheta, sinTheta * cosPhi);
  return(d);
}

vec3 InfiniteAreaLight::random(const vec3& o, Sampler* sampler) {
  vec2 u = sampler->Get2D();
  Float mapPdf;
  vec2 uv = distribution->SampleContinuous(u, &mapPdf);
  if (mapPdf == 0) {
    return(vec3(0.f,0.f,0.f));
  }
  //theta vertical, phi horizontal
  Float theta = (1-uv[1]) * M_PI, phi = (1-uv[0]) * 2.0f * M_PI;
  Float cosTheta = std::cos(theta), sinTheta = std::sin(theta);
  Float sinPhi = std::sin(phi), cosPhi = std::cos(phi);
  vec3 d(sinTheta * sinPhi, cosTheta, sinTheta * cosPhi);
  return(d);
}


bool InfiniteAreaLight::bounding_box(Float t0, Float t1, aabb& box) const {
  box = aabb(center - vec3(radius,radius,radius), center + vec3(radius,radius,radius));
  return(true);
}


#endif
