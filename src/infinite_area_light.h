#ifndef INFINITEAREALIGHTH
#define INFINITEAREALIGHTH

#include "hitable.h"
#include "distributions.h"
#include "mathinline.h"
#include "texture.h"
#include "onbh.h"
#include "vec3.h"
#include "material.h"

class InfiniteAreaLight: public hitable {
public:
  InfiniteAreaLight() {}
  ~InfiniteAreaLight() {
    delete distribution;
  }
  InfiniteAreaLight(int width, int height, Float r, point3f center, 
                    std::shared_ptr<texture> image,  std::shared_ptr<material> mat,
                    Transform* ObjectToWorld, Transform* WorldToObject, bool reverseOrientation);
  virtual const bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng) const;
  virtual const bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, Sampler* sampler) const;
  virtual bool HitP(const ray &r, Float t_min, Float t_max, random_gen& rng) const;
  virtual bool HitP(const ray &r, Float t_min, Float t_max, Sampler* sampler) const;

  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  virtual Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0);
  virtual Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0);
  
  virtual vec3f random(const point3f& o, random_gen& rng, Float time = 0);
  virtual vec3f random(const point3f& o, Sampler* sampler, Float time = 0);
  virtual std::string GetName() const {
    return(std::string("EnvironmentLight"));
  }
  virtual void hitable_info_bounds(Float t0, Float t1) const {
    aabb box;
    bounding_box(t0, t1, box);
    Rcpp::Rcout << GetName() << ": " <<  box.min() << "-" << box.max() << "\n";
  }
  size_t GetSize();
  int width, height;
  Float radius;
  point3f center;
  Distribution2D *distribution;
};



#endif
