#ifndef INFINITEAREALIGHTH
#define INFINITEAREALIGHTH

#include "hitable.h"
#include "distributions.h"
#include "mathinline.h"
#include "texture.h"
#include "onbh.h"
#include "vec3.h"

class InfiniteAreaLight: public hitable {
public:
  InfiniteAreaLight() {}
  ~InfiniteAreaLight() {
    delete distribution;
    // delete mat_ptr;
  }
  InfiniteAreaLight(int width, int height, Float r, vec3 center, 
                    std::shared_ptr<texture> image,  std::shared_ptr<material> mat);
  virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng);
  virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, Sampler* sampler);
  
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng, Float time = 0);
  virtual Float pdf_value(const vec3& o, const vec3& v, Sampler* sampler, Float time = 0);
  
  virtual vec3 random(const vec3& o, random_gen& rng, Float time = 0);
  virtual vec3 random(const vec3& o, Sampler* sampler, Float time = 0);
  
  int width, height;
  Float radius;
  vec3 center;
  std::shared_ptr<material> mat_ptr;
  Distribution2D *distribution;
};



#endif
