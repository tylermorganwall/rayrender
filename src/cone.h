#ifndef CONEH
#define CONEH

#include "hitable.h"
#include "material.h"
#include "mathinline.h"

class cone: public hitable {
public:
  cone() {}
  ~cone() {}
  cone(Float r, Float h, std::shared_ptr<material>  mat, 
       std::shared_ptr<alpha_texture> alpha_mask, std::shared_ptr<bump_texture> bump_tex,
       std::shared_ptr<Transform> ObjectToWorld, std::shared_ptr<Transform> WorldToObject, bool reverseOrientation) : 
    hitable(ObjectToWorld, WorldToObject, reverseOrientation), 
     radius(r), height(h), mat_ptr(mat), alpha_mask(alpha_mask), bump_tex(bump_tex) {};
  virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng);
  virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, Sampler* sampler);
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  virtual Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0);
  virtual Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0);
  
  virtual vec3f random(const point3f& o, random_gen& rng, Float time = 0);
  virtual vec3f random(const point3f& o, Sampler* sampler, Float time = 0);
  virtual std::string GetName() const {
    return(std::string("Cone"));
  }
  Float radius;
  Float height;
  //Float phi2; //unwrapped cone in a circle, angle around origin (radians)
  std::shared_ptr<material> mat_ptr;
  std::shared_ptr<alpha_texture> alpha_mask;
  std::shared_ptr<bump_texture> bump_tex;
};

#endif
