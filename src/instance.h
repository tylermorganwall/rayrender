#ifndef INSTANCEH
#define INSTANCEH


#include <Rcpp.h>
#include "rng.h"
#include "hitablelist.h"

class instance : public hitable {
public:
  instance() {}
  instance(hitable* scene, 
           std::shared_ptr<Transform> ObjectToWorld, 
           std::shared_ptr<Transform> WorldToObject,
           hitable_list* imp_list);
  
  virtual const bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const;
  virtual const bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const;
  Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0);
  Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0);
  virtual bool HitP(const ray &r, Float t_min, Float t_max, random_gen& rng) const;
  virtual bool HitP(const ray &r, Float t_min, Float t_max, Sampler* sampler) const;

  vec3f random(const point3f& o, random_gen& rng, Float time = 0);
  vec3f random(const point3f& o, Sampler* sampler, Float time = 0);
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  
  virtual std::string GetName() const {
    return(std::string("Instance"));
  }
  size_t GetSize();
  
  //Embedded scene
  hitable* original_scene;
  hitable_list* importance_sampled_objects;
};

#endif

