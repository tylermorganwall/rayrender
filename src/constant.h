#ifndef CONSTANTH
#define CONSTANTH

#include "hitable.h"
#include "material.h"
#include <float.h>
#include <memory>

class material;

class constant_medium : public hitable {
public:
  constant_medium(std::shared_ptr<hitable> b, Float d, std::shared_ptr<texture> a ) : boundary(b), density(d) {
    phase_function = std::make_shared<isotropic>(a);
  }
  virtual const bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const;
  virtual const bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const;
  // virtual bool HitP(const ray &r, Float t_min, Float t_max, random_gen& rng) const;
  // virtual bool HitP(const ray &r, Float t_min, Float t_max, Sampler* sampler) const;

  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0);
  Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0);
  vec3f random(const point3f& o, random_gen& rng, Float time = 0);
  vec3f random(const point3f& o, Sampler* sampler, Float time = 0);
  std::string GetName() const {
    return(std::string("Constant"));
  }
  size_t GetSize()  {
    return(sizeof(*this) + boundary->GetSize());
  }
  virtual void hitable_info_bounds(Float t0, Float t1) const {
    aabb box;
    bounding_box(t0, t1, box);
    Rcpp::Rcout << GetName() << ": " <<  box.min() << "-" << box.max() << "\n";
  }
  
  std::shared_ptr<hitable> boundary;
  Float density;
  std::shared_ptr<material> phase_function;
};

#endif
