#ifndef BOXH
#define BOXH

#include "hitablelist.h"
#include "rectangle.h"
#include <memory>

class box : public hitable {
public:
  box() {}
  box(const vec3f& p0, const vec3f& p1, std::shared_ptr<material> ptr, 
      std::shared_ptr<alpha_texture> alpha_mask, std::shared_ptr<bump_texture> bump_tex,
      Transform* ObjectToWorld, Transform* WorldToObject, bool reverseOrientation);
  virtual const bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const;
  virtual const bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const;
  
  virtual bool HitP(const ray &r, Float t_min, Float t_max, random_gen& rng) const;
  virtual bool HitP(const ray &r, Float t_min, Float t_max, Sampler* sampler) const;

  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  virtual Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0);
  virtual Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0);
  virtual vec3f random(const point3f& o, random_gen& rng, Float time = 0);
  virtual vec3f random(const point3f& o, Sampler* sampler, Float time = 0);
  virtual std::string GetName() const;
  size_t GetSize()  {
    return(sizeof(*this));
  }
  virtual void hitable_info_bounds(Float t0, Float t1) const {
    aabb box;
    bounding_box(t0, t1, box);
    Rcpp::Rcout << GetName() << ": " <<  box.min() << "-" << box.max() << "\n";
  }
  vec3f pmin, pmax;
  hitable_list list;
};


#endif
