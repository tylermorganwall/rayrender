#ifndef DISKH
#define DISKH

#include "hitable.h"
#include "onbh.h"
#include "material.h"

class disk : public hitable {
public:
  disk() {}
  disk(vec3f cen, Float r, Float i_r, std::shared_ptr<material> mat, std::shared_ptr<alpha_texture> alpha_mask, 
       std::shared_ptr<bump_texture> bump_tex,
       Transform* ObjectToWorld, Transform* WorldToObject, bool reverseOrientation) : 
      hitable(ObjectToWorld, WorldToObject, mat, reverseOrientation), 
      center(cen), radius(r), 
       inner_radius(i_r), alpha_mask(alpha_mask), bump_tex(bump_tex) {};
  ~disk() {}
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
    return(std::string("Disk"));
  }
  size_t GetSize()  {
    return(sizeof(*this));
  }
  virtual void hitable_info_bounds(Float t0, Float t1) const {
    aabb box;
    bounding_box(t0, t1, box);
    Rcpp::Rcout << GetName() << ": " <<  box.min() << "-" << box.max() << "\n";
  }
  vec3f center;
  Float radius;
  Float inner_radius;
  std::shared_ptr<alpha_texture> alpha_mask;
  std::shared_ptr<bump_texture> bump_tex;
};


#endif
