#ifndef RECTH
#define RECTH

#include "hitable.h"
#include "material.h"

class xy_rect : public hitable {
public:
  xy_rect() {}
  xy_rect(Float _x0, Float _x1, Float _y0, Float _y1, Float _k, 
          std::shared_ptr<material> mat, 
          std::shared_ptr<alpha_texture> alpha_mask, 
          std::shared_ptr<bump_texture> bump_tex, 
            Transform* ObjectToWorld, Transform* WorldToObject, bool reverseOrientation) : 
    hitable(ObjectToWorld, WorldToObject, mat, reverseOrientation), 
    x0(_x0), x1(_x1), y0(_y0), y1(_y1), k(_k), alpha_mask(alpha_mask),
    bump_tex(bump_tex) {};
  ~xy_rect() {}
  virtual const bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const;
  virtual const bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const;
  virtual bool HitP(const ray &r, Float t_min, Float t_max, random_gen& rng) const;
  virtual bool HitP(const ray &r, Float t_min, Float t_max, Sampler* sampler) const;

  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  virtual Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0);
  virtual Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0);
  
  virtual vec3f random(const point3f& o, random_gen& rng, Float time = 0);
  virtual vec3f random(const point3f& o, Sampler* sampler, Float time = 0);
  
  virtual std::string GetName() const {
    return(std::string("XY Rectangle"));
  }
  size_t GetSize() {
    return(sizeof(*this));
  }
  virtual void hitable_info_bounds(Float t0, Float t1) const {
    aabb box;
    bounding_box(t0, t1, box);
    Rcpp::Rcout << GetName() << ": " <<  box.min() << "-" << box.max() << "\n";
  }
  Float x0, x1, y0, y1, k;
  std::shared_ptr<alpha_texture> alpha_mask;
  std::shared_ptr<bump_texture> bump_tex;
};

class xz_rect : public hitable {
public:
  xz_rect() {}
  xz_rect(Float _x0, Float _x1, Float _z0, Float _z1, Float _k, 
          std::shared_ptr<material> mat, 
          std::shared_ptr<alpha_texture> alpha_mask, 
          std::shared_ptr<bump_texture> bump_tex, 
          Transform* ObjectToWorld, Transform* WorldToObject, bool reverseOrientation) : 
    hitable(ObjectToWorld, WorldToObject, mat, reverseOrientation), 
  x0(_x0), x1(_x1), z0(_z0), z1(_z1), k(_k), alpha_mask(alpha_mask), 
  bump_tex(bump_tex) {};
  ~xz_rect() {}
  virtual const bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const;
  virtual const bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const;
  virtual bool HitP(const ray &r, Float t_min, Float t_max, random_gen& rng) const;
  virtual bool HitP(const ray &r, Float t_min, Float t_max, Sampler* sampler) const;

  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  virtual Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0);
  virtual Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0);
  virtual vec3f random(const point3f& o, random_gen& rng, Float time = 0);
  virtual vec3f random(const point3f& o, Sampler* sampler, Float time = 0);
  virtual std::string GetName() const {
    return(std::string("XZ Rectangle"));
  }
  size_t GetSize() {
    return(sizeof(*this));
  }
  virtual void hitable_info_bounds(Float t0, Float t1) const {
    aabb box;
    bounding_box(t0, t1, box);
    Rcpp::Rcout << GetName() << ": " <<  box.min() << "-" << box.max() << "\n";
  }
  Float x0, x1, z0, z1, k;
  std::shared_ptr<alpha_texture> alpha_mask;
  std::shared_ptr<bump_texture> bump_tex;
};

class yz_rect : public hitable {
public:
  yz_rect() {}
  yz_rect(Float _y0, Float _y1, Float _z0, Float _z1, Float _k, 
          std::shared_ptr<material> mat, 
          std::shared_ptr<alpha_texture> alpha_mask, std::shared_ptr<bump_texture> bump_tex, 
          Transform* ObjectToWorld, Transform* WorldToObject, bool reverseOrientation) : 
    hitable(ObjectToWorld, WorldToObject, mat, reverseOrientation), 
  y0(_y0), y1(_y1), z0(_z0), z1(_z1), k(_k), alpha_mask(alpha_mask), 
  bump_tex(bump_tex) {};
  ~yz_rect() {}
  virtual const bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const;
  virtual const bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const;
  virtual bool HitP(const ray &r, Float t_min, Float t_max, random_gen& rng) const;
  virtual bool HitP(const ray &r, Float t_min, Float t_max, Sampler* sampler) const;

  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  virtual Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0);
  virtual Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0);
  virtual vec3f random(const point3f& o, random_gen& rng, Float time = 0);
  virtual vec3f random(const point3f& o, Sampler* sampler, Float time = 0);
  virtual std::string GetName() const {
    return(std::string("YZ Rectangle"));
  }
  size_t GetSize() {
    return(sizeof(*this));
  }
  virtual void hitable_info_bounds(Float t0, Float t1) const {
    aabb box;
    bounding_box(t0, t1, box);
    Rcpp::Rcout << GetName() << ": " <<  box.min() << "-" << box.max() << "\n";
  }
  Float y0, y1, z0, z1, k;
  std::shared_ptr<alpha_texture> alpha_mask;
  std::shared_ptr<bump_texture> bump_tex;
};

#endif
