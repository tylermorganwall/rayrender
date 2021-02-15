#ifndef SPHEREH
#define SPHEREH

#include "hitable.h"
#include "material.h"
#include "mathinline.h"

class sphere: public hitable {
  public:
    sphere() {}
    ~sphere() {}
    sphere(vec3 cen, Float r, std::shared_ptr<material> mat, 
           std::shared_ptr<alpha_texture> alpha_mask, std::shared_ptr<bump_texture> bump_tex) : center(cen), radius(r), 
           mat_ptr(mat), alpha_mask(alpha_mask), bump_tex(bump_tex) {};
    virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng);
    virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, Sampler* sampler);
    
    virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
    virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng, Float time = 0);
    virtual Float pdf_value(const vec3& o, const vec3& v, Sampler* sampler, Float time = 0);
    
    virtual vec3 random(const vec3& o, random_gen& rng, Float time = 0);
    virtual vec3 random(const vec3& o, Sampler* sampler, Float time = 0);
    vec3 center;
    Float radius;
    std::shared_ptr<material> mat_ptr;
    std::shared_ptr<alpha_texture> alpha_mask;
    std::shared_ptr<bump_texture> bump_tex;
};

class moving_sphere: public hitable {
  public:
    moving_sphere() {}
    moving_sphere(vec3 cen0, vec3 cen1, Float t0, Float t1, Float r, 
                  std::shared_ptr<material> mat,
                  std::shared_ptr<alpha_texture> alpha_mask, std::shared_ptr<bump_texture> bump_tex) : 
                  center0(cen0), center1(cen1), time0(t0), time1(t1), radius(r), 
                  mat_ptr(mat), alpha_mask(alpha_mask), bump_tex(bump_tex) {};
    virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng);
    virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, Sampler* sampler);
    
    virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
    virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng, Float time = 0);
    virtual Float pdf_value(const vec3& o, const vec3& v, Sampler* sampler, Float time = 0);
    
    virtual vec3 random(const vec3& o, random_gen& rng, Float time = 0);
    virtual vec3 random(const vec3& o, Sampler* sampler, Float time = 0);
    vec3 center(Float time) const;
    vec3 center0, center1;
    Float time0, time1;
    Float radius;
    std::shared_ptr<material> mat_ptr;
    std::shared_ptr<alpha_texture> alpha_mask;
    std::shared_ptr<bump_texture> bump_tex;
};

#endif
