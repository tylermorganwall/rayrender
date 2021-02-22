#ifndef ELLIPSOIDH
#define ELLIPSOIDH

#include "sphere.h"
#include "material.h"
#include "mathinline.h"

class ellipsoid: public hitable {
  public:
    ellipsoid() {}
    ellipsoid(vec3 cen, Float r, vec3 axes, std::shared_ptr<material> mat, 
              std::shared_ptr<alpha_texture> alpha_mask, std::shared_ptr<bump_texture> bump_tex) : 
      center(cen), radius(r), axes(axes),  mat_ptr(mat), alpha_mask(alpha_mask), bump_tex(bump_tex) {
      inv_axes = vec3(1.0f/axes.x(), 1.0f/axes.y(), 1.0f/axes.z());
      largest_proj_axis = axes.x() * axes.y() * axes.z() / ffmin(axes.x(), ffmin(axes.y(), axes.z()));
    };
    virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng);
    virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, Sampler* sampler);
    
    virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
    virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng, Float time = 0);
    virtual Float pdf_value(const vec3& o, const vec3& v, Sampler* sampler, Float time = 0);
    
    virtual vec3 random(const vec3& o, random_gen& rng, Float time = 0);
    virtual vec3 random(const vec3& o, Sampler* sampler, Float time = 0);
    
    vec3 center;
    Float radius;
    vec3 axes;
    vec3 inv_axes;
    Float largest_proj_axis;
    std::shared_ptr<material> mat_ptr;
    std::shared_ptr<alpha_texture> alpha_mask;
    std::shared_ptr<bump_texture> bump_tex;
};


#endif
