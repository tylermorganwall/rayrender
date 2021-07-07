#ifndef TRIANGLEH
#define TRIANGLEH

#include "hitable.h"
#include "material.h"
#include "onbh.h"

class triangle : public hitable {
public:
  triangle() {}
  triangle(vec3f _a, vec3f _b, vec3f _c, bool _single, std::shared_ptr<material> mat, 
           std::shared_ptr<alpha_texture> alpha_mask, std::shared_ptr<bump_texture> bump_tex) :
  a(_a), b(_b), c(_c), single(_single), mp(mat), alpha_mask(alpha_mask), bump_tex(bump_tex) {
    edge1 = b-a;
    edge2 = c-a;
    normal = cross(edge1, edge2);
    area = normal.length()/2;
    normal.make_unit_vector();
    normals_provided = false;
  };
  triangle(vec3f _a, vec3f _b, vec3f _c, vec3f _na, vec3f _nb, vec3f _nc, bool _single, 
           std::shared_ptr<material> mat, std::shared_ptr<alpha_texture> alpha_mask, std::shared_ptr<bump_texture> bump_tex) :
    a(_a), b(_b), c(_c), na(_na), nb(_nb), nc(_nc), single(_single), mp(mat), 
    alpha_mask(alpha_mask), bump_tex(bump_tex) {
    edge1 = b-a;
    edge2 = c-a;
    normal = cross(edge1, edge2);
    area = normal.length()/2;
    normals_provided = true;
  };
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler);
  
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  virtual Float pdf_value(const vec3f& o, const vec3f& v, random_gen& rng, Float time = 0);
  virtual Float pdf_value(const vec3f& o, const vec3f& v, Sampler* sampler, Float time = 0);
  
  virtual vec3f random(const vec3f& origin, random_gen& rng, Float time = 0);
  virtual vec3f random(const vec3f& origin, Sampler* sampler, Float time = 0);
  
  vec3f normal;
  vec3f a, b, c, na, nb, nc;
  vec3f edge1, edge2;
  Float area;
  bool normals_provided;
  bool single;
  std::shared_ptr<material> mp;
  std::shared_ptr<alpha_texture> alpha_mask;
  std::shared_ptr<bump_texture> bump_tex;
};

#endif
