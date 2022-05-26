#ifndef TRIANGLEH
#define TRIANGLEH

#include "hitable.h"
#include "material.h"
#include "onbh.h"

class triangle : public hitable {
public:
  triangle() {}
  triangle(point3f _a, point3f _b, point3f _c, bool _single, std::shared_ptr<material> mat, 
           std::shared_ptr<alpha_texture> alpha_mask, std::shared_ptr<bump_texture> bump_tex,
             std::shared_ptr<Transform> ObjectToWorld, std::shared_ptr<Transform> WorldToObject, bool reverseOrientation) :
    hitable(ObjectToWorld, WorldToObject, reverseOrientation), 
  a((*ObjectToWorld)(_a)), b((*ObjectToWorld)(_b)), c((*ObjectToWorld)(_c)), 
  single(_single), mp(mat), alpha_mask(alpha_mask), bump_tex(bump_tex) {
    edge1 = b-a;
    edge2 = c-a;
    normal = cross(edge1, edge2);
    area = normal.length()/2;
    normal.make_unit_vector();
    normals_provided = false;
    uv_provided = false;
  };
  triangle(point3f _a, point3f _b, point3f _c, normal3f _na, normal3f _nb, normal3f _nc, bool _single, 
           std::shared_ptr<material> mat, std::shared_ptr<alpha_texture> alpha_mask, std::shared_ptr<bump_texture> bump_tex,
           std::shared_ptr<Transform> ObjectToWorld, std::shared_ptr<Transform> WorldToObject, bool reverseOrientation) :
    hitable(ObjectToWorld, WorldToObject, reverseOrientation), 
    a((*ObjectToWorld)(_a)), b((*ObjectToWorld)(_b)), c((*ObjectToWorld)(_c)), 
    na((*ObjectToWorld)(_na)), nb((*ObjectToWorld)(_nb)), nc((*ObjectToWorld)(_nc)), single(_single), mp(mat), 
    alpha_mask(alpha_mask), bump_tex(bump_tex) {
    edge1 = b-a;
    edge2 = c-a;
    normal = cross(edge1, edge2);
    area = normal.length()/2;
    normals_provided = true;
    uv_provided = false;
  };
  triangle(point3f _a, point3f _b, point3f _c, 
           point2f uva, point2f uvb, point2f uvc,
           bool _single, 
           std::shared_ptr<material> mat, std::shared_ptr<alpha_texture> alpha_mask, std::shared_ptr<bump_texture> bump_tex,
           std::shared_ptr<Transform> ObjectToWorld, std::shared_ptr<Transform> WorldToObject, bool reverseOrientation) :
    hitable(ObjectToWorld, WorldToObject, reverseOrientation), 
    a((*ObjectToWorld)(_a)), b((*ObjectToWorld)(_b)), c((*ObjectToWorld)(_c)), 
    uv_a(uva), uv_b(uvb), uv_c(uvc), 
    single(_single), mp(mat), 
    alpha_mask(alpha_mask), bump_tex(bump_tex) {
    edge1 = b-a;
    edge2 = c-a;
    normal = cross(edge1, edge2);
    area = normal.length()/2;
    normals_provided = false;
    uv_provided = true;
  };
  triangle(point3f _a, point3f _b, point3f _c, 
           normal3f _na, normal3f _nb, normal3f _nc, 
           point2f uva, point2f uvb, point2f uvc,
           bool _single, 
           std::shared_ptr<material> mat, std::shared_ptr<alpha_texture> alpha_mask, std::shared_ptr<bump_texture> bump_tex,
           std::shared_ptr<Transform> ObjectToWorld, std::shared_ptr<Transform> WorldToObject, bool reverseOrientation) :
    hitable(ObjectToWorld, WorldToObject, reverseOrientation), 
    a((*ObjectToWorld)(_a)), b((*ObjectToWorld)(_b)), c((*ObjectToWorld)(_c)), 
    uv_a(uva), uv_b(uvb), uv_c(uvc), 
    na((*ObjectToWorld)(_na)), nb((*ObjectToWorld)(_nb)), nc((*ObjectToWorld)(_nc)), single(_single), mp(mat), 
    alpha_mask(alpha_mask), bump_tex(bump_tex) {
    edge1 = b-a;
    edge2 = c-a;
    normal = cross(edge1, edge2);
    area = normal.length()/2;
    normals_provided = true;
    uv_provided = true;
  };
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler);
  
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  virtual Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0);
  virtual Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0);
  
  virtual vec3f random(const point3f& origin, random_gen& rng, Float time = 0);
  virtual vec3f random(const point3f& origin, Sampler* sampler, Float time = 0);
  virtual std::string GetName() const {
    return(std::string("Triangle"));
  }
  vec3f normal;
  vec3f a, b, c;
  point2f uv_a, uv_b, uv_c;
  normal3f na, nb, nc;
  vec3f edge1, edge2;
  Float area;
  bool normals_provided;
  bool uv_provided;
  bool single;
  std::shared_ptr<material> mp;
  std::shared_ptr<alpha_texture> alpha_mask;
  std::shared_ptr<bump_texture> bump_tex;
};

#endif
