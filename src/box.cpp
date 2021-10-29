#include "box.h"

box::box(const vec3f& p0, const vec3f& p1, std::shared_ptr<material> ptr, 
         std::shared_ptr<alpha_texture> alpha_mask, std::shared_ptr<bump_texture> bump_tex,
         std::shared_ptr<Transform> ObjectToWorld, std::shared_ptr<Transform> WorldToObject, bool reverseOrientation) :
  hitable(ObjectToWorld, WorldToObject, reverseOrientation) {
  pmin = (p0);
  pmax = (p1);
  list.add(std::make_shared<xy_rect>(p0.x(), p1.x(), p0.y(), p1.y(), p1.z(), ptr, alpha_mask, bump_tex, 
                                     ObjectToWorld, WorldToObject, false));
  list.add(std::make_shared<xy_rect>(p0.x(), p1.x(), p0.y(), p1.y(), p0.z(), ptr, alpha_mask, bump_tex, 
                                     ObjectToWorld, WorldToObject, true));
  list.add(std::make_shared<xz_rect>(p0.x(), p1.x(), p0.z(), p1.z(), p1.y(), ptr, alpha_mask, bump_tex, 
                                     ObjectToWorld, WorldToObject, false));
  list.add(std::make_shared<xz_rect>(p0.x(), p1.x(), p0.z(), p1.z(), p0.y(), ptr, alpha_mask, bump_tex, 
                                     ObjectToWorld, WorldToObject, true));
  list.add(std::make_shared<yz_rect>(p0.y(), p1.y(), p0.z(), p1.z(), p1.x(), ptr, alpha_mask, bump_tex, 
                                     ObjectToWorld, WorldToObject, false));
  list.add(std::make_shared<yz_rect>(p0.y(), p1.y(), p0.z(), p1.z(), p0.x(), ptr, alpha_mask, bump_tex, 
                                     ObjectToWorld, WorldToObject, true));
}

bool box::bounding_box(Float t0, Float t1, aabb& box) const {
  box = (*ObjectToWorld)(aabb(pmin, pmax));
  return(true);
}
Float box::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) {
  return(list.pdf_value(o,v, rng, time));
}
Float box::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time ) {
  return(list.pdf_value(o,v, sampler, time));
}
vec3f box::random(const point3f& o, random_gen& rng, Float time) {
  return(list.random(o, rng, time));
}
vec3f box::random(const point3f& o, Sampler* sampler, Float time) {
  return(list.random(o, sampler, time));
}
std::string box::GetName() const {
  return(std::string("Box"));
}

bool box::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  return(list.hit(r,t_min,t_max,rec, rng));
}

bool box::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
  return(list.hit(r,t_min,t_max,rec, sampler));
}
