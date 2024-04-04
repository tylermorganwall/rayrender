#include "instance.h"

instance::instance(hitable* scene, 
                   std::shared_ptr<Transform> ObjectToWorld, 
                   std::shared_ptr<Transform> WorldToObject) : 
  hitable(ObjectToWorld, WorldToObject, false), original_scene(scene) {
}

bool instance::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  ray r2 = (*WorldToObject)(r);
  if(original_scene->hit(r2, t_min, t_max, rec, rng)) {
    rec = (*ObjectToWorld)(rec);
    return(true);
  }
  return(original_scene->hit(r2, t_min, t_max, rec, rng));
}

bool instance::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
  ray r2 = (*WorldToObject)(r);
  if(original_scene->hit(r2, t_min, t_max, rec, sampler)) {
    rec = (*ObjectToWorld)(rec);
    return(true);
  }
  return(false);
}

bool instance::bounding_box(Float t0, Float t1, aabb& box) const {
  bool bbox_bool = original_scene->bounding_box(t0,t1,box);
  box = (*ObjectToWorld)(box);
  return(bbox_bool);
}


Float instance::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) {
  return(original_scene->pdf_value(o,v, rng, time));
}

Float instance::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) {
  return(original_scene->pdf_value(o,v, sampler, time));
  
}

vec3f instance::random(const point3f& o, random_gen& rng, Float time) {
  return(original_scene->random(o, rng, time));
}

vec3f instance::random(const point3f& o, Sampler* sampler, Float time) {
  return(original_scene->random(o, sampler, time));
}

size_t instance::GetSize() {
  size_t total_size = original_scene->GetSize() + sizeof(*this);
  return(total_size);
}
