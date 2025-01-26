#include "instance.h"
#include "raylog.h"

instance::instance(hitable* scene, 
                   Transform* ObjectToWorld, 
                   Transform* WorldToObject,
                   hitable_list* imp_list) : 
  hitable(ObjectToWorld, WorldToObject, nullptr, false), 
  original_scene(scene), importance_sampled_objects(imp_list) {
}

const bool instance::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const {
  SCOPED_CONTEXT("MultiHit");
  SCOPED_TIMER_COUNTER("Instance");
  
  ray r2 = (*WorldToObject)(r);
  if(original_scene->hit(r2, t_min, t_max, rec, rng)) {
    rec = (*ObjectToWorld)(rec);
    if(rec.alpha_miss) {
      return(false);
    }
    return(true);
  }
  return(false);
}

const bool instance::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const {
  SCOPED_CONTEXT("MultiHit");
  SCOPED_TIMER_COUNTER("Instance");
  
  ray r2 = (*WorldToObject)(r);
  if(original_scene->hit(r2, t_min, t_max, rec, sampler)) {
    rec = (*ObjectToWorld)(rec);
    return(true);
  }
  return(false);
}

bool instance::HitP(const ray& r, Float t_min, Float t_max, random_gen& rng) const {
  SCOPED_CONTEXT("MultiHit");
  SCOPED_TIMER_COUNTER("Instance");
  
  ray r2 = (*WorldToObject)(r);
  return(original_scene->HitP(r2, t_min, t_max, rng));
}

bool instance::HitP(const ray& r, Float t_min, Float t_max, Sampler* sampler) const {
  SCOPED_CONTEXT("MultiHit");
  SCOPED_TIMER_COUNTER("Instance");
  
  ray r2 = (*WorldToObject)(r);
  return(original_scene->HitP(r2, t_min, t_max, sampler));
}

bool instance::bounding_box(Float t0, Float t1, aabb& box) const {
  bool bbox_bool = original_scene->bounding_box(t0,t1,box);
  box = (*ObjectToWorld)(box);
  return(bbox_bool);
}


Float instance::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) {
  ray r2 = (*WorldToObject)(ray(o,v));
  return(importance_sampled_objects->pdf_value(r2.origin(),r2.direction(), rng, time));
}

Float instance::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) {
  ray r2 = (*WorldToObject)(ray(o,v));
  return(importance_sampled_objects->pdf_value(r2.origin(),r2.direction(), sampler, time));
}

vec3f instance::random(const point3f& o, random_gen& rng, Float time) {
  return((*ObjectToWorld)(importance_sampled_objects->random((*WorldToObject)(o), rng, time)));
}

vec3f instance::random(const point3f& o, Sampler* sampler, Float time) {
  return((*ObjectToWorld)(importance_sampled_objects->random((*WorldToObject)(o), sampler, time)));
}

size_t instance::GetSize() {
  size_t total_size = original_scene->GetSize() + sizeof(*this);
  return(total_size);
}
