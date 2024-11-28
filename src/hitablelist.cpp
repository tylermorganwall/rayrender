#include "hitablelist.h"
#include "raylog.h"

const bool hitable_list::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const {
  SCOPED_CONTEXT("MultiHit");
  SCOPED_TIMER_COUNTER("Hitable List");
  hit_record temp_rec;
#ifdef DEBUGBVH
  temp_rec.bvh_nodes = rec.bvh_nodes;
#endif
  bool hit_anything = false;
  Float closest_so_far = t_max;
  for (const auto& object : objects) {
    if (object->hit(r, t_min, closest_so_far, temp_rec, rng)) {
      hit_anything = true;
      closest_so_far = temp_rec.t;
      rec = temp_rec;
    }
  }
  return(hit_anything);
}

const bool hitable_list::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const {
  SCOPED_CONTEXT("MultiHit");
  SCOPED_TIMER_COUNTER("Hitable List");
  hit_record temp_rec;
#ifdef DEBUGBVH
  temp_rec.bvh_nodes = rec.bvh_nodes;
#endif
  bool hit_anything = false;
  Float closest_so_far = t_max;
  for (const auto& object : objects) {
    if (object->hit(r, t_min, closest_so_far, temp_rec, sampler)) {
      hit_anything = true;
      closest_so_far = temp_rec.t;
      rec = temp_rec;
    }
  }
  return(hit_anything);
}

bool hitable_list::HitP(const ray& r, Float t_min, Float t_max, random_gen& rng) const {
  SCOPED_CONTEXT("MultiHit");
  SCOPED_TIMER_COUNTER("Hitable List");
  for (const auto& object : objects) {
    if (object->HitP(r, t_min, t_max, rng)) {
      return(true);
    }
  }
  return(false);
}

bool hitable_list::HitP(const ray& r, Float t_min, Float t_max, Sampler* sampler) const {
  SCOPED_CONTEXT("MultiHit");
  SCOPED_TIMER_COUNTER("Hitable List");
  for (const auto& object : objects) {
    if (object->HitP(r, t_min, t_max, sampler)) {
      return(true);
    }
  }
  return(false);
}

bool hitable_list::bounding_box(Float t0, Float t1, aabb& box) const {
  SCOPED_CONTEXT("Bounding Box");
  SCOPED_TIMER_COUNTER("Hitable List");
  if(objects.empty()) {
    return(false);
  }
  aabb temp_box;
  bool first_true = objects[0]->bounding_box(t0,t1,temp_box);
  if(!first_true) {
    return(false);
  } else {
    box = temp_box;
  }
  for (const auto& object : objects) {
    if(object->bounding_box(t0,t1, temp_box)) {
      box = surrounding_box(box, temp_box);
    } else {
      return(false);
    }
  }
  return(true);
}

Float hitable_list::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) {
  Float weight = 1.0 / objects.size();
  Float sum = 0;
  for (const auto& object : objects) {
    sum += weight*object->pdf_value(o,v, rng, time);
  }
  return(sum);
}

Float hitable_list::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) {
  Float weight = 1.0 / objects.size();
  Float sum = 0;
  for (const auto& object : objects) {
    sum += weight*object->pdf_value(o,v, sampler, time);
  }
  return(sum);
}

vec3f hitable_list::random(const point3f& o, random_gen& rng, Float time) {
  int index = int(rng.unif_rand() * objects.size() * 0.99999999);
  return(objects[index]->random(o, rng, time));
}

vec3f hitable_list::random(const point3f& o, Sampler* sampler, Float time) {
  int index = int(sampler->Get1D() * objects.size() * 0.99999999);
  return(objects[index]->random(o, sampler, time));
}

std::string hitable_list::GetName() const {
  return(std::string("Hitable List"));
}

size_t hitable_list::GetSize() {
  size_t total_size = sizeof(*this) + objects.size() * sizeof(std::shared_ptr<hitable>);
  for(size_t i = 0; i < objects.size(); i++) {
    total_size += objects[i]->GetSize();;
  }
  return(total_size);
}

void hitable_list::validate() const {
  // Objects Vector Validation
  for (const auto& object : objects) {
    if (!object) {
      throw std::runtime_error("Detected a null or uninitialized object in hitable_list.");
    }
  }
}
