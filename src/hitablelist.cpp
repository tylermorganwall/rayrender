#include "hitablelist.h"


bool hitable_list::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  hit_record temp_rec;
#ifdef DEBUGBVH
  temp_rec.bvh_nodes = rec.bvh_nodes;
#endif
  bool hit_anything = false;
  double closest_so_far = t_max;
  for (const auto& object : objects) {
    if (object->hit(r, t_min, closest_so_far, temp_rec, rng)) {
      hit_anything = true;
      closest_so_far = temp_rec.t;
      rec = temp_rec;
    }
  }
  return(hit_anything);
}

bool hitable_list::bounding_box(Float t0, Float t1, aabb& box) const {
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

Float hitable_list::pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
  Float weight = 1.0 / objects.size();
  Float sum = 0;
  for (const auto& object : objects) {
    sum += weight*object->pdf_value(o,v, rng);
  }
  return(sum);
}

vec3 hitable_list::random(const vec3& o, random_gen& rng) {
  int index = int(rng.unif_rand() * objects.size() * 0.99999999);
  return(objects[index]->random(o, rng));
}

vec3 hitable_list::random(const vec3& o, Sampler* sampler) {
  int index = int(sampler->Get1D() * objects.size() * 0.99999999);
  return(objects[index]->random(o, sampler));
}