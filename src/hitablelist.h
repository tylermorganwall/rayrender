#ifndef HITABLELISTH
#define HITABLELISTH

#include "hitable.h"

class hitable_list: public hitable {
  public:
    hitable_list() {}
    hitable_list(hitable **l, int n) {list = l; list_size = n;}
    virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng);
    virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
    virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng);
    virtual vec3 random(const vec3& o, random_gen& rng);
    hitable **list;
    int list_size;
};

bool hitable_list::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  hit_record temp_rec;
  bool hit_anything = false;
  double closest_so_far = t_max;
  for(int i = 0; i < list_size; i++) {
    if(list[i]->hit(r, t_min, closest_so_far, temp_rec, rng)) {
      hit_anything = true;
      closest_so_far = temp_rec.t;
      rec = temp_rec;
    }
  }
  return(hit_anything);
}

bool hitable_list::bounding_box(Float t0, Float t1, aabb& box) const {
  if(list_size < 1) {
    return(false);
  }
  aabb temp_box;
  bool first_true = list[0]->bounding_box(t0,t1,temp_box);
  if(!first_true) {
    return(false);
  } else {
    box = temp_box;
  }
  for(int i = 1; i < list_size; i++) {
    if(list[0]->bounding_box(t0,t1, temp_box)) {
      box = surrounding_box(box, temp_box);
    } else {
      return(false);
    }
  }
  return(true);
}

Float hitable_list::pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
  Float weight = 1.0 / list_size;
  Float sum = 0;
  for (int i = 0; i < list_size; i++) {
    sum += weight*list[i]->pdf_value(o,v, rng);
  }
  return(sum);
}

vec3 hitable_list::random(const vec3& o, random_gen& rng) {
  int index = int(rng.unif_rand() * list_size * 0.99999999);
  return(list[index]->random(o, rng));
}

#endif
