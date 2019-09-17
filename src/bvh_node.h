#ifndef BVHNODEH
#define BVHNODEH

#include "hitable.h"
#include "aabb.h"
#include <Rcpp.h>


class bvh_node : public hitable {
  public:
    bvh_node() {}
    bvh_node(hitable **l, int n, Float time0, Float time1, random_gen rng);
    virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
    virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
    hitable *left;
    hitable *right;
    aabb box;
};

bool bvh_node::bounding_box(Float t0, Float t1, aabb& b) const {
  b = box;
  return(true);
}

bool bvh_node::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  if(box.hit(r, t_min, t_max, rng)) {
    rec.bvh_nodes++;
    if(left->hit(r,t_min,t_max,rec, rng)) {
      right->hit(r,t_min,rec.t,rec, rng);
      return(true);
    } else {
      return(right->hit(r,t_min,t_max,rec, rng));
    }
  }
  return(false);
}
int box_x_compare(const void * a, const void * b) {
  aabb box_left, box_right;
  hitable *ah = *(hitable**)a;
  hitable *bh = *(hitable**)b;
  if(!ah->bounding_box(0,0,box_left) || !bh->bounding_box(0,0,box_right)) {
  }
  if(box_left.centroid.x() - box_right.centroid.x() < 0.0) {
    return(-1);
  } else {
    return(1);
  }
}

int box_y_compare(const void * a, const void * b) {
  aabb box_left, box_right;
  hitable *ah = *(hitable**)a;
  hitable *bh = *(hitable**)b;
  if(!ah->bounding_box(0,0,box_left) || !bh->bounding_box(0,0,box_right)) {
  }
  if(box_left.centroid.y() - box_right.centroid.y() < 0.0) {
    return(-1);
  } else {
    return(1);
  }
}

int box_z_compare(const void * a, const void * b) {
  aabb box_left, box_right;
  hitable *ah = *(hitable**)a;
  hitable *bh = *(hitable**)b;
  if(!ah->bounding_box(0,0,box_left) || !bh->bounding_box(0,0,box_right)) {
  }
  if(box_left.centroid.z() - box_right.centroid.z() < 0.0) {
    return(-1);
  } else {
    return(1);
  }
}

bvh_node::bvh_node(hitable **l, int n, Float time0, Float time1, random_gen rng) {
  aabb centroid_bounds;
  l[0]->bounding_box(time0, time1, centroid_bounds);
  for (int i = 0; i < n; ++i) {
    aabb tempbox;
    if(l[i]->bounding_box(time0,time1,tempbox)) {
      centroid_bounds = surrounding_box(centroid_bounds, tempbox);
    }
  }
  vec3 centroid_bounds_values = centroid_bounds.max() - centroid_bounds.min();
  int axis = centroid_bounds_values.x() > centroid_bounds_values.y() ? 0 : 1;
  if(axis == 0) {
    axis = centroid_bounds_values.x() > centroid_bounds_values.z() ? 0 : 2;
  } else {
    axis = centroid_bounds_values.y() > centroid_bounds_values.z() ? 1 : 2;
  }
  if(axis == 0) {
    std::qsort(l, n, sizeof(hitable *), box_x_compare);
  } else if (axis == 1) {
    std::qsort(l, n, sizeof(hitable *), box_y_compare);
  } else {
    std::qsort(l, n, sizeof(hitable *), box_z_compare);
  }
  if(n == 1) {
    left = right = l[0];
  } else if (n == 2) {
    left = l[0];
    right = l[1];
  } else {
    left = new bvh_node(l, n/2, time0, time1, rng);
    right = new bvh_node(l + n/2, n - n/2, time0, time1, rng);
  }
  aabb box_left, box_right;
  if(!left->bounding_box(time0,time1,box_left) || !right->bounding_box(time0,time1,box_right)) {
  }
  box = surrounding_box(box_left,box_right);
}

#endif
