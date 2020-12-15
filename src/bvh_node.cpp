#include "bvh_node.h"

bool bvh_node::bounding_box(Float t0, Float t1, aabb& b) const {
  b = box;
  return(true);
}

bool bvh_node::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
#ifndef DEBUGBVH
  if(box.hit(r, t_min, t_max, rng)) {
    if(left->hit(r,t_min,t_max,rec, rng)) {
      right->hit(r,t_min,rec.t,rec, rng);
      return(true);
    } else {
      return(right->hit(r,t_min,t_max,rec, rng));
    }
  }
  return(false);
#endif
#ifdef DEBUGBVH
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
#endif
}

// bool bvh_node::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng)  {
//   if (!box.hit(r, t_min, t_max, rng))
//     return false;
//   
//   bool hit_left = left->hit(r, t_min, t_max, rec, rng);
//   bool hit_right = right->hit(r, t_min, hit_left ? rec.t : t_max, rec, rng);
//   
//   return hit_left || hit_right;
// }

inline bool box_compare(const std::shared_ptr<hitable> a, const std::shared_ptr<hitable> b, int axis) {
  aabb box_a;
  aabb box_b;
  
  if (!a->bounding_box(0,0, box_a) || !b->bounding_box(0,0, box_b)){}

  return(box_a.min().e[axis] < box_b.min().e[axis]);
}


bool box_x_compare (const std::shared_ptr<hitable> a, const std::shared_ptr<hitable> b) {
  return box_compare(a, b, 0);
}

bool box_y_compare (const std::shared_ptr<hitable> a, const std::shared_ptr<hitable> b) {
  return box_compare(a, b, 1);
}

bool box_z_compare (const std::shared_ptr<hitable> a, const std::shared_ptr<hitable> b) {
  return box_compare(a, b, 2);
}

bvh_node::bvh_node(std::vector<std::shared_ptr<hitable> >& l, 
                   size_t start, size_t end,
                   Float time0, Float time1, int bvh_type, random_gen &rng) {
  aabb centroid_bounds;
  if(bvh_type == 1) {
    sah = true;
  } else {
    sah = false;
  }
  constexpr int nBuckets = 12;
  size_t n = end - start;

  std::vector<aabb> primitiveBounds(n);

  l[start]->bounding_box(time0, time1, centroid_bounds);
  aabb central_bounds(centroid_bounds.centroid);
  
  struct BucketInfo {
    int count = 0;
    aabb bounds;
  };
  
  BucketInfo buckets[nBuckets];

  for (int i = start; i < end; ++i) {
    aabb tempbox;
    if(l[i]->bounding_box(time0,time1,tempbox)) {
      centroid_bounds = surrounding_box(centroid_bounds, tempbox);
      central_bounds = surrounding_box(central_bounds, tempbox.centroid);
      primitiveBounds[i-start] = tempbox;
    }
  }
  
  vec3 centroid_bounds_values = central_bounds.max() - central_bounds.min();
  int axis = centroid_bounds_values.x() > centroid_bounds_values.y() ? 0 : 1;
  if(axis == 0) {
    axis = centroid_bounds_values.x() > centroid_bounds_values.z() ? 0 : 2;
  } else {
    axis = centroid_bounds_values.y() > centroid_bounds_values.z() ? 1 : 2;
  }
  auto comparator = (axis == 0) ? box_x_compare
    : (axis == 1) ? box_y_compare
    : box_z_compare;
  
  if (n == 1) {
    left = right = l[start];
  } else if (n == 2) {
    if (comparator(l[start], l[start+1])) {
      left = l[start];
      right = l[start+1];
    } else {
      left = l[start+1];
      right = l[start];
    }
  } else {
    std::sort(l.begin() + start, l.begin() + end, comparator);
    //SAH 
    if(sah) {
      for (int i = 0; i < n; ++i) {
        int b = nBuckets * central_bounds.offset(primitiveBounds[i].centroid)[axis];
        if (b == nBuckets) {
          b = nBuckets - 1;
        }
        buckets[b].count++;
        buckets[b].bounds = surrounding_box(buckets[b].bounds, primitiveBounds[i]);
      }
      int nSplits = nBuckets - 1;
      int countBelow[nSplits], countAbove[nSplits];
      aabb boundsBelow[nSplits], boundsAbove[nSplits];
      countBelow[0] = buckets[0].count;
      boundsBelow[0] = buckets[0].bounds;
      for (int i = 1; i < nSplits; ++i) {
        countBelow[i] = countBelow[i - 1] + buckets[i].count;
        if(buckets[i].count != 0) {
          boundsBelow[i] = surrounding_box(boundsBelow[i - 1], buckets[i].bounds);
        }
      }
      
      countAbove[nSplits - 1] = buckets[nBuckets - 1].count;
      boundsAbove[nSplits - 1] = buckets[nBuckets - 1].bounds;
      for (int i = nSplits - 2; i >= 0; --i) {
        countAbove[i] = countAbove[i + 1] + buckets[i + 1].count;
        if(buckets[i + 1].count != 0) {
          boundsAbove[i] = surrounding_box(boundsAbove[i + 1], buckets[i + 1].bounds);
        }
      }
      
      int minCostSplitBucket = -1;
      Float minCost = INFINITY;
      for (int i = 0; i < nSplits; ++i) {
        if (countBelow[i] == 0 || countAbove[i] == 0) {
          continue;
        }
        Float cost = (countBelow[i] * boundsBelow[i].surface_area() +
          countAbove[i] * boundsAbove[i].surface_area()) / centroid_bounds.surface_area();
        if (cost < minCost) {
          minCost = cost;
          minCostSplitBucket = i;
        }
      }
      int halfCount = countBelow[minCostSplitBucket];
      //End SAH
      left = std::make_shared<bvh_node>(l, start, start + halfCount, time0, time1, bvh_type, rng);
      right = std::make_shared<bvh_node>(l, start + halfCount, end, time0, time1, bvh_type, rng);
    } else {
      auto mid = start + n/2;
      left = std::make_shared<bvh_node>(l, start, mid, time0, time1, bvh_type, rng);
      right = std::make_shared<bvh_node>(l, mid, end, time0, time1, bvh_type, rng);
    }
  }

  aabb box_left, box_right;
  if(!left->bounding_box(time0,time1,box_left) || !right->bounding_box(time0,time1,box_right)) {
  }
  box = surrounding_box(box_left,box_right);
}
