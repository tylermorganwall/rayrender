#ifndef BVHNODEH
#define BVHNODEH

#include "hitable.h"
#include "aabb.h"
#include <Rcpp.h>


class bvh_node : public hitable {
  public:
    bvh_node() {}
    bvh_node(hitable **l, int n, float time0, float time1, random_gen rng);
    virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec, random_gen& rng);
    virtual bool bounding_box(float t0, float t1, aabb& box) const;
    hitable *left;
    hitable *right;
    aabb box;
};

struct BucketInfo {
  int count = 0;
  bool empty = true;
  aabb bounds;
};

struct BVHPrimitiveInfo {
  BVHPrimitiveInfo() {}
  BVHPrimitiveInfo(int primitive_number, const aabb &bounds)
    : primitive_number(primitive_number),
      bounds(bounds),
      centroid(0.5f * bounds._min + 0.5f * bounds._max) {}
  int primitive_number;
  aabb bounds;
  vec3 centroid;
};

bool bvh_node::bounding_box(float t0, float t1, aabb& b) const {
  b = box;
  return(true);
}

bool bvh_node::hit(const ray& r, float t_min, float t_max, hit_record& rec, random_gen& rng) {
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

bvh_node::bvh_node(hitable **l, int n, float time0, float time1, random_gen rng) {
  aabb centroid_bounds;
  l[0]->bounding_box(time0, time1, centroid_bounds);
  std::vector<BVHPrimitiveInfo> primitiveInfo(n);
  for (int i = 0; i < n; ++i) {
    aabb tempbox;
    if(l[i]->bounding_box(time0,time1,tempbox)) {
      primitiveInfo[i] = { i, tempbox };
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
  if(n <= 4) {
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
  } else {
    constexpr int nBuckets = 12;
    BucketInfo buckets[nBuckets];
    aabb tempbox;
    for (int i = 0; i < n; ++i) {
      if(l[i]->bounding_box(time0, time1, tempbox)) {
        int b = nBuckets * centroid_bounds.offset(tempbox.centroid)[axis];
        if (b == nBuckets) {
          b = nBuckets - 1;
        }
        buckets[b].count++;
        if(!buckets[b].empty) {
          buckets[b].bounds = surrounding_box(buckets[b].bounds, primitiveInfo[i].bounds);
        } else {
          buckets[b].bounds = primitiveInfo[i].bounds;
          buckets[b].empty = false;
        }
      }
    }
    int boundaries[nBuckets-1];
    int totalcount = buckets[0].count;
    boundaries[0] = totalcount;
    box = buckets[0].bounds;
    for(int i = 1; i < nBuckets; i++) {
      box = surrounding_box(box, buckets[i].bounds);
      totalcount = totalcount + buckets[i].count;
      boundaries[i] = totalcount;
    }
    
    //Compute cost
    float cost[nBuckets - 1];
    aabb box_left, box_right;
    
    for (int i = 0; i < nBuckets - 1; ++i) {
      box_left = buckets[0].bounds;
      int count0 = 0, count1 = 0;
      for (int j = 0; j <= i; ++j) {
        box_left = surrounding_box(box_left, buckets[j].bounds);
        count0 += buckets[j].count;
      }
      box_right = buckets[i+1].bounds;
      for (int j = i+1; j < nBuckets; ++j) {
        box_right = surrounding_box(box_right, buckets[j].bounds);
        count1 += buckets[j].count;
      }
      cost[i] = 0.125f + (count0 * box_left.surface_area() + 
        count1 * box_right.surface_area()) / box.surface_area();
    }
     
    
    float minCost = cost[0];
    int minCostSplitBucket = 0;
    int splitpoint = boundaries[0];
    for (int i = 1; i < nBuckets - 1; ++i) { 
      if (cost[i] < minCost) {
        minCost = cost[i];
        minCostSplitBucket = i;
        splitpoint = boundaries[i];
      }
    }
    if(splitpoint == n) {
      splitpoint--;
    }
    if(splitpoint == 0) {
      splitpoint++;
    }
    if(splitpoint == 1) {
      left = l[0];
    } else {
      left = new bvh_node(l, splitpoint , time0, time1, rng);
    }
    if(splitpoint == n - 1) {
      right = l[n];
    } else {
      right = new bvh_node(l + splitpoint + 1, n - splitpoint - 1, time0, time1, rng);
    }
  }
}

#endif
