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


bvh_node::bvh_node(hitable **l, int n, Float time0, Float time1, random_gen &rng) {
  aabb centroid_bounds;
  constexpr int nBuckets = 12;
  std::vector<aabb> primitiveBounds(n);
  
  l[0]->bounding_box(time0, time1, centroid_bounds);
  struct BucketInfo {
    int count = 0;
    aabb bounds;
  };
  
  BucketInfo buckets[nBuckets];
  Rcpp::Rcout << "Surface Areas: ";
  for (int i = 0; i < n; ++i) {
    aabb tempbox;
    if(l[i]->bounding_box(time0,time1,tempbox)) {
      centroid_bounds = surrounding_box(centroid_bounds, tempbox);
      primitiveBounds[i] = tempbox;
      // Rcpp::Rcout << tempbox.centroid << " i/n: " << i <<"/"<< n <<"\n";
    }
    Rcpp::Rcout << tempbox.surface_area() << " ";
  }
  Rcpp::Rcout <<  "\n";
  
  vec3 centroid_bounds_values = centroid_bounds.max() - centroid_bounds.min();
  int axis = centroid_bounds_values.x() > centroid_bounds_values.y() ? 0 : 1;
  if(axis == 0) {
    axis = centroid_bounds_values.x() > centroid_bounds_values.z() ? 0 : 2;
  } else {
    axis = centroid_bounds_values.y() > centroid_bounds_values.z() ? 1 : 2;
  }
  
  if(n == 1) {
    Rcpp::Rcout << "Single" << "\n";
    left = right = l[0];
  } else if (n == 2) {
    Rcpp::Rcout << "Double" << "\n";
    left = l[0];
    right = l[1];
  } else if (n == 3) {
    Rcpp::Rcout << "Triple" << "\n";
    
    if(axis == 0) {
      std::qsort(l, n, sizeof(hitable *), box_x_compare);
    } else if (axis == 1) {
      std::qsort(l, n, sizeof(hitable *), box_y_compare);
    } else {
      std::qsort(l, n, sizeof(hitable *), box_z_compare);
    }
    left = new bvh_node(l, 1, time0, time1, rng);
    right = new bvh_node(l + 1, n - 1, time0, time1, rng);
  } else {
    Rcpp::Rcout << "Multi" << "\n";
    //SAH
    // int axis = rng.unif_rand() * 2.99999;
    if(axis == 0) {
      std::qsort(l, n, sizeof(hitable *), box_x_compare);
    } else if (axis == 1) {
      std::qsort(l, n, sizeof(hitable *), box_y_compare);
    } else {
      std::qsort(l, n, sizeof(hitable *), box_z_compare);
    }
    for (int i = 0; i < n; ++i) {
      int b = nBuckets * centroid_bounds.offset(primitiveBounds[i].centroid)[axis];
      if (b == nBuckets) {
        b = nBuckets - 1;
      }
      buckets[b].count++;
      buckets[b].bounds = surrounding_box(buckets[b].bounds, primitiveBounds[i]);
    }
    // Rcpp::Rcout << "Centroid bounds: " << centroid_bounds.bounds[0] 
    //             << "--" << centroid_bounds.bounds[1] << "\n";
    // Rcpp::Rcout << "Centroid area: " << centroid_bounds.surface_area() << "\n";
    // Rcpp::Rcout << "Buckets: ";
    Float cost[nBuckets - 1];
    for (int i = 0; i < nBuckets - 1; ++i) {
      aabb b0, b1;
      int count0 = 0, count1 = 0;
      for (int j = 0; j <= i; ++j) {
        if(buckets[j].count > 0) {
          b0 = surrounding_box(b0, buckets[j].bounds);
          count0 += buckets[j].count;
        }
      }
      for (int j = i+1; j < nBuckets; ++j) {
        if(buckets[j].count > 0) {
          b1 = surrounding_box(b1, buckets[j].bounds);
          count1 += buckets[j].count;
        }
      }
      cost[i] = 0.125f + (count0 * b0.surface_area() +
        count1 *  b1.surface_area()) / centroid_bounds.surface_area();
      cost[i] = !std::isnan(cost[i]) ? cost[i] : INFINITY;
      Rcpp::Rcout << buckets[i].count << "/" << cost[i] << " " ;
    }
    if(n == 0) {
      throw std::runtime_error("bucket size");
    }
    Rcpp::Rcout << "\n";
    Float minCost = cost[0];
    int minCostSplitBucket = 0;
    for (int i = 1; i < nBuckets - 1; ++i) {
      // Rcpp::Rcout << buckets[i].count << " ";
      if (cost[i] < minCost) {
        minCost = cost[i];
        minCostSplitBucket = i;
      }
    }
    Rcpp::Rcout << "\n" ;//<< "costs: ";
    
    // for (int i = 0; i < nBuckets - 1; ++i) {
    //   Rcpp::Rcout << cost[i] << " ";
    // }
    // Rcpp::Rcout << "\n";
    int halfCount = 0;
    for (int i = 0; i < minCostSplitBucket+1; ++i) {
      halfCount += buckets[i].count;
    }
    Rcpp::Rcout << "halfCount: " << halfCount << " minBucket: " << minCostSplitBucket << " N: " << n << "\n";
    // Rcpp::Rcout << "done buckets: " << buckets[minCostSplitBucket].count  << " " << n << "\n" ;
    left = new bvh_node(l, halfCount, time0, time1, rng);
    right = new bvh_node(l + halfCount, n - halfCount, time0, time1, rng);
  }
  aabb box_left, box_right;
  if(!left->bounding_box(time0,time1,box_left) || !right->bounding_box(time0,time1,box_right)) {
  }
  box = surrounding_box(box_left,box_right);
}
