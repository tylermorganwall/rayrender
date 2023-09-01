#include "bvh_node.h"


#ifdef DEBUGBBOX
#include <iostream>
#include <fstream>
using namespace std;
#endif


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
    rec.bvh_nodes += 1.0;
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

bool bvh_node::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
#ifndef DEBUGBVH
  if(box.hit(r, t_min, t_max, sampler)) {
    if(left->hit(r,t_min,t_max,rec, sampler)) {
      right->hit(r,t_min,rec.t,rec, sampler);
      return(true);
    } else {
      return(right->hit(r,t_min,t_max,rec, sampler));
    }
  }
  return(false);
#endif
#ifdef DEBUGBVH
  if(box.hit(r, t_min, t_max, sampler)) {
    rec.bvh_nodes += 1.0;
    if(left->hit(r,t_min,t_max,rec, sampler)) {
      right->hit(r,t_min,rec.t,rec, sampler);
      return(true);
    } else {
      return(right->hit(r,t_min,t_max,rec, sampler));
    }
  }
  return(false);
#endif
}

inline bool box_compare(const std::shared_ptr<hitable> a, const std::shared_ptr<hitable> b, int axis) {
  aabb box_a;
  aabb box_b;
  
  if (!a->bounding_box(0,0, box_a) || !b->bounding_box(0,0, box_b)){}

  return(box_a.Centroid().e[axis] < box_b.Centroid().e[axis]);
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

#ifdef DEBUGBBOX
bvh_node::bvh_node(std::vector<std::shared_ptr<hitable> >& l, 
                   size_t start, size_t end,
                   Float time0, Float time1, int bvh_type, int depth, random_gen &rng) {
  if(start == end) {
    throw std::runtime_error("start node must not equal end node");
  }
#endif
#ifndef DEBUGBBOX
  bvh_node::bvh_node(std::vector<std::shared_ptr<hitable> >& l, 
                     size_t start, size_t end,
                     Float time0, Float time1, int bvh_type, random_gen &rng) {
#endif
  aabb centroid_bounds;
  if(bvh_type == 1) {
    sah = true;
  } else {
    sah = false;
  }
  constexpr int nBuckets = 12;
  size_t n = end - start;
  
  //Don't need sorted, just used to count bins and generate bin bounds
  //Contains AABB of each primitve
  std::vector<aabb> primitiveBounds(n); 

  l[start]->bounding_box(time0, time1, centroid_bounds);
  aabb central_bounds;


  for (unsigned int i = start; i < end; ++i) {
    aabb tempbox;
    if(l[i]->bounding_box(time0,time1,tempbox)) {
      centroid_bounds = surrounding_box(centroid_bounds, tempbox);
      central_bounds = surrounding_box(central_bounds, tempbox.Centroid());
      primitiveBounds[i-start] = tempbox;
    }
  }
  
  vec3f centroid_bounds_values = central_bounds.max() - central_bounds.min();

#ifdef DEBUGBBOX
  if(centroid_bounds_values.x() < 0 || centroid_bounds_values.y() < 0 || centroid_bounds_values.z() < 0) {
    throw std::runtime_error("centroid extent less than 0");
  }
      ofstream myfile;
      myfile.open("bbox.txt", ios::app | ios::out);
      myfile << "Min: " << central_bounds.min() << ", Max: " << central_bounds.max() << 
        ", Diag: " << central_bounds.diag << ", Vol: " << central_bounds.Volume() << ", N: " << n << ", Depth:" << depth << "\n";
      myfile.close();
  #endif
  
  int axis = centroid_bounds_values.x() > centroid_bounds_values.y() ? 0 : 1;
  if(axis == 0) {
    axis = centroid_bounds_values.x() > centroid_bounds_values.z() ? 0 : 2;
  } else {
    axis = centroid_bounds_values.y() > centroid_bounds_values.z() ? 1 : 2;
  }
  // int axis = 0;
  auto comparator = (axis == 0) ? box_x_compare
    : (axis == 1) ? box_y_compare
    : box_z_compare;

  if(central_bounds.Volume() < 1e-6) {
    sah = false;
    bvh_type = 2;
  }
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
    //Handle case where all shapes share the same centroid
    if(central_bounds.Diag().e[axis] == 0 ) {
      sah = false;
      bvh_type = 2;
    }      
    
    if(n <= 4) {
      sah = false;
      bvh_type = 2;
    }
    //SAH 
    if(sah) {
      struct BucketInfo {
        int count = 0;
        aabb bounds;
      };
      
      BucketInfo buckets[nBuckets];
      
      //Count number of objects in each bin and calculate bounding box for each bin.
      for (unsigned int i = 0; i < n; ++i) {
        int b = nBuckets * central_bounds.offset(primitiveBounds[i].Centroid())[axis];
        if (b == nBuckets) {
          b = nBuckets - 1;
        }
#ifdef DEBUGBBOX
        if(b < 0 || b > nBuckets - 1) {
          throw std::runtime_error("SAH bucket out of bounds");
        }
#endif
        buckets[b].count++;
        buckets[b].bounds = surrounding_box(buckets[b].bounds, primitiveBounds[i]);
      }
      constexpr int nSplits = nBuckets - 1;
      int countBelow[nSplits], countAbove[nSplits];
      aabb boundsBelow[nSplits], boundsAbove[nSplits];
      countBelow[0] = buckets[0].count;
      boundsBelow[0] = buckets[0].bounds;
      for (int i = 1; i < nSplits; ++i) {
        countBelow[i] = countBelow[i - 1] + buckets[i].count;
        boundsBelow[i] = surrounding_box(boundsBelow[i - 1], buckets[i].bounds);
      }
      
      countAbove[nSplits - 1] = buckets[nBuckets - 1].count;
      boundsAbove[nSplits - 1] = buckets[nBuckets - 1].bounds;
      for (int i = nSplits - 2; i >= 0; --i) {
        countAbove[i] = countAbove[i + 1] + buckets[i + 1].count;
        boundsAbove[i] = surrounding_box(boundsAbove[i + 1], buckets[i + 1].bounds);
      }

      int minCostSplitBucket = -1;
      Float minCost = INFINITY;
      Float costs[nBuckets];
      
      for (unsigned int i = 0; i < nSplits; ++i) {
        if (countBelow[i] == 0 || countAbove[i] == 0) {
          continue;
        }
        Float cost = (countBelow[i] * boundsBelow[i].surface_area() +
          countAbove[i] * boundsAbove[i].surface_area());
        costs[i] = cost;
        
        if (cost < minCost) {
          minCost = cost;
          minCostSplitBucket = i;
        }
      }
      int halfCount = countBelow[minCostSplitBucket];
      //End SAH
#ifdef DEBUGBBOX
      depth++;
      left = std::make_shared<bvh_node>(l, start, start + halfCount, time0, time1, bvh_type, depth, rng);
      right = std::make_shared<bvh_node>(l, start + halfCount, end, time0, time1, bvh_type, depth, rng);
#endif
#ifndef DEBUGBBOX
      left = std::make_shared<bvh_node>(l, start, start + halfCount, time0, time1, bvh_type, rng);
      right = std::make_shared<bvh_node>(l, start + halfCount, end, time0, time1, bvh_type, rng);
#endif
    } else {
#ifdef DEBUGBBOX
      depth++;
      auto mid = start + n/2;
      left = std::make_shared<bvh_node>(l, start, mid, time0, time1, bvh_type, depth, rng);
      right = std::make_shared<bvh_node>(l, mid, end, time0, time1, bvh_type,  depth, rng);
#endif
#ifndef DEBUGBBOX
      auto mid = start + n/2;
      left = std::make_shared<bvh_node>(l, start, mid, time0, time1, bvh_type, rng);
      right = std::make_shared<bvh_node>(l, mid, end, time0, time1, bvh_type, rng);
#endif
    }
  }

  aabb box_left, box_right;
  if(!left->bounding_box(time0,time1,box_left) || !right->bounding_box(time0,time1,box_right)) {
  }
  box = surrounding_box(box_left,box_right);
}
  
  
Float bvh_node::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) {
  return(0.5*left->pdf_value(o,v, rng, time) + 0.5*right->pdf_value(o,v, rng, time));
}

Float bvh_node::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) {
  return(0.5*left->pdf_value(o,v, sampler, time) + 0.5*right->pdf_value(o,v, sampler, time));
  
}

vec3f bvh_node::random(const point3f& o, random_gen& rng, Float time) {
  return(rng.unif_rand() > 0.5 ? left->random(o, rng, time) : right->random(o, rng, time));
}

vec3f bvh_node::random(const point3f& o, Sampler* sampler, Float time) {
  return(sampler->Get1D() > 0.5 ? left->random(o, sampler, time) : right->random(o, sampler, time));
  
}
  
size_t bvh_node::GetSize()  {
  return(left != right ? sizeof(*this) + left->GetSize() + right->GetSize() :
           sizeof(*this) + left->GetSize());
}
  
std::pair<size_t,size_t> operator+(const std::pair<size_t,size_t> & l,const std::pair<size_t,size_t> & r) {   
  return {l.first+r.first,l.second+r.second};                                    
}   

std::pair<size_t,size_t> bvh_node::CountNodeLeaf()  {
  std::pair<size_t,size_t> count(1,0);
  if(left != right) {
    count = count + (left->CountNodeLeaf() + right->CountNodeLeaf());
  } else {
    count = count + left->CountNodeLeaf();
  }
  return count;
}
  
void bvh_node::validate_bvh() {
  validate_bvh_node(this);
}

void bvh_node::validate_bvh_node(const bvh_node* node) {
  if (!node) {
    throw std::runtime_error("Encountered a nullptr node in BVH.");
  }
  
  // Retrieve bounding boxes
  aabb box_left, box_right;
  bool has_left_bbox = node->left->bounding_box(0, 0, box_left);  // Using 0 for time as a placeholder
  bool has_right_bbox = node->right->bounding_box(0, 0, box_right); // Using 0 for time as a placeholder
  
  // Ensure both children have valid bounding boxes
  if (!has_left_bbox || !has_right_bbox) {
    throw std::runtime_error("A child node doesn't have a valid bounding box.");
  }
  
  // If a node has only one child, both left and right should be the same
  if ((node->left == node->right) && !node->left) {
    throw std::runtime_error("Node with single child doesn't set both left and right pointers to the same child.");
  }
  
  // If it's not a leaf node, continue validating the left and right children
  if (node->left != node->right) {
    if (auto left_bvh = dynamic_cast<const bvh_node*>(node->left.get())) {
      validate_bvh_node(left_bvh);
    }
    if (auto right_bvh = dynamic_cast<const bvh_node*>(node->right.get())) {
      validate_bvh_node(right_bvh);
    }
  }
}
  
  
