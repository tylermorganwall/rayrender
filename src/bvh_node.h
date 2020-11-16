#ifndef BVHNODEH
#define BVHNODEH

#include "hitable.h"
#include "aabb.h"
#include <Rcpp.h>
#include "material.h"

class bvh_node : public hitable {
  public:
    bvh_node() {}
    ~bvh_node() {
      if(left != right) {
        delete left;
        delete right;
      } else {
        delete left;
      }
    }
    
    bvh_node(hitable **l, int n, Float time0, Float time1, random_gen &rng);
    virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
    virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
    hitable *left;
    hitable *right;
    aabb box;
};

#endif
