#ifndef BVHNODEH
#define BVHNODEH

#include "hitable.h"
#include "aabb.h"
#include <Rcpp.h>
#include "material.h"

class bvh_node : public hitable {
  public:
    bvh_node() {}
    bvh_node(hitable_list& l,  
             Float time0, Float time1, int bvh_type, random_gen &rng) :
      bvh_node(l.objects, 0 ,l.objects.size(), time0, time1, bvh_type, rng) {};
    bvh_node(std::vector<std::shared_ptr<hitable> >& l, 
             size_t start, size_t end, 
             Float time0, Float time1, int bvh_type, random_gen &rng);
    virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
    virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
    std::shared_ptr<hitable> left;
    std::shared_ptr<hitable> right;
    aabb box;
    bool sah;
};

#endif
