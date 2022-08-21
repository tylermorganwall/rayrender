#ifndef BVHNODEH
#define BVHNODEH

#include "hitable.h"
#include "aabb.h"
#include <Rcpp.h>
#include "material.h"

class bvh_node : public hitable {
  public:
    bvh_node() {}
#ifndef DEBUGBBOX 
    bvh_node(hitable_list& l,  
             Float time0, Float time1, int bvh_type, random_gen &rng) :
      bvh_node(l.objects, 0 ,l.objects.size(), time0, time1, bvh_type, rng) {};
    bvh_node(std::vector<std::shared_ptr<hitable> >& l, 
             size_t start, size_t end, 
             Float time0, Float time1, int bvh_type, random_gen &rng);
#endif
#ifdef DEBUGBBOX
    bvh_node(hitable_list& l,  
             Float time0, Float time1, int bvh_type, random_gen &rng) :
      bvh_node(l.objects, 0 ,l.objects.size(), time0, time1, bvh_type, 0, rng) {};
    bvh_node(std::vector<std::shared_ptr<hitable> >& l, 
             size_t start, size_t end, 
             Float time0, Float time1, int bvh_type, random_gen &rng);
    bvh_node(hitable_list& l,
             Float time0, Float time1, int bvh_type, int depth, random_gen &rng) :
      bvh_node(l.objects, 0 ,l.objects.size(), time0, time1, bvh_type, depth, rng) {};
    bvh_node(std::vector<std::shared_ptr<hitable> >& l,
             size_t start, size_t end,
             Float time0, Float time1, int bvh_type, int depth, random_gen &rng);
#endif

    virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
    virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler);

    virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
    
    Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0);
    Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0);
    vec3f random(const point3f& o, random_gen& rng, Float time = 0);
    vec3f random(const point3f& o, Sampler* sampler, Float time = 0);
    
    std::string GetName() const {
      return(std::string("BVH Node"));
    }
    size_t GetSize();
    std::pair<size_t,size_t> CountNodeLeaf();
    
    std::shared_ptr<hitable> left;
    std::shared_ptr<hitable> right;
    aabb box;
    bool sah;
#ifdef DEBUGBBOX
    int depth;
#endif
};

#endif
