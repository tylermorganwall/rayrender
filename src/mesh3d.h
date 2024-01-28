#ifndef MESH3DH
#define MESH3DH

#include "triangle.h"
#include "bvh_node.h"
#ifndef STBIMAGEH
#define STBIMAGEH
#include "stb_image.h"
#endif
#include <Rcpp.h>

class mesh3d : public hitable {
  public:
    mesh3d() {}
    ~mesh3d() {}
    mesh3d(Rcpp::List mesh_info, std::shared_ptr<material>  mat, 
           Float shutteropen, Float shutterclose, int bvh_type, random_gen rng,
           std::shared_ptr<Transform> ObjectToWorld, std::shared_ptr<Transform> WorldToObject, bool reverseOrientation);
    virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
    virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler);
    
    virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
    virtual std::string GetName() const {
      return(std::string("Mesh3d"));
    }
    size_t GetSize()  {
      return(sizeof(*this) + mesh_bvh->GetSize());
    }
    std::pair<size_t,size_t> CountNodeLeaf();
    
    
    //Data Members
    std::unique_ptr<TriangleMesh> mesh;
    hitable_list triangles;
    
    //Hitable extra
    std::shared_ptr<bvh_node> mesh_bvh;
};


#endif

