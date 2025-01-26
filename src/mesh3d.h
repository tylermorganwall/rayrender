#ifndef MESH3DH
#define MESH3DH

#include "triangle.h"
#include "bvh.h"


#include <Rcpp.h>

class mesh3d : public hitable {
  public:
    mesh3d() {}
    ~mesh3d() {}
    mesh3d(Rcpp::List mesh_info, std::shared_ptr<material>  mat, 
           std::string displacement_texture, Float displacement, bool displacement_vector, 
           TextureCache &texCache, bool recalculate_normals,
           bool verbose,
           Float shutteropen, Float shutterclose, int bvh_type, random_gen rng,
           Transform* ObjectToWorld, Transform* WorldToObject, bool reverseOrientation);
    virtual const bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const;
    virtual const bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const;
    virtual bool HitP(const ray &r, Float t_min, Float t_max, random_gen& rng) const;
    virtual bool HitP(const ray &r, Float t_min, Float t_max, Sampler* sampler) const;

    virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
    virtual std::string GetName() const {
      return(std::string("Mesh3d"));
    }
    size_t GetSize()  {
      return(sizeof(*this) + mesh_bvh->GetSize());
    }
    virtual void hitable_info_bounds(Float t0, Float t1) const {
      aabb box;
      bounding_box(t0, t1, box);
      Rcpp::Rcout << GetName() << ": " <<  box.min() << "-" << box.max() << "\n";
    }
    std::pair<size_t,size_t> CountNodeLeaf();
    
    
    //Data Members
    std::unique_ptr<TriangleMesh> mesh;
    hitable_list triangles;
    
    //Hitable extra
    std::shared_ptr<BVHAggregate> mesh_bvh;
};


#endif

