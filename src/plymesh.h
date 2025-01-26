#ifndef PLYMESHH
#define PLYMESHH

#include "triangle.h"

#include "bvh.h"
#include <Rcpp.h>


enum class Topology;
struct TriMesh;

class plymesh : public hitable {
  public:
    plymesh() {}
   ~plymesh() {}
  plymesh(std::string inputfile, std::string basedir, std::shared_ptr<material> mat, 
          std::shared_ptr<alpha_texture> alpha, std::shared_ptr<bump_texture> bump,
          Float scale, int subdivision_levels, bool recalculate_normals,
          bool verbose,
          Float shutteropen, Float shutterclose, int bvh_type, random_gen rng,
          Transform* ObjectToWorld, Transform* WorldToObject, bool reverseOrientation);
  virtual const bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const;
  virtual const bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const;
  virtual bool HitP(const ray &r, Float t_min, Float t_max, random_gen& rng) const;
  virtual bool HitP(const ray &r, Float t_min, Float t_max, Sampler* sampler) const;

  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  virtual std::string GetName() const {
    return(std::string("Plymesh"));
  }
  size_t GetSize()  {
    return(ply_mesh_bvh->GetSize() + sizeof(triangles));
  }
  virtual void hitable_info_bounds(Float t0, Float t1) const {
    aabb box;
    bounding_box(t0, t1, box);
    Rcpp::Rcout << GetName() << ": " <<  box.min() << "-" << box.max() << "\n";
  }
  std::pair<size_t,size_t> CountNodeLeaf();
  
  std::unique_ptr<TriangleMesh> mesh;
  std::shared_ptr<BVHAggregate> ply_mesh_bvh;
  hitable_list triangles;
};


#endif

