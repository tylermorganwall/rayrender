#ifndef PLYMESHH
#define PLYMESHH

#include "triangle.h"
#include "bvh_node.h"
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
          std::shared_ptr<Transform> ObjectToWorld, std::shared_ptr<Transform> WorldToObject, bool reverseOrientation);
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler);
  
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  virtual std::string GetName() const {
    return(std::string("Plymesh"));
  }
  size_t GetSize()  {
    return(ply_mesh_bvh->GetSize() + sizeof(triangles));
  }
  std::pair<size_t,size_t> CountNodeLeaf();
  
  std::unique_ptr<TriangleMesh> mesh;
  std::shared_ptr<bvh_node> ply_mesh_bvh;
  hitable_list triangles;
};


#endif

