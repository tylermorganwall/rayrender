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
   ~plymesh() {
    delete ply_mesh_bvh;
    delete mat_ptr;
  }
  plymesh(std::string inputfile, std::string basedir, material *mat, 
          Float scale, Float shutteropen, Float shutterclose, random_gen rng);
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  bvh_node* ply_mesh_bvh;
  material *mat_ptr;
  std::vector<hitable* > triangles;
};


#endif

