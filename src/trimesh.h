#ifndef TRIMESHH
#define TRIMESHH

#include "trianglemesh.h"
#include "triangle.h"
#include "bvh_node.h"
#include "rng.h"
#ifndef STBIMAGEH
#define STBIMAGEH
#include "stb_image.h"
#endif
#include <Rcpp.h>


class trimesh : public hitable {
public:
  trimesh() {}
  trimesh(std::string inputfile, std::string basedir, Float scale, Float sigma,
          std::shared_ptr<material> default_material, bool load_materials, bool load_textures, bool load_vertex_colors,
          bool importance_sample_lights, bool load_normals, bool calculate_consistent_normals,
          hitable_list& imp_sample_objects,
          Float shutteropen, Float shutterclose, int bvh_type, random_gen rng, bool verbose,
          std::shared_ptr<Transform> ObjectToWorld, std::shared_ptr<Transform> WorldToObject, bool reverseOrientation);
  
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler);
  
  Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0);
  Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0);
  vec3f random(const point3f& o, random_gen& rng, Float time = 0);
  vec3f random(const point3f& o, Sampler* sampler, Float time = 0);
  
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  virtual std::string GetName() const {
    return(std::string("TriangleMesh"));
  }
  size_t GetSize();
  std::pair<size_t,size_t> CountNodeLeaf();
  
  //Data Members
  std::unique_ptr<TriangleMesh> mesh;
  hitable_list triangles;

  //Hitable extras
  std::shared_ptr<material> mat_ptr;
  std::shared_ptr<bvh_node> tri_mesh_bvh;
};


#endif
