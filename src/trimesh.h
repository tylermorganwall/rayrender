#ifndef TRIMESHH
#define TRIMESHH

#include "triangle.h"
#include "bvh_node.h"
#include "rng.h"
#ifndef STBIMAGEH
#define STBIMAGEH
#include "stb_image.h"
#endif
#include <Rcpp.h>

inline char separator() {
  #if defined _WIN32 || defined __CYGWIN__
    return '\\';
  #else
    return '/';
  #endif
}

class trimesh : public hitable {
public:
  trimesh() {}
  ~trimesh() {
    delete tri_mesh_bvh;
    for(auto mat : obj_materials) {
      if(mat) stbi_image_free(mat);
    }
    for(auto bump : bump_materials) {
      if(bump) stbi_image_free(bump);
    }
    delete mat_ptr;
  }
  trimesh(std::string inputfile, std::string basedir, Float scale, 
          Float shutteropen, Float shutterclose, random_gen rng);
  trimesh(std::string inputfile, std::string basedir, Float scale, Float sigma,
          Float shutteropen, Float shutterclose, random_gen rng);
  trimesh(std::string inputfile, std::string basedir, material *mat, 
          Float scale, Float shutteropen, Float shutterclose, random_gen rng);
  trimesh(std::string inputfile, std::string basedir, float vertex_color_sigma,
          Float scale, bool is_vertex_color, Float shutteropen, Float shutterclose, random_gen rng);
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  
  bvh_node* tri_mesh_bvh;
  material *mat_ptr;
  std::vector<Float* > obj_materials;
  std::vector<Float* > bump_materials;
  std::vector<bump_texture* > bump_textures;
  std::vector<alpha_texture* > alpha_materials;
  std::vector<hitable* > triangles;
};

#endif
