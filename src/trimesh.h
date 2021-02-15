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
    for(auto mat : obj_materials) {
      if(mat) stbi_image_free(mat);
    }
    for(auto bump : bump_materials) {
      if(bump) stbi_image_free(bump);
    }
  }
  trimesh(std::string inputfile, std::string basedir, Float scale, 
          Float shutteropen, Float shutterclose, int bvh_type, random_gen rng);
  trimesh(std::string inputfile, std::string basedir, Float scale, Float sigma,
          Float shutteropen, Float shutterclose, int bvh_type, random_gen rng);
  trimesh(std::string inputfile, std::string basedir, std::shared_ptr<material> mat, 
          Float scale, Float shutteropen, Float shutterclose, int bvh_type, random_gen rng);
  trimesh(std::string inputfile, std::string basedir, float vertex_color_sigma,
          Float scale, bool is_vertex_color, Float shutteropen, Float shutterclose, int bvh_type, 
          random_gen rng);
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler);
  
  Float pdf_value(const vec3& o, const vec3& v, random_gen& rng, Float time = 0);
  Float pdf_value(const vec3& o, const vec3& v, Sampler* sampler, Float time = 0);
  vec3 random(const vec3& o, random_gen& rng, Float time = 0);
  vec3 random(const vec3& o, Sampler* sampler, Float time = 0);
  
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  
  std::shared_ptr<bvh_node> tri_mesh_bvh;
  std::shared_ptr<material> mat_ptr;
  std::vector<Float* > obj_materials;
  std::vector<Float* > bump_materials;
  std::vector<std::shared_ptr<bump_texture> > bump_textures;
  std::vector<std::shared_ptr<alpha_texture> > alpha_materials;
  hitable_list triangles;
};

#endif
