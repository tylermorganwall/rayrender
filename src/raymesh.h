#ifndef RAYMESHH
#define RAYMESHH

#include <Rcpp.h>
#include "float.h"
#include "hitable.h"
#include "hitablelist.h"
#include "bvh.h"

class TextureCache;
struct TriangleMesh;
class random_gen;
class Transform;

class raymesh : public hitable {
public:
  raymesh() {}
  raymesh(Rcpp::List raymesh_list, 
          std::shared_ptr<material> default_material, 
          std::shared_ptr<alpha_texture> alpha_mask, std::shared_ptr<bump_texture> bump_tex,
          bool importance_sample_lights, bool calculate_consistent_normals, bool override_material,
          bool flip_transmittance, int subdivision_levels,
          std::string displacement_texture, Float displacement, bool displacement_vector,
          TextureCache &texCache, bool recalculate_normals,
          hitable_list& imp_sample_objects, bool verbose,
          Float shutteropen, Float shutterclose, int bvh_type, random_gen rng, 
          Transform* ObjectToWorld, Transform* WorldToObject, bool reverseOrientation);
  
  virtual const bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const;
  virtual const bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const;
  virtual bool HitP(const ray &r, Float t_min, Float t_max, random_gen& rng) const;
  virtual bool HitP(const ray &r, Float t_min, Float t_max, Sampler* sampler) const;

  Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0);
  Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0);
  vec3f random(const point3f& o, random_gen& rng, Float time = 0);
  vec3f random(const point3f& o, Sampler* sampler, Float time = 0);
  
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  virtual std::string GetName() const {
    return(std::string("RayMesh"));
  }
  size_t GetSize();
  virtual void hitable_info_bounds(Float t0, Float t1) const {
    aabb box;
    bounding_box(t0, t1, box);
    Rcpp::Rcout << GetName() << ": " <<  box.min() << "-" << box.max() << "\n";
  }
  std::pair<size_t,size_t> CountNodeLeaf();
  
  //Data Members
  std::unique_ptr<TriangleMesh> mesh;
  hitable_list triangles;
  
  //Hitable extras
  std::shared_ptr<BVHAggregate> tri_mesh_bvh;
};

#endif

