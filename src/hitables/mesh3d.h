#ifndef MESH3DH
#define MESH3DH

#include "../hitables/triangle.h"
#include "../core/bvh.h"
#include "../math/model_animation.h"


#include <Rcpp.h>

class mesh3d : public hitable {
  public:
    mesh3d() {}
    ~mesh3d() {}
    mesh3d(Rcpp::List mesh_info, std::shared_ptr<material>  mat, 
           std::string displacement_texture, Float displacement, bool displacement_vector, 
           TextureCache &texCache, bool recalculate_normals,
           bool verbose,
           Float shutteropen, Float shutterclose, const BVHBuildOptions& bvh_type, random_gen rng,
           Transform* ObjectToWorld, Transform* WorldToObject, bool reverseOrientation);
    virtual const bool hit(const Ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const override;
    virtual const bool hit(const Ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const override;
    virtual bool HitP(const Ray &r, Float t_min, Float t_max, random_gen& rng) const override;
    virtual bool HitP(const Ray &r, Float t_min, Float t_max, Sampler* sampler) const override;

  OpaqueShadowType ShadowType() const override {
    return mesh_bvh->ShadowType();
  }
  bool OpaqueHit(const Ray& r, Float t_min, Float t_max, random_gen& rng) const override {
    return mesh_bvh->OpaqueHit(r, t_min, t_max, rng);
  }

    virtual bool bounding_box(Float t0, Float t1, aabb& box) const override;
    virtual std::string GetName() const override {
      return(std::string("Mesh3d"));
    }
    size_t GetSize() override {
      return(sizeof(*this) + mesh_bvh->GetSize());
    }
    virtual void hitable_info_bounds(Float t0, Float t1) const override {
      aabb box;
      bounding_box(t0, t1, box);
      Rcpp::Rcout << GetName() << ": " <<  box.min() << "-" << box.max() << "\n";
    }
    std::pair<size_t,size_t> CountNodeLeaf() override;
    void BindDeformation(std::shared_ptr<const animation::ModelInstance>, size_t node, size_t primitive);
    std::function<void()> PrepareDeformation(double time);
    bool HasDeformation() const { return bool(deformation); }
    uint64_t DeformationRevision() const { return deformation_revision; }
    Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0) override;
    Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0) override;
    vec3f random(const point3f& o, random_gen& rng, Float time = 0) override { return mesh_bvh->random(o,rng,time); }
    vec3f random(const point3f& o, Sampler* sampler, Float time = 0) override { return mesh_bvh->random(o,sampler,time); }
    
    
    //Data Members
    std::unique_ptr<TriangleMesh> mesh;
    hitable_list triangles;
    
    //Hitable extra
    std::shared_ptr<BVHAggregate> mesh_bvh;
private:
    std::shared_ptr<const animation::ModelInstance> deformation;
    size_t deformation_node = 0, deformation_primitive = 0;
    std::optional<animation::ModelPose> last_pose;
    uint64_t deformation_revision = 0;
};


#endif
