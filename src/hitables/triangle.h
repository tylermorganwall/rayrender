#ifndef TRIANGLEH
#define TRIANGLEH

#include "../hitables/hitable.h"
#include "../materials/material.h"
#include "../math/onbh.h"
#include "../hitables/trianglemesh.h"

class triangle : public hitable {
public:
  triangle() : face_number(0) {}
  triangle(TriangleMesh* mesh, const int *v, const int *n, const int *t, const int face_number,
           Transform* ObjectToWorld, Transform* WorldToObject, bool reverseOrientation) : 
    hitable(ObjectToWorld, WorldToObject, nullptr, reverseOrientation), mesh(mesh), v(v), n(n), t(t), face_number(face_number) {
    
  }
  virtual const bool hit(const Ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const override;
  virtual const bool hit(const Ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const override;
  virtual bool HitP(const Ray &r, Float t_min, Float t_max, random_gen& rng) const override;
  virtual bool HitP(const Ray &r, Float t_min, Float t_max, Sampler* sampler) const override;

  OpaqueShadowType ShadowType() const override;
  int TextureFaceIndex() const override {
    // PBRT defaults to face zero when faceIndices is absent.
    return mesh->ptex_face_indices.empty() ? 0 : mesh->ptex_face_indices[face_number];
  }
  bool OpaqueHit(const Ray&, Float, Float, random_gen&) const override;

  virtual bool bounding_box(Float t0, Float t1, aabb& box) const override;
  virtual Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0) override;
  virtual Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0) override;
  
  virtual vec3f random(const point3f& origin, random_gen& rng, Float time = 0) override;
  virtual vec3f random(const point3f& origin, Sampler* sampler, Float time = 0) override;
  void GetUVs(point2f uv[3]) const;
  Float Area() const;
  Float SolidAngle(point3f p) const;
  
  virtual std::string GetName() const override {
    return(std::string("Triangle"));
  }
  size_t GetSize() override {
    return(sizeof(*this));
  }
  virtual void hitable_info_bounds(Float t0, Float t1) const override {
    aabb box;
    bounding_box(t0, t1, box);
    Rcpp::Rcout << GetName() << ": " <<  box.min() << "-" << box.max() << "\n";
  }
  TriangleMesh* mesh;
  const int *v;
  const int *n;
  const int *t;
  const int face_number;
};

#endif
