#ifndef TRIANGLEH
#define TRIANGLEH

#include "hitable.h"
#include "material.h"
#include "onbh.h"
#include "trianglemesh.h"

class triangle : public hitable {
public:
  triangle() : face_number(0) {}
  triangle(TriangleMesh* mesh, const int *v, const int *n, const int *t, const int face_number,
           std::shared_ptr<Transform> ObjectToWorld, std::shared_ptr<Transform> WorldToObject, bool reverseOrientation) : 
    hitable(ObjectToWorld, WorldToObject, reverseOrientation), mesh(mesh), v(v), n(n), t(t), face_number(face_number) {
    
  }
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler);
  
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  virtual Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0);
  virtual Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0);
  
  virtual vec3f random(const point3f& origin, random_gen& rng, Float time = 0);
  virtual vec3f random(const point3f& origin, Sampler* sampler, Float time = 0);
  void GetUVs(point2f uv[3]) const;
  Float Area() const;
  Float SolidAngle(point3f p) const;
  
  virtual std::string GetName() const {
    return(std::string("Triangle"));
  }
  size_t GetSize()  {
    return(sizeof(*this));
  }
  TriangleMesh* mesh;
  const int *v;
  const int *n;
  const int *t;
  const int face_number;
};

#endif
