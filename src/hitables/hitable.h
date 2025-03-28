#ifndef HITABLEH
#define HITABLEH

#include "../math/aabb.h"
#include "../core/ray.h"
#include "../math/vectypes.h"
#include "../materials/texture.h"
#include "../math/rng.h"
#include "../math/sampler.h"
#include <Rcpp.h>
#include "../math/transform.h"
#include "../math/animatedtransform.h"
#include <memory>
#include <cfloat>
#include <tuple>

class material;
class hitable;

void get_sphere_uv(const vec3f& p, Float& u, Float& v);
void get_sphere_uv(const normal3f& p, Float& u, Float& v);



struct alignas(16) hit_record {
  hit_record() : has_bump(false), alpha_miss(false), infinite_area_hit(false) {};
  point3f p; //PBRT: In Interaction
  Float t; //PBRT: In Interaction
  normal3f normal; //PBRT: In interaction

#ifdef DEBUGBVH
  Float bvh_nodes;
#endif
  vec3f dpdu, dpdv; //PBRT: In SurfaceInteraction
  vec3f pError; //PBRT: In Interaction
  Float u; //PBRT: In SurfaceInteraction
  Float v; //PBRT: In SurfaceInteraction
  bool has_bump; 
  bool alpha_miss;
  bool infinite_area_hit;
  normal3f bump_normal; 

  // vec3f wo; //PBRT: In Interaction, negative ray direction
  const hitable* shape = nullptr; //PBRT: In SurfaceInteraction, const Shape *shape
  material* mat_ptr; //PBRT: In SurfaceInteraction as bsdf or bssrdf

  // mutable vec3f dpdx, dpdy;
  // mutable normal3f dndu, dndv;
  // mutable Float dudx, dvdx, dudy, dvdy;
  //const Shape *shape (recording the shape)
  //const Primitive *primitive (recording the primitive)
  //int faceIndex (for ptex lookups)
};

class hitable {
  public:
    hitable() : reverseOrientation(false), transformSwapsHandedness(false) {}
    hitable(Transform* ObjectToWorld, 
            Transform* WorldToObject, 
            std::shared_ptr<material> mat_ptr,
            bool reverseOrientation) : 
      ObjectToWorld(ObjectToWorld), 
      WorldToObject(WorldToObject), 
      mat_ptr(mat_ptr), 
      reverseOrientation(reverseOrientation),
      transformSwapsHandedness(ObjectToWorld->SwapsHandedness()) {}
    virtual const bool hit(const Ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const = 0;
    virtual const bool hit(const Ray& r, Float tmin, Float tmax, hit_record& rec, Sampler* sampler) const = 0;
    virtual bool HitP(const Ray &r, Float t_min, Float t_max, random_gen& rng) const {
      hit_record tmp;
      return hit(r, t_min, t_max, tmp, rng);
    }
    virtual bool HitP(const Ray &r, Float t_min, Float t_max, Sampler* sampler) const {
      hit_record tmp;
      return hit(r, t_min, t_max, tmp, sampler);
    }


    // virtual const bool hit(const CompactRay& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const = 0;
    // virtual const bool hit(const CompactRay& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const = 0;

    virtual bool bounding_box(Float t0, Float t1, aabb& box) const = 0;
    virtual Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0) {
      return(0.0);
    }
    virtual Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0) {
      return(0.0);
    }
    virtual vec3f random(const point3f& o, random_gen& rng, Float time = 0) {
      return(vec3f(0,1,0));
    }
    virtual vec3f random(const point3f& o, Sampler* sampler, Float time = 0) {
      return(vec3f(0,1,0));
    }
    virtual std::string GetName() const = 0;

    virtual size_t GetSize() = 0;
    virtual std::pair<size_t,size_t> CountNodeLeaf() {
      return(std::pair<size_t,size_t>(0,1));
    }

    virtual void hitable_info_bounds(Float t0, Float t1) const = 0;
    
    virtual ~hitable() {}
    Transform* ObjectToWorld;
    Transform* WorldToObject;
    const std::shared_ptr<material> mat_ptr;
    const bool reverseOrientation;
    const bool transformSwapsHandedness;
};


class AnimatedHitable: public hitable {
public:
  AnimatedHitable(std::shared_ptr<hitable> &primitive,
                  const AnimatedTransform &PrimitiveToWorld)
    : primitive(primitive), PrimitiveToWorld(PrimitiveToWorld) {} 
  ~AnimatedHitable() {}
  const bool hit(const Ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const;
  const bool hit(const Ray& r, Float tmin, Float tmax, hit_record& rec, Sampler* sampler) const;
  virtual bool HitP(const Ray &r, Float t_min, Float t_max, random_gen& rng) const {
    hit_record tmp;
    return hit(r, t_min, t_max, tmp, rng);
  }
  virtual bool HitP(const Ray &r, Float t_min, Float t_max, Sampler* sampler) const {
    hit_record tmp;
    return hit(r, t_min, t_max, tmp, sampler);
  }

  Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0);
  Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0);
  vec3f random(const point3f& o, random_gen& rng, Float time = 0);
  vec3f random(const point3f& o, Sampler* sampler, Float time = 0);
  std::string GetName() const;
  size_t GetSize() {
    return(sizeof(*this));
  }
  virtual void hitable_info_bounds(Float t0, Float t1) const {
    aabb box;
    bounding_box(t0, t1, box);
    Rcpp::Rcout << GetName() << ": " <<  box.min() << "-" << box.max() << "\n";
  }
  bool bounding_box(Float t0, Float t1, aabb& box) const;
  std::shared_ptr<hitable> primitive;
  const AnimatedTransform PrimitiveToWorld;
};



#endif
