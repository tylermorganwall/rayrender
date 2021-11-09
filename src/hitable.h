#ifndef HITABLEH
#define HITABLEH

#include "aabb.h"
#include "vec3.h"
#include "texture.h"
#include "rng.h"
#include "sampler.h"
#include <Rcpp.h>
#include "transform.h"
#include "animatedtransform.h"
#include <memory>
#include <cfloat>

class material;

void get_sphere_uv(const vec3f& p, Float& u, Float& v);
void get_sphere_uv(const normal3f& p, Float& u, Float& v);

struct hit_record;

class hitable {
  public:
    hitable() : reverseOrientation(false), transformSwapsHandedness(false) {}
    hitable(std::shared_ptr<Transform> ObjectToWorld, std::shared_ptr<Transform> WorldToObject, bool reverseOrientation) : 
      ObjectToWorld(ObjectToWorld), WorldToObject(WorldToObject), reverseOrientation(reverseOrientation),
      transformSwapsHandedness(ObjectToWorld->SwapsHandedness()) {}
    virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) = 0;
    virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, Sampler* sampler) = 0;
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
    virtual std::string GetName() const {
      return(std::string("Hitable"));
    }
    virtual ~hitable() {}
    const std::shared_ptr<Transform> ObjectToWorld, WorldToObject;
    const bool reverseOrientation;
    const bool transformSwapsHandedness;
};


struct hit_record {
  hit_record() : has_bump(false), alpha_miss(false) {};
  
  point3f p; //PBRT: In Interaction
  Float t; //PBRT: In Interaction
  Float u; //PBRT: In SurfaceInteraction
  Float v; //PBRT: In SurfaceInteraction
#ifdef DEBUGBVH
  Float bvh_nodes;
#endif
  normal3f normal; //PBRT: In interaction
  vec3f dpdu, dpdv; //PBRT: In SurfaceInteraction
  vec3f pError; //PBRT: In Interaction
  vec3f wo; //PBRT: In Interaction, negative ray direction
  normal3f bump_normal; 
  bool has_bump; 
  const hitable* shape = nullptr; //PBRT: In SurfaceInteraction, const Shape *shape
  material* mat_ptr; //PBRT: In SurfaceInteraction as bsdf or bssrdf
  bool alpha_miss;
  //Missing from PBRT: 
  //mutable Float dudx, dvdx, dudy, dvdy (for texture sampling)
  //mutable vec3 dpdx, dpdy
  //const Shape *shape (recording the shape)
  //const Primitive *primitive (recording the primitive)
  //int faceIndex (for ptex lookups)
};


class AnimatedHitable: public hitable {
public:
  AnimatedHitable(std::shared_ptr<hitable> &primitive,
                  const AnimatedTransform &PrimitiveToWorld)
    : primitive(primitive), PrimitiveToWorld(PrimitiveToWorld) {} 
  ~AnimatedHitable() {}
  bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, Sampler* sampler);
  Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0);
  Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0);
  vec3f random(const point3f& o, random_gen& rng, Float time = 0);
  vec3f random(const point3f& o, Sampler* sampler, Float time = 0);
  std::string GetName() const;
  bool bounding_box(Float t0, Float t1, aabb& box) const;
  std::shared_ptr<hitable> primitive;
  const AnimatedTransform PrimitiveToWorld;
};



#endif
