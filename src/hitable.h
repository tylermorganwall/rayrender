#ifndef HITABLEH
#define HITABLEH

#include "aabb.h"
#include "vec3.h"
#include "texture.h"
#include "rng.h"
#include "sampler.h"
#include <Rcpp.h>
#include "transform.h"
#include <memory>

class material;
class hitable;


void get_sphere_uv(const vec3f& p, Float& u, Float& v);

struct hit_record {
  point3f p; //PBRT: In Interaction
  Float t; //PBRT: In Interaction
  Float u; //PBRT: In SurfaceInteraction
  Float v; //PBRT: In SurfaceInteraction
#ifdef DEBUGBVH
  Float bvh_nodes;
#endif
  vec3f normal; //PBRT: In interaction
  vec3f dpdu, dpdv; //PBRT: In SurfaceInteraction
  vec3f pError; //PBRT: In Interaction
  vec3f wo; //PBRT: In Interaction, negative ray direction
  vec3f bump_normal; 
  bool has_bump; 
  const hitable* hitable = nullptr; //PBRT: In SurfaceInteraction, const Shape *shape
  material* mat_ptr; //PBRT: In SurfaceInteraction as bsdf or bssrdf
  //Missing from PBRT: 
  //mutable Float dudx, dvdx, dudy, dvdy (for texture sampling)
  //mutable vec3 dpdx, dpdy
  //const Shape *shape (recording the shape)
  //const Primitive *primitive (recording the primitive)
  //int faceIndex (for ptex lookups)
  
};

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
    virtual ~hitable() {}
    const std::shared_ptr<Transform> ObjectToWorld, WorldToObject;
    const bool reverseOrientation;
    const bool transformSwapsHandedness;
};

// class flip_normals : public hitable {
//   public:
//     flip_normals(std::shared_ptr<hitable> p) : ptr(p) {}
//     ~flip_normals() {}
//     virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
//       if(ptr->hit(r, t_min, t_max, rec, rng)) {
//         return(true);
//       } else {
//         return(false);
//       }
//     }
//     virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
//       if(ptr->hit(r, t_min, t_max, rec, sampler)) {
//         return(true);
//       } else {
//         return(false);
//       }
//     }
//     virtual bool bounding_box(Float t0, Float t1, aabb& box) const {
//       return(ptr->bounding_box(t0,t1,box));
//     }
//     Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0) {
//       return(ptr->pdf_value(o,v, rng, time));
//     }
//     Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0) {
//       return(ptr->pdf_value(o,v, sampler, time));
//     }
//     vec3f random(const point3f& o, random_gen& rng, Float time = 0) {
//       return(ptr->random(o, rng, time));
//     }
//     vec3f random(const point3f& o, Sampler* sampler, Float time = 0) {
//       return(ptr->random(o, sampler, time));
//     }
//     
//     std::shared_ptr<hitable> ptr;
// };
// 
// 
// 
// class translate : public hitable {
// public: 
//   translate(std::shared_ptr<hitable> p, const vec3f& displacement) : ptr(p), offset(displacement) {}
//   ~translate() {}
//   virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
//   virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler);
//   
//   virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
//   Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0) {
//     return(ptr->pdf_value(o - offset, v, rng, time));
//   }
//   Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0) {
//     return(ptr->pdf_value(o - offset, v, sampler, time));
//   }
//   vec3f random(const point3f& o, random_gen& rng, Float time = 0) {
//     return(ptr->random(o - offset, rng, time));
//   }
//   vec3f random(const point3f& o, Sampler* sampler, Float time = 0) {
//     return(ptr->random(o - offset, sampler, time));
//   }
//   std::shared_ptr<hitable> ptr;
//   vec3f offset;
// };
// 
// 
// class scale : public hitable {
// public: 
//   scale(std::shared_ptr<hitable> p, const vec3f& scale_factor) : ptr(p), scale_factor(scale_factor) {
//     inv_scale.e[0] = 1 / scale_factor.x();
//     inv_scale.e[1] = 1 / scale_factor.y();
//     inv_scale.e[2] = 1 / scale_factor.z();
//   }
//   ~scale() {}
//   virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
//   virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler);
//   
//   virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
//   Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0) {
//     return(ptr->pdf_value(o * inv_scale, v, rng, time));
//   }
//   Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0) {
//     return(ptr->pdf_value(o * inv_scale, v, sampler, time));
//   }
//   vec3f random(const point3f& o, random_gen& rng, Float time = 0) {
//     return(ptr->random(o * inv_scale, rng, time));
//   }
//   vec3f random(const point3f& o, Sampler* sampler, Float time = 0) {
//     return(ptr->random(o * inv_scale, sampler, time));
//   }
//   std::shared_ptr<hitable> ptr;
//   vec3f scale_factor;
//   vec3f inv_scale;
// };
// 
// class rotate_y : public hitable {
// public:
//   rotate_y(std::shared_ptr<hitable> p, Float angle);
//   ~rotate_y() {}
//   virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
//   virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler);
//   virtual bool bounding_box(Float t0, Float t1, aabb& box) const {
//     box = bbox; 
//     return(hasbox);
//   }
//   Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0) {
//     vec3f v2 = v;
//     v2.e[0] = cos_theta*v.x() - sin_theta*v.z();
//     v2.e[2] = sin_theta*v.x() + cos_theta*v.z();
//     vec3f o2 = o;
//     o2.e[0] = cos_theta*o.x() - sin_theta*o.z();
//     o2.e[2] = sin_theta*o.x() + cos_theta*o.z();
//     return(ptr->pdf_value(o2,v2, rng, time));
//   }
//   Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0) {
//     vec3f v2 = v;
//     v2.e[0] = cos_theta*v.x() - sin_theta*v.z();
//     v2.e[2] = sin_theta*v.x() + cos_theta*v.z();
//     vec3f o2 = o;
//     o2.e[0] = cos_theta*o.x() - sin_theta*o.z();
//     o2.e[2] = sin_theta*o.x() + cos_theta*o.z();
//     return(ptr->pdf_value(o2,v2, sampler, time));
//   }
//   vec3f random(const point3f& o, random_gen& rng, Float time = 0) {
//     vec3f o2 = o;
//     o2.e[0] = cos_theta*o.x() - sin_theta*o.z();
//     o2.e[2] = sin_theta*o.x() + cos_theta*o.z();
//     vec3f temp = ptr->random(o2, rng, time);
//     vec3f temp2 = temp;
//     temp2.e[0] = cos_theta*temp.x() + sin_theta*temp.z();
//     temp2.e[2] = -sin_theta*temp.x() + cos_theta*temp.z(); 
//     return(temp2);
//   }
//   vec3f random(const point3f& o, Sampler* sampler, Float time = 0) {
//     vec3f o2 = o;
//     o2.e[0] = cos_theta*o.x() - sin_theta*o.z();
//     o2.e[2] = sin_theta*o.x() + cos_theta*o.z();
//     vec3f temp = ptr->random(o2, sampler, time);
//     vec3f temp2 = temp;
//     temp2.e[0] = cos_theta*temp.x() + sin_theta*temp.z();
//     temp2.e[2] = -sin_theta*temp.x() + cos_theta*temp.z(); 
//     return(temp2);
//   }
//   std::shared_ptr<hitable> ptr;
//   Float sin_theta;
//   Float cos_theta;
//   bool hasbox;
//   aabb bbox;
// };
// 
// class rotate_x : public hitable {
// public:
//   rotate_x(std::shared_ptr<hitable> p, Float angle);
//   ~rotate_x() {}
//   virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
//   virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler);
//   
//   virtual bool bounding_box(Float t0, Float t1, aabb& box) const {
//     box = bbox; 
//     return(hasbox);
//   }
//   Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0) {
//     vec3f v2 = v;
//     v2.e[1] = cos_theta*v.y() - sin_theta*v.z();
//     v2.e[2] = sin_theta*v.y() + cos_theta*v.z();
//     vec3f o2 = o;
//     o2.e[1] = cos_theta*o.y() - sin_theta*o.z();
//     o2.e[2] = sin_theta*o.y() + cos_theta*o.z();
//     return(ptr->pdf_value(o2,v2, rng, time));
//   }
//   Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0) {
//     vec3f v2 = v;
//     v2.e[1] = cos_theta*v.y() - sin_theta*v.z();
//     v2.e[2] = sin_theta*v.y() + cos_theta*v.z();
//     vec3f o2 = o;
//     o2.e[1] = cos_theta*o.y() - sin_theta*o.z();
//     o2.e[2] = sin_theta*o.y() + cos_theta*o.z();
//     return(ptr->pdf_value(o2,v2, sampler, time));
//   }
//   vec3f random(const point3f& o, random_gen& rng, Float time = 0) {
//     vec3f o2 = o;
//     o2.e[1] = cos_theta*o.y() - sin_theta*o.z();
//     o2.e[2] = sin_theta*o.y() + cos_theta*o.z();
//     vec3f temp = ptr->random(o2, rng, time);
//     vec3f temp2 = temp;
//     temp2.e[1] = cos_theta*temp.y() + sin_theta*temp.z();
//     temp2.e[2] = -sin_theta*temp.y() + cos_theta*temp.z(); 
//     return(temp2);
//   }
//   vec3f random(const point3f& o, Sampler* sampler, Float time = 0) {
//     vec3f o2 = o;
//     o2.e[1] = cos_theta*o.y() - sin_theta*o.z();
//     o2.e[2] = sin_theta*o.y() + cos_theta*o.z();
//     vec3f temp = ptr->random(o2, sampler, time);
//     vec3f temp2 = temp;
//     temp2.e[1] = cos_theta*temp.y() + sin_theta*temp.z();
//     temp2.e[2] = -sin_theta*temp.y() + cos_theta*temp.z(); 
//     return(temp2);
//   }
//   std::shared_ptr<hitable> ptr;
//   Float sin_theta;
//   Float cos_theta;
//   bool hasbox;
//   aabb bbox;
// };
// 
// class rotate_z : public hitable {
// public:
//   rotate_z(std::shared_ptr<hitable> p, Float angle);
//   ~rotate_z() {}
//   virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
//   virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler);
//   
//   virtual bool bounding_box(Float t0, Float t1, aabb& box) const {
//     box = bbox; 
//     return(hasbox);
//   }
//   Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0) {
//     vec3f v2 = v;
//     v2.e[0] = cos_theta*v.x() - sin_theta*v.y();
//     v2.e[1] = sin_theta*v.x() + cos_theta*v.y();
//     vec3f o2 = o;
//     o2.e[0] = cos_theta*o.x() - sin_theta*o.y();
//     o2.e[1] = sin_theta*o.x() + cos_theta*o.y();
//     return(ptr->pdf_value(o2,v2, rng, time));
//   }
//   Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0) {
//     vec3f v2 = v;
//     v2.e[0] = cos_theta*v.x() - sin_theta*v.y();
//     v2.e[1] = sin_theta*v.x() + cos_theta*v.y();
//     vec3f o2 = o;
//     o2.e[0] = cos_theta*o.x() - sin_theta*o.y();
//     o2.e[1] = sin_theta*o.x() + cos_theta*o.y();
//     return(ptr->pdf_value(o2,v2, sampler, time));
//   }
//   vec3f random(const point3f& o, random_gen& rng, Float time = 0) {
//     vec3f o2 = o;
//     o2.e[0] = cos_theta*o.x() - sin_theta*o.y();
//     o2.e[1] = sin_theta*o.x() + cos_theta*o.y();
//     vec3f temp = ptr->random(o2, rng, time);
//     vec3f temp2 = temp;
//     temp2.e[0] = cos_theta*temp.x() + sin_theta*temp.y();
//     temp2.e[1] = -sin_theta*temp.x() + cos_theta*temp.y(); 
//     return(temp2);
//   }
//   vec3f random(const point3f& o, Sampler* sampler, Float time = 0) {
//     vec3f o2 = o;
//     o2.e[0] = cos_theta*o.x() - sin_theta*o.y();
//     o2.e[1] = sin_theta*o.x() + cos_theta*o.y();
//     vec3f temp = ptr->random(o2, sampler, time);
//     vec3f temp2 = temp;
//     temp2.e[0] = cos_theta*temp.x() + sin_theta*temp.y();
//     temp2.e[1] = -sin_theta*temp.x() + cos_theta*temp.y(); 
//     return(temp2);
//   }
//   std::shared_ptr<hitable> ptr;
//   Float sin_theta;
//   Float cos_theta;
//   bool hasbox;
//   aabb bbox;
// };

#endif
