#ifndef CURVEH
#define CURVEH

#include "hitable.h"
#include "material.h"
#include "mathinline.h"

enum class CurveType { Flat, Cylinder, Ribbon };

struct CurveCommon {
  CurveCommon(const vec3f c[4], Float w0, Float w1, CurveType type,
              const vec3f *norm);
  const CurveType type;
  point3f cpObj[4];
  Float width[2];
  normal3f n[2];
  Float normalAngle, invSinNormalAngle;
};

class curve: public hitable {
  public:
    curve() : uMin(0), uMax(0) {}
    ~curve() {
      // delete mat_ptr;
    }
    curve(Float uMin, Float uMax, 
          const std::shared_ptr<CurveCommon> common, std::shared_ptr<material> mat,
          std::shared_ptr<Transform> ObjectToWorld, std::shared_ptr<Transform> WorldToObject, bool reverseOrientation) : 
      hitable(ObjectToWorld, WorldToObject, reverseOrientation), 
      mat_ptr(mat), common(common), uMin(uMin), uMax(uMax) {};
    virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng);
    virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, Sampler* sampler);
    
    virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
    virtual Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0);
    virtual Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0);
    
    virtual vec3f random(const point3f& o, random_gen& rng, Float time = 0);
    virtual vec3f random(const point3f& o, Sampler* sampler, Float time = 0);
    virtual std::string GetName() const {
      return(std::string("Curve"));
    }
    std::shared_ptr<material> mat_ptr;
  private:
    bool recursiveIntersect(const ray& r, Float tmin, Float tmax, hit_record& rec,
                            const point3f cp[4], Float u0, Float u1, int depth,
                            const Transform &rayToObject) const;

    const std::shared_ptr<CurveCommon> common;
    Float uMin, uMax;
};


#endif
