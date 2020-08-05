#ifndef CURVEH
#define CURVEH

#include "hitable.h"
#include "material.h"
#include "mathinline.h"

enum class CurveType { Flat, Cylinder, Ribbon };

struct CurveCommon {
  CurveCommon(const vec3 c[4], Float w0, Float w1, CurveType type,
              const vec3 *norm);
  const CurveType type;
  vec3 cpObj[4];
  Float width[2];
  vec3 n[2];
  Float normalAngle, invSinNormalAngle;
};

class curve: public hitable {
  public:
    curve() : uMin(0), uMax(0) {}
    ~curve() {
      delete mat_ptr;
      delete common;
    }
    curve(Float uMin, Float uMax, 
          const CurveCommon* common, material *mat) : 
      mat_ptr(mat), common(common), uMin(uMin), uMax(uMax) {};
    virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng);
    virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
    virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng);
    virtual vec3 random(const vec3& o, random_gen& rng);
    virtual vec3 random(const vec3& o, Sampler* sampler);
    material *mat_ptr;
  private:
    bool recursiveIntersect(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng,
                            const vec3 cp[4], Float u0, Float u1, int depth) const;

    const CurveCommon* common;
    Float uMin, uMax;
};


#endif
