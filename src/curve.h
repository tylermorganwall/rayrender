#ifndef CURVEH
#define CURVEH

#include "hitable.h"
#include "material.h"
#include "mathinline.h"

enum class CurveType { 
  Flat = 1, 
  Cylinder = 2, 
  Ribbon = 3 
};

struct CurveCommon {
  CurveCommon(const point3f c[4], Float w0, Float w1, CurveType type,
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
          Transform* ObjectToWorld, Transform* WorldToObject, bool reverseOrientation) : 
      hitable(ObjectToWorld, WorldToObject, mat, reverseOrientation), 
      common(common), uMin(uMin), uMax(uMax) {};
    virtual const bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng) const;
    virtual const bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, Sampler* sampler) const;
    // virtual bool HitP(const ray &r, Float t_min, Float t_max, random_gen& rng) const;
    // virtual bool HitP(const ray &r, Float t_min, Float t_max, Sampler* sampler) const;

    virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
    virtual Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0);
    virtual Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0);
    
    virtual vec3f random(const point3f& o, random_gen& rng, Float time = 0);
    virtual vec3f random(const point3f& o, Sampler* sampler, Float time = 0);
    virtual std::string GetName() const {
      return(std::string("Curve"));
    }
    size_t GetSize()  {
      return(sizeof(*this) + sizeof(*common));
    }
    virtual void hitable_info_bounds(Float t0, Float t1) const {
      aabb box;
      bounding_box(t0, t1, box);
      Rcpp::Rcout << GetName() << ": " <<  box.min() << "-" << box.max() << "\n";
    }
  private:
    bool recursiveIntersect(const ray& r, Float tmin, Float tmax, hit_record& rec,
                            const point3f cp[4], Float u0, Float u1, int depth,
                            const Transform &rayToObject) const;

    const std::shared_ptr<CurveCommon> common;
    Float uMin, uMax;
};


#endif
