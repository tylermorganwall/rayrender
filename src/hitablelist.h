#ifndef HITABLELISTH
#define HITABLELISTH

#include "hitable.h"
#include "sampler.h"
#include <memory>

class hitable_list: public hitable {
  public:
    hitable_list() {}
    hitable_list(std::shared_ptr<hitable> object) {add(object);}
    const bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng) const;
    const bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, Sampler* sampler) const;
    virtual bool HitP(const ray &r, Float t_min, Float t_max, random_gen& rng) const;
    virtual bool HitP(const ray &r, Float t_min, Float t_max, Sampler* sampler) const;

    bool bounding_box(Float t0, Float t1, aabb& box) const;
    Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0);
    Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0);
    vec3f random(const point3f& o, random_gen& rng, Float time = 0);
    vec3f random(const point3f& o, Sampler* sampler, Float time = 0);
    void add(std::shared_ptr<hitable> object) { objects.push_back(object); }
    std::shared_ptr<hitable> back() {return(objects.back());}
    void validate() const;
    virtual void hitable_info_bounds(Float t0, Float t1) const {
      aabb box_top;
      bounding_box(t0, t1, box_top);
      Rcpp::Rcout << GetName() << ": " <<  box_top.min() << "-" << box_top.max() << "\n";
      for(size_t i = 0; i < objects.size(); i++) {
        aabb box;
        objects[i]->bounding_box(t0, t1, box);
        Rcpp::Rcout << "   " << objects[i]->GetName() << ": " <<  box.min() << "-" << box.max() << "\n";
      }
    }
    
    size_t size() {return(objects.size());}
    std::vector<std::shared_ptr<hitable>> objects;
    std::string GetName() const;
    size_t GetSize();
};

#endif
