#ifndef PDFH
#define PDFH

#include "onbh.h"
#include "hitable.h"
#include "rng.h"
#include "vec2.h"
#include "microfacetdist.h"
#include "hitablelist.h"
#include "sampler.h"
#include "mathinline.h"
#include <array>

class pdf {
public: 
  virtual Float value(const vec3& direction, random_gen& rng, Float time = 0) = 0;
  virtual Float value(const vec3& direction, Sampler* sampler, Float time = 0) = 0;
  
  virtual vec3 generate(random_gen& rng, bool& diffuse_bounce, Float time = 0) = 0;
  virtual vec3 generate(Sampler* sampler, bool& diffuse_bounce, Float time = 0) = 0;
  virtual ~pdf(){};
};

class cosine_pdf : public pdf {
public:
  cosine_pdf(const vec3& w) {
    uvw.build_from_w(w);
  }
  virtual Float value(const vec3& direction, random_gen& rng, Float time = 0) {
    Float cosine = dot(unit_vector(direction), uvw.w());
    if(cosine > 0) {
      return(cosine/M_PI);
    } else {
      return(0);
    }
  } 
  virtual Float value(const vec3& direction, Sampler* sampler, Float time = 0) {
    Float cosine = dot(unit_vector(direction), uvw.w());
    if(cosine > 0) {
      return(cosine/M_PI);
    } else {
      return(0);
    }
  } 
  virtual vec3 generate(random_gen& rng, bool& diffuse_bounce, Float time = 0) {
    diffuse_bounce = true;
    return(uvw.local_to_world(rng.random_cosine_direction()));
  }
  virtual vec3 generate(Sampler* sampler, bool& diffuse_bounce, Float time = 0) {
    diffuse_bounce = true;
    return(uvw.local_to_world(rand_cosine_direction(sampler->Get2D())));
  }
  onb uvw;
};

class micro_pdf : public pdf {
public:
  micro_pdf(const vec3& w, const vec3& wi_, MicrofacetDistribution* distribution) : distribution(distribution) {
    uvw.build_from_w(w);
    wi = -unit_vector(uvw.world_to_local(wi_));;
  }
  virtual Float value(const vec3& direction, random_gen& rng, Float time = 0) {
    vec3 wo = unit_vector(uvw.world_to_local(direction));
    vec3 wh = unit_vector(wi + wo);
    return(distribution->Pdf(wo, wi, wh) / ( 4 * dot(wo, wh) ));
  }
  virtual Float value(const vec3& direction, Sampler* sampler, Float time = 0) {
    vec3 wo = unit_vector(uvw.world_to_local(direction));
    vec3 wh = unit_vector(wi + wo);
    return(distribution->Pdf(wo, wi, wh) / ( 4 * dot(wo, wh) ));
  }
  virtual vec3 generate(random_gen& rng, bool& diffuse_bounce, Float time = 0) {
    vec3 wh = distribution->Sample_wh(wi, rng.unif_rand(), rng.unif_rand());
    return(uvw.local_to_world(Reflect(wi, wh)));
  }
  virtual vec3 generate(Sampler* sampler, bool& diffuse_bounce, Float time = 0) {
    vec2 u = sampler->Get2D();
    vec3 wh = distribution->Sample_wh(wi, u.x(), u.y());
    return(uvw.local_to_world(Reflect(wi, wh)));
  }
  onb uvw;
  vec3 wi;
  MicrofacetDistribution *distribution;
};

class glossy_pdf : public pdf {
public:
  glossy_pdf(const vec3& w, const vec3& wi_, MicrofacetDistribution* distribution) : distribution(distribution) {
    uvw.build_from_w(w);
    wi = -unit_vector(uvw.world_to_local(wi_));;
  }
  virtual Float value(const vec3& direction, random_gen& rng, Float time = 0) {
    vec3 wo = unit_vector(uvw.world_to_local(direction));
    if(wo.z() * wi.z() < 0) {
      return(INFINITY);
    }
    vec3 wh = unit_vector(wi + wo);
    return(0.5f * (AbsCosTheta(wi) * M_1_PI + distribution->Pdf(wo, wi, wh) / (4 * dot(wo, wh))));
  }
  virtual Float value(const vec3& direction, Sampler* sampler, Float time = 0) {
    vec3 wo = unit_vector(uvw.world_to_local(direction));
    if(wo.z() * wi.z() < 0) {
      return(INFINITY);
    }
    vec3 wh = unit_vector(wi + wo);
    return(0.5f * (AbsCosTheta(wi) * M_1_PI + distribution->Pdf(wo, wi, wh) / (4 * dot(wo, wh))));
  }
  virtual vec3 generate(random_gen& rng, bool& diffuse_bounce, Float time = 0) {
    if(rng.unif_rand() < 0.5) {
      vec3 wh = distribution->Sample_wh(wi, rng.unif_rand(), rng.unif_rand());
      return(uvw.local_to_world(Reflect(wi, wh)));
    } else {
      diffuse_bounce = true;
      return(uvw.local_to_world(rng.random_cosine_direction()));
    }
  }
  virtual vec3 generate(Sampler* sampler, bool& diffuse_bounce, Float time = 0) {
    if(sampler->Get1D() < 0.5) {
      vec2 u = sampler->Get2D();
      vec3 wh = distribution->Sample_wh(wi, u.x(), u.y());
      return(uvw.local_to_world(Reflect(wi, wh)));
    } else {
      diffuse_bounce = true;
      return(uvw.local_to_world(rand_cosine_direction(sampler->Get2D())));
    }
  }
  onb uvw;
  vec3 wi;
  MicrofacetDistribution *distribution;
};

class hair_pdf : public pdf {
  public:
    hair_pdf(const onb uvw_, const vec3& wi_, const vec3& wo_, 
             Float eta_, Float h_, Float gammaO_, Float s_, vec3 sigma_a_,
             const Float cos2kAlpha_[3], const Float sin2kAlpha_[3], const Float v_[pMax + 1]) {
      uvw = uvw_;
      wi = wi_;
      wo = wo_;
      for (int i = 0; i < 3; ++i) {
        sin2kAlpha[i] = sin2kAlpha_[i];
        cos2kAlpha[i] = cos2kAlpha_[i];
      }
      for (int p = 0; p <= pMax; ++p) {
        v[p] = v_[p];
      }
      eta = eta_;
      h = h_;
      gammaO = gammaO_;
      s = s_;
      sigma_a = sigma_a_;
    }
    Float value(const vec3& direction, random_gen& rng, Float time = 0);
    Float value(const vec3& direction, Sampler* sampler, Float time = 0);
    
    virtual vec3 generate(random_gen& rng, bool& diffuse_bounce, Float time = 0);
    virtual vec3 generate(Sampler* sampler, bool& diffuse_bounce, Float time = 0);
    onb uvw;
    vec3 wi;
    vec3 wo;
    Float eta, h, gammaO, s;
    vec3 sigma_a;
    Float sin2kAlpha[3], cos2kAlpha[3];
  private:
    std::array<Float, pMax + 1> ComputeApPdf(Float cosThetaO) const;
    Float v[pMax + 1];
    
};

class hitable_pdf : public pdf {
public:
  hitable_pdf(hitable_list *p, const vec3& origin) : ptr(p), o(origin) {}
  virtual Float value(const vec3& direction, random_gen& rng, Float time = 0) {
    return(ptr->pdf_value(o, direction, rng, time));
  }
  virtual Float value(const vec3& direction, Sampler* sampler, Float time = 0) {
    return(ptr->pdf_value(o, direction, sampler, time));
  }
  virtual vec3 generate(random_gen& rng, bool& diffuse_bounce, Float time = 0) {
    diffuse_bounce = true;
    return(ptr->random(o, rng, time)); 
  }
  virtual vec3 generate(Sampler* sampler, bool& diffuse_bounce, Float time = 0) {
    diffuse_bounce = true;
    return(ptr->random(o, sampler, time)); 
  }
  hitable_list *ptr;
  vec3 o;
};

class mixture_pdf : public pdf {
public:
  mixture_pdf(pdf *p0, pdf *p1) {
    p[0] = p0; //Importance Sampling List
    p[1] = p1; //Surface PDF
  }
  virtual Float value(const vec3& direction, random_gen& rng, Float time = 0) {
    return(0.5 * p[0]->value(direction, rng, time) + 0.5 * p[1]->value(direction, rng, time));
  }
  virtual Float value(const vec3& direction, Sampler* sampler, Float time = 0) {
    return(0.5 * p[0]->value(direction, sampler, time) + 0.5 * p[1]->value(direction, sampler, time));
  }
  virtual vec3 generate(random_gen& rng, bool& diffuse_bounce, Float time = 0) {
    if(rng.unif_rand() < 0.5) {
      return(p[0]->generate(rng, diffuse_bounce, time));
    } else {
      return(p[1]->generate(rng, diffuse_bounce, time));
    } 
  }
  virtual vec3 generate(Sampler* sampler, bool& diffuse_bounce, Float time = 0) {
    if(sampler->Get1D() < 0.5) {
      return(p[0]->generate(sampler, diffuse_bounce, time));
    } else {
      return(p[1]->generate(sampler, diffuse_bounce, time));
    } 
  }
  pdf *p[2];
};

#endif
