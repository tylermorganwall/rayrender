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
  virtual Float value(const vec3f& direction, random_gen& rng, Float time = 0) = 0;
  virtual Float value(const vec3f& direction, Sampler* sampler, Float time = 0) = 0;
  
  virtual vec3f generate(random_gen& rng, bool& diffuse_bounce, Float time = 0) = 0;
  virtual vec3f generate(Sampler* sampler, bool& diffuse_bounce, Float time = 0) = 0;
  virtual ~pdf(){};
};

class cosine_pdf : public pdf {
public:
  cosine_pdf(const normal3f& w) {
    uvw.build_from_w(w);
  }
  virtual Float value(const vec3f& direction, random_gen& rng, Float time = 0);
  virtual Float value(const vec3f& direction, Sampler* sampler, Float time = 0);
  virtual vec3f generate(random_gen& rng, bool& diffuse_bounce, Float time = 0);
  virtual vec3f generate(Sampler* sampler, bool& diffuse_bounce, Float time = 0);
  onb uvw;
};

class micro_pdf : public pdf {
public:
  micro_pdf(const normal3f& w, const vec3f& wi_, MicrofacetDistribution* distribution, 
            Float uu, Float vv) : distribution(distribution),  u(uu), v(vv) {
    uvw.build_from_w(w);
    wi = -unit_vector(uvw.world_to_local(wi_));;
  }
  virtual Float value(const vec3f& direction, random_gen& rng, Float time = 0);
  virtual Float value(const vec3f& direction, Sampler* sampler, Float time = 0);
  virtual vec3f generate(random_gen& rng, bool& diffuse_bounce, Float time = 0);
  virtual vec3f generate(Sampler* sampler, bool& diffuse_bounce, Float time = 0);
  onb uvw;
  vec3f wi;
  MicrofacetDistribution *distribution;
  Float u,v;
  
};

class micro_transmission_pdf : public pdf {
public:
  micro_transmission_pdf(const normal3f& w, const vec3f& wi_, MicrofacetDistribution* distribution,
                         Float eta, Float uu, Float vv);
  virtual Float value(const vec3f& direction, random_gen& rng, Float time = 0);
  virtual Float value(const vec3f& direction, Sampler* sampler, Float time = 0);
  virtual vec3f generate(random_gen& rng, bool& diffuse_bounce, Float time = 0);
  virtual vec3f generate(Sampler* sampler, bool& diffuse_bounce, Float time = 0);
  onb uvw;
  vec3f wi;
  Float eta;
  MicrofacetDistribution *distribution;
  Float u,v;
};

class glossy_pdf : public pdf {
public:
  glossy_pdf(const normal3f& w, const vec3f& wi_, MicrofacetDistribution* distribution, 
             Float uu, Float vv) : distribution(distribution), u(uu), v(vv) {
    uvw.build_from_w(w);
    wi = -unit_vector(uvw.world_to_local(wi_));;
  }
  virtual Float value(const vec3f& direction, random_gen& rng, Float time = 0);
  virtual Float value(const vec3f& direction, Sampler* sampler, Float time = 0);
  virtual vec3f generate(random_gen& rng, bool& diffuse_bounce, Float time = 0);
  virtual vec3f generate(Sampler* sampler, bool& diffuse_bounce, Float time = 0);
  onb uvw;
  vec3f wi;
  MicrofacetDistribution *distribution;
  Float u,v;
  
};

class hair_pdf : public pdf {
  public:
    hair_pdf(const onb uvw_, const vec3f& wi_, const vec3f& wo_, 
             Float eta_, Float h_, Float gammaO_, Float s_, point3f sigma_a_,
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
    Float value(const vec3f& direction, random_gen& rng, Float time = 0);
    Float value(const vec3f& direction, Sampler* sampler, Float time = 0);
    
    virtual vec3f generate(random_gen& rng, bool& diffuse_bounce, Float time = 0);
    virtual vec3f generate(Sampler* sampler, bool& diffuse_bounce, Float time = 0);
    onb uvw;
    vec3f wi;
    vec3f wo;
    Float eta, h, gammaO, s;
    point3f sigma_a;
    Float sin2kAlpha[3], cos2kAlpha[3];
  private:
    std::array<Float, pMax + 1> ComputeApPdf(Float cosThetaO) const;
    Float v[pMax + 1];
    
};

class hitable_pdf : public pdf {
public:
  hitable_pdf(hitable_list *p, const point3f& origin) : ptr(p), o(origin) {}
  virtual Float value(const vec3f& direction, random_gen& rng, Float time = 0);
  virtual Float value(const vec3f& direction, Sampler* sampler, Float time = 0);
  virtual vec3f generate(random_gen& rng, bool& diffuse_bounce, Float time = 0);
  virtual vec3f generate(Sampler* sampler, bool& diffuse_bounce, Float time = 0);
  hitable_list *ptr;
  point3f o;
};

class mixture_pdf : public pdf {
public:
  mixture_pdf(pdf *p0, pdf *p1) {
    p[0] = p0; //Importance Sampling List
    p[1] = p1; //Surface PDF
  }
  virtual Float value(const vec3f& direction, random_gen& rng, Float time = 0);
  virtual Float value(const vec3f& direction, Sampler* sampler, Float time = 0);
  virtual vec3f generate(random_gen& rng, bool& diffuse_bounce, Float time = 0);
  virtual vec3f generate(Sampler* sampler, bool& diffuse_bounce, Float time = 0);
  pdf *p[2];
};

#endif
