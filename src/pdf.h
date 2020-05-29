#ifndef PDFH
#define PDFH

#include "onbh.h"
#include "hitable.h"
#include "rng.h"
#include "vec2.h"
#include "microfacetdist.h"
#include "hitablelist.h"

class pdf {
public: 
  virtual Float value(const vec3& direction, random_gen& rng) = 0;
  virtual vec3 generate(random_gen& rng) = 0;
  virtual ~pdf(){};
};

class cosine_pdf : public pdf {
public:
  cosine_pdf(const vec3& w) {
    uvw.build_from_w(w);
  }
  virtual Float value(const vec3& direction, random_gen& rng) {
    Float cosine = dot(unit_vector(direction), uvw.w());
    if(cosine > 0) {
      return(cosine/M_PI);
    } else {
      return(0);
    }
  } 
  virtual vec3 generate(random_gen& rng) {
    return(uvw.local_to_world(rng.random_cosine_direction()));
  }
  onb uvw;
};

class micro_pdf : public pdf {
public:
  micro_pdf(const vec3& w, const vec3& wi_, MicrofacetDistribution* distribution) : distribution(distribution) {
    uvw.build_from_w(w);
    wi = -unit_vector(uvw.world_to_local(wi_));;
  }
  virtual Float value(const vec3& direction, random_gen& rng) {
    vec3 wo = unit_vector(uvw.world_to_local(direction));
    vec3 wh = unit_vector(wi + wo);
    return(distribution->Pdf(wo, wi, wh) / ( 4 * dot(wo, wh) ));
  }
  virtual vec3 generate(random_gen& rng) {
    vec3 wh = distribution->Sample_wh(wi, rng.unif_rand(), rng.unif_rand());
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
  virtual Float value(const vec3& direction, random_gen& rng) {
    vec3 wo = unit_vector(uvw.world_to_local(direction));
    vec3 wh = unit_vector(wi + wo);
    return(0.5f  * M_1_PI);// + distribution->D(wh) * distribution->G(wo,wi,wh) * AbsDot(wo, wh));
  }
  virtual vec3 generate(random_gen& rng) {
    vec3 wh = distribution->Sample_wh(wi, rng.unif_rand(), rng.unif_rand());
    return(uvw.local_to_world(Reflect(wi, wh)));
  }
  onb uvw;
  vec3 wi;
  MicrofacetDistribution *distribution;
};

class hitable_pdf : public pdf {
public:
  hitable_pdf(hitable_list *p, const vec3& origin) : ptr(p), o(origin) {}
  virtual Float value(const vec3& direction, random_gen& rng) {
    return(ptr->pdf_value(o, direction, rng));
  }
  virtual vec3 generate(random_gen& rng) {
    return(ptr->random(o, rng)); 
  }
  hitable_list *ptr;
  vec3 o;
};

class mixture_pdf : public pdf {
public:
  mixture_pdf(pdf *p0, pdf *p1) {
    p[0] = p0;
    p[1] = p1;
  }
  virtual Float value(const vec3& direction, random_gen& rng) {
    return(0.5 * p[0]->value(direction, rng) + 0.5 * p[1]->value(direction, rng));
  }
  virtual vec3 generate(random_gen& rng) {
    if(rng.unif_rand() < 0.5) {
      return(p[0]->generate(rng));
    } else {
      return(p[1]->generate(rng));
    } 
  }
  pdf *p[2];
};

#endif
