#ifndef PDFH
#define PDFH

#include "onbh.h"
#include "hitable.h"
#include "rng.h"
#include "vec2.h"
#include "microfacetdist.h"

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

class micro_trow_pdf : public pdf {
public:
  micro_trow_pdf(const vec3& w, const vec3& wi_, MicrofacetDistribution* distribution) : distribution(distribution) {
    alphas = distribution->GetAlphas();
    uvw.build_from_w(w);
    wi = -unit_vector(uvw.world_to_local(wi_));;
  }
  virtual Float value(const vec3& direction, random_gen& rng) {
    vec3 wo = unit_vector(uvw.world_to_local(direction));
    if (!SameHemisphere(wo, wi)) {
      return(INFINITY);
    }
    vec3 wh = unit_vector(wi + wo);
    return(distribution->D(wh) * distribution->G(wo, wi, wh) * std::fabs(dot(wo, wh)) / AbsCosTheta(wo));
  }
  virtual vec3 generate(random_gen& rng) {
    vec3 wh = distribution->Sample_wh(wi, rng.unif_rand(), rng.unif_rand());
    // Float theta = atan(sqrt(-distribution->GetAlpha() * distribution->GetAlpha() * log(1.0-rng.unif_rand())));
    // Float phi = 2 * M_PI * rng.unif_rand();
    // return(uvw.local_to_world(vec3(sin(theta) * cos(phi),sin(theta) * sin(phi),cos(theta))));
    
    return(uvw.local_to_world(Reflect(wi, wh)));
  }
  onb uvw;
  vec3 wi;
  vec2 alphas;
  MicrofacetDistribution *distribution;
};

class micro_beck_pdf : public pdf {
public:
  micro_beck_pdf(const vec3& w, Float roughness, vec3 wo_) {
    wo_.make_unit_vector();
    wo = wo_;
    a2 = roughness * roughness * roughness * roughness;
    uvw.build_from_w(w);
  }
  virtual Float value(const vec3& direction, random_gen& rng) {
    Float cosTheta = dot(unit_vector(direction), uvw.w());;
    Float expval = (a2 - 1.0f) * cosTheta + 1;
    Float D = a2 / (M_PI * expval * expval);
    return(D/4);
  }
  virtual vec3 generate(random_gen& rng) {
    Float theta = std::acos(std::sqrtf(1.0f / (-a2 * std::log(1.0f - rng.unif_rand()))));
    Float phi = 2 * M_PI * rng.unif_rand();
    vec3 wm = vec3(std::sin(theta) * std::cos(phi),
                   std::sin(theta) * std::sin(phi),
                   std::cos(theta));
    wm.make_unit_vector();
    vec3 wi = 2.0f * dot(wo, wm) * wm - wo;
    return(uvw.local_to_world(wi));
  }
  onb uvw;
  Float a2;
  vec3 wo;
};

class hitable_pdf : public pdf {
public:
  hitable_pdf(hitable *p, const vec3& origin) : ptr(p), o(origin) {}
  virtual Float value(const vec3& direction, random_gen& rng) {
    return(ptr->pdf_value(o, direction, rng));
  }
  virtual vec3 generate(random_gen& rng) {
    return(ptr->random(o, rng)); 
  }
  hitable *ptr;
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
