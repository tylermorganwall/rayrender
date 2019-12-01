#ifndef PDFH
#define PDFH

#include "onbh.h"
#include "hitable.h"
#include "rng.h"

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
    return(uvw.local(rng.random_cosine_direction()));
  }
  onb uvw;
};

class micro_pdf : public pdf {
public:
  micro_pdf(const vec3& w, Float roughness, vec3 wo_) {
    wo_.make_unit_vector();
    wo = wo_;
    a2 = roughness * roughness * roughness * roughness;
    uvw.build_from_w(w);
  }
  virtual Float value(const vec3& direction, random_gen& rng) {
    Float cosTheta = dot(unit_vector(direction), uvw.w());;
    Float expval = (a2 - 1.0f) * cosTheta + 1;
    Float D = a2 / (M_PI * expval * expval);
    // vec3 wi = wo - 2.0f * dot(uvw.w(), direction) * uvw.w();
    return(D/4);
  }
  virtual vec3 generate(random_gen& rng) {
    Float r0 = rng.unif_rand();
    Float theta = std::acos(std::sqrt((1 - r0) / ((a2-1)*r0 + 1)));
    Float phi = 2 * M_PI * rng.unif_rand();
    vec3 wm = vec3(std::sin(theta) * std::cos(phi),
                   std::cos(theta),
                   std::sin(theta) * std::sin(phi));
    wm.make_unit_vector();
    vec3 wi = 2.0f * dot(wo, wm) * wm - wo;
    return(uvw.local(wi));
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
