#ifndef PDFH
#define PDFH

#include "onbh.h"
#include "hitable.h"
#include "rng.h"
#include "vec2.h"
#include "microfacetdist.h"

#include <chrono>
#include <thread>

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

class micro_trow_pdf : public pdf {
public:
  micro_trow_pdf(const vec3& w, const vec3& wi_, MicrofacetDistribution* distribution) : distribution(distribution) {
    wi = wi_;
    wi.make_unit_vector();
    alphas = distribution->GetAlphas();
    uvw.build_from_w(w);
  }
  virtual Float value(const vec3& direction, random_gen& rng) {
    vec3 wo = unit_vector(direction);
    if (!SameHemisphere(wo, wi)) {
      return(0);
    }
    // vec3 wh = (wo + wi);
    // wh.make_unit_vector();
    // std::this_thread::sleep_for(std::chrono::milliseconds(10));
    // Rcpp::Rcout << distribution->D(wh) *  AbsCosTheta(wh) << "\n";
    // vec3 colorval = hrec.mat_ptr->f(r, hrec, scattered);
    // Rcpp::Rcout << colorval.x() << " " << colorval.y() << " " << colorval.z() << "\n";
    // vec3 wi = wo - 2.0f * dot(uvw.v(), wh) * uvw.v();
    // Rcpp::Rcout << distribution->D(wh) *  AbsCosTheta(wh) << "\n";
    return(distribution->D(wo) *  AbsCosTheta(wo) / (4 * dot(wo,wi)));
    // Float cosTheta = dot(unit_vector(direction), uvw.w());;
    // Float expval = (a2 - 1.0f) * cosTheta * cosTheta + 1;
    // // Float D = a2 * cosTheta * std::sqrtf(1-cosTheta * cosTheta)/ (M_PI * expval * expval);
    // // vec3 wi = wo - 2.0f * dot(uvw.w(), direction) * uvw.w();
    // return(a2 * cosTheta * std::sqrtf(1-cosTheta * cosTheta)/ (M_PI * expval * expval));
  }
  virtual vec3 generate(random_gen& rng) {
    Float tan2Theta, phi;
    Float u0 = rng.unif_rand();
    Float u1 = rng.unif_rand();
    if (alphas.x() == alphas.y()) {
      Float logSample = std::log(1 - u0);
      if (std::isinf(logSample)) {
        logSample = 0;
      }
      tan2Theta = -alphas.x() * alphas.x() * logSample;
      phi = u1 * 2 * M_PI;
    } else {
      Float logSample = std::log(u0);
      phi = std::atan(alphas.y() / alphas.x() * std::tan(2 * M_PI * u1 + 0.5f * M_PI));
      if (u1 > 0.5f) {
        phi += M_PI;
      }
      Float sinPhi = std::sin(phi), cosPhi = std::cos(phi);
      Float alphax2 = alphas.x() * alphas.x();
      Float alphay2 = alphas.y() * alphas.y();
      tan2Theta = -logSample /(cosPhi * cosPhi / alphax2 + sinPhi * sinPhi / alphay2);
    }
    Float cosTheta = 1 / std::sqrt(1 + tan2Theta);
    Float sinTheta = std::sqrt(std::max((Float)0, 1 - cosTheta * cosTheta));
    vec3 wh = SphericalDirection(sinTheta, cosTheta, phi);
    if (!SameHemisphere(wi, wh)) {
      wh = -wh;
    }
    vec3 wo = -wi + 2 * dot(wi, wh) * wh;
    return(uvw.local(wo));
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
    // Rcpp::Rcout << "DUM" << "\n";
    
    // vec3 wi = wo - 2.0f * dot(uvw.w(), direction) * uvw.w();
    return(D/4);
  }
  virtual vec3 generate(random_gen& rng) {
    Float theta = std::acos(std::sqrtf(1.0f / (-a2 * std::log(1.0f - rng.unif_rand()))));
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
