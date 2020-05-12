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
    return(AbsCosTheta(wo) * M_1_PI);
    // return(distribution->D(wo) *  AbsCosTheta(wo) / (4 * dot(wo,wi)));
  }
  static void TrowbridgeReitzSample11(Float cosTheta, Float U1, Float U2,
                                      Float *slope_x, Float *slope_y);
  virtual vec3 generate(random_gen& rng) {
    vec3 wiStretched = unit_vector(vec3(alphas.x() * wi.x(),  wi.y(), alphas.y() *wi.z()));

    // 2. simulate P22_{wi}(x_slope, y_slope, 1, 1)
    Float slope_x, slope_y;
    TrowbridgeReitzSample11(CosTheta(wiStretched), rng.unif_rand(), rng.unif_rand(), 
                            &slope_x, &slope_y);

    // 3. rotate
    Float tmp = CosPhi(wiStretched) * slope_x - SinPhi(wiStretched) * slope_y;
    slope_y = SinPhi(wiStretched) * slope_x + CosPhi(wiStretched) * slope_y;
    slope_x = tmp;

    // 4. unstretch
    slope_x = alphas.x() * slope_x;
    slope_y = alphas.y() * slope_y;
    // 5. compute normal
    return(uvw.local(unit_vector(vec3(-slope_x, 1.0f, -slope_y))));
    // Float tan2Theta, phi;
    // Float u0 = rng.unif_rand();
    // Float u1 = rng.unif_rand();
    // if (alphas.x() == alphas.y()) {
    //   Float logSample = std::log(1 - u0);
    //   if (std::isinf(logSample)) {
    //     logSample = 0;
    //   }
    //   tan2Theta = -alphas.x() * alphas.x() * logSample;
    //   phi = u1 * 2 * M_PI;
    // } else {
    //   Float logSample = std::log(u0);
    //   phi = std::atan(alphas.y() / alphas.x() * std::tan(2 * M_PI * u1 + 0.5f * M_PI));
    //   if (u1 > 0.5f) {
    //     phi += M_PI;
    //   }
    //   Float sinPhi = std::sin(phi), cosPhi = std::cos(phi);
    //   Float alphax2 = alphas.x() * alphas.x();
    //   Float alphay2 = alphas.y() * alphas.y();
    //   tan2Theta = -logSample /(cosPhi * cosPhi / alphax2 + sinPhi * sinPhi / alphay2);
    // }
    // Float cosTheta = 1 / std::sqrt(1 + tan2Theta);
    // Float sinTheta = std::sqrt(std::max((Float)0, 1 - cosTheta * cosTheta));
    // vec3 wh = SphericalDirection(sinTheta, cosTheta, phi);
    // if (!SameHemisphere(wi, wh)) {
    //   wh = -wh;
    // }
    // vec3 wo = -wi + 2 * dot(wi, wh) * wh;
    // return(uvw.local(wo));
  }
  onb uvw;
  vec3 wi;
  vec2 alphas;
  MicrofacetDistribution *distribution;
};

void micro_trow_pdf::TrowbridgeReitzSample11(Float cosTheta, Float U1, Float U2,
                                    Float *slope_x, Float *slope_y) {
  // special case (normal incidence)
  if (cosTheta > .9999) {
    Float r = sqrt(U1 / (1 - U1));
    Float phi = 2 * M_PI * U2;
    *slope_x = r * cos(phi);
    *slope_y = r * sin(phi);
    return;
  }
  
  Float sinTheta = std::sqrt(std::max((Float)0, (Float)1 - cosTheta * cosTheta));
  Float tanTheta = sinTheta / cosTheta;
  Float a = 1 / tanTheta;
  Float G1 = 2 / (1 + std::sqrt(1.f + 1.f / (a * a)));
  
  // sample slope_x
  Float A = 2 * U1 / G1 - 1;
  Float tmp = 1.f / (A * A - 1.f);
  if (tmp > 1e10) tmp = 1e10;
  Float B = tanTheta;
  Float D = std::sqrt(std::max(Float(B * B * tmp * tmp - (A * A - B * B) * tmp), Float(0)));
  Float slope_x_1 = B * tmp - D;
  Float slope_x_2 = B * tmp + D;
  *slope_x = (A < 0 || slope_x_2 > 1.f / tanTheta) ? slope_x_1 : slope_x_2;
  
  // sample slope_y
  Float S;
  if (U2 > 0.5f) {
    S = 1.f;
    U2 = 2.f * (U2 - .5f);
  } else {
    S = -1.f;
    U2 = 2.f * (.5f - U2);
  }
  Float z =
    (U2 * (U2 * (U2 * 0.27385f - 0.73369f) + 0.46341f)) /
      (U2 * (U2 * (U2 * 0.093073f + 0.309420f) - 1.000000f) + 0.597999f);
  *slope_y = S * z * std::sqrt(1.f + *slope_x * *slope_x);
}

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
