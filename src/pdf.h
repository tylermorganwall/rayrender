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

class pdf {
public: 
  virtual Float value(const vec3& direction, random_gen& rng) = 0;
  virtual vec3 generate(random_gen& rng) = 0;
  virtual vec3 generate(Sampler* sampler) = 0;
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
  virtual vec3 generate(Sampler* sampler) {
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
  virtual Float value(const vec3& direction, random_gen& rng) {
    vec3 wo = unit_vector(uvw.world_to_local(direction));
    vec3 wh = unit_vector(wi + wo);
    return(distribution->Pdf(wo, wi, wh) / ( 4 * dot(wo, wh) ));
  }
  virtual vec3 generate(random_gen& rng) {
    vec3 wh = distribution->Sample_wh(wi, rng.unif_rand(), rng.unif_rand());
    return(uvw.local_to_world(Reflect(wi, wh)));
  }
  virtual vec3 generate(Sampler* sampler) {
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
  virtual Float value(const vec3& direction, random_gen& rng) {
    vec3 wo = unit_vector(uvw.world_to_local(direction));
    if(wo.z() * wi.z() < 0) {
      return(INFINITY);
    }
    vec3 wh = unit_vector(wi + wo);
    return(0.5f * (AbsCosTheta(wi) * M_1_PI + distribution->Pdf(wo, wi, wh) / (4 * dot(wo, wh))));
  }
  virtual vec3 generate(random_gen& rng) {
    if(rng.unif_rand() < 0.5) {
      vec3 wh = distribution->Sample_wh(wi, rng.unif_rand(), rng.unif_rand());
      return(uvw.local_to_world(Reflect(wi, wh)));
    } else {
      return(uvw.local_to_world(rng.random_cosine_direction()));
    }
  }
  virtual vec3 generate(Sampler* sampler) {
    if(sampler->Get1D() < 0.5) {
      vec2 u = sampler->Get2D();
      vec3 wh = distribution->Sample_wh(wi, u.x(), u.y());
      return(uvw.local_to_world(Reflect(wi, wh)));
    } else {
      return(uvw.local_to_world(rand_cosine_direction(sampler->Get2D())));
    }
  }
  onb uvw;
  vec3 wi;
  MicrofacetDistribution *distribution;
};


static uint32_t Compact1By1(uint32_t x) {
  // TODO: as of Haswell, the PEXT instruction could do all this in a
  // single instruction.
  // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
  x &= 0x55555555;
  // x = --fe --dc --ba --98 --76 --54 --32 --10
  x = (x ^ (x >> 1)) & 0x33333333;
  // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
  x = (x ^ (x >> 2)) & 0x0f0f0f0f;
  // x = ---- ---- fedc ba98 ---- ---- 7654 3210
  x = (x ^ (x >> 4)) & 0x00ff00ff;
  // x = ---- ---- ---- ---- fedc ba98 7654 3210
  x = (x ^ (x >> 8)) & 0x0000ffff;
  return(x);
}

static vec2 DemuxFloat(Float f) {
  uint64_t v = f * (1ull << 32);
  uint32_t bits[2] = {Compact1By1(v), Compact1By1(v >> 1)};
  return(vec2(bits[0] / Float(1 << 16), bits[1] / Float(1 << 16)));
}

static const int pMax = 3;

class hair_pdf : public pdf {
  public:
    hair_pdf(const vec3& w, const vec3& wi_, const vec3& direction, Float eta, Float h,
            Float **cos2kAlpha_, Float **sin2kAlpha_) {
      uvw.build_from_w(w);
      wi = -unit_vector(uvw.world_to_local(wi_));
      wo = unit_vector(uvw.world_to_local(direction));
      sin2kAlpha = *sin2kAlpha_;
      cos2kAlpha = *cos2kAlpha_;
    }
    virtual Float value(const vec3& direction, random_gen& rng) {
      onb uvw;
      
      Float sinThetaO = wo.x();
      Float cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
      Float phiO = std::atan2(wo.z(), wo.y());
      
      // Compute hair coordinate system terms related to _wi_
      Float sinThetaI = wi.x();
      Float cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));
      Float phiI = std::atan2(wi.z(), wi.y());
      
      // Compute $\cos \thetat$ for refracted ray
      Float sinThetaT = sinThetaO / eta;
      Float cosThetaT = SafeSqrt(1 - Sqr(sinThetaT));
      
      // Compute $\gammat$ for refracted ray
      Float etap = std::sqrt(eta * eta - Sqr(sinThetaO)) / cosThetaO;
      Float sinGammaT = h / eta;
      Float cosGammaT = SafeSqrt(1 - Sqr(sinGammaT));
      Float gammaT = SafeASin(sinGammaT);
      vec3 wo = unit_vector(uvw.world_to_local(direction));
      
      for (int p = 0; p < pMax; ++p) {
        // Compute $\sin \thetao$ and $\cos \thetao$ terms accounting for scales
        Float sinThetaOp, cosThetaOp;
        if (p == 0) {
          sinThetaOp = sinThetaO * cos2kAlpha[1] - cosThetaO * sin2kAlpha[1];
          cosThetaOp = cosThetaO * cos2kAlpha[1] + sinThetaO * sin2kAlpha[1];
        }
        
        // Handle remainder of $p$ values for hair scale tilt
        else if (p == 1) {
          sinThetaOp = sinThetaO * cos2kAlpha[0] + cosThetaO * sin2kAlpha[0];
          cosThetaOp = cosThetaO * cos2kAlpha[0] - sinThetaO * sin2kAlpha[0];
        } else if (p == 2) {
          sinThetaOp = sinThetaO * cos2kAlpha[2] + cosThetaO * sin2kAlpha[2];
          cosThetaOp = cosThetaO * cos2kAlpha[2] - sinThetaO * sin2kAlpha[2];
        } else {
          sinThetaOp = sinThetaO;
          cosThetaOp = cosThetaO;
        }
        
        // Handle out-of-range $\cos \thetao$ from scale adjustment
        cosThetaOp = std::abs(cosThetaOp);
        *pdf += Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, v[p]) *
          apPdf[p] * Np(dphi, p, s, gammaO, gammaT);
      }
      *pdf += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, v[pMax]) *
        apPdf[pMax] * (1 / (2 * Pi));
      // if (std::abs(wi->x) < .9999) CHECK_NEAR(*pdf, Pdf(wo, *wi), .01);
      return f(wo, *wi);
      // if(wo.z() * wi.z() < 0) {
      //   return(INFINITY);
      // }
      // vec3 wh = unit_vector(wi + wo);
      // return(0.5f);
    }
    virtual vec3 generate(random_gen& rng) {
      // Float sinThetaO = wo.x();
      // Float cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
      // Float phiO = std::atan2(wo.z(), wo.y());
      // vec2 u[2] = {vec2(rng.unif_rand(),rng.unif_rand()), vec2(rng.unif_rand(),rng.unif_rand())};
      // // Determine which term $p$ to sample for hair scattering
      // std::array<Float, pMax + 1> apPdf = ComputeApPdf(cosThetaO);
      // int p;
      // for (p = 0; p < pMax; ++p) {
      //   if (u[0][0] < apPdf[p]) break;
      //   u[0][0] -= apPdf[p];
      // }
      // 
      // // Rotate $\sin \thetao$ and $\cos \thetao$ to account for hair scale tilt
      // Float sinThetaOp, cosThetaOp;
      // if (p == 0) {
      //   sinThetaOp = sinThetaO * cos2kAlpha[1] - cosThetaO * sin2kAlpha[1];
      //   cosThetaOp = cosThetaO * cos2kAlpha[1] + sinThetaO * sin2kAlpha[1];
      // }
      // else if (p == 1) {
      //   sinThetaOp = sinThetaO * cos2kAlpha[0] + cosThetaO * sin2kAlpha[0];
      //   cosThetaOp = cosThetaO * cos2kAlpha[0] - sinThetaO * sin2kAlpha[0];
      // } else if (p == 2) {
      //   sinThetaOp = sinThetaO * cos2kAlpha[2] + cosThetaO * sin2kAlpha[2];
      //   cosThetaOp = cosThetaO * cos2kAlpha[2] - sinThetaO * sin2kAlpha[2];
      // } else {
      //   sinThetaOp = sinThetaO;
      //   cosThetaOp = cosThetaO;
      // }
      // 
      // // Sample $M_p$ to compute $\thetai$
      // u[1][0] = std::max(u[1][0], Float(1e-5));
      // Float cosTheta =
      //   1 + v[p] * std::log(u[1][0] + (1 - u[1][0]) * std::exp(-2 / v[p]));
      // Float sinTheta = SafeSqrt(1 - Sqr(cosTheta));
      // Float cosPhi = std::cos(2 * Pi * u[1][1]);
      // Float sinThetaI = -cosTheta * sinThetaOp + sinTheta * cosPhi * cosThetaOp;
      // Float cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));
      // 
      // // Sample $N_p$ to compute $\Delta\phi$
      // 
      // // Compute $\gammat$ for refracted ray
      // Float etap = std::sqrt(eta * eta - Sqr(sinThetaO)) / cosThetaO;
      // Float sinGammaT = h / etap;
      // Float gammaT = SafeASin(sinGammaT);
      // Float dphi;
      // if (p < pMax)
      //   dphi =
      //     Phi(p, gammaO, gammaT) + SampleTrimmedLogistic(u[0][1], s, -Pi, Pi);
      // else
      //   dphi = 2 * Pi * u[0][1];
      // 
      // // Compute _wi_ from sampled hair scattering angles
      // Float phiI = phiO + dphi;
      // *wi = Vector3f(sinThetaI, cosThetaI * std::cos(phiI),
      //                cosThetaI * std::sin(phiI));
      // // if(rng.unif_rand() < 0.5) {
      //   // vec3 wh = distribution->Sample_wh(wi, rng.unif_rand(), rng.unif_rand());
      //   // return(uvw.local_to_world(Reflect(wi, wh)));
      // } else {
      //   // return(uvw.local_to_world(rng.random_cosine_direction()));
      // }
    }
    virtual vec3 generate(Sampler* sampler) {
      Float sinThetaO = wo.x();
      Float cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
      Float phiO = std::atan2(wo.z(), wo.y());
      
      // Derive four random samples from _u2_
        vec2 u2 = sampler->Get2D();
        vec2 u[2] = {DemuxFloat(u2.e[0]), DemuxFloat(u2.e[1])};
        std::array<Float, pMax + 1> apPdf = ComputeApPdf(cosThetaO);
        int p;
        for (p = 0; p < pMax; ++p) {
          if (u[0][0] < apPdf[p]) break;
          u[0][0] -= apPdf[p];
        }
        
        // Rotate $\sin \thetao$ and $\cos \thetao$ to account for hair scale tilt
        Float sinThetaOp, cosThetaOp;
        if (p == 0) {
          sinThetaOp = sinThetaO * cos2kAlpha[1] - cosThetaO * sin2kAlpha[1];
          cosThetaOp = cosThetaO * cos2kAlpha[1] + sinThetaO * sin2kAlpha[1];
        }
        else if (p == 1) {
          sinThetaOp = sinThetaO * cos2kAlpha[0] + cosThetaO * sin2kAlpha[0];
          cosThetaOp = cosThetaO * cos2kAlpha[0] - sinThetaO * sin2kAlpha[0];
        } else if (p == 2) {
          sinThetaOp = sinThetaO * cos2kAlpha[2] + cosThetaO * sin2kAlpha[2];
          cosThetaOp = cosThetaO * cos2kAlpha[2] - sinThetaO * sin2kAlpha[2];
        } else {
          sinThetaOp = sinThetaO;
          cosThetaOp = cosThetaO;
        }
        
        // Sample $M_p$ to compute $\thetai$
        u[1][0] = std::max(u[1][0], Float(1e-5));
        Float cosTheta =
          1 + v[p] * std::log(u[1][0] + (1 - u[1][0]) * std::exp(-2 / v[p]));
        Float sinTheta = SafeSqrt(1 - Sqr(cosTheta));
        Float cosPhi = std::cos(2 * Pi * u[1][1]);
        Float sinThetaI = -cosTheta * sinThetaOp + sinTheta * cosPhi * cosThetaOp;
        Float cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));
        
        // Sample $N_p$ to compute $\Delta\phi$
        
        // Compute $\gammat$ for refracted ray
        Float etap = std::sqrt(eta * eta - Sqr(sinThetaO)) / cosThetaO;
        Float sinGammaT = h / etap;
        Float gammaT = SafeASin(sinGammaT);
        Float dphi;
        if (p < pMax)
          dphi =
            Phi(p, gammaO, gammaT) + SampleTrimmedLogistic(u[0][1], s, -Pi, Pi);
        else
          dphi = 2 * Pi * u[0][1];
        
        // Compute _wi_ from sampled hair scattering angles
        Float phiI = phiO + dphi;
        *wi = vec3(sinThetaI, cosThetaI * std::cos(phiI),
                       cosThetaI * std::sin(phiI));
        // vec3 wh = distribution->Sample_wh(wi, u.x(), u.y());
        // return(uvw.local_to_world(Reflect(wi, wh)));
    }
    onb uvw;
    vec3 wi;
    vec3 wo;
    Float eta, h;
    Float *sin2kAlpha, *cos2kAlpha;
  private:
    std::array<Float, pMax + 1> ComputeApPdf(Float cosThetaO) const;
};

std::array<Float, pMax + 1> hair_pdf::ComputeApPdf(Float cosThetaO) const {
  // Compute array of $A_p$ values for _cosThetaO_
  Float sinThetaO = SafeSqrt(1 - cosThetaO * cosThetaO);
  
  // Compute $\cos \thetat$ for refracted ray
  Float sinThetaT = sinThetaO / eta;
  Float cosThetaT = SafeSqrt(1 - Sqr(sinThetaT));
  
  // Compute $\gammat$ for refracted ray
  Float etap = std::sqrt(eta * eta - Sqr(sinThetaO)) / cosThetaO;
  Float sinGammaT = h / etap;
  Float cosGammaT = SafeSqrt(1 - Sqr(sinGammaT));
  
  // Compute the transmittance _T_ of a single path through the cylinder
  vec3 T = Exp(-sigma_a * (2 * cosGammaT / cosThetaT));
  std::array<vec3, pMax + 1> ap = Ap(cosThetaO, eta, h, T);
  
  // Compute $A_p$ PDF from individual $A_p$ terms
  std::array<Float, pMax + 1> apPdf;
  Float sumY = std::accumulate(ap.begin(), ap.end(), Float(0),
                               [](Float s, const vec3 &ap) { return s + ap.y(); });
  for (int i = 0; i <= pMax; ++i) {
    apPdf[i] = ap[i].y() / sumY;
  }
  return apPdf;
}

class hitable_pdf : public pdf {
public:
  hitable_pdf(hitable_list *p, const vec3& origin) : ptr(p), o(origin) {}
  virtual Float value(const vec3& direction, random_gen& rng) {
    return(ptr->pdf_value(o, direction, rng));
  }
  virtual vec3 generate(random_gen& rng) {
    return(ptr->random(o, rng)); 
  }
  virtual vec3 generate(Sampler* sampler) {
    return(ptr->random(o, sampler)); 
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
  virtual vec3 generate(Sampler* sampler) {
    if(sampler->Get1D() < 0.5) {
      return(p[0]->generate(sampler));
    } else {
      return(p[1]->generate(sampler));
    } 
  }
  pdf *p[2];
};

#endif
