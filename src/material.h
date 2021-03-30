#ifndef MATERIALH
#define MATERIALH 

#include "ray.h"
#include "hitable.h"
#include "onbh.h"
#include "pdf.h"
#include "rng.h"
#include "mathinline.h"
#include "microfacetdist.h"
#include <array>

// #define DEBUG2

// #ifdef DEBUG2
// #include <iostream>
// #include <fstream>
// using namespace std;
// #endif

// #ifdef DEBUG2
//     ray debugrayi(rec.p, wi);
//     ray debugrayo(rec.p, wo);
//     ofstream myfile;
//     myfile.open("onb.txt", ios::app | ios::out);
//     myfile << rec.p << ", " << uvw.world_to_local(wi) << ", " << uvw.world_to_local(wo) << ", " << 
//       uvw.w() << ", " << uvw.v() << ", " << uvw.u() << "\n";
//     myfile.close();
// #endif

#include <utility>

struct hit_record;

inline vec3 reflect(const vec3& v, const vec3& n) {
  return(v - 2*dot(v,n) * n);
}

inline Float schlick(Float cosine, Float ref_idx, Float ref_idx2) {
  Float r0 = (ref_idx2 - ref_idx) / (ref_idx2 + ref_idx);
  r0 = r0 * r0;
  return(r0 + (1-r0) * pow((1-cosine),5));
}

inline Float schlick_reflection(Float cosine, Float r0) {
  Float r02 = (1 - r0) / (1 + r0);
  r02 = r02 * r02;
  return(r02 + (1-r02) * pow((1-cosine),5));
}

inline bool refract(const vec3& v, const vec3& n, Float ni_over_nt, vec3& refracted) {
  vec3 uv = unit_vector(v);
  Float dt = dot(uv, n);
  Float discriminant = 1.0 - ni_over_nt * ni_over_nt * (1 - dt * dt);
  if(discriminant > 0) {
    refracted = ni_over_nt * (uv - n * dt) - n * sqrt(discriminant);
    return(true);
  } else {
    return(false);
  }
}

inline vec3 refract(const vec3& uv, const vec3& n, Float ni_over_nt) {
  Float cos_theta = dot(-uv, n);
  vec3 r_out_parallel =  ni_over_nt * (uv + cos_theta*n);
  vec3 r_out_perp = -sqrt(1.0 - r_out_parallel.squared_length()) * n;
  return(r_out_parallel + r_out_perp);
}

struct scatter_record {
  ray specular_ray;
  bool is_specular;
  vec3 attenuation;
  pdf *pdf_ptr = nullptr;
  ~scatter_record() { if(pdf_ptr) delete pdf_ptr; }
};

inline vec3 FrCond(Float cosi, const vec3 &eta, const vec3 &k) {
  vec3 tmp = (eta*eta + k*k) * cosi*cosi;
  vec3 Rparl2 = (tmp - (2.0f * eta * cosi) + vec3(1.0f)) /
    (tmp + (2.0f * eta * cosi) + vec3(1.0f));
  vec3 tmp_f = eta*eta + k*k;
  vec3 Rperp2 = (tmp_f - (2.0f * eta * cosi) + cosi*cosi) /
    (tmp_f + (2.0f * eta * cosi) + cosi*cosi);
  return((Rparl2 + Rperp2) / 2.0f);
}

class material {
  public:
    virtual bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng) {
      return(false);
    };
    virtual bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, Sampler* sampler) {
      return(false);
    };
    virtual vec3 f(const ray& r_in, const hit_record& rec, const ray& scattered) const {
      return(vec3(0,0,0));
    }
    virtual vec3 emitted(const ray& r_in, const hit_record& rec, Float u, Float v, const vec3& p, bool& is_invisible) {
      return(vec3(0,0,0));
    }
    virtual vec3 get_albedo(const ray& r_in, const hit_record& rec) const {
      return(vec3(0,0,0));
    }
    virtual ~material() {};
};


class lambertian : public material {
  public: 
    lambertian(std::shared_ptr<texture> a) : albedo(a) {}
    vec3 f(const ray& r_in, const hit_record& rec, const ray& scattered) const;
    bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng);
    bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, Sampler* sampler);
    vec3 get_albedo(const ray& r_in, const hit_record& rec) const;
    
    std::shared_ptr<texture> albedo;
};

class metal : public material {
  public:
    metal(std::shared_ptr<texture> a, Float f, vec3 eta, vec3 k) : albedo(a), eta(eta), k(k) { 
      if (f < 1) fuzz = f; else fuzz = 1;
    }
    virtual bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng);
    virtual bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, Sampler* sampler);
    vec3 get_albedo(const ray& r_in, const hit_record& rec) const;
    std::shared_ptr<texture> albedo;
    vec3 eta, k;
    Float fuzz;
};
// 

class dielectric : public material {
  public:
    dielectric(const vec3& a, Float ri, const vec3& atten, int priority2, random_gen& rng) : ref_idx(ri), 
               albedo(a), attenuation(atten), priority(priority2), rng(rng) {};
    bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng);
    bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, Sampler* sampler);
    
    vec3 get_albedo(const ray& r_in, const hit_record& rec) const {
      return(vec3(1,1,1));
    }
    Float ref_idx;
    vec3 albedo;
    vec3 attenuation;
    size_t priority;
    random_gen rng;
};

class diffuse_light : public material {
public:
  diffuse_light(std::shared_ptr<texture>  a, Float intensity, bool invisible) : 
    emit(a), intensity(intensity), invisible(invisible) {}
  ~diffuse_light() {}
  bool scatter(const ray& r_in, const hit_record& rec, scatter_record& srec, random_gen& rng) {
    return(false);
  }
  bool scatter(const ray& r_in, const hit_record& rec, scatter_record& srec, Sampler* sampler) {
    return(false);
  }
  vec3 emitted(const ray& r_in, const hit_record& rec, Float u, Float v, const vec3& p, bool& is_invisible);
  vec3 get_albedo(const ray& r_in, const hit_record& rec) const;
  std::shared_ptr<texture>  emit;
  Float intensity;
  bool invisible;
};

class spot_light : public material {
  public:
    spot_light(std::shared_ptr<texture>  a, vec3 dir, Float cosTotalWidth, Float cosFalloffStart, bool invisible) : 
      emit(a), spot_direction(unit_vector(dir)), cosTotalWidth(cosTotalWidth), 
      cosFalloffStart(cosFalloffStart), invisible(invisible) {}
    ~spot_light() {}
    bool scatter(const ray& r_in, const hit_record& rec, scatter_record& srec, random_gen& rng) {
      return(false);
    }
    bool scatter(const ray& r_in, const hit_record& rec, scatter_record& srec, Sampler* sampler) {
      return(false);
    }
    vec3 emitted(const ray& r_in, const hit_record& rec, Float u, Float v, const vec3& p, bool& is_invisible);
    Float falloff(const vec3 &w) const;
    vec3 get_albedo(const ray& r_in, const hit_record& rec) const;
    std::shared_ptr<texture>  emit;
    vec3 spot_direction;
    const Float cosTotalWidth, cosFalloffStart;
    bool invisible;
};

class isotropic : public material {
public:
  isotropic(std::shared_ptr<texture> a) : albedo(a) {}
  ~isotropic() {}
  virtual bool scatter(const ray& r_in, const hit_record& rec, scatter_record& srec, random_gen& rng);
  virtual bool scatter(const ray& r_in, const hit_record& rec, scatter_record& srec, Sampler* sampler);
  vec3 f(const ray& r_in, const hit_record& rec, const ray& scattered) const;
  vec3 get_albedo(const ray& r_in, const hit_record& rec) const;
  std::shared_ptr<texture> albedo;
};

class orennayar : public material {
public:
  orennayar(std::shared_ptr<texture>  a, Float sigma) : albedo(a)  {
    Float sigma2 = sigma*sigma;
    A = 1.0f - (sigma2 / (2.0f * (sigma2 + 0.33f)));
    B = 0.45f * sigma2 / (sigma2 + 0.09f);
  }
  ~orennayar() {}
  bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng);
  bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, Sampler* sampler);
  vec3 f(const ray& r_in, const hit_record& rec, const ray& scattered) const;
  vec3 get_albedo(const ray& r_in, const hit_record& rec) const;
  Float A, B;
  std::shared_ptr<texture> albedo;
};

class MicrofacetReflection : public material {
public:
  MicrofacetReflection(std::shared_ptr<texture> a, MicrofacetDistribution *distribution, 
                       vec3 eta, vec3 k)
    : albedo(a), distribution(distribution), eta(eta), k(k) {}
  ~MicrofacetReflection() {
    if(distribution) delete distribution;
  }
  
  bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng);
  bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, Sampler* sampler);
  vec3 f(const ray& r_in, const hit_record& rec, const ray& scattered) const;
  vec3 get_albedo(const ray& r_in, const hit_record& rec) const;
private:
  std::shared_ptr<texture> albedo;
  MicrofacetDistribution *distribution;
  vec3 eta;
  vec3 k;
};

class glossy : public material {
public:
  glossy(std::shared_ptr<texture> a, MicrofacetDistribution *distribution, 
         vec3 Rs, vec3 Rd2)
    : albedo(a), distribution(distribution), Rs(Rs) {}
  ~glossy() {
    if(distribution) delete distribution;
  }
  bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng);
  bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, Sampler* sampler);
  vec3 f(const ray& r_in, const hit_record& rec, const ray& scattered) const;
  vec3 SchlickFresnel(Float cosTheta) const;
  vec3 get_albedo(const ray& r_in, const hit_record& rec) const;
  
private:
  std::shared_ptr<texture>  albedo;
  MicrofacetDistribution *distribution;
  vec3 Rs;
};

// Hair Local Functions

class hair : public material {
  public:
    hair(vec3 sigma_a, 
         Float eta, Float beta_m, Float beta_n, Float alpha) :
      sigma_a(sigma_a), 
      eta(eta), 
      beta_m(beta_m),  //longitudinal roughness
      beta_n(beta_n), 
      alpha(alpha) {
      v[0] = Sqr(0.726f * beta_m + 0.812f * Sqr(beta_m) + 3.7f * Pow<20>(beta_m));
      v[1] = 0.25 * v[0];
      v[2] = 4 * v[0];
      for (int p = 3; p <= pMax; ++p) {
        v[p] = v[2];
      }
      // Compute azimuthal logistic scale factor from $\beta_n$
      s = SqrtPiOver8 *(0.265f * beta_n + 1.194f * Sqr(beta_n) + 5.372f * Pow<22>(beta_n));

      // Compute $\alpha$ terms for hair scales
      sin2kAlpha[0] = std::sin(mpi_over_180 * alpha);
      cos2kAlpha[0] = SafeSqrt(1 - Sqr(sin2kAlpha[0]));
      for (int i = 1; i < 3; ++i) {
        sin2kAlpha[i] = 2 * cos2kAlpha[i - 1] * sin2kAlpha[i - 1];
        cos2kAlpha[i] = Sqr(cos2kAlpha[i - 1]) - Sqr(sin2kAlpha[i - 1]);
      }
    }
    ~hair() {}
    bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng);
    bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, Sampler* sampler);
    
    vec3 f(const ray& r_in, const hit_record& rec, const ray& scattered) const;
    
    vec3 sigma_a;
    Float eta;
    Float beta_m, beta_n;
    Float alpha;
  private:
    std::array<Float, pMax + 1> ComputeApPdf(Float cosThetaO, Float h) const;
    Float v[pMax + 1];
    Float s;
    Float sin2kAlpha[3], cos2kAlpha[3];
};

#endif
