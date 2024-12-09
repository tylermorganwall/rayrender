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

inline vec3f reflect(const vec3f& v, const vec3f& n) {
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


struct scatter_record {
  ray specular_ray;
  bool is_specular;
  point3f attenuation;
  pdf *pdf_ptr = nullptr;
  ~scatter_record() { if(pdf_ptr) delete pdf_ptr; }
};

inline point3f FrCond(Float cosi, const point3f &eta, const point3f &k) {
  point3f tmp = (eta*eta + k*k) * cosi*cosi;
  point3f Rparl2 = (tmp + (-2.0f * eta * cosi) + point3f(1.0f)) /
    (tmp + (2.0f * eta * cosi) + vec3f(1.0f));
  point3f tmp_f = eta*eta + k*k;
  point3f Rperp2 = (tmp_f + (-2.0f * eta * cosi) + point3f(1)*cosi*cosi) /
    (tmp_f + (2.0f * eta * cosi) + point3f(1)*cosi*cosi);
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
    virtual point3f f(const ray& r_in, const hit_record& rec, const vec3f& scattered) const {
      return(point3f(0,0,0));
    }
    virtual point3f emitted(const ray& r_in, const hit_record& rec, Float u, Float v, const point3f& p, bool& is_invisible) {
      return(point3f(0,0,0));
    }
    virtual point3f get_albedo(const hit_record& rec) const {
      return(point3f(0,0,0));
    }
    virtual ~material() {};
    virtual const std::string GetName() = 0;
    virtual size_t GetSize() = 0;
};


class lambertian : public material {
  public: 
    lambertian(std::shared_ptr<texture> a) : albedo(a) {}
    point3f f(const ray& r_in, const hit_record& rec, const vec3f& scattered) const;
    bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng);
    bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, Sampler* sampler);
    point3f get_albedo(const hit_record& rec) const;
    size_t GetSize();
    const std::string GetName() {
      return(std::string("diffuse"));
    };
    std::shared_ptr<texture> albedo;
};

class metal : public material {
  public:
    metal(std::shared_ptr<texture> a, Float f, point3f eta, point3f k) : albedo(a), eta(eta), k(k) { 
      if (f < 1) fuzz = f; else fuzz = 1;
    }
    virtual bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng);
    virtual bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, Sampler* sampler);
    point3f get_albedo(const hit_record& rec) const;
    size_t GetSize();
    const std::string GetName() {
      return(std::string("metal"));
    };
    std::shared_ptr<texture> albedo;
    point3f eta, k;
    Float fuzz;
};
// 

class dielectric : public material {
  public:
    dielectric(const point3f& a, Float ri, const point3f& atten, int priority2) : ref_idx(ri), 
               albedo(a), attenuation(atten), priority(priority2) {};
    bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng);
    bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, Sampler* sampler);
    
    point3f get_albedo(const hit_record& rec) const {
      return(point3f(1,1,1));
    }
    size_t GetSize();
    const std::string GetName() {
      return(std::string("dielectric"));
    };
    Float ref_idx;
    point3f albedo;
    point3f attenuation;
    size_t priority;
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
  point3f emitted(const ray& r_in, const hit_record& rec, Float u, Float v, const point3f& p, bool& is_invisible);
  point3f get_albedo(const hit_record& rec) const;
  size_t GetSize();
  const std::string GetName() {
    return(std::string("diffuse_light"));
  };
  std::shared_ptr<texture>  emit;
  Float intensity;
  bool invisible;
};

class spot_light : public material {
  public:
    spot_light(std::shared_ptr<texture>  a, vec3f dir, Float cosTotalWidth, Float cosFalloffStart, Float intensity, bool invisible) : 
      emit(a), spot_direction(unit_vector(dir)), intensity(intensity), cosTotalWidth(cosTotalWidth), 
      cosFalloffStart(cosFalloffStart), invisible(invisible) {}
    ~spot_light() {}
    bool scatter(const ray& r_in, const hit_record& rec, scatter_record& srec, random_gen& rng) {
      return(false);
    }
    bool scatter(const ray& r_in, const hit_record& rec, scatter_record& srec, Sampler* sampler) {
      return(false);
    }
    point3f emitted(const ray& r_in, const hit_record& rec, Float u, Float v, const point3f& p, bool& is_invisible);
    Float falloff(const vec3f &w) const;
    point3f get_albedo(const hit_record& rec) const;
    size_t GetSize();
    const std::string GetName() {
      return(std::string("spot_light"));
    };
    std::shared_ptr<texture>  emit;
    vec3f spot_direction;
    Float intensity;
    const Float cosTotalWidth, cosFalloffStart;
    bool invisible;
};

class isotropic : public material {
public:
  isotropic(std::shared_ptr<texture> a) : albedo(a) {}
  ~isotropic() {}
  virtual bool scatter(const ray& r_in, const hit_record& rec, scatter_record& srec, random_gen& rng);
  virtual bool scatter(const ray& r_in, const hit_record& rec, scatter_record& srec, Sampler* sampler);
  point3f f(const ray& r_in, const hit_record& rec, const vec3f& scattered) const;
  point3f get_albedo(const hit_record& rec) const;
  size_t GetSize();
  const std::string GetName() {
    return(std::string("isotropic"));
  };
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
  point3f f(const ray& r_in, const hit_record& rec, const vec3f& scattered) const;
  point3f get_albedo(const hit_record& rec) const;
  size_t GetSize();
  const std::string GetName() {
    return(std::string("orennayer"));
  };
  Float A, B;
  std::shared_ptr<texture> albedo;
};

class MicrofacetReflection : public material {
public:
  MicrofacetReflection(std::shared_ptr<texture> a, MicrofacetDistribution *distribution, 
                       point3f eta, point3f k)
    : albedo(a), distribution(distribution), eta(eta), k(k) {}
  ~MicrofacetReflection() {
    if(distribution) delete distribution;
  }
  
  bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng);
  bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, Sampler* sampler);
  point3f f(const ray& r_in, const hit_record& rec, const vec3f& scattered) const;
  point3f get_albedo(const hit_record& rec) const;
  size_t GetSize();
  const std::string GetName() {
    return(std::string("MicrofacetReflection"));
  };
private:
  std::shared_ptr<texture> albedo;
  MicrofacetDistribution *distribution;
  point3f eta;
  point3f k;
};

class MicrofacetTransmission : public material {
public:
  MicrofacetTransmission(std::shared_ptr<texture> a, MicrofacetDistribution *distribution, 
                         point3f eta2, point3f k) : 
    albedo(a), distribution(distribution), eta(eta2[0]), k(k) {}
  ~MicrofacetTransmission() {
    if(distribution) delete distribution;
  }
  
  bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng);
  bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, Sampler* sampler);
  point3f f(const ray& r_in, const hit_record& rec, const vec3f& scattered) const;
  point3f get_albedo(const hit_record& rec) const;
  point3f SchlickFresnel(Float cosTheta) const;
  size_t GetSize();
  const std::string GetName() {
    return(std::string("MicrofacetTransmission"));
  };
  
private:
  std::shared_ptr<texture> albedo;
  MicrofacetDistribution *distribution;
  Float eta;
  point3f k;
};

class glossy : public material {
public:
  glossy(std::shared_ptr<texture> a, MicrofacetDistribution *distribution, 
         point3f Rs, point3f Rd2)
    : albedo(a), distribution(distribution), Rs(Rs), Rd(Rd2) {}
  ~glossy() {
    if(distribution) delete distribution;
  }
  bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng);
  bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, Sampler* sampler);
  point3f f(const ray& r_in, const hit_record& rec, const vec3f& scattered) const;
  point3f SchlickFresnel(Float cosTheta) const;
  point3f get_albedo(const hit_record& rec) const;
  size_t GetSize();
  const std::string GetName() {
    return(std::string("glossy"));
  };
  
private:
  std::shared_ptr<texture>  albedo;
  MicrofacetDistribution *distribution;
  point3f Rs;
  point3f Rd;
};

// Hair Local Functions

class hair : public material {
  public:
    hair(point3f sigma_a, 
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
    
    point3f f(const ray& r_in, const hit_record& rec, const vec3f& scattered) const;
    size_t GetSize();
    const std::string GetName() {
      return(std::string("hair"));
    };
    
    point3f sigma_a;
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
