#ifndef MATERIALH
#define MATERIALH 

#include "ray.h"
#include "hitable.h"
#include "onbh.h"
#include "pdf.h"
#include "rng.h"
#include "mathinline.h"
#include "microfacetdist.h"

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

vec3 reflect(const vec3& v, const vec3& n) {
  return(v - 2*dot(v,n) * n);
}

Float schlick(Float cosine, Float ref_idx, Float ref_idx2) {
  Float r0 = (ref_idx2 - ref_idx) / (ref_idx2 + ref_idx);
  r0 = r0 * r0;
  return(r0 + (1-r0) * pow((1-cosine),5));
}

Float schlick_reflection(Float cosine, Float r0) {
  Float r02 = (1 - r0) / (1 + r0);
  r02 = r02 * r02;
  return(r02 + (1-r02) * pow((1-cosine),5));
}

bool refract(const vec3& v, const vec3& n, Float ni_over_nt, vec3& refracted) {
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

vec3 refract(const vec3& uv, const vec3& n, Float ni_over_nt) {
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

class material {
  public:
    virtual bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng) {
      return(false);
    };
    virtual vec3 f(const ray& r_in, const hit_record& rec, const ray& scattered) const {
      return(vec3(0,0,0));
    }
    virtual vec3 emitted(const ray& r_in, const hit_record& rec, Float u, Float v, const vec3& p) const {
      return(vec3(0,0,0));
    }
    virtual ~material() {};
};


class lambertian : public material {
  public: 
    lambertian(texture *a) : albedo(a) {}
    ~lambertian() {
      if(albedo) delete albedo;
    }
    vec3 f(const ray& r_in, const hit_record& rec, const ray& scattered) const {
      //unit_vector(scattered.direction()) == wo
      //r_in.direction() == wi
      vec3 wi = unit_vector(scattered.direction());
      Float cosine = dot(rec.normal, wi);
      //Shadow terminator if bump map
      Float G = 1.0f;
      if(rec.has_bump) {
        Float NsNlight = dot(rec.bump_normal, wi);
        Float NsNg = dot(rec.bump_normal, rec.normal);
        G = NsNlight > 0.0 && NsNg > 0.0 ? std::fmin(1.0, dot(wi, rec.normal) / (NsNlight * NsNg)) : 0;
        G = G > 0.0f ? -G * G * G + G * G + G : 0.0f;
        cosine = dot(rec.bump_normal, wi);
      }
      if(cosine < 0) {
        cosine = 0;
      }
      return(G * albedo->value(rec.u, rec.v, rec.p) * cosine * M_1_PI);
    }
    bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng) {
      srec.is_specular = false;
      srec.attenuation = albedo->value(hrec.u, hrec.v, hrec.p);
      srec.pdf_ptr = new cosine_pdf(hrec.normal);
      return(true);
    }
    
    texture *albedo;
};

class metal : public material {
  public:
    metal(texture* a, Float f, vec3 eta, vec3 k) : albedo(a), eta(eta), k(k) { 
      if (f < 1) fuzz = f; else fuzz = 1;
    }
    ~metal() {
      if(albedo) delete albedo;
    }
    virtual bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng) {
      vec3 normal = !hrec.has_bump ? hrec.normal : hrec.bump_normal;
      vec3 wi = -unit_vector(r_in.direction());
      vec3 reflected = Reflect(wi,unit_vector(normal));
      Float cosine = AbsDot(unit_vector(reflected), unit_vector(normal));
      if(cosine < 0) {
        cosine = 0;
      }
      vec3 offset_p = offset_ray(hrec.p-r_in.A, hrec.normal) + r_in.A;
      
      srec.specular_ray = ray(offset_p, reflected + fuzz * rng.random_in_unit_sphere(), r_in.pri_stack, r_in.time());
      srec.attenuation = albedo->value(hrec.u, hrec.v, hrec.p) * FrCond(cosine, eta, k);
      srec.is_specular = true;
      srec.pdf_ptr = 0;
      return(true);
    }
    texture *albedo;
    vec3 eta, k;
    Float fuzz;
};
// 

class dielectric : public material {
  public:
    dielectric(const vec3& a, Float ri, const vec3& atten, int priority2, random_gen& rng) : ref_idx(ri), 
               albedo(a), attenuation(atten), priority(priority2), rng(rng) {};
    virtual bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng) {
      srec.is_specular = true;
      vec3 outward_normal;
      vec3 normal = !hrec.has_bump ? hrec.normal : hrec.bump_normal;
      
      vec3 reflected = reflect(r_in.direction(), normal);
      Float ni_over_nt;
      srec.attenuation = albedo;
      Float current_ref_idx = 1.0;
      
      size_t active_priority_value = priority;
      size_t next_down_priority = 100000;
      int current_layer = -1; //keeping track of index of current material
      int prev_active = -1;  //keeping track of index of active material (higher priority)

      bool entering = dot(hrec.normal, r_in.direction()) < 0;
      bool skip = false;
      
      for(size_t i = 0; i < r_in.pri_stack->size(); i++) {
        if(r_in.pri_stack->at(i) == this) {
          current_layer = i;
          continue;
        }
        if(r_in.pri_stack->at(i)->priority < active_priority_value) {
          active_priority_value = r_in.pri_stack->at(i)->priority;
          skip = true;
        }
        if(r_in.pri_stack->at(i)->priority < next_down_priority && r_in.pri_stack->at(i) != this) {
          prev_active = i;
          next_down_priority = r_in.pri_stack->at(i)->priority;
        }
      }
      if(entering) {
        r_in.pri_stack->push_back(this);
      }
      current_ref_idx = prev_active != -1 ? r_in.pri_stack->at(prev_active)->ref_idx : 1;
      if(skip) {
        srec.specular_ray = ray(hrec.p, r_in.direction(), r_in.pri_stack, r_in.time());
        Float distance = (hrec.p-r_in.point_at_parameter(0)).length();
        vec3 prev_atten = r_in.pri_stack->at(prev_active)->attenuation;
        srec.attenuation = vec3(std::exp(-distance * prev_atten.x()),
                                std::exp(-distance * prev_atten.y()),
                                std::exp(-distance * prev_atten.z()));
        if(!entering && current_layer != -1) {
          r_in.pri_stack->erase(r_in.pri_stack->begin() + static_cast<size_t>(current_layer));
        }
        return(true);
      }
      
      Float reflect_prob;
      if(!entering) {
        outward_normal = -normal;
        ni_over_nt = ref_idx / current_ref_idx ;
      } else {
        outward_normal = normal;
        ni_over_nt = current_ref_idx / ref_idx;
      }
      vec3 unit_direction = unit_vector(r_in.direction());
      Float cos_theta = ffmin(dot(-unit_direction, outward_normal), (Float)1.0);
      Float sin_theta = sqrt(1.0 - cos_theta*cos_theta);
      if(ni_over_nt * sin_theta <= 1.0 ) {
        reflect_prob = schlick(cos_theta, ref_idx, current_ref_idx);
        reflect_prob = ni_over_nt == 1 ? 0 : reflect_prob;
      } else {
        reflect_prob = 1.0;
      }
      if(!entering) {
        Float distance = (hrec.p-r_in.point_at_parameter(0)).length();
        srec.attenuation = vec3(std::exp(-distance * attenuation.x()),
                                std::exp(-distance * attenuation.y()),
                                std::exp(-distance * attenuation.z()));
      } else {
        if(prev_active != -1) {
          Float distance = (hrec.p-r_in.point_at_parameter(0)).length();
          vec3 prev_atten = r_in.pri_stack->at(prev_active)->attenuation;

          srec.attenuation = albedo * vec3(std::exp(-distance * prev_atten.x()),
                                  std::exp(-distance * prev_atten.y()),
                                  std::exp(-distance * prev_atten.z()));
        } else {
          srec.attenuation = albedo;
        }
      }
      if(rng.unif_rand() < reflect_prob) {
        if(entering) {
          r_in.pri_stack->pop_back();
        }
        srec.specular_ray = ray(hrec.p, reflected, r_in.pri_stack, r_in.time());
      } else {
        if(!entering && current_layer != -1) {
          r_in.pri_stack->erase(r_in.pri_stack->begin() + current_layer);
        }
        vec3 refracted = refract(unit_direction, outward_normal, ni_over_nt);
        srec.specular_ray = ray(hrec.p, refracted, r_in.pri_stack, r_in.time());
      }
      return(true);
    }
    Float ref_idx;
    vec3 albedo;
    vec3 attenuation;
    size_t priority;
    random_gen rng;
};

class diffuse_light : public material {
public:
  diffuse_light(texture *a) : emit(a) {}
  ~diffuse_light() {
    if(emit) delete emit;
  }
  virtual bool scatter(const ray& r_in, const hit_record& rec, scatter_record& srec, random_gen& rng) {
    return(false);
  }
  virtual vec3 emitted(const ray& r_in, const hit_record& rec, Float u, Float v, const vec3& p) const {
    if(dot(rec.normal, r_in.direction()) < 0.0) {
      return(emit->value(u,v,p));
    } else {
      return(vec3(0,0,0));
    }
  }
  texture *emit;
};

class isotropic : public material {
public:
  isotropic(texture *a) : albedo(a) {}
  ~isotropic() {
    if(albedo) delete albedo;
  }
  virtual bool scatter(const ray& r_in, const hit_record& rec, scatter_record& srec, random_gen& rng) {
    srec.is_specular = true;
    srec.specular_ray = ray(rec.p, rng.random_in_unit_sphere(), r_in.pri_stack);
    srec.attenuation = albedo->value(rec.u,rec.v,rec.p);
    return(true);
  }
  vec3 f(const ray& r_in, const hit_record& rec, const ray& scattered) const {
    return(albedo->value(rec.u,rec.v,rec.p) * 0.25 * M_1_PI);
  }
  texture *albedo;
};

class orennayar : public material {
public:
  orennayar(texture *a, Float sigma) : albedo(a)  {
    Float sigma2 = sigma*sigma;
    A = 1.0f - (sigma2 / (2.0f * (sigma2 + 0.33f)));
    B = 0.45f * sigma2 / (sigma2 + 0.09f);
  }
  ~orennayar() {
    if(albedo) delete albedo;
  }
  bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng) {
    srec.is_specular = false;
    srec.attenuation = albedo->value(hrec.u, hrec.v, hrec.p);
    srec.pdf_ptr = new cosine_pdf(hrec.normal);
    return(true);
  }
vec3 f(const ray& r_in, const hit_record& rec, const ray& scattered) const {
    onb uvw;
    if(!rec.has_bump) {
      uvw.build_from_w(rec.normal);
    } else {
      uvw.build_from_w(rec.bump_normal);
    }
    vec3 wi = -unit_vector(uvw.world_to_local(r_in.direction()));
    vec3 wo = unit_vector(uvw.world_to_local(scattered.direction()));

    Float cosine = wo.z();
    
    if(cosine < 0) {
      cosine = 0;
    }
    
    Float sinThetaI = SinTheta(wi);
    Float sinThetaO = SinTheta(wo);
    Float maxCos = 0;
    if(sinThetaI > 1e-4 && sinThetaO > 1e-4) {
      Float sinPhiI = SinPhi(wi);
      Float cosPhiI = CosPhi(wi);
      Float sinPhiO = SinPhi(wo);
      Float cosPhiO = CosPhi(wo);
      Float dCos = cosPhiI * cosPhiO + sinPhiI * sinPhiO;
      maxCos = std::fmax((Float)0, dCos);
    }
    Float sinAlpha, tanBeta;
    if(AbsCosTheta(wi) > AbsCosTheta(wo)) {
      sinAlpha = sinThetaO;
      tanBeta = sinThetaI / AbsCosTheta(wi);
    } else {
      sinAlpha = sinThetaI;
      tanBeta = sinThetaO / AbsCosTheta(wo);
    }
    return(albedo->value(rec.u, rec.v, rec.p) * (A + B * maxCos * sinAlpha * tanBeta ) * cosine * M_1_PI );
  }
  Float A, B;
  texture *albedo;
};

class MicrofacetReflection : public material {
public:
  MicrofacetReflection(texture* a, MicrofacetDistribution *distribution, 
                       vec3 eta, vec3 k)
    : albedo(a), distribution(distribution), eta(eta), k(k) {}
  ~MicrofacetReflection() {
    if(albedo) delete albedo;
    if(distribution) delete distribution;
  }
  
  bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng) {
    srec.is_specular = false;
    srec.attenuation = albedo->value(hrec.u, hrec.v, hrec.p);
    if(!hrec.has_bump) {
      srec.pdf_ptr = new micro_pdf(hrec.normal, r_in.direction(), distribution);
    } else {
      srec.pdf_ptr = new micro_pdf(hrec.bump_normal, r_in.direction(), distribution);
    }
    return(true);
  }

  vec3 f(const ray& r_in, const hit_record& rec, const ray& scattered) const {
    onb uvw;
    if(!rec.has_bump) {
      uvw.build_from_w(rec.normal);
    } else {
      uvw.build_from_w(rec.bump_normal);
    }
    vec3 wi = -unit_vector(uvw.world_to_local(r_in.direction()));
    vec3 wo = unit_vector(uvw.world_to_local(scattered.direction()));
    
    Float cosThetaO = AbsCosTheta(wo);
    Float cosThetaI = AbsCosTheta(wi);
    vec3 normal = unit_vector(wi + wo);
    if (cosThetaI == 0 || cosThetaO == 0) {
      return(vec3(0,0,0));
    }
    if (normal.x() == 0 && normal.y() == 0 && normal.z() == 0) {
      return(vec3(0,0,0));
    }
    vec3 F = FrCond(cosThetaO, eta, k);
    Float G = distribution->G(wo,wi,normal);
    Float D = distribution->D(normal);
    return(albedo->value(rec.u, rec.v, rec.p) * F * G * D  * cosThetaI / (4 * CosTheta(wo) * CosTheta(wi) ));
  }
private:
  texture *albedo;
  MicrofacetDistribution *distribution;
  vec3 eta;
  vec3 k;
};

class glossy : public material {
public:
  glossy(texture* a, MicrofacetDistribution *distribution, 
         vec3 Rs, vec3 Rd2)
    : albedo(a), distribution(distribution), Rs(Rs) {}
  ~glossy() {
    if(albedo) delete albedo;
    if(distribution) delete distribution;
  }
  
  bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng) {
    srec.is_specular = false;
    srec.attenuation = albedo->value(hrec.u, hrec.v, hrec.p);
    srec.pdf_ptr = new glossy_pdf(hrec.normal, r_in.direction(), distribution);
    return(true);
  }
  
  vec3 f(const ray& r_in, const hit_record& rec, const ray& scattered) const {
    onb uvw;
    if(!rec.has_bump) {
      uvw.build_from_w(rec.normal);
    } else {
      uvw.build_from_w(rec.bump_normal);
    }
    vec3 wi = -unit_vector(uvw.world_to_local(r_in.direction()));
    vec3 wo = unit_vector(uvw.world_to_local(scattered.direction()));
    
    auto pow5 = [](Float v) { return (v * v) * (v * v) * v; };
    vec3 diffuse = (28.0f/(23.0f*M_PI)) * albedo->value(rec.u, rec.v, rec.p) *
      (vec3(1.0f) - Rs) *
      (1.0 - pow5(1 - 0.5f * AbsCosTheta(wi))) *
      (1.0 - pow5(1 - 0.5f * AbsCosTheta(wo)));
    vec3 wh = unit_vector(wi + wo);
    if (wh.x() == 0 && wh.y() == 0 && wh.z() == 0) {
      return(vec3(0.0f));
    }
    Float cosine  = dot(wh,wi);
    if(cosine < 0) {
      cosine = 0;
    }
    vec3 specular = distribution->D(wh) /
        (4 * AbsDot(wi, wh) *
        std::fmax(AbsCosTheta(wi), AbsCosTheta(wo))) *
        SchlickFresnel(dot(wo, wh));
    return((diffuse + specular) * cosine);
  }
  vec3 SchlickFresnel(Float cosTheta) const {
    auto pow5 = [](Float v) { return (v * v) * (v * v) * v; };
    return Rs + pow5(1 - cosTheta) * (vec3(1.0f) - Rs);
  }
  
private:
  texture *albedo;
  MicrofacetDistribution *distribution;
  vec3 Rs;
};

#endif
