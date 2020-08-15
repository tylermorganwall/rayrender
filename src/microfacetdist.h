#ifndef MICROFACETDIST
#define MICROFACETDIST

#include "mathinline.h"
#include "rng.h"
#include "vec2.h"

class MicrofacetDistribution {
public:
  virtual ~MicrofacetDistribution() {}
  virtual Float D(const vec3 &wh) const = 0;
  virtual Float Lambda(const vec3 &w) const = 0;
  Float G1(const vec3 &w) const {
    return(1 / (1 + Lambda(w)));
  }
  Float G(const vec3 &wo, const vec3 &wi, const vec3 &wh) const {
    Float NdotWh = AbsCosTheta(wh);
    Float NdotWo = AbsCosTheta(wo);
    Float NdotWi = AbsCosTheta(wi);
    Float WOdotWh = std::fabs(dot(wo, wh));
    return(std::fmin(1.f, std::fmin((2.0f * NdotWh * NdotWo / WOdotWh), (2.0f * NdotWh * NdotWi / WOdotWh))));
  }
  virtual Float GetAlpha() const = 0;
  virtual vec2 GetAlphas() const = 0;
  virtual bool GetType() const = 0;
  virtual vec3 Sample_wh(const vec3 &wi, const Float u1, const Float u2) const = 0;
  Float Pdf(const vec3 &wo,const vec3 &wi, const vec3 &wh) {
    return(D(wh) * G(wo, wi, wh) * AbsDot(wo, wh) / AbsCosTheta(wo));
  }
protected:
  MicrofacetDistribution(bool sampleVisibleArea) 
      : sampleVisibleArea(sampleVisibleArea) {}
  const bool sampleVisibleArea;
  bool type;
};

class BeckmannDistribution : public MicrofacetDistribution {
public:
  static Float RoughnessToAlpha(Float roughness) {
    roughness = std::fmax(roughness, (Float)0.0001550155);
    Float x = std::log(roughness);
    return(1.62142f + 0.819955f * x + 0.1734f * x * x +
           0.0171201f * x * x * x + 0.000640711f * x * x * x * x );
  }
  BeckmannDistribution(Float alphax_, Float alphay_, bool type, bool samplevis = true)
    : MicrofacetDistribution(samplevis), type(type) {
    alphax = RoughnessToAlpha(alphax_);
    alphay = RoughnessToAlpha(alphay_);
    alphax *= alphax;
    alphay *= alphay;
  }
  ~BeckmannDistribution() {}
  Float D(const vec3 &wh) const;
  Float GetAlpha() const {
    return(std::sqrt(alphax * alphax + alphay * alphay));
  }
  vec2 GetAlphas() const {
    return(vec2(alphax, alphay));
  }
  bool GetType() const { 
    return(type);
  };
  vec3 Sample_wh(const vec3 &wi, const Float u1, const Float u2) const;
private:
  Float Lambda(const vec3 &w) const;
  Float alphax, alphay;
  bool type;
};

class TrowbridgeReitzDistribution : public MicrofacetDistribution {
public:
  TrowbridgeReitzDistribution(const Float alphax_, const Float alphay_, const bool type_, bool samplevis = true)
    : MicrofacetDistribution(samplevis) {
    alphax = RoughnessToAlpha(alphax_);
    alphay = RoughnessToAlpha(alphay_);
    alphax *= alphax;
    alphay *= alphay;
    type = type_;
  }
  ~TrowbridgeReitzDistribution() {}
  static Float RoughnessToAlpha(Float roughness) {
    roughness = std::fmax(roughness, (Float)0.0001550155);
    Float x = std::log(roughness);
    return(1.62142f + 0.819955f * x + 0.1734f * x * x +
           0.0171201f * x * x * x + 0.000640711f * x * x * x * x );
  }
  Float D(const vec3 &w) const;
  Float GetAlpha() const {
    return(std::sqrt(alphax * alphax + alphay * alphay));
  }
  vec2 GetAlphas() const {
    return(vec2(alphax,alphay));
  }
  bool GetType() const { 
    return(type);
  };
  vec3 Sample_wh(const vec3 &wi, const Float u1, const Float u2) const;
private:
  Float Lambda(const vec3 &w) const;
  Float alphax, alphay;
  bool type;
};

#endif
