#ifndef MICROFACETDIST
#define MICROFACETDIST

#include "mathinline.h"

class MicrofacetDistribution {
public:
  virtual Float D(const vec3 &wh) const = 0;
  virtual Float Lambda(const vec3 &w) const = 0;
  Float G1(const vec3 &w) const {
    return(1 / (1 + Lambda(w)));
  }
  Float G(const vec3 &wo, const vec3 &wi) const {
    return(1 / (1 + Lambda(wo) + Lambda(wi)));
  }
  Float Pdf(const vec3 &wo, const vec3 &wi) const;
  Float GetAlpha() const;
protected:
  MicrofacetDistribution(bool sampleVisibleArea) 
      : sampleVisibleArea(sampleVisibleArea) {}
  const bool sampleVisibleArea;
};

class BeckmannDistribution : public MicrofacetDistribution {
public:
  static Float RoughnessToAlpha(Float roughness) {
    roughness = std::max(roughness, (Float)1e-3);
    Float x = std::log(roughness);
    return(1.62142f + 0.819955f * x + 0.1734f * x * x +
           0.0171201f * x * x * x + 0.000640711f * x * x * x * x);
  }
  BeckmannDistribution(Float alphax, Float alphay, bool samplevis = true)
    : MicrofacetDistribution(samplevis), alphax(alphax), alphay(alphay) {}
  Float D(const vec3 &wh) const;
  Float GetAlpha() const {
    return(std::sqrt(alphax * alphax + alphay * alphay));
  }
private:
  Float Lambda(const vec3 &w) const;
  const Float alphax, alphay;
};

Float BeckmannDistribution::D(const vec3 &wh) const {
  Float tan2Theta = Tan2Theta(wh);
  if(std::isinf(tan2Theta)) {
    return(0.0);
  }
  Float cos4Theta = Cos2Theta(wh) * Cos2Theta(wh);
  return(std::exp(-tan2Theta * (Cos2Phi(wh) / (alphax * alphax) + 
         Sin2Phi(wh) / (alphay * alphay))) / (M_PI * alphax * alphay * cos4Theta));
}

Float BeckmannDistribution::Lambda(const vec3 &w) const {
  Float absTanTheta = std::abs(TanTheta(w));
  if (std::isinf(absTanTheta)) {
    return(0.0);
  }
  Float alpha = std::sqrt(Cos2Phi(w) * alphax * alphax + Sin2Phi(w) * alphay * alphay);
  Float a = 1 / (alpha * absTanTheta);
  if (a >= 1.6f) {
    return(0.0);
  }
  return((1 - 1.259f * a + 0.396f * a * a) / (3.535f * a + 2.181f * a * a));
}

class TrowbridgeReitzDistribution : public MicrofacetDistribution {
public:
  static Float RoughnessToAlpha(Float roughness) {
    roughness = std::max(roughness, (Float)1e-3);
    Float x = std::log(roughness);
    return(1.62142f + 0.819955f * x + 0.1734f * x * x +
           0.0171201f * x * x * x + 0.000640711f * x * x * x * x);
  }
  TrowbridgeReitzDistribution(Float alphax, Float alphay, bool samplevis = true)
    : MicrofacetDistribution(samplevis), alphax(alphax), alphay(alphay) {}
  Float D(const vec3 &w) const;
  Float GetAlpha() const {
    return(std::sqrt(alphax * alphax + alphay * alphay));
  }
private:
  Float Lambda(const vec3 &w) const;
  const Float alphax, alphay;
};

Float TrowbridgeReitzDistribution::D(const vec3 &wh) const {
  Float tan2Theta = Tan2Theta(wh);
  if(std::isinf(tan2Theta)) {
    return(0.0);
  }
  const Float cos4Theta = Cos2Theta(wh) * Cos2Theta(wh);
  Float e = tan2Theta *(Cos2Phi(wh) / (alphax * alphax) + 
                        Sin2Phi(wh) / (alphay * alphay));
  return(1/(M_PI * alphax * alphay * cos4Theta * (1+e) * (1+e)));
}

Float TrowbridgeReitzDistribution::Lambda(const vec3 &w) const {
  Float absTanTheta = std::abs(TanTheta(w));
  if (std::isinf(absTanTheta)) {
    return(0.0);
  }
  Float alpha = std::sqrt(Cos2Phi(w) * alphax * alphax + Sin2Phi(w) * alphay * alphay);
  Float alpha2Tan2Theta = (alpha * absTanTheta) * (alpha * absTanTheta);
  return((-1 + std::sqrt(1.0f + alpha2Tan2Theta)) / 2);
}

#endif
