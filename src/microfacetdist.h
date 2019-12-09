#ifndef MICROFACETDIST
#define MICROFACETDIST

#include "mathinline.h"
#include "rng.h"
#include "vec2.h"

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
  virtual Float GetAlpha() const = 0;
  virtual vec2 GetAlphas() const = 0;
  virtual bool GetType() const = 0;
protected:
  MicrofacetDistribution(bool sampleVisibleArea) 
      : sampleVisibleArea(sampleVisibleArea) {}
  const bool sampleVisibleArea;
  bool type;
};

class BeckmannDistribution : public MicrofacetDistribution {
public:
  static Float RoughnessToAlpha(Float roughness) {
    roughness = std::max(roughness, (Float)1e-3);
    Float x = std::log(roughness);
    return(1.62142f + 0.819955f * x + 0.1734f * x * x +
           0.0171201f * x * x * x + 0.000640711f * x * x * x * x);
  }
  BeckmannDistribution(Float alphax, Float alphay, bool type, bool samplevis = true)
    : MicrofacetDistribution(samplevis), alphax(alphax), alphay(alphay),type(type) {
  }
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
private:
  Float Lambda(const vec3 &w) const;
  const Float alphax, alphay;
  bool type;
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

// vec3 BeckmannDistribution::scatter(const vec3 &wo, random_gen &rng) const {
//   Float tan2Theta, phi;
//   Float u0 = rng.unif_rand();
//   Float u1 = rng.unif_rand();
//   if (alphax == alphay) {
//     Float logSample = std::log(1 - u0);
//     if (std::isinf(logSample)) {
//       logSample = 0;
//     }
//     tan2Theta = -alphax * alphax * logSample;
//     phi = u1 * 2 * M_PI;
//   } else {
//     Float logSample = std::log(u0);
//     phi = std::atan(alphay / alphax * std::tan(2 * M_PI * u1 + 0.5f * M_PI));
//     if (u1 > 0.5f) {
//       phi += M_PI;
//     }
//     Float sinPhi = std::sin(phi), cosPhi = std::cos(phi);
//     Float alphax2 = alphax * alphax, alphay2 = alphay * alphay;
//     tan2Theta = -logSample /
//       (cosPhi * cosPhi / alphax2 + sinPhi * sinPhi / alphay2);
//     
//   }
//   
//   Float cosTheta = 1 / std::sqrt(1 + tan2Theta);
//   Float sinTheta = std::sqrt(std::max((Float)0, 1 - cosTheta * cosTheta));
//   vec3 wh = SphericalDirection(sinTheta, cosTheta, phi);
//   if (!SameHemisphere(wo, wh)) wh = -wh;
//   return(wh);
// }

class TrowbridgeReitzDistribution : public MicrofacetDistribution {
public:
  TrowbridgeReitzDistribution(const Float alphax_, const Float alphay_, const bool type_, bool samplevis = true)
    : MicrofacetDistribution(samplevis) {
    // type(type), alphax(alphax), alphay(alphay)
    alphax = alphax_;
    alphay = alphay_;
    type = type_;
    // Rcpp::Rcout << alphax << " " << alphay << " " << type << "\n";
  }
  static Float RoughnessToAlpha(Float roughness) {
    roughness = std::max(roughness, (Float)1e-3);
    Float x = std::log(roughness);
    return(1.62142f + 0.819955f * x + 0.1734f * x * x +
           0.0171201f * x * x * x + 0.000640711f * x * x * x * x);
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
private:
  Float Lambda(const vec3 &w) const;
  Float alphax, alphay;
  bool type;
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
