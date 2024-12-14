#include "microfacetdist.h"


static void BeckmannSample11(Float cosThetaI, Float U1, Float U2,
                             Float *slope_x, Float *slope_y) {
  /* Special case (normal incidence) */
  if (cosThetaI > .9999) {
    Float r = std::sqrt(-std::log(static_cast<Float>(1) - U1));
    Float sinPhi = std::sin(2 * static_cast<Float>(M_PI) * U2);
    Float cosPhi = std::cos(2 * static_cast<Float>(M_PI) * U2);
    *slope_x = r * cosPhi;
    *slope_y = r * sinPhi;
    return;
  }
  
  /* The original inversion routine from the paper contained
   discontinuities, which causes issues for QMC integration
   and techniques like Kelemen-style MLT. The following code
   performs a numerical inversion with better behavior */
  Float sinThetaI = std::sqrt(std::fmax((Float)0, (Float)1 - cosThetaI * cosThetaI));
  Float tanThetaI = sinThetaI / cosThetaI;
  Float cotThetaI = 1 / tanThetaI;
  
  /* Search interval -- everything is parameterized
   in the Erf() domain */
  Float a = -1, c = Erf(cotThetaI);
  Float sample_x = std::fmax(U1, (Float)1e-6f);
  
  /* Start with a good initial guess */
  // Float b = (1-sample_x) * a + sample_x * c;
  
  /* We can do better (inverse of an approximation computed in
   * Mathematica) */
  Float thetaI = std::acos(cosThetaI);
  Float fit = 1 + thetaI * (-0.876f + thetaI * (0.4265f - 0.0594f * thetaI));
  Float b = c - (1 + c) * std::pow(1 - sample_x, fit);
  
  /* Normalization factor for the CDF */
  static const Float SQRT_PI_INV = static_cast<Float>(1) / std::sqrt(static_cast<Float>(M_PI));
  Float normalization = 1 / (1 + c + SQRT_PI_INV * tanThetaI * std::exp(-cotThetaI * cotThetaI));
  
  int it = 0;
  while (++it < 10) {
    /* Bisection criterion -- the oddly-looking
     Boolean expression are intentional to check
     for NaNs at little additional cost */
    if (!(b >= a && b <= c)) b = 0.5f * (a + c);
    
    /* Evaluate the CDF and its derivative
     (i.e. the density function) */
    Float invErf = ErfInv(b);
    Float value =
      normalization *
      (1 + b + SQRT_PI_INV * tanThetaI * std::exp(-invErf * invErf)) -
      sample_x;
    Float derivative = normalization * (1 - invErf * tanThetaI);
    
    if (std::fabs(value) < 1e-5f) break;
    
    /* Update bisection intervals */
    if (value > 0) {
      c = b;
    } else {
      a = b;
    }
    
    b -= value / derivative;
  }
  
  /* Now convert back into a slope value */
  *slope_x = ErfInv(b);
  
  /* Simulate Y component */
  *slope_y = ErfInv(2.0f * std::fmax(U2, (Float)1e-6f) - 1.0f);
}


static vec3f BeckmannSample(const vec3f &wi, Float alpha_x, Float alpha_y,
                           Float U1, Float U2) {
  // 1. stretch wi
  vec3f wiStretched = unit_vector(vec3f(alpha_x * wi.x(), alpha_y * wi.y(), wi.z()));
  
  // 2. simulate P22_{wi}(x_slope, y_slope, 1, 1)
  Float slope_x, slope_y;
  BeckmannSample11(CosTheta(wiStretched), U1, U2, &slope_x, &slope_y);
  
  // 3. rotate
  Float tmp = CosPhi(wiStretched) * slope_x - SinPhi(wiStretched) * slope_y;
  slope_y = SinPhi(wiStretched) * slope_x + CosPhi(wiStretched) * slope_y;
  slope_x = tmp;
  
  // 4. unstretch
  slope_x = alpha_x * slope_x;
  slope_y = alpha_y * slope_y;
  
  // 5. compute normal
  return(unit_vector(vec3f(-slope_x, -slope_y, 1.f)));
}


static void TrowbridgeReitzSample11(Float cosTheta, Float U1, Float U2,
                                    Float *slope_x, Float *slope_y) {
  // special case (normal incidence)
  if (cosTheta > 0.9999f) {
    Float r = sqrt(U1 / (1 - U1));
    Float phi = 2 * static_cast<Float>(M_PI) * U2;
    *slope_x = r * cos(phi);
    *slope_y = r * sin(phi);
    return;
  }
  
  Float sinTheta = std::sqrt(std::fmax((Float)0, (Float)1 - cosTheta * cosTheta));
  Float tanTheta = sinTheta / cosTheta;
  Float a = 1 / tanTheta;
  Float G1 = 2 / (1 + std::sqrt(static_cast<Float>(1) + static_cast<Float>(1) / (a * a)));
  
  // sample slope_x
  Float A = 2 * U1 / G1 - 1;
  Float tmp = 1.f / (A * A - 1.f);
  if (tmp > 1e10) {
    tmp = 1e10;
  }
  Float B = tanTheta;
  Float D = std::sqrt(std::fmax(Float(B * B * tmp * tmp - (A * A - B * B) * tmp), Float(0)));
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
  Float z = (U2 * (U2 * (U2 * 0.27385f - 0.73369f) + 0.46341f)) /
    (U2 * (U2 * (U2 * 0.093073f + 0.309420f) - 1.000000f) + 0.597999f);
  *slope_y = S * z * std::sqrt(static_cast<Float>(1) + *slope_x * *slope_x);
}

static vec3f TrowbridgeReitzSample(const vec3f &wi, Float alpha_x,
                                  Float alpha_y, Float U1, Float U2) {
  // 1. stretch wi
  vec3f wiStretched = unit_vector(vec3f(alpha_x * wi.x(), alpha_y * wi.y(), wi.z()));
  
  // 2. simulate P22_{wi}(x_slope, y_slope, 1, 1)
  Float slope_x, slope_y;
  TrowbridgeReitzSample11(CosTheta(wiStretched), U1, U2, &slope_x, &slope_y);
  
  // 3. rotate
  Float tmp = CosPhi(wiStretched) * slope_x - SinPhi(wiStretched) * slope_y;
  slope_y = SinPhi(wiStretched) * slope_x + CosPhi(wiStretched) * slope_y;
  slope_x = tmp;
  
  // 4. unstretch
  slope_x = alpha_x * slope_x;
  slope_y = alpha_y * slope_y;
  
  // 5. compute normal
  return(unit_vector(vec3f(-slope_x, -slope_y, 1.)));
}

Float TrowbridgeReitzDistribution::D(const vec3f &normal) const {
  const Float tan2Theta = Tan2Theta(normal);
  if(std::isinf(tan2Theta)) {
    return(0.0);
  }
  const Float cos4Theta = Cos2Theta(normal) * Cos2Theta(normal);
  Float e = tan2Theta * (Cos2Phi(normal) / (alphax_constant * alphax_constant) + 
    Sin2Phi(normal) / (alphay_constant * alphay_constant));
  return(1.0/(static_cast<Float>(M_PI) * alphax_constant * alphay_constant * cos4Theta * (1.0+e) * (1.0+e)));
}

Float TrowbridgeReitzDistribution::Lambda(const vec3f &w) const {
  Float absTanTheta = std::fabs(TanTheta(w));
  if (std::isinf(absTanTheta)) {
    return(0.0);
  }
  Float alpha = std::sqrt(Cos2Phi(w) * alphax_constant * alphax_constant + Sin2Phi(w) * alphay_constant * alphay_constant);
  Float alpha2Tan2Theta = (alpha * absTanTheta) * (alpha * absTanTheta);
  return((-1 + std::sqrt(static_cast<Float>(1) + alpha2Tan2Theta)) / 2);
}


vec3f TrowbridgeReitzDistribution::Sample_wh(const vec3f &wi,
                                            const Float u1, const Float u2) const {
  bool flip = wi.z() < 0;
  vec3f wh = TrowbridgeReitzSample(flip ? -wi : wi, alphax_constant, alphay_constant, u1, u2);
  if (flip) {
    wh = -wh;
  }
  return(wh);
}


Float BeckmannDistribution::D(const vec3f &wh) const {
  Float tan2Theta = Tan2Theta(wh);
  if(std::isinf(tan2Theta)) {
    return(0.0);
  }
  Float cos4Theta = Cos2Theta(wh) * Cos2Theta(wh);
  return(std::exp(-tan2Theta * (Cos2Phi(wh) / (alphax_constant * alphax_constant) + 
         Sin2Phi(wh) / (alphay_constant * alphay_constant))) / (static_cast<Float>(M_PI) * alphax_constant * alphay_constant * cos4Theta));
}

Float BeckmannDistribution::Lambda(const vec3f &w) const {
  Float absTanTheta = std::fabs(TanTheta(w));
  if (std::isinf(absTanTheta)) {
    return(0.0);
  }
  Float alpha = std::sqrt(Cos2Phi(w) * alphax_constant * alphax_constant + Sin2Phi(w) * alphay_constant * alphay_constant);
  Float a = 1 / (alpha * absTanTheta);
  if (a >= 1.6f) {
    return(0.0);
  }
  return((1 - 1.259f * a + 0.396f * a * a) / (3.535f * a + 2.181f * a * a));
}


vec3f BeckmannDistribution::Sample_wh(const vec3f &wi, const Float u1, const Float u2) const {
  // Sample visible area of normals for Beckmann distribution
  bool flip = wi.z() < 0;
  vec3f wh = BeckmannSample(flip ? -wi : wi, alphax_constant, alphay_constant, u1, u2);
  if (flip) {
    wh = -wh;
  }
  return(wh);
}

Float BeckmannDistribution::GetAlpha(Float u, Float v) const {
  point2f alphas = GetAlphas(u,v);
  return(std::sqrt(alphas.x() * alphas.x() + alphas.y() * alphas.y()));
}
point2f BeckmannDistribution::GetAlphas(Float u, Float v) const {
  return(!has_roughness ? point2f(alphax_constant,alphay_constant) : roughness->value(u,v));
}

Float TrowbridgeReitzDistribution::GetAlpha(Float u, Float v) const {
  point2f alphas = GetAlphas(u,v);
  return(std::sqrt(alphas.x() * alphas.x() + alphas.y() * alphas.y()));
}
point2f TrowbridgeReitzDistribution::GetAlphas(Float u, Float v) const {
  return(!has_roughness ? point2f(alphax_constant,alphay_constant) : roughness->value(u,v));
}

//Texture roughness

Float TrowbridgeReitzDistribution::D(const vec3f &normal, Float u, Float v) const {
  const Float tan2Theta = Tan2Theta(normal);
  if(std::isinf(tan2Theta)) {
    return(0.0);
  }
  const Float cos4Theta = Cos2Theta(normal) * Cos2Theta(normal);
  point2f alphas = GetAlphas(u,v);
  Float e = tan2Theta * (Cos2Phi(normal) / (alphas.x() * alphas.x()) + 
    Sin2Phi(normal) / (alphas.y() * alphas.y()));
  return(1.0/(static_cast<Float>(M_PI) * alphas.x() * alphas.y() * cos4Theta * (1.0+e) * (1.0+e)));
}

Float TrowbridgeReitzDistribution::Lambda(const vec3f &w, Float u, Float v) const {
  Float absTanTheta = std::fabs(TanTheta(w));
  if (std::isinf(absTanTheta)) {
    return(0.0);
  }
  point2f alphas = GetAlphas(u,v);
  
  Float alpha = std::sqrt(Cos2Phi(w) * alphas.x() * alphas.x() + Sin2Phi(w) * alphas.y() * alphas.y());
  Float alpha2Tan2Theta = (alpha * absTanTheta) * (alpha * absTanTheta);
  return((-1 + std::sqrt(static_cast<Float>(1) + alpha2Tan2Theta)) / 2);
}


vec3f TrowbridgeReitzDistribution::Sample_wh(const vec3f &wi,
                                             const Float u1, const Float u2, Float u, Float v) const {
  point2f alphas = GetAlphas(u,v);
  bool flip = wi.z() < 0;
  vec3f wh = TrowbridgeReitzSample(flip ? -wi : wi, alphas.x(), alphas.y(), u1, u2);
  if (flip) {
    wh = -wh;
  }
  return(wh);
}


Float BeckmannDistribution::D(const vec3f &wh, Float u, Float v) const {
  point2f alphas = GetAlphas(u,v);
  
  Float tan2Theta = Tan2Theta(wh);
  if(std::isinf(tan2Theta)) {
    return(0.0);
  }
  Float cos4Theta = Cos2Theta(wh) * Cos2Theta(wh);
  return(std::exp(-tan2Theta * (Cos2Phi(wh) / (alphas.x() * alphas.x()) + 
         Sin2Phi(wh) / (alphas.y() * alphas.y()))) / (static_cast<Float>(M_PI) * alphas.x() * alphas.y() * cos4Theta));
}

Float BeckmannDistribution::Lambda(const vec3f &w, Float u, Float v) const {
  point2f alphas = GetAlphas(u,v);
  
  Float absTanTheta = std::fabs(TanTheta(w));
  if (std::isinf(absTanTheta)) {
    return(0.0);
  }
  Float alpha = std::sqrt(Cos2Phi(w) * alphas.x() * alphas.x() + Sin2Phi(w) * alphas.y() * alphas.y());
  Float a = 1 / (alpha * absTanTheta);
  if (a >= 1.6f) {
    return(0.0);
  }
  return((1 - 1.259f * a + 0.396f * a * a) / (3.535f * a + 2.181f * a * a));
}


vec3f BeckmannDistribution::Sample_wh(const vec3f &wi, const Float u1, const Float u2, Float u, Float v) const {
  point2f alphas = GetAlphas(u,v);
  
  // Sample visible area of normals for Beckmann distribution
  bool flip = wi.z() < 0;
  vec3f wh = BeckmannSample(flip ? -wi : wi, alphas.x(), alphas.y(), u1, u2);
  if (flip) {
    wh = -wh;
  }
  return(wh);
}
