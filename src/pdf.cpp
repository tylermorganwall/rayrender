#include "pdf.h"
#include "mathinline.h"


inline bool Refract(const vec3f &wi, const normal3f &n, Float eta, vec3f *wt) {
  // Compute $\cos \theta_\roman{t}$ using Snell's law
  Float cosThetaI = dot(n, wi);
  Float sin2ThetaI = std::fmax(Float(0), Float(1 - cosThetaI * cosThetaI));
  Float sin2ThetaT = eta * eta * sin2ThetaI;
  
  // Handle total internal reflection for transmission
  if (sin2ThetaT >= 1) return(false);
  Float cosThetaT = std::sqrt(1 - sin2ThetaT);
  *wt = eta * -wi + (eta * cosThetaI - cosThetaT) * vec3f(n.x(),n.y(),n.z());
  return(true);
}

Float cosine_pdf::value(const vec3f& direction, random_gen& rng, Float time) {
  Float cosine = dot(unit_vector(direction), uvw.w());
  if(cosine > 0) {
    return(cosine/M_PI);
  } else {
    return(0);
  }
}

Float cosine_pdf::value(const vec3f& direction, Sampler* sampler, Float time) {
  Float cosine = dot(unit_vector(direction), uvw.w());
  if(cosine > 0) {
    return(cosine/M_PI);
  } else {
    return(0);
  }
} 
vec3f cosine_pdf::generate(random_gen& rng, bool& diffuse_bounce, Float time) {
  diffuse_bounce = true;
  return(uvw.local_to_world(rng.random_cosine_direction()));
}

vec3f cosine_pdf::generate(Sampler* sampler, bool& diffuse_bounce, Float time) {
  diffuse_bounce = true;
  return(uvw.local_to_world(rand_cosine_direction(sampler->Get2D())));
}


Float micro_pdf::value(const vec3f& direction, random_gen& rng, Float time) {
  vec3f wo = unit_vector(uvw.world_to_local(direction));
  vec3f wh = unit_vector(wi + wo);
  return(distribution->Pdf(wo, wi, wh, u, v) / ( 4 * dot(wo, wh) ));
}

Float micro_pdf::value(const vec3f& direction, Sampler* sampler, Float time) {
  vec3f wo = unit_vector(uvw.world_to_local(direction));
  vec3f wh = unit_vector(wi + wo);
  return(distribution->Pdf(wo, wi, wh, u, v) / ( 4 * dot(wo, wh) ));
}

vec3f micro_pdf::generate(random_gen& rng, bool& diffuse_bounce, Float time) {
  vec3f wh = distribution->Sample_wh(wi, rng.unif_rand(), rng.unif_rand(), u, v);
  return(uvw.local_to_world(Reflect(wi, wh)));
}

vec3f micro_pdf::generate(Sampler* sampler, bool& diffuse_bounce, Float time) {
  vec2f u_r = sampler->Get2D();
  vec3f wh = distribution->Sample_wh(wi, u_r.x(), u_r.y(), u, v);
  return(uvw.local_to_world(Reflect(wi, wh)));
}

micro_transmission_pdf::micro_transmission_pdf(const normal3f& w, const vec3f& wi_, MicrofacetDistribution* distribution,
                       Float eta, Float uu, Float vv) : eta(eta), distribution(distribution), 
                       u(uu), v(vv) {
  uvw.build_from_w(w);
  wi = -unit_vector(uvw.world_to_local(wi_));
}

Float micro_transmission_pdf::value(const vec3f& direction, random_gen& rng, Float time) {
  vec3f wo = unit_vector(uvw.world_to_local(direction));
  Float cosTheta_o = CosTheta(wo), cosTheta_i = CosTheta(wi);
  bool reflect = cosTheta_i * cosTheta_o > 0;
  bool entering = CosTheta(wi) > 0;
  // Compute $\wh$ from $\wo$ and $\wi$ for microfacet transmission
  Float eta2 = 1;
  if (!reflect) {
    eta2 = entering ? (1.0/eta) : (eta);
  }
  vec3f wh = unit_vector(wi * eta2 + wo);

  wh = Faceforward(wh, normal3f(0, 0, 1));

  Float R = FrDielectric(-dot(wo,wh), eta);

  if (dot(wo, wh) * dot(wi, wh) > 0)  {
    return distribution->Pdf(wo, wi, wh, u, v) / ( 4 * AbsDot(wo,wh) ) * R;
  } 

  // Compute change of variables _dwh\_dwi_ for microfacet transmission
  Float sqrtDenom = eta2 * dot(wo, wh) + dot(wi, wh);
  Float dwh_dwi = std::fabs(eta2 * eta2 * dot(wo, wh) / (sqrtDenom * sqrtDenom));
  return distribution->Pdf(wo, wi, wh, u, v) * dwh_dwi * (1.0-R);
}

Float micro_transmission_pdf::value(const vec3f& direction, Sampler* sampler, Float time) {
  vec3f wo = unit_vector(uvw.world_to_local(direction));
  Float cosTheta_o = CosTheta(wo), cosTheta_i = CosTheta(wi);
  bool reflect = cosTheta_i * cosTheta_o > 0;
  bool entering = CosTheta(wi) > 0;
  // Compute $\wh$ from $\wo$ and $\wi$ for microfacet transmission
  Float eta2 = 1;
  if (!reflect) {
    eta2 = entering ? (1.0/eta) : (eta);
  }
  vec3f wh = unit_vector(wi  * eta2 + wo);
  wh = Faceforward(wh, normal3f(0, 0, 1));
  Float R = FrDielectric(-dot(wo,wh), eta);

  if (dot(wo, wh) * dot(wi, wh) > 0)  {
    return distribution->Pdf(wo, wi, wh, u, v) / ( 4 * AbsDot(wo,wh) ) * R;
  }

  // Compute change of variables _dwh\_dwi_ for microfacet transmission
  Float sqrtDenom = eta2 * dot(wo, wh) + dot(wi, wh);
  Float dwh_dwi = std::fabs(eta2 * eta2 * dot(wo, wh) / (sqrtDenom * sqrtDenom));
  return distribution->Pdf(wo, wi, wh, u, v) * dwh_dwi * (1.0-R);
}

inline Float schlick(Float cosine, Float ref_idx, Float ref_idx2) {
  Float r0 = (ref_idx2 - ref_idx) / (ref_idx2 + ref_idx);
  r0 = r0 * r0;
  return(r0 + (1-r0) * pow((1-cosine),5));
}

vec3f micro_transmission_pdf::generate(random_gen& rng, bool& diffuse_bounce, Float time) {
  vec3f wh = distribution->Sample_wh(wi, rng.unif_rand(),rng.unif_rand(), u, v);

  bool entering = CosTheta(wi) > 0;
  Float eta2 = entering ? (1.0/eta) : (eta);  
  Float R = FrDielectric(dot(wi,wh), eta2);
  vec3f dir;
  if (!Refract(wi, wh, eta2, &dir) || rng.unif_rand() < R) {
    dir = Reflect(wi,wh);
  }
  return(uvw.local_to_world(dir));
}

vec3f micro_transmission_pdf::generate(Sampler* sampler, bool& diffuse_bounce, Float time) {
  vec2f u_r = sampler->Get2D();
  normal3f wh = distribution->Sample_wh(wi, u_r.x(), u_r.y(), u, v);

  bool entering = CosTheta(wi) > 0;
  Float eta2 = entering ? (1.0/eta) : (eta);  
  Float R = FrDielectric(dot(wi,wh), eta2);
  vec3f dir;
  if (!Refract(wi, wh, eta2, &dir) || sampler->Get1D() < R) {
    dir = Reflect(wi,wh);
  }
  return(uvw.local_to_world(dir));
}


Float glossy_pdf::value(const vec3f& direction, random_gen& rng, Float time) {
  vec3f wo = unit_vector(uvw.world_to_local(direction));
  if(!SameHemisphere(wo,wi)) {
    return(0);
  }
  vec3f wh = unit_vector(wi + wo);
  return(0.5f * (AbsCosTheta(wi) * M_1_PI + distribution->Pdf(wo, wi, wh, u, v) / (4 * dot(wo, wh))));
}

Float glossy_pdf::value(const vec3f& direction, Sampler* sampler, Float time) {
  vec3f wo = unit_vector(uvw.world_to_local(direction));
  if(!SameHemisphere(wo,wi)) {
    return(0);
  }
  vec3f wh = unit_vector(wi + wo);
  return(0.5f * (AbsCosTheta(wi) * M_1_PI + distribution->Pdf(wo, wi, wh, u, v) / (4 * dot(wo, wh))) );
}

vec3f glossy_pdf::generate(random_gen& rng, bool& diffuse_bounce, Float time) {
  if(rng.unif_rand() < 0.5) {
    vec3f wh = distribution->Sample_wh(wi, rng.unif_rand(), rng.unif_rand(), u, v);
    return(uvw.local_to_world(Reflect(wi, wh)));
  } else {
    diffuse_bounce = true;
    vec3f wo = uvw.local_to_world(rng.random_cosine_direction());
    if (wi.z() < 0) wo.e[2] *= -1;
    return(wo);
  }
}
vec3f glossy_pdf::generate(Sampler* sampler, bool& diffuse_bounce, Float time) {
  if(sampler->Get1D() < 0.5) {
    vec2f u_r = sampler->Get2D();
    vec3f wh = distribution->Sample_wh(wi, u_r.x(), u_r.y(), u, v);
    return(uvw.local_to_world(Reflect(wi, wh)));
  } else {
    diffuse_bounce = true;
    vec3f wo = uvw.local_to_world(rand_cosine_direction(sampler->Get2D()));
    if (wi.z() < 0) wo.e[2] *= -1;
    return(wo);
  }
}

Float hitable_pdf::value(const vec3f& direction, random_gen& rng, Float time) {
  return(ptr->pdf_value(o, direction, rng, time));
}
Float hitable_pdf::value(const vec3f& direction, Sampler* sampler, Float time) {
  return(ptr->pdf_value(o, direction, sampler, time));
}
vec3f hitable_pdf::generate(random_gen& rng, bool& diffuse_bounce, Float time) {
  diffuse_bounce = true;
  return(ptr->random(o, rng, time)); 
}
vec3f hitable_pdf::generate(Sampler* sampler, bool& diffuse_bounce, Float time) {
  diffuse_bounce = true;
  return(ptr->random(o, sampler, time)); 
}


Float mixture_pdf::value(const vec3f& direction, random_gen& rng, Float time) {
  return(0.5 * p[0]->value(direction, rng, time) + 0.5 * p[1]->value(direction, rng, time));
}

Float mixture_pdf::value(const vec3f& direction, Sampler* sampler, Float time) {
  return(0.5 * p[0]->value(direction, sampler, time) + 0.5 * p[1]->value(direction, sampler, time));
}

vec3f mixture_pdf::generate(random_gen& rng, bool& diffuse_bounce, Float time) {
  if(rng.unif_rand() < 0.5) {
    return(p[0]->generate(rng, diffuse_bounce, time));
  } else {
    return(p[1]->generate(rng, diffuse_bounce, time));
  } 
}

vec3f mixture_pdf::generate(Sampler* sampler, bool& diffuse_bounce, Float time) {
  if(sampler->Get1D() < 0.5) {
    return(p[0]->generate(sampler, diffuse_bounce, time));
  } else {
    return(p[1]->generate(sampler, diffuse_bounce, time));
  } 
}

Float hair_pdf::value(const vec3f& direction, random_gen& rng, Float time) {
  Float sinThetaO = wo.x();
  Float cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
  Float phiO = std::atan2(wo.z(), wo.y());
  
  // Compute hair coordinate system terms related to _wi_
  Float sinThetaI = wi.x();
  Float cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));
  Float phiI = std::atan2(wi.z(), wi.y());
  
  // Compute $\gammat$ for refracted ray
  Float etap = std::sqrt(eta * eta - Sqr(sinThetaO)) / cosThetaO;
  Float sinGammaT = h / etap;
  Float gammaT = SafeASin(sinGammaT);
  
  // Compute PDF for $A_p$ terms
  std::array<Float, pMax + 1> apPdf = ComputeApPdf(cosThetaO);
  
  // Compute PDF sum for hair scattering events
  Float phi = phiI - phiO;
  Float pdf = 0.;
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
    cosThetaOp = std::fabs(cosThetaOp);
    pdf += Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, v[p]) *
      apPdf[p] * Np(phi, p, s, gammaO, gammaT);
  }
  pdf += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, v[pMax]) * 
    apPdf[pMax] * ONE_OVER_2_PI;
  return(pdf);
}


Float hair_pdf::value(const vec3f& direction, Sampler* sampler, Float time) {
  Float sinThetaO = wo.x();
  Float cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
  Float phiO = std::atan2(wo.z(), wo.y());
  
  // Compute hair coordinate system terms related to _wi_
  Float sinThetaI = wi.x();
  Float cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));
  Float phiI = std::atan2(wi.z(), wi.y());
  
  // Compute $\gammat$ for refracted ray
  Float etap = std::sqrt(eta * eta - Sqr(sinThetaO)) / cosThetaO;
  Float sinGammaT = h / etap;
  Float gammaT = SafeASin(sinGammaT);
  
  // Compute PDF for $A_p$ terms
  std::array<Float, pMax + 1> apPdf = ComputeApPdf(cosThetaO);
  
  // Compute PDF sum for hair scattering events
  Float phi = phiI - phiO;
  Float pdf = 0.;
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
    cosThetaOp = std::fabs(cosThetaOp);
    pdf += Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, v[p]) *
      apPdf[p] * Np(phi, p, s, gammaO, gammaT);
  }
  pdf += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, v[pMax]) * 
    apPdf[pMax] * ONE_OVER_2_PI;
  return(pdf);
}


vec3f hair_pdf::generate(random_gen& rng, bool& diffuse_bounce, Float time) {
  diffuse_bounce = true;
  Float sinThetaO = wo.x();
  Float cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
  Float phiO = std::atan2(wo.z(), wo.y());
  
  // Derive four random samples from _u2_
  vec2f u2 = vec2f(rng.unif_rand(),rng.unif_rand());
  vec2f u[2] = {DemuxFloat(u2.e[0]), DemuxFloat(u2.e[1])};
  std::array<Float, pMax + 1> apPdf = ComputeApPdf(cosThetaO);
  int p = 0;
  for (; p < pMax; ++p) {
    if (u[0].e[0] < apPdf[p]) break;
    u[0].e[0] -= apPdf[p];
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
  u[1].e[0] = std::fmax(u[1].e[0], Float(1e-5));
  Float cosTheta = 1 + v[p] * std::log(u[1].e[0] + (1 - u[1].e[0]) * std::exp(-2 / v[p]));
  Float sinTheta = SafeSqrt(1 - Sqr(cosTheta));
  Float cosPhi = std::cos(2 * M_PI * u[1].e[1]);
  Float sinThetaI = -cosTheta * sinThetaOp + sinTheta * cosPhi * cosThetaOp;
  Float cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));
  
  // Sample $N_p$ to compute $\Delta\phi$
  
  // Compute $\gammat$ for refracted ray
  Float etap = std::sqrt(eta * eta - Sqr(sinThetaO)) / cosThetaO;
  Float sinGammaT = h / etap;
  Float gammaT = SafeASin(sinGammaT);
  Float dphi;
  if (p < pMax) {
    dphi = Phi(p, gammaO, gammaT) + SampleTrimmedLogistic(u[0].e[1], s, -M_PI, M_PI);
  } else {
    dphi = 2 * M_PI * u[0].e[1];
  }
  // Compute _wi_ from sampled hair scattering angles
  Float phiI = phiO + dphi;
  return(uvw.local_to_world(vec3f(sinThetaI, cosThetaI * std::cos(phiI),
                                 cosThetaI * std::sin(phiI))));
}

vec3f hair_pdf::generate(Sampler* sampler, bool& diffuse_bounce, Float time) {
  diffuse_bounce = true;
  Float sinThetaO = wo.x();
  Float cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
  Float phiO = std::atan2(wo.z(), wo.y());
  
  // Derive four random samples from _u2_
  vec2f u2 = sampler->Get2D();
  vec2f u[2] = {DemuxFloat(u2.e[0]), DemuxFloat(u2.e[1])};
  std::array<Float, pMax + 1> apPdf = ComputeApPdf(cosThetaO);
  int p = 0;
  for (; p < pMax; ++p) {
    if (u[0].e[0] < apPdf[p]) break;
    u[0].e[0] -= apPdf[p];
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
  u[1].e[0] = std::fmax(u[1].e[0], Float(1e-5));
  Float cosTheta = 1 + v[p] * std::log(u[1].e[0] + (1 - u[1].e[0]) * std::exp(-2 / v[p]));
  Float sinTheta = SafeSqrt(1 - Sqr(cosTheta));
  Float cosPhi = std::cos(2 * M_PI * u[1].e[1]);
  Float sinThetaI = -cosTheta * sinThetaOp + sinTheta * cosPhi * cosThetaOp;
  Float cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));
  
  // Sample $N_p$ to compute $\Delta\phi$
  
  // Compute $\gammat$ for refracted ray
  Float etap = std::sqrt(eta * eta - Sqr(sinThetaO)) / cosThetaO;
  Float sinGammaT = h / etap;
  Float gammaT = SafeASin(sinGammaT);
  Float dphi;
  if (p < pMax) {
    dphi = Phi(p, gammaO, gammaT) + SampleTrimmedLogistic(u[0].e[1], s, -M_PI, M_PI);
  } else {
    dphi = 2 * M_PI * u[0].e[1];
  }
  // Compute _wi_ from sampled hair scattering angles
  Float phiI = phiO + dphi;
  return(uvw.local_to_world(vec3f(sinThetaI, cosThetaI * std::cos(phiI),
                                 cosThetaI * std::sin(phiI))));
}

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
  point3f T = Exp(-sigma_a * (2 * cosGammaT / cosThetaT));
  std::array<point3f, pMax + 1> ap = Ap(cosThetaO, eta, h, T);
  
  // Compute $A_p$ PDF from individual $A_p$ terms
  std::array<Float, pMax + 1> apPdf;
  Float sumY = std::accumulate(ap.begin(), ap.end(), Float(0),
                               [](Float s, const point3f &ap) { return s + ap.y(); });
  for (int i = 0; i <= pMax; ++i) {
    apPdf[i] = ap[i].y() / sumY;
  }
  return apPdf;
}
