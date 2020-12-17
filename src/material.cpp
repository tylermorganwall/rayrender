#include "material.h"
#include "mathinline.h"

bool dielectric::scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng) {
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
  vec3 offset_p = hrec.p;
  
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
    srec.specular_ray = ray(offset_p, r_in.direction(), r_in.pri_stack, r_in.time());
    Float distance = (offset_p-r_in.point_at_parameter(0)).length();
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
    Float distance = (offset_p-r_in.point_at_parameter(0)).length();
    srec.attenuation = vec3(std::exp(-distance * attenuation.x()),
                            std::exp(-distance * attenuation.y()),
                            std::exp(-distance * attenuation.z()));
  } else {
    if(prev_active != -1) {
      Float distance = (offset_p-r_in.point_at_parameter(0)).length();
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
    srec.specular_ray = ray(offset_p, reflected, r_in.pri_stack, r_in.time());
  } else {
    if(!entering && current_layer != -1) {
      r_in.pri_stack->erase(r_in.pri_stack->begin() + current_layer);
    }
    vec3 refracted = refract(unit_direction, outward_normal, ni_over_nt);
    srec.specular_ray = ray(offset_p, refracted, r_in.pri_stack, r_in.time());
  }
  return(true);
}

std::array<Float, pMax + 1> hair::ComputeApPdf(Float cosThetaO, Float h) const {
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
  return(apPdf);
}

bool hair::scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng) {
  onb uvw(hrec.dpdu, hrec.dpdv, hrec.normal);
  vec3 wo = unit_vector(uvw.world_to_local(r_in.direction()));

  Float sinThetaO = wo.x();
  Float cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
  Float phiO = std::atan2(wo.z(), wo.y());
  Float h = -1 + 2 * hrec.v;
  Float gammaO = SafeASin(h);
  
  // Derive four random samples from _u2_
  vec2 u2 = vec2(rng.unif_rand(),rng.unif_rand());
  vec2 u[2] = {DemuxFloat(u2.e[0]), DemuxFloat(u2.e[1])};
  
  // Determine which term $p$ to sample for hair scattering
  std::array<Float, pMax + 1> apPdf = ComputeApPdf(cosThetaO, h);
  int p;
  for (p = 0; p < pMax; ++p) {
    if (u[0].e[0] < apPdf[p]) break;
    u[0].e[0] -= apPdf[p];
  }
  
  // Rotate $\sin \thetao$ and $\cos \thetao$ to account for hair scale tilt
  Float sinThetaOp, cosThetaOp;
  if (p == 0) {
    sinThetaOp = sinThetaO * cos2kAlpha[1] - cosThetaO * sin2kAlpha[1];
    cosThetaOp = cosThetaO * cos2kAlpha[1] + sinThetaO * sin2kAlpha[1];
  } else if (p == 1) {
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
  u[1].e[0] = std::max(u[1].e[0], Float(1e-5));
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
  vec3 wi(sinThetaI, cosThetaI * std::cos(phiI),
          cosThetaI * std::sin(phiI));
  srec.is_specular = false;
  srec.attenuation = vec3(1,1,1);
  
  srec.pdf_ptr = new hair_pdf(uvw, wi, wo, 
                              eta, h, gammaO,  s, sigma_a,
                              cos2kAlpha, sin2kAlpha, v);
  return(true);
}

vec3 hair::f(const ray& r_in, const hit_record& rec, const ray& scattered) const {
  onb uvw(rec.dpdu, rec.dpdv, rec.normal);
  // vec3 wo = -unit_vector(uvw.world_to_local(r_in.direction()));
  vec3 wo = -unit_vector(uvw.world_to_local(r_in.direction()));
  
  vec3 wi = unit_vector(uvw.world_to_local(scattered.direction()));
  
  Float h = -1 + 2 * rec.v;
  Float gammaO = SafeASin(h);
  
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
  Float sinGammaT = h / etap;
  Float cosGammaT = SafeSqrt(1 - Sqr(sinGammaT));
  Float gammaT = SafeASin(sinGammaT);
  
  // Compute the transmittance _T_ of a single path through the cylinder
  vec3 T = Exp(-sigma_a * (2 * cosGammaT / cosThetaT));
  
  // Evaluate hair BSDF
  Float phi = phiI - phiO;
  std::array<vec3, pMax + 1> ap = Ap(cosThetaO, eta, h, T);
  vec3 fsum(0.0);
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
    fsum += Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, v[p]) * ap[p] *
      Np(phi, p, s, gammaO, gammaT);
  }
  
  // Compute contribution of remaining terms after _pMax_
  fsum += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, v[pMax]) * ap[pMax] * ONE_OVER_2_PI;
  if (AbsCosTheta(wi) > 0) {
    fsum /= AbsCosTheta(wi);
  }
  return(fsum);
}

// vec3 hair::SigmaAFromConcentration(Float ce, Float cp) {
//   vec3 sigma_a;
//   Float eumelaninSigmaA[3] = {0.419f, 0.697f, 1.37f};
//   Float pheomelaninSigmaA[3] = {0.187f, 0.4f, 1.05f};
//   for (int i = 0; i < 3; ++i)
//     sigma_a.e[i] = (ce * eumelaninSigmaA[i] + cp * pheomelaninSigmaA[i]);
//   return(sigma_a);
// }
// 
// vec3 hair::SigmaAFromReflectance(const vec3 &c, Float beta_n) {
//   vec3 sigma_a;
//   for (int i = 0; i < 3; ++i)
//     sigma_a.e[i] = Sqr(std::log(c.e[i]) /
//       (5.969f - 0.215f * beta_n + 2.532f * Sqr(beta_n) -
//         10.73f * Pow<3>(beta_n) + 5.574f * Pow<4>(beta_n) +
//         0.245f * Pow<5>(beta_n)));
//   return(sigma_a);
// }