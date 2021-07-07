#include "material.h"
#include "mathinline.h"

//
//Lambertian
//

vec3f lambertian::f(const ray& r_in, const hit_record& rec, const ray& scattered) const {
  //unit_vector(scattered.direction()) == wo
  //r_in.direction() == wi
  vec3f wi = unit_vector(scattered.direction());
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
bool lambertian::scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng) {
  srec.is_specular = false;
  srec.attenuation = albedo->value(hrec.u, hrec.v, hrec.p);
  srec.pdf_ptr = new cosine_pdf(hrec.normal);
  return(true);
}
bool lambertian::scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, Sampler* sampler) {
  srec.is_specular = false;
  srec.attenuation = albedo->value(hrec.u, hrec.v, hrec.p);
  srec.pdf_ptr = new cosine_pdf(hrec.normal);
  return(true);
}
vec3f lambertian::get_albedo(const ray& r_in, const hit_record& rec) const {
  return(albedo->value(rec.u, rec.v, rec.p));
}

//
//Metal
//

bool metal::scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng) {
  vec3f normal = !hrec.has_bump ? hrec.normal : hrec.bump_normal;
  vec3f wi = -unit_vector(r_in.direction());
  vec3f reflected = Reflect(wi,unit_vector(normal));
  Float cosine = AbsDot(unit_vector(reflected), unit_vector(normal));
  if(cosine < 0) {
    cosine = 0;
  }
  vec3f offset_p = offset_ray(hrec.p-r_in.A, hrec.normal) + r_in.A;
  
  srec.specular_ray = ray(offset_p, reflected + fuzz * rng.random_in_unit_sphere(), r_in.pri_stack, r_in.time());
  srec.attenuation = albedo->value(hrec.u, hrec.v, hrec.p) * FrCond(cosine, eta, k);
  srec.is_specular = true;
  srec.pdf_ptr = 0;
  return(true);
}
bool metal::scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, Sampler* sampler) {
  vec3f normal = !hrec.has_bump ? hrec.normal : hrec.bump_normal;
  vec3f wi = -unit_vector(r_in.direction());
  vec3f reflected = Reflect(wi,unit_vector(normal));
  Float cosine = AbsDot(unit_vector(reflected), unit_vector(normal));
  if(cosine < 0) {
    cosine = 0;
  }
  vec3f offset_p = offset_ray(hrec.p-r_in.A, hrec.normal) + r_in.A;
  
  srec.specular_ray = ray(offset_p, reflected + fuzz * rand_to_unit(sampler->Get2D()), r_in.pri_stack, r_in.time());
  srec.attenuation = albedo->value(hrec.u, hrec.v, hrec.p) * FrCond(cosine, eta, k);
  srec.is_specular = true;
  srec.pdf_ptr = 0;
  return(true);
}
vec3f metal::get_albedo(const ray& r_in, const hit_record& rec) const {
  return(albedo->value(rec.u, rec.v, rec.p));
}

//
//Dielectric
//

bool dielectric::scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng) {
  srec.is_specular = true;
  vec3f outward_normal;
  vec3f normal = !hrec.has_bump ? hrec.normal : hrec.bump_normal;
  
  vec3f reflected = reflect(r_in.direction(), normal);
  Float ni_over_nt;
  srec.attenuation = albedo;
  Float current_ref_idx = 1.0;
  
  size_t active_priority_value = priority;
  size_t next_down_priority = 100000;
  int current_layer = -1; //keeping track of index of current material
  int prev_active = -1;  //keeping track of index of active material (higher priority)
  
  bool entering = dot(hrec.normal, r_in.direction()) < 0;
  bool skip = false;
  vec3f offset_p = hrec.p;
  
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
    vec3f prev_atten = r_in.pri_stack->at(prev_active)->attenuation;
    srec.attenuation = vec3f(std::exp(-distance * prev_atten.x()),
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
  vec3f unit_direction = unit_vector(r_in.direction());
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
    srec.attenuation = vec3f(std::exp(-distance * attenuation.x()),
                            std::exp(-distance * attenuation.y()),
                            std::exp(-distance * attenuation.z()));
  } else {
    if(prev_active != -1) {
      Float distance = (offset_p-r_in.point_at_parameter(0)).length();
      vec3f prev_atten = r_in.pri_stack->at(prev_active)->attenuation;
      
      srec.attenuation = albedo * vec3f(std::exp(-distance * prev_atten.x()),
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
    vec3f refracted = refract(unit_direction, outward_normal, ni_over_nt);
    srec.specular_ray = ray(offset_p, refracted, r_in.pri_stack, r_in.time());
  }
  return(true);
}


bool dielectric::scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, Sampler* sampler) {
  srec.is_specular = true;
  vec3f outward_normal;
  vec3f normal = !hrec.has_bump ? hrec.normal : hrec.bump_normal;
  
  vec3f reflected = reflect(r_in.direction(), normal);
  Float ni_over_nt;
  srec.attenuation = albedo;
  Float current_ref_idx = 1.0;
  
  size_t active_priority_value = priority;
  size_t next_down_priority = 100000;
  int current_layer = -1; //keeping track of index of current material
  int prev_active = -1;  //keeping track of index of active material (higher priority)
  
  bool entering = dot(hrec.normal, r_in.direction()) < 0;
  bool skip = false;
  vec3f offset_p = hrec.p;
  
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
    vec3f prev_atten = r_in.pri_stack->at(prev_active)->attenuation;
    srec.attenuation = vec3f(std::exp(-distance * prev_atten.x()),
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
  vec3f unit_direction = unit_vector(r_in.direction());
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
    srec.attenuation = vec3f(std::exp(-distance * attenuation.x()),
                            std::exp(-distance * attenuation.y()),
                            std::exp(-distance * attenuation.z()));
  } else {
    if(prev_active != -1) {
      Float distance = (offset_p-r_in.point_at_parameter(0)).length();
      vec3f prev_atten = r_in.pri_stack->at(prev_active)->attenuation;
      
      srec.attenuation = albedo * vec3f(std::exp(-distance * prev_atten.x()),
                                       std::exp(-distance * prev_atten.y()),
                                       std::exp(-distance * prev_atten.z()));
    } else {
      srec.attenuation = albedo;
    }
  }
  if(sampler->Get1D() < reflect_prob) {
    if(entering) {
      r_in.pri_stack->pop_back();
    }
    srec.specular_ray = ray(offset_p, reflected, r_in.pri_stack, r_in.time());
  } else {
    if(!entering && current_layer != -1) {
      r_in.pri_stack->erase(r_in.pri_stack->begin() + current_layer);
    }
    vec3f refracted = refract(unit_direction, outward_normal, ni_over_nt);
    srec.specular_ray = ray(offset_p, refracted, r_in.pri_stack, r_in.time());
  }
  return(true);
}

//
//Diffuse Light
//

vec3f diffuse_light::emitted(const ray& r_in, const hit_record& rec, Float u, Float v, const vec3f& p, bool& is_invisible) {
  is_invisible = invisible;
  if(dot(rec.normal, r_in.direction()) < 0.0) {
    return(emit->value(u,v,p) * intensity);
  } else {
    return(vec3f(0,0,0));
  }
}

vec3f diffuse_light::get_albedo(const ray& r_in, const hit_record& rec) const {
  return(emit->value(rec.u, rec.v, rec.p));
}

//
//Spot Light
//

vec3f spot_light::emitted(const ray& r_in, const hit_record& rec, Float u, Float v, const vec3f& p, bool& is_invisible) {
  is_invisible = invisible;
  if(dot(rec.normal, r_in.direction()) < 0.0) {
    return(falloff(r_in.origin() - rec.p) * emit->value(u,v,p) );
  } else {
    return(vec3f(0,0,0));
  }
}

Float spot_light::falloff(const vec3f &w) const {
  Float cosTheta = dot(spot_direction, unit_vector(w));
  if (cosTheta < cosTotalWidth) {
    return(0);
  }
  if (cosTheta > cosFalloffStart) {
    return(1);
  }
  Float delta = (cosTheta - cosTotalWidth) /(cosFalloffStart - cosTotalWidth);
  return((delta * delta) * (delta * delta));
}

vec3f spot_light::get_albedo(const ray& r_in, const hit_record& rec) const {
  return(emit->value(rec.u, rec.v, rec.p));
}

//
//Isotropic
//

bool isotropic::scatter(const ray& r_in, const hit_record& rec, scatter_record& srec, random_gen& rng) {
  srec.is_specular = true;
  srec.specular_ray = ray(rec.p, rng.random_in_unit_sphere(), r_in.pri_stack);
  srec.attenuation = albedo->value(rec.u,rec.v,rec.p);
  return(true);
}
bool isotropic::scatter(const ray& r_in, const hit_record& rec, scatter_record& srec, Sampler* sampler) {
  srec.is_specular = true;
  srec.specular_ray = ray(rec.p, rand_to_sphere(1, 1, sampler->Get2D()), r_in.pri_stack);
  srec.attenuation = albedo->value(rec.u,rec.v,rec.p);
  return(true);
}
vec3f isotropic::f(const ray& r_in, const hit_record& rec, const ray& scattered) const {
  return(albedo->value(rec.u,rec.v,rec.p) * 0.25 * M_1_PI);
}
vec3f isotropic::get_albedo(const ray& r_in, const hit_record& rec) const {
  return(albedo->value(rec.u, rec.v, rec.p));
}

//
//Oren Nayar
//

bool orennayar::scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng) {
  srec.is_specular = false;
  srec.attenuation = albedo->value(hrec.u, hrec.v, hrec.p);
  srec.pdf_ptr = new cosine_pdf(hrec.normal);
  return(true);
}

bool orennayar::scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, Sampler* sampler) {
  srec.is_specular = false;
  srec.attenuation = albedo->value(hrec.u, hrec.v, hrec.p);
  srec.pdf_ptr = new cosine_pdf(hrec.normal);
  return(true);
}

vec3f orennayar::f(const ray& r_in, const hit_record& rec, const ray& scattered) const {
  onb uvw;
  if(!rec.has_bump) {
    uvw.build_from_w(rec.normal);
  } else {
    uvw.build_from_w(rec.bump_normal);
  }
  vec3f wi = -unit_vector(uvw.world_to_local(r_in.direction()));
  vec3f wo = unit_vector(uvw.world_to_local(scattered.direction()));
  
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

vec3f orennayar::get_albedo(const ray& r_in, const hit_record& rec) const {
  return(albedo->value(rec.u, rec.v, rec.p));
}

//
//MicrofacetReflection
//


bool MicrofacetReflection::scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng) {
  srec.is_specular = false;
  srec.attenuation = albedo->value(hrec.u, hrec.v, hrec.p);
  if(!hrec.has_bump) {
    srec.pdf_ptr = new micro_pdf(hrec.normal, r_in.direction(), distribution);
  } else {
    srec.pdf_ptr = new micro_pdf(hrec.bump_normal, r_in.direction(), distribution);
  }
  return(true);
}

bool MicrofacetReflection::scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, Sampler* sampler) {
  srec.is_specular = false;
  srec.attenuation = albedo->value(hrec.u, hrec.v, hrec.p);
  if(!hrec.has_bump) {
    srec.pdf_ptr = new micro_pdf(hrec.normal, r_in.direction(), distribution);
  } else {
    srec.pdf_ptr = new micro_pdf(hrec.bump_normal, r_in.direction(), distribution);
  }
  return(true);
}

vec3f MicrofacetReflection::f(const ray& r_in, const hit_record& rec, const ray& scattered) const {
  onb uvw;
  if(!rec.has_bump) {
    uvw.build_from_w(rec.normal);
  } else {
    uvw.build_from_w(rec.bump_normal);
  }
  vec3f wi = -unit_vector(uvw.world_to_local(r_in.direction()));
  vec3f wo = unit_vector(uvw.world_to_local(scattered.direction()));
  
  Float cosThetaO = AbsCosTheta(wo);
  Float cosThetaI = AbsCosTheta(wi);
  vec3f normal = unit_vector(wi + wo);
  if (cosThetaI == 0 || cosThetaO == 0) {
    return(vec3f(0,0,0));
  }
  if (normal.x() == 0 && normal.y() == 0 && normal.z() == 0) {
    return(vec3f(0,0,0));
  }
  vec3f F = FrCond(cosThetaO, eta, k);
  Float G = distribution->G(wo,wi,normal);
  Float D = distribution->D(normal);
  return(albedo->value(rec.u, rec.v, rec.p) * F * G * D  * cosThetaI / (4 * CosTheta(wo) * CosTheta(wi) ));
}

vec3f MicrofacetReflection::get_albedo(const ray& r_in, const hit_record& rec) const {
  return(albedo->value(rec.u, rec.v, rec.p));
}

//
//Glossy
//


bool glossy::scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng) {
  srec.is_specular = false;
  srec.attenuation = albedo->value(hrec.u, hrec.v, hrec.p);
  srec.pdf_ptr = new glossy_pdf(hrec.normal, r_in.direction(), distribution);
  return(true);
}

bool glossy::scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, Sampler* sampler) {
  srec.is_specular = false;
  srec.attenuation = albedo->value(hrec.u, hrec.v, hrec.p);
  srec.pdf_ptr = new glossy_pdf(hrec.normal, r_in.direction(), distribution);
  return(true);
}

vec3f glossy::f(const ray& r_in, const hit_record& rec, const ray& scattered) const {
  onb uvw;
  if(!rec.has_bump) {
    uvw.build_from_w(rec.normal);
  } else {
    uvw.build_from_w(rec.bump_normal);
  }
  vec3f wi = -unit_vector(uvw.world_to_local(r_in.direction()));
  vec3f wo = unit_vector(uvw.world_to_local(scattered.direction()));
  
  auto pow5 = [](Float v) { return (v * v) * (v * v) * v; };
  vec3f diffuse = (28.0f/(23.0f*M_PI)) * albedo->value(rec.u, rec.v, rec.p) *
    (vec3f(1.0f) - Rs) *
    (1.0 - pow5(1 - 0.5f * AbsCosTheta(wi))) *
    (1.0 - pow5(1 - 0.5f * AbsCosTheta(wo)));
  vec3f wh = unit_vector(wi + wo);
  if (wh.x() == 0 && wh.y() == 0 && wh.z() == 0) {
    return(vec3f(0.0f));
  }
  Float cosine  = dot(wh,wi);
  if(cosine < 0) {
    cosine = 0;
  }
  vec3f specular = distribution->D(wh) /
    (4 * AbsDot(wi, wh) *
      std::fmax(AbsCosTheta(wi), AbsCosTheta(wo))) *
      SchlickFresnel(dot(wo, wh));
  return((diffuse + specular) * cosine);
}

vec3f glossy::SchlickFresnel(Float cosTheta) const {
  auto pow5 = [](Float v) { return (v * v) * (v * v) * v; };
  return Rs + pow5(1 - cosTheta) * (vec3f(1.0f) - Rs);
}

vec3f glossy::get_albedo(const ray& r_in, const hit_record& rec) const {
  return(albedo->value(rec.u, rec.v, rec.p));
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
  vec3f T = Exp(-sigma_a * (2 * cosGammaT / cosThetaT));
  std::array<vec3f, pMax + 1> ap = Ap(cosThetaO, eta, h, T);
  
  // Compute $A_p$ PDF from individual $A_p$ terms
  std::array<Float, pMax + 1> apPdf;
  Float sumY = std::accumulate(ap.begin(), ap.end(), Float(0),
                               [](Float s, const vec3f &ap) { return s + ap.y(); });
  for (int i = 0; i <= pMax; ++i) {
    apPdf[i] = ap[i].y() / sumY;
  }
  return(apPdf);
}

bool hair::scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng) {
  onb uvw(hrec.dpdu, hrec.dpdv, hrec.normal);
  vec3f wo = unit_vector(uvw.world_to_local(r_in.direction()));

  Float sinThetaO = wo.x();
  Float cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
  Float phiO = std::atan2(wo.z(), wo.y());
  Float h = -1 + 2 * hrec.v;
  Float gammaO = SafeASin(h);
  
  // Derive four random samples from _u2_
  vec2f u2 = vec2f(rng.unif_rand(),rng.unif_rand());
  vec2f u[2] = {DemuxFloat(u2.e[0]), DemuxFloat(u2.e[1])};
  
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
  vec3f wi(sinThetaI, cosThetaI * std::cos(phiI),
          cosThetaI * std::sin(phiI));
  srec.is_specular = false;
  srec.attenuation = vec3f(1,1,1);
  
  srec.pdf_ptr = new hair_pdf(uvw, wi, wo, 
                              eta, h, gammaO,  s, sigma_a,
                              cos2kAlpha, sin2kAlpha, v);
  return(true);
}


bool hair::scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec, Sampler* sampler) {
  onb uvw(hrec.dpdu, hrec.dpdv, hrec.normal);
  vec3f wo = unit_vector(uvw.world_to_local(r_in.direction()));
  
  Float sinThetaO = wo.x();
  Float cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
  Float phiO = std::atan2(wo.z(), wo.y());
  Float h = -1 + 2 * hrec.v;
  Float gammaO = SafeASin(h);
  
  // Derive four random samples from _u2_
  vec2f u2 = vec2f(sampler->Get1D(),sampler->Get1D());
  vec2f u[2] = {DemuxFloat(u2.e[0]), DemuxFloat(u2.e[1])};
  
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
  vec3f wi(sinThetaI, cosThetaI * std::cos(phiI),
          cosThetaI * std::sin(phiI));
  srec.is_specular = false;
  srec.attenuation = vec3f(1,1,1);
  
  srec.pdf_ptr = new hair_pdf(uvw, wi, wo, 
                              eta, h, gammaO,  s, sigma_a,
                              cos2kAlpha, sin2kAlpha, v);
  return(true);
}


vec3f hair::f(const ray& r_in, const hit_record& rec, const ray& scattered) const {
  onb uvw(rec.dpdu, rec.dpdv, rec.normal);
  // vec3f wo = -unit_vector(uvw.world_to_local(r_in.direction()));
  vec3f wo = -unit_vector(uvw.world_to_local(r_in.direction()));
  
  vec3f wi = unit_vector(uvw.world_to_local(scattered.direction()));
  
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
  vec3f T = Exp(-sigma_a * (2 * cosGammaT / cosThetaT));
  
  // Evaluate hair BSDF
  Float phi = phiI - phiO;
  std::array<vec3f, pMax + 1> ap = Ap(cosThetaO, eta, h, T);
  vec3f fsum(0.0);
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

// vec3f hair::SigmaAFromConcentration(Float ce, Float cp) {
//   vec3f sigma_a;
//   Float eumelaninSigmaA[3] = {0.419f, 0.697f, 1.37f};
//   Float pheomelaninSigmaA[3] = {0.187f, 0.4f, 1.05f};
//   for (int i = 0; i < 3; ++i)
//     sigma_a.e[i] = (ce * eumelaninSigmaA[i] + cp * pheomelaninSigmaA[i]);
//   return(sigma_a);
// }
// 
// vec3f hair::SigmaAFromReflectance(const vec3f &c, Float beta_n) {
//   vec3f sigma_a;
//   for (int i = 0; i < 3; ++i)
//     sigma_a.e[i] = Sqr(std::log(c.e[i]) /
//       (5.969f - 0.215f * beta_n + 2.532f * Sqr(beta_n) -
//         10.73f * Pow<3>(beta_n) + 5.574f * Pow<4>(beta_n) +
//         0.245f * Pow<5>(beta_n)));
//   return(sigma_a);
// }
