#include "rough_dielectric.h"
#include "../materials/material.h"
#include "../math/mathinline.h"
#include "../utils/raylog.h"
#include "../utils/assert.h"

//Lambertian
//

inline bool Refract(const vec3f &wi, const normal3f &n, Float eta, vec3f *wt) {
  // Compute $\cos \theta_\roman{t}$ using Snell's law
  if(eta == 1) {
    *wt = -wi;
    return(true);
  }
  Float cosThetaI = dot(n, wi);
  Float sin2ThetaI = std::fmax(Float(0), Float(1 - cosThetaI * cosThetaI));
  Float sin2ThetaT = eta * eta * sin2ThetaI;
  
  // Handle total internal reflection for transmission
  if (sin2ThetaT >= 1) return(false);
  Float cosThetaT = std::sqrt(static_cast<Float>(1) - sin2ThetaT);
  *wt = eta * -wi + (eta * cosThetaI - cosThetaT) * convert_to_vec3(n);
  return(true);
}


//Metal
//

bool metal::scatter(const Ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng) {
  SCOPED_CONTEXT("Material");
  SCOPED_TIMER_COUNTER("Metal Scatter");
  
  
  normal3f normal = !hrec.has_bump ? hrec.normal : hrec.bump_normal;
  vec3f wi = -unit_vector(r_in.direction());
  vec3f reflected = Reflect(wi,normal);
  Float cosine = AbsDot(unit_vector(reflected), normal);
  if(cosine < 0) {
    cosine = 0;
  }
  point3f offset_p = offset_ray(hrec.p-r_in.o, hrec.normal) + r_in.o;
  
  srec.specular_ray = Ray(offset_p, reflected + fuzz * rng.random_in_unit_sphere(), r_in.pri_stack, r_in.time());
  srec.attenuation = albedo->value(hrec) * FrCond(cosine, eta, k);
  srec.is_specular = true;
  srec.pdf_ptr = 0;

  return(true);
}

bool metal::scatter(const Ray& r_in, const hit_record& hrec, scatter_record& srec, Sampler* sampler) {
  SCOPED_CONTEXT("Material");
  SCOPED_TIMER_COUNTER("Metal Scatter");
  
  normal3f normal = !hrec.has_bump ? hrec.normal : hrec.bump_normal;
  vec3f wi = -unit_vector(r_in.direction());
  vec3f reflected = Reflect(wi,normal);
  Float cosine = AbsDot(unit_vector(reflected), normal);
  if(cosine < 0) {
    cosine = 0;
  }
  point3f offset_p = offset_ray(hrec.p-r_in.o, hrec.normal) + r_in.o;
  
  srec.specular_ray = Ray(offset_p, reflected + fuzz * rand_to_unit(sampler->Get2D()), r_in.pri_stack, r_in.time());
  srec.attenuation = albedo->value(hrec) * FrCond(cosine, eta, k);
  srec.is_specular = true;
  srec.pdf_ptr = 0;

  return(true);
}

point3f metal::get_albedo(const hit_record& rec) const {
  return(albedo->value(rec));
}

size_t metal::GetSize()  {
  return(sizeof(*this));
}

//
//Dielectric
//

bool dielectric::scatter(const Ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng) {
  SCOPED_CONTEXT("Material");
  SCOPED_TIMER_COUNTER("Dielectric Scatter");
  
  srec.is_specular = true;
  normal3f outward_normal;
  normal3f wh = !hrec.has_bump ? hrec.normal : hrec.bump_normal;
  
  vec3f wi = -unit_vector(r_in.direction());
  
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
    //Determine current layer and continue to iterate through stack
    if(r_in.pri_stack->at(i) == this) {
      current_layer = i;
      continue;
    }
    //If layer's priority value less than active priority value, skip the layer
    if(r_in.pri_stack->at(i)->priority < active_priority_value) {
      active_priority_value = r_in.pri_stack->at(i)->priority;
      skip = true;
    }
    //If layer's priority value less than next_down_priority and the current layer isn't this material,
    //set the previous active layer to the current iterator on the stack
    if(r_in.pri_stack->at(i)->priority < next_down_priority && r_in.pri_stack->at(i) != this) {
      prev_active = i;
      next_down_priority = r_in.pri_stack->at(i)->priority;
    }
  }
  current_ref_idx = prev_active != -1 ? r_in.pri_stack->at(prev_active)->ref_idx : 1;
  
  outward_normal = entering ? wh : -wh;
  ni_over_nt = entering ? current_ref_idx / ref_idx : ref_idx / current_ref_idx ;
  
  //Never reflect if ni_over_nt == 1
  Float reflect_prob = FrDielectric(dot(wi, outward_normal), ni_over_nt);
  
  //If entering, push current layer to stack
  if(entering) {
    r_in.pri_stack->push_back(this);
  }
  // point3f offset_p = OffsetRayOrigin(hrec.p, hrec.pError, hrec.normal, r_in.direction());
  point3f offset_p = hrec.p;
  
  if(skip) {
    srec.specular_ray = Ray(offset_p, r_in.direction(), r_in.pri_stack, r_in.time());
    Float distance = (offset_p-r_in(0)).length();
    point3f prev_atten = r_in.pri_stack->at(prev_active)->attenuation;
    srec.attenuation = point3f(std::exp(-distance * prev_atten.xyz.x),
                               std::exp(-distance * prev_atten.xyz.y),
                               std::exp(-distance * prev_atten.xyz.z));
    if(r_in.segment_absorption) srec.attenuation = point3f(1);
    srec.is_passthrough = true;
    srec.is_transmission = true;
    if(!entering && current_layer != -1) {
      r_in.pri_stack->erase(r_in.pri_stack->begin() + static_cast<size_t>(current_layer));
    }
    return(true);
  }
  
  
  //Calculate attenuation color
  if(!entering) {
    Float distance = (offset_p-r_in(0)).length();
    srec.attenuation = point3f(std::exp(-distance * attenuation.xyz.x),
                               std::exp(-distance * attenuation.xyz.y),
                               std::exp(-distance * attenuation.xyz.z));
  } else {
    if(prev_active != -1) {
      Float distance = (offset_p-r_in(0)).length();
      point3f prev_atten = r_in.pri_stack->at(prev_active)->attenuation;
      
      srec.attenuation = albedo * point3f(std::exp(-distance * prev_atten.xyz.x),
                                          std::exp(-distance * prev_atten.xyz.y),
                                          std::exp(-distance * prev_atten.xyz.z));
    } else {
      srec.attenuation = albedo;
    }
  }
  if(r_in.segment_absorption) srec.attenuation = entering ? albedo : point3f(1);
  if(rng.unif_rand() < reflect_prob) {
    if(entering) {
      r_in.pri_stack->pop_back();
    }
    vec3f reflected = Reflect(wi, outward_normal);
    srec.specular_ray = Ray(offset_p, reflected, r_in.pri_stack, r_in.time());
  } else {
    if(!entering && current_layer != -1) {
      r_in.pri_stack->erase(r_in.pri_stack->begin() + current_layer);
    }
    vec3f refracted(0,0,0);
    bool valid_refraction = Refract(wi, outward_normal, ni_over_nt, &refracted);
    if(!valid_refraction && r_in.segment_absorption) return false;
    srec.is_transmission = true;
    srec.eta = 1/ni_over_nt;
    srec.specular_ray = Ray(offset_p, refracted, r_in.pri_stack, r_in.time());
  }
  return(true);
}


bool dielectric::scatter(const Ray& r_in, const hit_record& hrec, scatter_record& srec, Sampler* sampler) {
  SCOPED_CONTEXT("Material");
  SCOPED_TIMER_COUNTER("Dielectric Scatter");
  
  srec.is_specular = true;
  normal3f outward_normal;
  normal3f wh = !hrec.has_bump ? hrec.normal : hrec.bump_normal;
  
  vec3f wi = -unit_vector(r_in.direction());
  
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
    //Determine current layer and continue to iterate through stack
    if(r_in.pri_stack->at(i) == this) {
      current_layer = i;
      continue;
    }
    //If layer's priority value less than active priority value, skip the layer
    if(r_in.pri_stack->at(i)->priority < active_priority_value) {
      active_priority_value = r_in.pri_stack->at(i)->priority;
      skip = true;
    }
    //If layer's priority value less than next_down_priority and the current layer isn't this material,
    //set the previous active layer to the current iterator on the stack
    if(r_in.pri_stack->at(i)->priority < next_down_priority && r_in.pri_stack->at(i) != this) {
      prev_active = i;
      next_down_priority = r_in.pri_stack->at(i)->priority;
    }
  }
  point3f offset_p = OffsetRayOrigin(hrec.p, hrec.pError, hrec.normal, r_in.direction());
  current_ref_idx = prev_active != -1 ? r_in.pri_stack->at(prev_active)->ref_idx : 1;
  
  //If entering, push current layer to stack
  if(entering) {
    r_in.pri_stack->push_back(this);
  }

  if(skip) {
    srec.specular_ray = Ray(offset_p, r_in.direction(), r_in.pri_stack, r_in.time());
    Float distance = (offset_p-r_in(0)).length();
    point3f prev_atten = r_in.pri_stack->at(prev_active)->attenuation;
    srec.attenuation = point3f(std::exp(-distance * prev_atten.xyz.x),
                               std::exp(-distance * prev_atten.xyz.y),
                               std::exp(-distance * prev_atten.xyz.z));
    if(r_in.segment_absorption) srec.attenuation = point3f(1);
    srec.is_passthrough = true;
    srec.is_transmission = true;
    if(!entering && current_layer != -1) {
      r_in.pri_stack->erase(r_in.pri_stack->begin() + static_cast<size_t>(current_layer));
    }
    return(true);
  }
  
  outward_normal = entering ? wh : -wh;
  ni_over_nt = entering ? current_ref_idx / ref_idx : ref_idx / current_ref_idx ;
  
  //Never reflect if ni_over_nt == 1
  // Float reflect_prob = ni_over_nt != 1 ? FrDielectric(dot(wi, wh), ni_over_nt) : 0.0;
  Float reflect_prob = FrDielectric(dot(wi, outward_normal), ni_over_nt);
  
  //Calculate attenuation color
  if(!entering) {
    Float distance = (offset_p-r_in(0)).length();
    srec.attenuation = point3f(std::exp(-distance * attenuation.xyz.x),
                               std::exp(-distance * attenuation.xyz.y),
                               std::exp(-distance * attenuation.xyz.z));
  } else {
    if(prev_active != -1) {
      Float distance = (offset_p-r_in(0)).length();
      point3f prev_atten = r_in.pri_stack->at(prev_active)->attenuation;
      
      srec.attenuation = albedo * point3f(std::exp(-distance * prev_atten.xyz.x),
                                          std::exp(-distance * prev_atten.xyz.y),
                                          std::exp(-distance * prev_atten.xyz.z));
    } else {
      srec.attenuation = albedo;
    }
  }
  if(r_in.segment_absorption) srec.attenuation = entering ? albedo : point3f(1);
  if(sampler->Get1D() < reflect_prob) {
    if(entering) {
      r_in.pri_stack->pop_back();
    }
    vec3f reflected = Reflect(wi, outward_normal);
    srec.specular_ray = Ray(offset_p, reflected, r_in.pri_stack, r_in.time());
  } else {
    if(!entering && current_layer != -1) {
      r_in.pri_stack->erase(r_in.pri_stack->begin() + current_layer);
    }
    vec3f refracted(-wi);    
    bool valid_refraction = Refract(wi, outward_normal, ni_over_nt, &refracted);
    if(!valid_refraction && r_in.segment_absorption) return false;
    srec.is_transmission = true;
    srec.eta = 1/ni_over_nt;
    srec.specular_ray = Ray(offset_p, refracted, r_in.pri_stack, r_in.time());
  }
  return(true);
}

size_t dielectric::GetSize()  {
  return(sizeof(*this));
}

//
//Diffuse Light
//

point3f diffuse_light::emitted(const Ray& r_in, const hit_record& rec, Float u, Float v, const point3f& p, bool& is_invisible) {
  SCOPED_CONTEXT("Material");
  SCOPED_TIMER_COUNTER("Light Emit");
  
  is_invisible = invisible;
  if(dot(rec.normal, r_in.direction()) < 0.0) {
    return(emit->value(u,v,p) * intensity);
  } else {
    return(point3f(0,0,0));
  }
}

point3f diffuse_light::get_albedo(const hit_record& rec) const {
  return(emit->value(rec));
}

size_t diffuse_light::GetSize()  {
  return(sizeof(*this));
}

//
//Spot Light
//

point3f spot_light::emitted(const Ray& r_in, const hit_record& rec, Float u, Float v, const point3f& p, bool& is_invisible) {
  SCOPED_CONTEXT("Material");
  SCOPED_TIMER_COUNTER("Spotlight Emit");
  
  is_invisible = invisible;
  if(dot(rec.normal, r_in.direction()) < 0.0) {
    return(falloff(r_in.origin() - rec.p) * emit->value(u,v,p) * intensity );
  } else {
    return(point3f(0,0,0));
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

point3f spot_light::get_albedo(const hit_record& rec) const {
  return(emit->value(rec) * intensity);
}

size_t spot_light::GetSize()  {
  return(sizeof(*this));
}

//
//Isotropic
//

bool isotropic::scatter(const Ray& r_in, const hit_record& rec, scatter_record& srec, random_gen& rng) {
  SCOPED_CONTEXT("Material");
  SCOPED_TIMER_COUNTER("Isotropic Scatter");
  
  srec.is_specular = true;
  srec.specular_ray = Ray(rec.p, rng.random_in_unit_sphere(), r_in.pri_stack);
  srec.attenuation = albedo->value(rec);
  return(true);
}
bool isotropic::scatter(const Ray& r_in, const hit_record& rec, scatter_record& srec, Sampler* sampler) {
  SCOPED_CONTEXT("Material");
  SCOPED_TIMER_COUNTER("Isotropic Scatter");
  
  srec.is_specular = true;
  srec.specular_ray = Ray(rec.p, rand_to_sphere(1, 1, sampler->Get2D()), r_in.pri_stack);
  srec.attenuation = albedo->value(rec);
  return(true);
}

point3f isotropic::f(const Ray& r_in, const hit_record& rec, const vec3f& scattered) const {
  SCOPED_CONTEXT("Material");
  SCOPED_TIMER_COUNTER("Isotropic F");
  
  return(albedo->value(rec) * static_cast<Float>(0.25) * static_cast<Float>(M_1_PI));
}
point3f isotropic::get_albedo(const hit_record& rec) const {
  return(albedo->value(rec));
}

size_t isotropic::GetSize()  {
  return(sizeof(*this));
}

//
//MicrofacetReflection
//


bool MicrofacetReflection::scatter(const Ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng) {
  SCOPED_CONTEXT("Material");
  SCOPED_TIMER_COUNTER("MicrofacetReflection Scatter");
  
  srec.is_specular = false;
  srec.attenuation = albedo->value(hrec);
  if(!hrec.has_bump) {
    srec.pdf_ptr = new micro_pdf(hrec.normal, r_in.direction(), distribution, hrec.u, hrec.v);
    static_cast<micro_pdf*>(srec.pdf_ptr)->alphas = distribution->Resolve(hrec);
  } else {
    srec.pdf_ptr = new micro_pdf(hrec.bump_normal, r_in.direction(), distribution, hrec.u, hrec.v);
    static_cast<micro_pdf*>(srec.pdf_ptr)->alphas = distribution->Resolve(hrec);
  }
  return(true);
}

bool MicrofacetReflection::scatter(const Ray& r_in, const hit_record& hrec, scatter_record& srec, Sampler* sampler) {
  SCOPED_CONTEXT("Material");
  SCOPED_TIMER_COUNTER("MicrofacetReflection Scatter");
  
  srec.is_specular = false;
  srec.attenuation = albedo->value(hrec);
  if(!hrec.has_bump) {
    srec.pdf_ptr = new micro_pdf(hrec.normal, r_in.direction(), distribution, hrec.u, hrec.v);
    static_cast<micro_pdf*>(srec.pdf_ptr)->alphas = distribution->Resolve(hrec);
  } else {
    srec.pdf_ptr = new micro_pdf(hrec.bump_normal, r_in.direction(), distribution, hrec.u, hrec.v);
    static_cast<micro_pdf*>(srec.pdf_ptr)->alphas = distribution->Resolve(hrec);
  }
  return(true);
}

point3f MicrofacetReflection::f(const Ray& r_in, const hit_record& rec, const vec3f& scattered) const {
  SCOPED_CONTEXT("Material");
  SCOPED_TIMER_COUNTER("MicrofacetReflection F");
  
  onb uvw;
  if(!rec.has_bump) {
    uvw.build_from_w_normalized(rec.normal);
  } else {
    uvw.build_from_w_normalized(rec.bump_normal);
  }
  vec3f wi = -unit_vector(uvw.world_to_local(r_in.direction()));
  vec3f wo = unit_vector(uvw.world_to_local(scattered));
  
  Float cosThetaO = AbsCosTheta(wo);
  Float cosThetaI = AbsCosTheta(wi);
  vec3f normal = unit_vector(wi + wo);
  if (cosThetaI == 0 || cosThetaO == 0) {
    return(point3f(0,0,0));
  }
  if (normal.xyz.x == 0 && normal.xyz.y == 0 && normal.xyz.z == 0) {
    return(point3f(0,0,0));
  }
  point3f F = FrCond(cosThetaO, eta, k);
  Float G = distribution->G(wo,wi,normal);
  Float D = distribution->D(normal, distribution->Resolve(rec));
  return(albedo->value(rec) * F * G * D  * cosThetaI / (4 * CosTheta(wo) * CosTheta(wi) ));
}

point3f MicrofacetReflection::get_albedo(const hit_record& rec) const {
  return(albedo->value(rec));
}

size_t MicrofacetReflection::GetSize()  {
  return(sizeof(*this));
}



//
//MicrofacetTransmission
//

bool MicrofacetTransmission::Scatter(const Ray& ray, const hit_record& hit,
                                     scatter_record& result) {
  SCOPED_CONTEXT("Material");
  SCOPED_TIMER_COUNTER("MicrofacetTransmission Scatter");
  const normal3f normal = hit.has_bump ? hit.bump_normal : hit.normal;
  result.attenuation = albedo->value(hit);
  if (eta == 1) {
    // Index-matched roughness cannot change direction. Handle this as a delta
    // event instead of attempting to reconstruct a zero transmission half-vector.
    result.is_specular = true;
    result.is_transmission = true;
    result.eta = 1;
    if (dot(ray.direction(), normal) > 0) {
      Float distance = (hit.p - ray(0)).length();
      result.attenuation *= point3f(std::exp(-distance * k[0]),
        std::exp(-distance * k[1]), std::exp(-distance * k[2]));
    }
    result.specular_ray = Ray(OffsetRayOrigin(hit.p, hit.pError, hit.normal,
      ray.direction()), ray.direction(), ray.pri_stack, ray.time());
  } else {
    result.is_specular = false;
    auto density = new micro_transmission_pdf(normal, ray.direction(), distribution,
      eta, hit.u, hit.v);
    density->alphas = distribution->Resolve(hit);
    result.pdf_ptr = density;
  }
  return true;
}

bool MicrofacetTransmission::scatter(const Ray& ray, const hit_record& hit,
                                    scatter_record& result, random_gen&) {
  return Scatter(ray, hit, result);
}

bool MicrofacetTransmission::scatter(const Ray& ray, const hit_record& hit,
                                    scatter_record& result, Sampler*) {
  return Scatter(ray, hit, result);
}

point3f MicrofacetTransmission::f(const Ray& ray, const hit_record& hit,
                                const vec3f& scattered) const {
  SCOPED_CONTEXT("Material");
  SCOPED_TIMER_COUNTER("MicrofacetTransmission F");
  if (!(scattered.squared_length() > 0)) return point3f(0);
  onb frame;
  frame.build_from_w_normalized(hit.has_bump ? hit.bump_normal : hit.normal);
  const vec3f view = -unit_vector(frame.world_to_local(ray.direction()));
  const vec3f outgoing = unit_vector(frame.world_to_local(scattered));
  auto result = EvaluateRoughDielectric(view, outgoing, eta, *distribution,
    distribution->Resolve(hit));
  if (!result.transmission) return point3f(result.f_cos);
  point3f attenuation(1);
  if (view[2] < 0) {
    Float distance = (hit.p - ray(0)).length();
    attenuation = point3f(std::exp(-distance * k[0]), std::exp(-distance * k[1]),
      std::exp(-distance * k[2]));
  }
  return result.f_cos * albedo->value(hit) * attenuation;
}

point3f MicrofacetTransmission::get_albedo(const hit_record& rec) const {
  return(albedo->value(rec));
}

point3f MicrofacetTransmission::SchlickFresnel(Float cosTheta) const {
  auto pow5 = [](Float v) { return (v * v) * (v * v) * v; };
  return pow5(1 - cosTheta) * (point3f(1.0f));
}

size_t MicrofacetTransmission::GetSize()  {
  return(sizeof(*this));
}

//
//Glossy
//


bool glossy::scatter(const Ray& r_in, const hit_record& hrec, scatter_record& srec, random_gen& rng) {
  SCOPED_CONTEXT("Material");
  SCOPED_TIMER_COUNTER("Glossy Scatter");
  
  srec.is_specular = false;
  srec.attenuation = albedo->value(hrec);
  srec.pdf_ptr = new glossy_pdf(hrec.normal, r_in.direction(), distribution, hrec.u, hrec.v);
    static_cast<glossy_pdf*>(srec.pdf_ptr)->alphas = distribution->Resolve(hrec);
  return(true);
}

bool glossy::scatter(const Ray& r_in, const hit_record& hrec, scatter_record& srec, Sampler* sampler) {
  SCOPED_CONTEXT("Material");
  SCOPED_TIMER_COUNTER("Glossy Scatter");
  
  srec.is_specular = false;
  srec.attenuation = albedo->value(hrec);
  srec.pdf_ptr = new glossy_pdf(hrec.normal, r_in.direction(), distribution, hrec.u, hrec.v);
    static_cast<glossy_pdf*>(srec.pdf_ptr)->alphas = distribution->Resolve(hrec);
  return(true);
}

point3f glossy::f(const Ray& r_in, const hit_record& rec, const vec3f& scattered) const {
  SCOPED_CONTEXT("Material");
  SCOPED_TIMER_COUNTER("Glossy F");
  
  onb uvw;
  if(!rec.has_bump) {
    uvw.build_from_w_normalized(rec.normal);
  } else {
    uvw.build_from_w_normalized(rec.bump_normal);
  }
  vec3f wi = -unit_vector(uvw.world_to_local(r_in.direction()));
  vec3f wo = unit_vector(uvw.world_to_local(scattered));
  
  auto pow5 = [](Float v) { return (v * v) * (v * v) * v; };
  point3f diffuse = static_cast<Float>(28.0 / (23.0 * static_cast<Float>(M_PI))) * Rd * albedo->value(rec) *
    (point3f(1.0) + -Rs) *
    static_cast<Float>(1.0 - pow5(1 - 0.5f * AbsCosTheta(wi))) *
    static_cast<Float>(1.0 - pow5(1 - 0.5f * AbsCosTheta(wo)));
  vec3f wh = unit_vector(wi + wo);
  if (wh.xyz.x == 0 && wh.xyz.y == 0 && wh.xyz.z == 0) {
    return(point3f(0.0f));
  }
  Float cosine  = dot(wh,wi);
  if(cosine < 0 || !SameHemisphere(wi,wo)) {
    return(point3f(0));
  }
  point3f specular = distribution->D(wh, distribution->Resolve(rec)) /
      (4 * AbsDot(wi, wh) *
      std::fmax(AbsCosTheta(wi), AbsCosTheta(wo))) *
      SchlickFresnel(dot(wo, wh));
  return((diffuse + specular) * cosine );
}

point3f glossy::SchlickFresnel(Float cosTheta) const {
  auto pow5 = [](Float v) { return (v * v) * (v * v) * v; };
  return Rs + pow5(1 - cosTheta) * (point3f(1.0f) + -Rs);
}

point3f glossy::get_albedo(const hit_record& rec) const {
  return(albedo->value(rec));
}


size_t glossy::GetSize()  {
  return(sizeof(*this));
}

//
//Hair
//


namespace {
// Hair scattering uses the longitudinal tangent as x. Curve derivatives carry
// width/parameter scale and flat-curve normals need not be perpendicular to it.
onb HairFrame(const hit_record& hit) {
  vec3f tangent = hit.dpdu;
  if (!(tangent.squared_length() > 0)) tangent = vec3f(1, 0, 0);
  tangent = unit_vector(tangent);
  vec3f normal = convert_to_vec3(hit.normal);
  normal = normal - dot(normal, tangent) * tangent;
  if (normal.squared_length() < 1e-12f) {
    onb fallback;
    fallback.build_from_w(tangent);
    normal = fallback.v();
  } else normal = unit_vector(normal);
  return onb::FromXZ(tangent, normal);
}
}

bool hair::Scatter(const Ray& incoming, const hit_record& hit, scatter_record& scattering) {
  SCOPED_CONTEXT("Material");
  SCOPED_TIMER_COUNTER("Hair Scatter");
  const onb frame = HairFrame(hit);
  // Both directions point away from the interaction, as in PBRT HairBxDF.
  const vec3f outgoing = -unit_vector(frame.world_to_local(incoming.direction()));
  const Float h = clamp(-1 + 2 * hit.v, Float(-1), Float(1));
  scattering.is_specular = false;
  scattering.attenuation = point3f(1);
  scattering.pdf_ptr = new hair_pdf(frame, outgoing, eta, h, SafeASin(h), s,
                                   sigma_a, cos2kAlpha, sin2kAlpha, v);
  return true;
}

bool hair::scatter(const Ray& incoming, const hit_record& hit,
                   scatter_record& scattering, random_gen&) {
  return Scatter(incoming, hit, scattering);
}

bool hair::scatter(const Ray& incoming, const hit_record& hit,
                   scatter_record& scattering, Sampler*) {
  return Scatter(incoming, hit, scattering);
}

point3f hair::get_albedo(const hit_record&) const {
  return albedo;
}

point3f hair::f(const Ray& r_in, const hit_record& rec, const vec3f& scattered) const {
  SCOPED_CONTEXT("Material");
  SCOPED_TIMER_COUNTER("Hair F");
  
  onb uvw = HairFrame(rec);
  vec3f wo = -unit_vector(uvw.world_to_local(r_in.direction()));
  
  vec3f wi = unit_vector(uvw.world_to_local(scattered));
  
  Float h = clamp(-1 + 2 * rec.v, Float(-1), Float(1));
  Float gammaO = SafeASin(h);
  
  Float sinThetaO = wo.xyz.x;
  Float cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
  Float phiO = std::atan2(wo.xyz.z, wo.xyz.y);
  
  // Compute hair coordinate system terms related to _wi_
  Float sinThetaI = wi.xyz.x;
  Float cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));
  Float phiI = std::atan2(wi.xyz.z, wi.xyz.y);
  
  // Compute $\cos \thetat$ for refracted ray
  Float sinThetaT = sinThetaO / eta;
  Float cosThetaT = SafeSqrt(1 - Sqr(sinThetaT));
  
  // Compute $\gammat$ for refracted ray
  Float etap = std::sqrt(eta * eta - Sqr(sinThetaO)) / cosThetaO;
  Float sinGammaT = h / etap;
  Float cosGammaT = SafeSqrt(1 - Sqr(sinGammaT));
  Float gammaT = SafeASin(sinGammaT);
  
  // Compute the transmittance _T_ of a single path through the cylinder
  point3f T = Exp(-sigma_a * (2 * cosGammaT / cosThetaT));
  
  // Evaluate hair BSDF
  Float phi = phiI - phiO;
  std::array<point3f, pMax + 1> ap = Ap(cosThetaO, eta, h, T);
  point3f fsum(0.0);
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
    fsum += Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, v[p]) * ap[p] *
      Np(phi, p, s, gammaO, gammaT);
  }
  
  // Compute contribution of remaining terms after _pMax_
  fsum += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, v[pMax]) * ap[pMax] * ONE_OVER_2_PI;
  // rayrender consumes projected scattering (f * abs(cos)), whereas PBRT
  // returns raw f and applies this cosine in its integrator. Do not divide.
  return(fsum);
}

size_t hair::GetSize()  {
  return(sizeof(*this));
}
