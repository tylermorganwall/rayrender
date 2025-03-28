#include "../hitables/infinite_area_light.h"
#include "../math/mathinline.h"
#include "../utils/raylog.h"
#include "../math/vectypes.h"

InfiniteAreaLight::InfiniteAreaLight(int width, int height, Float r, point3f center, 
                                     std::shared_ptr<texture> image, std::shared_ptr<material> mat,
                                     Transform* ObjectToWorld, Transform* WorldToObject, bool reverseOrientation)
                                     : hitable(ObjectToWorld, WorldToObject, mat, reverseOrientation), 
                                       width(width), height(height), radius(r), center(center) {
  //Set up distribution
  std::unique_ptr<Float[]>  img(new Float[width * height]);
  for (int v = 0; v < height; ++v) {
    Float vp = (Float)v / (Float)height;
    Float sinTheta = std::sin(static_cast<Float>(M_PI) * Float(v + .5f) / Float(height));
    for (int u = 0; u < width; ++u) {
      Float up = (Float)u / (Float)width;
      point3f rgb = image->value(up, vp, center);
      img[u + v * width] =  0.212671f*rgb.xyz.x + 0.715160f*rgb.xyz.y + 0.072169f*rgb.xyz.z;
      img[u + v * width] *= sinTheta;
    }
  }
  distribution = new Distribution2D(img.get(), width, height);
}

const bool InfiniteAreaLight::hit(const Ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("InfLight");
  
  Ray r2 = (*WorldToObject)(r);
  
  vec3f oc = r2.origin() - center;
  Float a = dot(r2.direction(), r2.direction());
  Float b = 2 * dot(oc, r2.direction()); 
  Float c = dot(oc,oc) - radius * radius;
  Float temp1, temp2;
  if (!quadratic(a, b, c, &temp1, &temp2)) {
    return(false);
  }

  if(temp1 < t_max && temp1 > t_min) {
    rec.t = temp1;
    rec.p = r2(rec.t);
    rec.p *= radius / rec.p.length(); 
    rec.normal = convert_to_normal3(-(r2.direction())); 
    
    vec3f v2(-r2.direction().xyz.z,r2.direction().xyz.y,r2.direction().xyz.x);
    get_sphere_uv(unit_vector(v2), rec.u, rec.v);
    rec.u = 1 - rec.u;
    
    //Interaction information
    // Float zRadius = std::sqrt(rec.p.xyz.x * rec.p.xyz.x  + rec.p.xyz.z  * rec.p.xyz.z );
    // Float invZRadius = 1 / zRadius;
    // Float cosPhi = rec.p.xyz.x * invZRadius;
    // Float sinPhi = rec.p.xyz.z * invZRadius;
    // Float theta = std::acos(clamp(rec.p.xyz.z / radius, -1, 1));
    // rec.dpdu = 2 * static_cast<Float>(M_PI) * vec3f(-rec.p.xyz.z, 0, rec.p.xyz.x);
    // rec.dpdv = 2 * static_cast<Float>(M_PI) * vec3f(rec.p.xyz.z * cosPhi, rec.p.xyz.z * sinPhi, -radius * std::sin(theta));
    
    rec = (*ObjectToWorld)(rec);
    rec.shape = this;
    rec.pError = vec3f(0,0,0);
    
    rec.mat_ptr = mat_ptr.get();
    rec.alpha_miss = false;
    rec.infinite_area_hit = true;
    return(true);
  }
  if(temp2 < t_max && temp2 > t_min) {
    rec.t = temp2;
    rec.p = r(rec.t);
    rec.p *= radius / rec.p.length(); 
    rec.normal = convert_to_normal3(-(r2.direction())); 
    
    vec3f v2(-r2.direction().xyz.z,r2.direction().xyz.y,r2.direction().xyz.x);
    get_sphere_uv(unit_vector(v2), rec.u, rec.v);
    rec.u = 1 - rec.u;
    
    
  //   //Interaction information
  //   Float zRadius = std::sqrt(rec.p.xyz.x * rec.p.xyz.x  + rec.p.xyz.z  * rec.p.xyz.z );
  //   Float invZRadius = 1 / zRadius;
  //   Float cosPhi = rec.p.xyz.x * invZRadius;
  //   Float sinPhi = rec.p.xyz.z * invZRadius;
  //   Float theta = std::acos(clamp(rec.p.xyz.z / radius, -1, 1));
  //   rec.dpdu = 2 * static_cast<Float>(M_PI) * vec3f(-rec.p.xyz.z, 0, rec.p.xyz.x);
  //   rec.dpdv = 2 * static_cast<Float>(M_PI) * vec3f(rec.p.xyz.z * cosPhi, rec.p.xyz.z * sinPhi, -radius * std::sin(theta));
    rec = (*ObjectToWorld)(rec);
    rec.shape = this;
    rec.pError = vec3f(0,0,0);

    rec.mat_ptr = mat_ptr.get();
    rec.alpha_miss = false;
    rec.infinite_area_hit = true;
    
    return(true);
  }
  return(false);
}


const bool InfiniteAreaLight::hit(const Ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("InfLight");
  
  Ray r2 = (*WorldToObject)(r);
  vec3f oc = r2.origin() - center;
  Float a = dot(r2.direction(), r2.direction());
  Float b = 2 * dot(oc, r2.direction()); 
  Float c = dot(oc,oc) - radius * radius;
  Float temp1, temp2;
  if (!quadratic(a, b, c, &temp1, &temp2)) {
    return(false);
  }

  if(temp1 < t_max && temp1 > t_min) {
    rec.t = temp1;
    rec.p = r2(rec.t);
  //   rec.p *= radius / rec.p.length(); 
  //   rec.normal = convert_to_normal3(-(r2.direction()));
    
    vec3f v2(-r2.direction().xyz.z,r2.direction().xyz.y,r2.direction().xyz.x);
    get_sphere_uv(unit_vector(v2), rec.u, rec.v);
    rec.u = 1 - rec.u;    
    
    
    // //Interaction information
    // Float zRadius = std::sqrt(rec.p.xyz.x * rec.p.xyz.x  + rec.p.xyz.z  * rec.p.xyz.z );
    // Float invZRadius = 1 / zRadius;
    // Float cosPhi = rec.p.xyz.x * invZRadius;
    // Float sinPhi = rec.p.xyz.z * invZRadius;
    // Float theta = std::acos(clamp(rec.p.xyz.z / radius, -1, 1));
    // rec.dpdu = 2 * static_cast<Float>(M_PI) * vec3f(-rec.p.xyz.z, 0, rec.p.xyz.x);
    // rec.dpdv = 2 * static_cast<Float>(M_PI) * vec3f(rec.p.xyz.z * cosPhi, rec.p.xyz.z * sinPhi, -radius * std::sin(theta));
    
    rec = (*ObjectToWorld)(rec);
    rec.shape = this;
    rec.pError = vec3f(0,0,0);
    rec.mat_ptr = mat_ptr.get();
    rec.alpha_miss = false;
    rec.infinite_area_hit = true;
    
    return(true);
  }
  if(temp2 < t_max && temp2 > t_min) {
    rec.t = temp2;
    rec.p = r2(rec.t);
  //   rec.p *= radius / rec.p.length(); 
  //   rec.normal = convert_to_normal3(-(r2.direction()));
    
    vec3f v2(-r2.direction().xyz.z,r2.direction().xyz.y,r2.direction().xyz.x);
    get_sphere_uv(unit_vector(v2), rec.u, rec.v);
    rec.u = 1 - rec.u;
    
  //   //Interaction information
  //   Float zRadius = std::sqrt(rec.p.xyz.x * rec.p.xyz.x  + rec.p.xyz.z  * rec.p.xyz.z );
  //   Float invZRadius = 1 / zRadius;
  //   Float cosPhi = rec.p.xyz.x * invZRadius;
  //   Float sinPhi = rec.p.xyz.z * invZRadius;
  //   Float theta = std::acos(clamp(rec.p.xyz.z / radius, -1, 1));
  //   rec.dpdu = 2 * static_cast<Float>(M_PI) * vec3f(-rec.p.xyz.z, 0, rec.p.xyz.x);
  //   rec.dpdv = 2 * static_cast<Float>(M_PI) * vec3f(rec.p.xyz.z * cosPhi, rec.p.xyz.z * sinPhi, -radius * std::sin(theta));
    
    rec = (*ObjectToWorld)(rec);
    rec.shape = this;
    rec.pError = vec3f(0,0,0);
    
    rec.mat_ptr = mat_ptr.get();
    rec.alpha_miss = false;
    rec.infinite_area_hit = true;
    
    return(true);
  }
  return(false);
}

bool InfiniteAreaLight::HitP(const Ray& r, Float t_min, Float t_max, random_gen& rng) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("InfLight");
  
  Ray r2 = (*WorldToObject)(r);
  
  vec3f oc = r2.origin() - center;
  Float a = dot(r2.direction(), r2.direction());
  Float b = 2 * dot(oc, r2.direction()); 
  Float c = dot(oc,oc) - radius * radius;
  Float temp1, temp2;
  if (!quadratic(a, b, c, &temp1, &temp2)) {
    return(false);
  }

  if(temp1 < t_max && temp1 > t_min) {
    return(true);
  }
  if(temp2 < t_max && temp2 > t_min) {
    return(true);
  }
  return(false);
}


bool InfiniteAreaLight::HitP(const Ray& r, Float t_min, Float t_max, Sampler* sampler) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("InfLight");
  
  Ray r2 = (*WorldToObject)(r);
  vec3f oc = r2.origin() - center;
  Float a = dot(r2.direction(), r2.direction());
  Float b = 2 * dot(oc, r2.direction()); 
  Float c = dot(oc,oc) - radius * radius;
  Float temp1, temp2;
  if (!quadratic(a, b, c, &temp1, &temp2)) {
    return(false);
  }

  if(temp1 < t_max && temp1 > t_min) {
    return(true);
  }
  if(temp2 < t_max && temp2 > t_min) {
    return(true);
  }
  return(false);
}

Float InfiniteAreaLight::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) {
  hit_record rec;
  // if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, rng)) {
    vec3f d = (*WorldToObject)(v);
    vec3f v2(-d.xyz.z,d.xyz.y,d.xyz.x);
    get_sphere_uv(unit_vector(v2), rec.u, rec.v);
    rec.u = 1 - rec.u;
    Float sinTheta = std::sin(rec.v * static_cast<Float>(M_PI));
    if (sinTheta == 0) {
      return(0);
    }
    //u = phi, v = theta
    return(distribution->Pdf(vec2f(rec.u, rec.v)) /
           (2 * static_cast<Float>(M_PI) * static_cast<Float>(M_PI) * sinTheta));
  // } else {
  //   return(0);
  // }
}


Float InfiniteAreaLight::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) {
  hit_record rec;

  // if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, sampler)) {
    vec3f d = (*WorldToObject)(v);
    vec3f v2(-d.xyz.z,d.xyz.y,d.xyz.x);
    get_sphere_uv(unit_vector(v2), rec.u, rec.v);
    rec.u = 1 - rec.u;
    Float sinTheta = std::sin(rec.v * static_cast<Float>(M_PI));
    if (sinTheta == 0) {
      return(0);
    }
    //u = phi, v = theta
    return(distribution->Pdf(vec2f(rec.u, rec.v)) /
           (2 * static_cast<Float>(M_PI) * static_cast<Float>(M_PI) * sinTheta));
  // } else {
  //   return(0);
  // }
}

vec3f InfiniteAreaLight::random(const point3f& o, random_gen& rng, Float time) {
  vec2f u(rng.unif_rand(), rng.unif_rand());
  Float mapPdf;
  vec2f uv = distribution->SampleContinuous(u, &mapPdf);
  if (mapPdf == 0) {
    return(vec3f(0.f,0.f,0.f));
  }
  //theta vertical, phi horizontal
  Float theta = (1-uv[1]) * static_cast<Float>(M_PI), phi = (1-uv[0]) * 2.0f * static_cast<Float>(M_PI);
  Float cosTheta = std::cos(theta), sinTheta = std::sin(theta);
  Float sinPhi = std::sin(phi), cosPhi = std::cos(phi);
  vec3f d(sinTheta * sinPhi, cosTheta, sinTheta * cosPhi);
  return((*ObjectToWorld)(d));
}

vec3f InfiniteAreaLight::random(const point3f& o, Sampler* sampler, Float time) {
  vec2f u = sampler->Get2D();
  Float mapPdf;
  vec2f uv = distribution->SampleContinuous(u, &mapPdf);
  if (mapPdf == 0) {
    return(vec3f(0.f,0.f,0.f));
  }
  //theta vertical, phi horizontal
  Float theta = (1-uv[1]) * static_cast<Float>(M_PI), phi = (1-uv[0]) * 2.0f * static_cast<Float>(M_PI);
  Float cosTheta = std::cos(theta), sinTheta = std::sin(theta);
  Float sinPhi = std::sin(phi), cosPhi = std::cos(phi);
  vec3f d(sinTheta * sinPhi, cosTheta, sinTheta * cosPhi);
  return((*ObjectToWorld)(d));
}


bool InfiniteAreaLight::bounding_box(Float t0, Float t1, aabb& box) const {
  box = (*ObjectToWorld)(aabb(-vec3f(radius,radius,radius), vec3f(radius,radius,radius)));
  return(true);
}

size_t InfiniteAreaLight::GetSize()  {
  return(mat_ptr ? 
           sizeof(*this) + distribution->GetSize() + sizeof(Float) * width * height + mat_ptr->GetSize() :
           sizeof(*this) + distribution->GetSize() + sizeof(Float) * width * height);
}
