#include "infinite_area_light.h"

InfiniteAreaLight::InfiniteAreaLight(int width, int height, Float r, vec3f center, 
                                     std::shared_ptr<texture> image, std::shared_ptr<material> mat,
                                     std::shared_ptr<Transform> ObjectToWorld, std::shared_ptr<Transform> WorldToObject, bool reverseOrientation)
                                     : hitable(ObjectToWorld, WorldToObject, reverseOrientation), 
                                       width(width), height(height), radius(r), center(center), mat_ptr(mat) {
  //Set up distribution
  std::unique_ptr<Float[]>  img(new Float[width * height]);
  for (int v = 0; v < height; ++v) {
    Float vp = (Float)v / (Float)height;
    Float sinTheta = std::sin(M_PI * Float(v + .5f) / Float(height));
    for (int u = 0; u < width; ++u) {
      Float up = (Float)u / (Float)width;
      point3f rgb = image->value(up, vp, center);
      img[u + v * width] =  0.212671f*rgb.r() + 0.715160f*rgb.g() + 0.072169f*rgb.b();
      img[u + v * width] *= sinTheta;
    }
  }
  distribution = new Distribution2D(img.get(), width, height);
}

bool InfiniteAreaLight::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  ray r2 = (*WorldToObject)(r);
  
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
    rec.p = r2.point_at_parameter(rec.t);
    rec.p *= radius / rec.p.length(); 
    rec.normal = -unit_vector(r2.direction());
    
    vec3f v2(-r2.direction().z(),r2.direction().y(),r2.direction().x());
    get_sphere_uv(unit_vector(v2), rec.u, rec.v);
    rec.u = 1 - rec.u;
    
    //Interaction information
    Float zRadius = std::sqrt(rec.p.x() * rec.p.x()  + rec.p.z()  * rec.p.z() );
    Float invZRadius = 1 / zRadius;
    Float cosPhi = rec.p.x() * invZRadius;
    Float sinPhi = rec.p.z() * invZRadius;
    Float theta = std::acos(clamp(rec.p.z() / radius, -1, 1));
    rec.dpdu = 2 * M_PI * vec3f(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = 2 * M_PI * vec3f(rec.p.z() * cosPhi, rec.p.z() * sinPhi, -radius * std::sin(theta));
    
    rec = (*ObjectToWorld)(rec);
    rec.shape = this;
    rec.pError = vec3f(0,0,0);
    
    rec.mat_ptr = mat_ptr.get();
    rec.alpha_miss = false;
    
    return(true);
  }
  if(temp2 < t_max && temp2 > t_min) {
    rec.t = temp2;
    rec.p = r.point_at_parameter(rec.t);
    rec.p *= radius / rec.p.length(); 
    rec.normal = -unit_vector(r2.direction());
    
    vec3f v2(-r2.direction().z(),r2.direction().y(),r2.direction().x());
    get_sphere_uv(unit_vector(v2), rec.u, rec.v);
    rec.u = 1 - rec.u;
    
    
    //Interaction information
    Float zRadius = std::sqrt(rec.p.x() * rec.p.x()  + rec.p.z()  * rec.p.z() );
    Float invZRadius = 1 / zRadius;
    Float cosPhi = rec.p.x() * invZRadius;
    Float sinPhi = rec.p.z() * invZRadius;
    Float theta = std::acos(clamp(rec.p.z() / radius, -1, 1));
    rec.dpdu = 2 * M_PI * vec3f(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = 2 * M_PI * vec3f(rec.p.z() * cosPhi, rec.p.z() * sinPhi, -radius * std::sin(theta));
    rec = (*ObjectToWorld)(rec);
    rec.shape = this;
    rec.pError = vec3f(0,0,0);

    rec.mat_ptr = mat_ptr.get();
    rec.alpha_miss = false;
    
    return(true);
  }
  return(false);
}


bool InfiniteAreaLight::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
  ray r2 = (*WorldToObject)(r);
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
    rec.p = r2.point_at_parameter(rec.t);
    rec.p *= radius / rec.p.length(); 
    rec.normal = -unit_vector(r2.direction());
    
    vec3f v2(-r2.direction().z(),r2.direction().y(),r2.direction().x());
    get_sphere_uv(unit_vector(v2), rec.u, rec.v);
    rec.u = 1 - rec.u;    
    
    
    //Interaction information
    Float zRadius = std::sqrt(rec.p.x() * rec.p.x()  + rec.p.z()  * rec.p.z() );
    Float invZRadius = 1 / zRadius;
    Float cosPhi = rec.p.x() * invZRadius;
    Float sinPhi = rec.p.z() * invZRadius;
    Float theta = std::acos(clamp(rec.p.z() / radius, -1, 1));
    rec.dpdu = 2 * M_PI * vec3f(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = 2 * M_PI * vec3f(rec.p.z() * cosPhi, rec.p.z() * sinPhi, -radius * std::sin(theta));
    
    rec = (*ObjectToWorld)(rec);
    rec.shape = this;
    rec.pError = vec3f(0,0,0);
    rec.mat_ptr = mat_ptr.get();
    rec.alpha_miss = false;
    
    return(true);
  }
  if(temp2 < t_max && temp2 > t_min) {
    rec.t = temp2;
    rec.p = r2.point_at_parameter(rec.t);
    rec.p *= radius / rec.p.length(); 
    rec.normal = -unit_vector(r2.direction());
    
    vec3f v2(-r2.direction().z(),r2.direction().y(),r2.direction().x());
    get_sphere_uv(unit_vector(v2), rec.u, rec.v);
    rec.u = 1 - rec.u;
    
    //Interaction information
    Float zRadius = std::sqrt(rec.p.x() * rec.p.x()  + rec.p.z()  * rec.p.z() );
    Float invZRadius = 1 / zRadius;
    Float cosPhi = rec.p.x() * invZRadius;
    Float sinPhi = rec.p.z() * invZRadius;
    Float theta = std::acos(clamp(rec.p.z() / radius, -1, 1));
    rec.dpdu = 2 * M_PI * vec3f(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = 2 * M_PI * vec3f(rec.p.z() * cosPhi, rec.p.z() * sinPhi, -radius * std::sin(theta));
    
    rec = (*ObjectToWorld)(rec);
    rec.shape = this;
    rec.pError = vec3f(0,0,0);
    
    rec.mat_ptr = mat_ptr.get();
    rec.alpha_miss = false;
    
    return(true);
  }
  return(false);
}

Float InfiniteAreaLight::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, rng)) {
    vec3f d = (*WorldToObject)(v);
    vec3f v2(-d.z(),d.y(),d.x());
    get_sphere_uv(unit_vector(v2), rec.u, rec.v);
    rec.u = 1 - rec.u;
    Float sinTheta = std::sin(rec.v * M_PI);
    if (sinTheta == 0) {
      return(0);
    }
    //u = phi, v = theta
    return(distribution->Pdf(vec2f(rec.u, rec.v)) /
           (2 * M_PI * M_PI * sinTheta));
  } else {
    return(0);
  }
}


Float InfiniteAreaLight::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, sampler)) {
    vec3f d = (*WorldToObject)(v);
    vec3f v2(-d.z(),d.y(),d.x());
    get_sphere_uv(unit_vector(v2), rec.u, rec.v);
    rec.u = 1 - rec.u;
    Float sinTheta = std::sin(rec.v * M_PI);
    if (sinTheta == 0) {
      return(0);
    }
    //u = phi, v = theta
    return(distribution->Pdf(vec2f(rec.u, rec.v)) /
           (2 * M_PI * M_PI * sinTheta));
  } else {
    return(0);
  }
}

vec3f InfiniteAreaLight::random(const point3f& o, random_gen& rng, Float time) {
  vec2f u(rng.unif_rand(), rng.unif_rand());
  Float mapPdf;
  vec2f uv = distribution->SampleContinuous(u, &mapPdf);
  if (mapPdf == 0) {
    return(vec3f(0.f,0.f,0.f));
  }
  //theta vertical, phi horizontal
  Float theta = (1-uv[1]) * M_PI, phi = (1-uv[0]) * 2.0f * M_PI;
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
  Float theta = (1-uv[1]) * M_PI, phi = (1-uv[0]) * 2.0f * M_PI;
  Float cosTheta = std::cos(theta), sinTheta = std::sin(theta);
  Float sinPhi = std::sin(phi), cosPhi = std::cos(phi);
  vec3f d(sinTheta * sinPhi, cosTheta, sinTheta * cosPhi);
  return((*ObjectToWorld)(d));
}


bool InfiniteAreaLight::bounding_box(Float t0, Float t1, aabb& box) const {
  box = (*ObjectToWorld)(aabb(-vec3f(radius,radius,radius), vec3f(radius,radius,radius)));
  return(true);
}
