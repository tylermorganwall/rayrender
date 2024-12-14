#include "sphere.h"
#include "raylog.h"
#include "vectypes.h"

// #include "RcppThread.h"

const bool sphere::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("Sphere");
  vec3f oErr, dErr;
  ray r2 = (*WorldToObject)(r, &oErr, &dErr);
  // Compute quadratic sphere coefficients
  
  // Initialize _EFloat_ ray coordinate values
  EFloat ox(r2.origin().x(), oErr.x()), 
         oy(r2.origin().y(), oErr.y()), 
         oz(r2.origin().z(), oErr.z());
  EFloat dx(r2.direction().x(), dErr.x()), 
         dy(r2.direction().y(), dErr.y()), 
         dz(r2.direction().z(), dErr.z());
         
  EFloat a = dx * dx + dy * dy + dz * dz;
  EFloat b = 2 * (dx * ox + dy * oy + dz * oz);
  EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);
  
  // Solve quadratic equation for _t_ values
  EFloat temp1, temp2;
  if (!Quadratic(a, b, c, &temp1, &temp2)) {
    return(false);
  }
  bool is_hit = true;
  bool second_is_hit = true;
  bool alpha_miss = false;

  if(alpha_mask) {
    Float u;
    Float v;
    if(temp1 < t_max && temp1 > t_min) {
      point3f p1 = r2.point_at_parameter((Float)temp1);
      p1 *= radius / p1.length(); 
      vec3f normal = convert_to_vec3(p1 / radius);
      get_sphere_uv(normal, u, v);
      if(alpha_mask->value(u, v, rec.p) < rng.unif_rand()) {
        is_hit = false;
      }
    }
    if(temp2 < t_max && temp2 > t_min) {
      point3f p2 = r2.point_at_parameter((Float)temp2);
      p2 *= radius / p2.length(); 
      vec3f normal = convert_to_vec3(p2 / radius);
      get_sphere_uv(normal, u, v);
      if(alpha_mask->value(u, v, rec.p) < rng.unif_rand()) {
        if(!is_hit) {
          return(false);
        }
        second_is_hit = false;
      } 
    }
  }
  if(temp1 < t_max && temp1 > t_min && is_hit) {
    rec.t = (Float)temp1;
    rec.p = r2.point_at_parameter(rec.t);
    rec.p *= radius / rec.p.length(); 
    rec.normal = convert_to_normal3(rec.p) / radius;
    
    //Interaction information
    Float zRadius = std::sqrt(rec.p.x() * rec.p.x()  + rec.p.z()  * rec.p.z() );
    Float invZRadius = 1 / zRadius;
    Float cosPhi = rec.p.x() * invZRadius;
    Float sinPhi = rec.p.z() * invZRadius;
    Float theta = std::acos(clamp(rec.p.z() / radius, -1, 1));
    rec.dpdu = 2 * static_cast<Float>(M_PI) * vec3f(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = 2 * static_cast<Float>(M_PI) * vec3f(rec.p.z() * cosPhi, rec.p.z() * sinPhi, -radius * std::sin(theta));
    get_sphere_uv(rec.normal, rec.u, rec.v);
    rec.has_bump = bump_tex ? true : false;
    
    if(bump_tex) {
      point3f bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      rec.bump_normal = rec.normal + convert_to_normal3(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv); 
      rec.bump_normal.make_unit_vector();
    }
    rec.pError = convert_to_vec3(gamma(5) * Abs(rec.p));
    rec = (*ObjectToWorld)(rec);
    rec.normal *= reverseOrientation  ? -1 : 1;
    rec.bump_normal *= reverseOrientation  ? -1 : 1;
    rec.normal.make_unit_vector();
    rec.shape = this;
    rec.alpha_miss = alpha_miss;
    
    rec.mat_ptr = mat_ptr.get();
    return(true);
  }
  if(temp2 < t_max && temp2 > t_min && second_is_hit) {
    rec.t = (Float)temp2;
    rec.p = r2.point_at_parameter(rec.t);
    rec.p *= radius / rec.p.length();
    rec.normal = convert_to_normal3(rec.p) / radius;
    
    //Interaction information
    Float zRadius = std::sqrt(rec.p.x() * rec.p.x()  + rec.p.z()  * rec.p.z() );
    Float invZRadius = 1 / zRadius;
    Float cosPhi = rec.p.x() * invZRadius;
    Float sinPhi = rec.p.z() * invZRadius;
    Float theta = std::acos(clamp(rec.p.z() / radius, -1, 1));
    rec.dpdu = 2 * static_cast<Float>(M_PI) * vec3f(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = 2 * static_cast<Float>(M_PI) * vec3f(rec.p.z() * cosPhi, rec.p.z() * sinPhi, -radius * std::sin(theta));
    get_sphere_uv(rec.normal, rec.u, rec.v);
    rec.has_bump = bump_tex ? true : false;
    
    if(bump_tex) {
      point3f bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      rec.bump_normal = rec.normal +  convert_to_normal3(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv); 
      rec.bump_normal.make_unit_vector();
    }
    
    if(alpha_mask) {
      rec.normal = -rec.normal;
      rec.bump_normal = -rec.bump_normal;
    }
    rec.pError = convert_to_vec3(gamma(5) * Abs(rec.p));
    rec = (*ObjectToWorld)(rec);
    rec.normal *= reverseOrientation  ? -1 : 1;
    rec.bump_normal *= reverseOrientation  ? -1 : 1;
    rec.normal.make_unit_vector();
    
    rec.shape = this;
    rec.alpha_miss = alpha_miss;
    
    rec.mat_ptr = mat_ptr.get();
    return(true);
  }
  return(false);
}


const bool sphere::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("Sphere");
  
  vec3f oErr, dErr;
  ray r2 = (*WorldToObject)(r, &oErr, &dErr);
  // Compute quadratic sphere coefficients
  
  // Initialize _EFloat_ ray coordinate values
  EFloat ox(r2.origin().x(), oErr.x()), oy(r2.origin().y(), oErr.y()), oz(r2.origin().z(), oErr.z());
  EFloat dx(r2.direction().x(), dErr.x()), dy(r2.direction().y(), dErr.y()), dz(r2.direction().z(), dErr.z());
  EFloat a = dx * dx + dy * dy + dz * dz;
  EFloat b = 2 * (dx * ox + dy * oy + dz * oz);
  EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);
  
  // Solve quadratic equation for _t_ values
  EFloat temp1, temp2;
  if (!Quadratic(a, b, c, &temp1, &temp2)) {
    return(false);
  }
  bool is_hit = true;
  bool second_is_hit = true;
  bool alpha_miss = false;
  if(alpha_mask) {
    Float u;
    Float v;
    if(temp1 < t_max && temp1 > t_min) {
      point3f p1 = r2.point_at_parameter((Float)temp1);
      p1 *= radius / p1.length(); 
      vec3f normal = convert_to_vec3(p1) / radius;
      get_sphere_uv(normal, u, v);
      if(alpha_mask->value(u, v, rec.p) < sampler->Get1D()) {
        is_hit = false;
      }
    }
    if(temp2 < t_max && temp2 > t_min) {
      point3f p2 = r2.point_at_parameter((Float)temp2);
      p2 *= radius / p2.length(); 
      vec3f normal = convert_to_vec3(p2) / radius;
      get_sphere_uv(normal, u, v);
      if(alpha_mask->value(u, v, rec.p) < sampler->Get1D()) {
        if(!is_hit) {
          alpha_miss = true;
        }
        second_is_hit = false;
      } 
    }
  }
  if(temp1 < t_max && temp1 > t_min && is_hit) {
    rec.t = (Float)temp1;
    rec.p = r2.point_at_parameter(rec.t);
    rec.p *= radius / rec.p.length(); 
    rec.normal = convert_to_normal3(rec.p) / radius;
    
    //Interaction information
    Float zRadius = std::sqrt(rec.p.x() * rec.p.x()  + rec.p.z()  * rec.p.z() );
    Float invZRadius = 1 / zRadius;
    Float cosPhi = rec.p.x() * invZRadius;
    Float sinPhi = rec.p.z() * invZRadius;
    Float theta = std::acos(clamp(rec.p.z() / radius, -1, 1));
    rec.dpdu = 2 * static_cast<Float>(M_PI) * vec3f(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = 2 * static_cast<Float>(M_PI) * vec3f(rec.p.z() * cosPhi, rec.p.z() * sinPhi, -radius * std::sin(theta));
    get_sphere_uv(rec.normal, rec.u, rec.v);
    rec.has_bump = bump_tex ? true : false;
    
    if(bump_tex) {
      point3f bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      rec.bump_normal = rec.normal + convert_to_normal3(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv); 
      rec.bump_normal.make_unit_vector();
    }
    rec.pError = convert_to_vec3(gamma(5) * Abs(rec.p));
    rec = (*ObjectToWorld)(rec);
    rec.normal *= reverseOrientation  ? -1 : 1;
    rec.bump_normal *= reverseOrientation  ? -1 : 1;
    rec.normal.make_unit_vector();
    
    rec.shape = this;
    rec.alpha_miss = alpha_miss;
    
    rec.mat_ptr = mat_ptr.get();
    return(true);
  }
  if(temp2 < t_max && temp2 > t_min && second_is_hit) {
    rec.t = (Float)temp2;
    rec.p = r2.point_at_parameter(rec.t);
    rec.p *= radius / rec.p.length();
    rec.normal = convert_to_normal3(rec.p) / radius;
    
    //Interaction information
    Float zRadius = std::sqrt(rec.p.x() * rec.p.x()  + rec.p.z()  * rec.p.z() );
    Float invZRadius = 1 / zRadius;
    Float cosPhi = rec.p.x() * invZRadius;
    Float sinPhi = rec.p.z() * invZRadius;
    Float theta = std::acos(clamp(rec.p.z() / radius, -1, 1));
    rec.dpdu = 2 * static_cast<Float>(M_PI) * vec3f(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = 2 * static_cast<Float>(M_PI) * vec3f(rec.p.z() * cosPhi, rec.p.z() * sinPhi, -radius * std::sin(theta));
    get_sphere_uv(rec.normal, rec.u, rec.v);
    rec.has_bump = bump_tex ? true : false;
    
    if(bump_tex) {
      point3f bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      rec.bump_normal = rec.normal +  convert_to_normal3(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv); 
      rec.bump_normal.make_unit_vector();
    }
    
    if(alpha_mask) {
      rec.normal = -rec.normal;
      rec.bump_normal = -rec.bump_normal;
    }
    rec.pError = convert_to_vec3(gamma(5) * Abs(rec.p));
    rec = (*ObjectToWorld)(rec);
    rec.normal *= reverseOrientation  ? -1 : 1;
    rec.bump_normal *= reverseOrientation  ? -1 : 1;
    rec.normal.make_unit_vector();
    
    rec.shape = this;
    rec.alpha_miss = alpha_miss;
    
    rec.mat_ptr = mat_ptr.get();
    return(true);
  }
  return(false);
}

bool sphere::HitP(const ray& r, Float t_min, Float t_max, random_gen& rng) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("Sphere");
  vec3f oErr, dErr;
  ray r2 = (*WorldToObject)(r, &oErr, &dErr);
  // Compute quadratic sphere coefficients
  
  // Initialize _EFloat_ ray coordinate values
  EFloat ox(r2.origin().x(), oErr.x()), 
         oy(r2.origin().y(), oErr.y()), 
         oz(r2.origin().z(), oErr.z());
  EFloat dx(r2.direction().x(), dErr.x()), 
         dy(r2.direction().y(), dErr.y()), 
         dz(r2.direction().z(), dErr.z());
         
  EFloat a = dx * dx + dy * dy + dz * dz;
  EFloat b = 2 * (dx * ox + dy * oy + dz * oz);
  EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);
  
  // Solve quadratic equation for _t_ values
  EFloat temp1, temp2;
  if (!Quadratic(a, b, c, &temp1, &temp2)) {
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

bool sphere::HitP(const ray& r, Float t_min, Float t_max, Sampler* sampler) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("Sphere");
  vec3f oErr, dErr;
  ray r2 = (*WorldToObject)(r, &oErr, &dErr);
  // Compute quadratic sphere coefficients
  
  // Initialize _EFloat_ ray coordinate values
  EFloat ox(r2.origin().x(), oErr.x()), 
         oy(r2.origin().y(), oErr.y()), 
         oz(r2.origin().z(), oErr.z());
  EFloat dx(r2.direction().x(), dErr.x()), 
         dy(r2.direction().y(), dErr.y()), 
         dz(r2.direction().z(), dErr.z());
         
  EFloat a = dx * dx + dy * dy + dz * dz;
  EFloat b = 2 * (dx * ox + dy * oy + dz * oz);
  EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);
  
  // Solve quadratic equation for _t_ values
  EFloat temp1, temp2;
  if (!Quadratic(a, b, c, &temp1, &temp2)) {
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


Float sphere::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) {
  if(!this->HitP(ray(o,v), 0.001, FLT_MAX, rng)) {
    return(0);
  }
  point3f pCenter = (*ObjectToWorld)(point3f(0.f, 0.f, 0.f));
  // Return uniform PDF if point is inside sphere
  // point3f pOrigin =
  //   OffsetRayOrigin(ref.p, ref.pError, ref.n, pCenter - ref.p);
  // if (DistanceSquared(pOrigin, pCenter) <= radius * radius)
  //   return Shape::Pdf(ref, wi);
  
  // Compute general sphere PDF
  Float sinThetaMax2 = radius * radius / DistanceSquared(o, pCenter);
  Float cosThetaMax = std::sqrt(std::fmax((Float)0, 1 - sinThetaMax2));
  return UniformConePdf(cosThetaMax);
}


Float sphere::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) {
  if(!this->HitP(ray(o,v), 0.001, FLT_MAX, sampler)) {
    return(0);
  }
  point3f pCenter = (*ObjectToWorld)(point3f(0.f, 0.f, 0.f));
  // Return uniform PDF if point is inside sphere
  // point3f pOrigin =
  //   OffsetRayOrigin(ref.p, ref.pError, ref.n, pCenter - ref.p);
  // if (DistanceSquared(pOrigin, pCenter) <= radius * radius)
  //   return Shape::Pdf(ref, wi);
  
  // Compute general sphere PDF
  Float sinThetaMax2 = radius * radius / DistanceSquared(o, pCenter);
  Float cosThetaMax = std::sqrt(std::fmax((Float)0, 1 - sinThetaMax2));
  return UniformConePdf(cosThetaMax);
}

vec3f sphere::random(const point3f& o, random_gen& rng, Float time) {
  point3f pCenter = (*ObjectToWorld)(point3f(0.f, 0.f, 0.f));
  vec3f wc = pCenter - o;
  Float dc = wc.length();
  Float invDc = 1 / dc;
  wc *= dc;

  onb uvw;
  uvw.build_from_w(wc);
  vec2f u = vec2f(rng.unif_rand(),rng.unif_rand());
  
  Float sinThetaMax = radius * invDc;
  Float sinThetaMax2 = sinThetaMax * sinThetaMax;
  Float invSinThetaMax = 1 / sinThetaMax;
  
  Float cosThetaMax = std::sqrt(std::fmax((Float)0.f, 1 - sinThetaMax2));
  
  Float cosTheta  = (cosThetaMax - 1) * u[0] + 1;
  Float sinTheta2 = 1 - cosTheta * cosTheta;
  
  if (sinThetaMax2 < 0.00068523f /* sin^2(1.5 deg) */) {
    /* Fall back to a Taylor series expansion for small angles, where
     the standard approach suffers from severe cancellation errors */
    sinTheta2 = sinThetaMax2 * u[0];
    cosTheta = std::sqrt(static_cast<Float>(1) - sinTheta2);
  }
  
  // Compute angle $\alpha$ from center of sphere to sampled point on surface
  Float cosAlpha = sinTheta2 * invSinThetaMax +
    cosTheta * std::sqrt(std::fmax((Float)0.f, 1.f - sinTheta2 * invSinThetaMax * invSinThetaMax));
  Float sinAlpha = std::sqrt(std::fmax((Float)0.f, 1.f - cosAlpha*cosAlpha));
  Float phi = u.e[1] * 2 * static_cast<Float>(M_PI);
  
  // Compute surface normal and sampled point on sphere
  vec3f nWorld = SphericalDirection(sinAlpha, cosAlpha, phi, -uvw.u(), -uvw.v(), -uvw.w());
  point3f pWorld = pCenter + radius * point3f(nWorld.x(), nWorld.y(), nWorld.z());
  return (pWorld-o);
}

//Missing some error handling stuff from PBRT
vec3f sphere::random(const point3f& o, Sampler* sampler, Float time) {
  point3f pCenter = (*ObjectToWorld)(point3f(0, 0, 0));
  vec3f wc = pCenter - o;
  Float dc = wc.length();
  Float invDc = 1 / dc;
  wc *= dc;
  
  onb uvw;
  uvw.build_from_w(wc);
  vec2f u = sampler->Get2D();
  Float sinThetaMax = radius * invDc;
  Float sinThetaMax2 = sinThetaMax * sinThetaMax;
  Float invSinThetaMax = 1 / sinThetaMax;
  
  Float cosThetaMax = std::sqrt(std::fmax((Float)0.f, 1 - sinThetaMax2));
  
  Float cosTheta  = (cosThetaMax - 1) * u[0] + 1;
  Float sinTheta2 = 1 - cosTheta * cosTheta;
  
  if (sinThetaMax2 < 0.00068523f /* sin^2(1.5 deg) */) {
    /* Fall back to a Taylor series expansion for small angles, where
     the standard approach suffers from severe cancellation errors */
    sinTheta2 = sinThetaMax2 * u[0];
    cosTheta = std::sqrt(1 - sinTheta2);
  }
  
  // Compute angle $\alpha$ from center of sphere to sampled point on surface
  Float cosAlpha = sinTheta2 * invSinThetaMax +
    cosTheta * std::sqrt(std::fmax((Float)0.f, 1.f - sinTheta2 * invSinThetaMax * invSinThetaMax));
  Float sinAlpha = std::sqrt(std::fmax((Float)0.f, 1.f - cosAlpha*cosAlpha));
  Float phi = u.e[1] * 2 * static_cast<Float>(M_PI);
  
  // Compute surface normal and sampled point on sphere
  vec3f nWorld = SphericalDirection(sinAlpha, cosAlpha, phi, -uvw.u(), -uvw.v(), -uvw.w());
  point3f pWorld = pCenter + radius * point3f(nWorld.x(), nWorld.y(), nWorld.z());
  return unit_vector(pWorld-o); 
}

bool sphere::bounding_box(Float t0, Float t1, aabb& box) const {
  box = (*ObjectToWorld)(aabb(-vec3f(radius,radius,radius), vec3f(radius,radius,radius)));
  return(true);
}

size_t sphere::GetSize()  {
  return(mat_ptr ? sizeof(*this) + mat_ptr->GetSize() : sizeof(*this));
}
