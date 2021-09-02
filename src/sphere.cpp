#include "sphere.h"

#include "RcppThread.h"

bool sphere::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  ray r2 = (*WorldToObject)(r);
  vec3f oc = r2.origin() - center;
  Float a = dot(r2.direction(), r2.direction());
  Float b = 2 * dot(oc, r2.direction()); 
  Float c = dot(oc,oc) - radius * radius;
  Float temp1, temp2;
  if (!quadratic(a, b, c, &temp1, &temp2)) {
    return(false);
  }
  bool is_hit = true;
  bool second_is_hit = true;
  if(alpha_mask) {
    Float u;
    Float v;
    if(temp1 < t_max && temp1 > t_min) {
      point3f p1 = r2.point_at_parameter(temp1);
      p1 *= radius / p1.length(); 
      vec3f normal = (p1 - center) / radius;
      get_sphere_uv(normal, u, v);
      if(alpha_mask->value(u, v, rec.p).x() < rng.unif_rand()) {
        is_hit = false;
      }
    }
    if(temp2 < t_max && temp2 > t_min) {
      point3f p2 = r2.point_at_parameter(temp2);
      p2 *= radius / p2.length(); 
      vec3f normal = (p2 - center) / radius;
      get_sphere_uv(normal, u, v);
      if(alpha_mask->value(u, v, rec.p).x() < rng.unif_rand()) {
        if(!is_hit) {
          return(false);
        }
        second_is_hit = false;
      } 
    }
  }
  if(temp1 < t_max && temp1 > t_min && is_hit) {
    rec.t = temp1;
    rec.p = r2.point_at_parameter(rec.t);
    rec.p *= radius / rec.p.length(); 
    rec.normal = (rec.p - center) / radius;
    
    //Interaction information
    Float zRadius = std::sqrt(rec.p.x() * rec.p.x()  + rec.p.z()  * rec.p.z() );
    Float invZRadius = 1 / zRadius;
    Float cosPhi = rec.p.x() * invZRadius;
    Float sinPhi = rec.p.z() * invZRadius;
    Float theta = std::acos(clamp(rec.p.z() / radius, -1, 1));
    rec.dpdu = 2 * M_PI * vec3f(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = 2 * M_PI * vec3f(rec.p.z() * cosPhi, rec.p.z() * sinPhi, -radius * std::sin(theta));
    get_sphere_uv(rec.normal, rec.u, rec.v);
    rec.has_bump = bump_tex ? true : false;
    
    if(bump_tex) {
      point3f bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      rec.bump_normal = rec.normal + normal3f(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv); 
      rec.bump_normal.make_unit_vector();
      rec.bump_normal = (*ObjectToWorld)(rec.bump_normal);
    }
    rec.p = (*ObjectToWorld)(rec.p);
    rec.normal = !reverseOrientation ? (*ObjectToWorld)(rec.normal) : -(*ObjectToWorld)(rec.normal);
    
    rec.mat_ptr = mat_ptr.get();
    return(true);
  }
  if(temp2 < t_max && temp2 > t_min && second_is_hit) {
    rec.t = temp2;
    rec.p = r2.point_at_parameter(rec.t);
    rec.p *= radius / rec.p.length();
    rec.normal = (rec.p - center) / radius;
    
    //Interaction information
    Float zRadius = std::sqrt(rec.p.x() * rec.p.x()  + rec.p.z()  * rec.p.z() );
    Float invZRadius = 1 / zRadius;
    Float cosPhi = rec.p.x() * invZRadius;
    Float sinPhi = rec.p.z() * invZRadius;
    Float theta = std::acos(clamp(rec.p.z() / radius, -1, 1));
    rec.dpdu = 2 * M_PI * vec3f(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = 2 * M_PI * vec3f(rec.p.z() * cosPhi, rec.p.z() * sinPhi, -radius * std::sin(theta));
    get_sphere_uv(rec.normal, rec.u, rec.v);
    rec.has_bump = bump_tex ? true : false;
    
    if(bump_tex) {
      point3f bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      rec.bump_normal = rec.normal +  normal3f(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv); 
      rec.bump_normal.make_unit_vector();
      rec.bump_normal = (*ObjectToWorld)(rec.bump_normal);
    }
    
    if(alpha_mask) {
      rec.normal = -rec.normal;
      rec.bump_normal = -rec.bump_normal;
    }
    rec.p = (*ObjectToWorld)(rec.p);
    rec.normal = !reverseOrientation ? (*ObjectToWorld)(rec.normal) : -(*ObjectToWorld)(rec.normal);
    rec.mat_ptr = mat_ptr.get();
    return(true);
  }
  return(false);
}


bool sphere::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
  ray r2 = (*WorldToObject)(r);
  vec3f oc = r2.origin() - center;
  Float a = dot(r2.direction(), r2.direction());
  Float b = 2 * dot(oc, r2.direction()); 
  Float c = dot(oc,oc) - radius * radius;
  Float temp1, temp2;
  if (!quadratic(a, b, c, &temp1, &temp2)) {
    return(false);
  }
  bool is_hit = true;
  bool second_is_hit = true;
  if(alpha_mask) {
    Float u;
    Float v;
    if(temp1 < t_max && temp1 > t_min) {
      point3f p1 = r2.point_at_parameter(temp1);
      p1 *= radius / p1.length(); 
      vec3f normal = (p1 - center) / radius;
      get_sphere_uv(normal, u, v);
      if(alpha_mask->value(u, v, rec.p).x() < sampler->Get1D()) {
        is_hit = false;
      }
    }
    if(temp2 < t_max && temp2 > t_min) {
      point3f p2 = r2.point_at_parameter(temp2);
      p2 *= radius / p2.length(); 
      vec3f normal = (p2 - center) / radius;
      get_sphere_uv(normal, u, v);
      if(alpha_mask->value(u, v, rec.p).x() < sampler->Get1D()) {
        if(!is_hit) {
          return(false);
        }
        second_is_hit = false;
      } 
    }
  }
  if(temp1 < t_max && temp1 > t_min && is_hit) {
    rec.t = temp1;
    rec.p = r2.point_at_parameter(rec.t);
    rec.p *= radius / rec.p.length(); 
    rec.normal = (rec.p - center) / radius;
    
    //Interaction information
    Float zRadius = std::sqrt(rec.p.x() * rec.p.x()  + rec.p.z()  * rec.p.z() );
    Float invZRadius = 1 / zRadius;
    Float cosPhi = rec.p.x() * invZRadius;
    Float sinPhi = rec.p.z() * invZRadius;
    Float theta = std::acos(clamp(rec.p.z() / radius, -1, 1));
    rec.dpdu = 2 * M_PI * vec3f(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = 2 * M_PI * vec3f(rec.p.z() * cosPhi, rec.p.z() * sinPhi, -radius * std::sin(theta));
    get_sphere_uv(rec.normal, rec.u, rec.v);
    rec.has_bump = bump_tex ? true : false;
    
    if(bump_tex) {
      point3f bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      rec.bump_normal = rec.normal + normal3f(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv); 
      rec.bump_normal.make_unit_vector();
      rec.bump_normal = (*ObjectToWorld)(rec.bump_normal);
    }
    rec.p = (*ObjectToWorld)(rec.p);
    rec.normal = !reverseOrientation ? (*ObjectToWorld)(rec.normal) : -(*ObjectToWorld)(rec.normal);
    
    rec.mat_ptr = mat_ptr.get();
    return(true);
  }
  if(temp2 < t_max && temp2 > t_min && second_is_hit) {
    rec.t = temp2;
    rec.p = r2.point_at_parameter(rec.t);
    rec.p *= radius / rec.p.length();
    rec.normal = (rec.p - center) / radius;
    
    //Interaction information
    Float zRadius = std::sqrt(rec.p.x() * rec.p.x()  + rec.p.z()  * rec.p.z() );
    Float invZRadius = 1 / zRadius;
    Float cosPhi = rec.p.x() * invZRadius;
    Float sinPhi = rec.p.z() * invZRadius;
    Float theta = std::acos(clamp(rec.p.z() / radius, -1, 1));
    rec.dpdu = 2 * M_PI * vec3f(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = 2 * M_PI * vec3f(rec.p.z() * cosPhi, rec.p.z() * sinPhi, -radius * std::sin(theta));
    get_sphere_uv(rec.normal, rec.u, rec.v);
    rec.has_bump = bump_tex ? true : false;
    
    if(bump_tex) {
      point3f bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      rec.bump_normal = rec.normal +  normal3f(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv); 
      rec.bump_normal.make_unit_vector();
      rec.bump_normal = (*ObjectToWorld)(rec.bump_normal);
    }
    
    if(alpha_mask) {
      rec.normal = -rec.normal;
      rec.bump_normal = -rec.bump_normal;
    }
    rec.p = (*ObjectToWorld)(rec.p);
    rec.normal = !reverseOrientation ? (*ObjectToWorld)(rec.normal) : -(*ObjectToWorld)(rec.normal);
    rec.mat_ptr = mat_ptr.get();
    return(true);
  }
  return(false);
}


Float sphere::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, rng)) {
    point3f o2 = (*WorldToObject)(o);
    Float cos_theta_max = sqrt(1 - radius * radius/(center - o2).squared_length());
    Float solid_angle = 2 * M_PI * (1-cos_theta_max);
    return(1/solid_angle);
  } else {
    return(0);
  }
}


Float sphere::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, sampler)) {
    point3f o2 = (*WorldToObject)(o);
    Float cos_theta_max = sqrt(1 - radius * radius/(center - o2).squared_length());
    Float solid_angle = 2 * M_PI * (1-cos_theta_max);
    return(1/solid_angle);
  } else {
    return(0);
  }
}

vec3f sphere::random(const point3f& o, random_gen& rng, Float time) {
  point3f pCenter = (*ObjectToWorld)(point3f(0, 0, 0));
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
  
  Float cosThetaMax = std::sqrt(std::max((Float)0.f, 1 - sinThetaMax2));
  
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
  Float phi = u.e[1] * 2 * M_PI;
  
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
  
  Float cosThetaMax = std::sqrt(std::max((Float)0.f, 1 - sinThetaMax2));
  
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
  Float phi = u.e[1] * 2 * M_PI;
  
  // Compute surface normal and sampled point on sphere
  vec3f nWorld = SphericalDirection(sinAlpha, cosAlpha, phi, -uvw.u(), -uvw.v(), -uvw.w());
  point3f pWorld = pCenter + radius * point3f(nWorld.x(), nWorld.y(), nWorld.z());
  return (pWorld-o);
}

bool sphere::bounding_box(Float t0, Float t1, aabb& box) const {
  box = (*ObjectToWorld)(aabb(center - vec3f(radius,radius,radius), center + vec3f(radius,radius,radius)));
  return(true);
}

//
// Point3f pOrigin = OffsetRayOrigin(ref.p, ref.pError, ref.n,
//                                   pCenter - ref.p);
// Float sinThetaMax2 = radius * radius / (ref.p - pCenter).squared)dista;
// Float cosThetaMax = std::sqrt(std::max((Float)0, 1 - sinThetaMax2));
// Float cosTheta = (1 - u[0]) + u[0] * cosThetaMax;
// Float sinTheta = std::sqrt(std::max((Float)0, 1 - cosTheta * cosTheta));
// Float phi = u[1] * 2 * Pi;
// 
// <<Compute angle  from center of sphere to sampled point on surface>> 
// Float dc = Distance(ref.p, pCenter);
// Float ds = dc * cosTheta -
//   std::sqrt(std::max((Float)0,
//                      radius * radius - dc * dc * sinTheta * sinTheta));
// Float cosAlpha = (dc * dc + radius * radius - ds * ds) /
//   (2 * dc * radius);
// Float sinAlpha = std::sqrt(std::max((Float)0, 1 - cosAlpha * cosAlpha));
// 
// <<Compute surface normal and sampled point on sphere>> 
// Vector3f nObj = SphericalDirection(sinAlpha, cosAlpha, phi,
//                                    -wcX, -wcY, -wc);
// Point3f pObj = radius * Point3f(nObj.x, nObj.y, nObj.z);
// 
// <<Return Interaction for sampled point on sphere>> 
// Interaction it;
// <<Reproject pObj to sphere surface and compute pObjError>> 
// it.p = (*ObjectToWorld)(pObj, pObjError, &it.pError);
// it.n = (*ObjectToWorld)(Normal3f(nObj));
// if (reverseOrientation) it.n *= -1;
// return it;


point3f moving_sphere::center(Float time) const {
  return(center0 + ((time - time0) / (time1 - time0)) * (center1 - center0));
}

bool moving_sphere::bounding_box(Float t0, Float t1, aabb& box) const {
  aabb box0(center(t0) - vec3f(radius,radius,radius), center(t0) + vec3f(radius,radius,radius));
  aabb box1(center(t1) - vec3f(radius,radius,radius), center(t1) + vec3f(radius,radius,radius));
  box = (*ObjectToWorld)(surrounding_box(box0, box1));
  return(true);
}

bool moving_sphere::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  ray r2 = (*WorldToObject)(r);
  vec3f oc = r2.origin() - center(r2.time());
  Float a = dot(r2.direction(), r2.direction());
  Float b = 2 * dot(oc, r2.direction()); 
  Float c = dot(oc,oc) - radius * radius;
  Float temp1, temp2;
  if (!quadratic(a, b, c, &temp1, &temp2)) {
    return(false);
  }
  bool is_hit = true;
  if(alpha_mask) {
    Float u;
    Float v;
    if(temp1 < t_max && temp1 > t_min) {
      point3f p1 = r2.point_at_parameter(temp1);
      p1 *= radius / p1.length(); 
      vec3f normal = (p1 - center(r2.time())) / radius;
      get_sphere_uv(normal, u, v);
      if(alpha_mask->value(u, v, rec.p).x() < rng.unif_rand()) {
        is_hit = false;
      }
    }
    if(temp2 < t_max && temp2 > t_min) {
      point3f p2 = r2.point_at_parameter(temp2);
      p2 *= radius / p2.length(); 
      vec3f normal = (p2 - center(r2.time())) / radius;
      get_sphere_uv(normal, u, v);
      if(alpha_mask->value(u, v, rec.p).x() < rng.unif_rand()) {
        if(!is_hit) {
          return(false);
        }
      } 
    }
  }
  if(temp1 < t_max && temp1 > t_min && is_hit) {
    rec.t = temp1;
    rec.p = r2.point_at_parameter(rec.t);
    rec.p *= radius / (rec.p - center(r2.time())).length();
    rec.normal = (rec.p - center(r2.time())) / radius;
    
    //Interaction information
    Float zRadius = std::sqrt(rec.p.x() * rec.p.x()  + rec.p.z()  * rec.p.z() );
    Float invZRadius = 1 / zRadius;
    Float cosPhi = rec.p.x() * invZRadius;
    Float sinPhi = rec.p.z() * invZRadius;
    Float theta = std::acos(clamp(rec.p.z() / radius, -1, 1));
    rec.dpdu = 2 * M_PI * vec3f(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = 2 * M_PI * vec3f(rec.p.z() * cosPhi, rec.p.z() * sinPhi, -radius * std::sin(theta));
    
    if(bump_tex) {
      point3f bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      rec.bump_normal = rec.normal + normal3f(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv); 
      rec.bump_normal.make_unit_vector();
      rec.bump_normal = (*ObjectToWorld)(rec.bump_normal);
      
    }
    
    get_sphere_uv(rec.normal, rec.u, rec.v);
    rec.p = (*ObjectToWorld)(rec.p);
    rec.normal = !reverseOrientation ? (*ObjectToWorld)(rec.normal) : -(*ObjectToWorld)(rec.normal);
    rec.mat_ptr = mat_ptr.get();
    return(true);
  }
  if(temp2 < t_max && temp2 > t_min) {
    rec.t = temp2;
    rec.p = r2.point_at_parameter(rec.t);
    rec.p *= radius / (rec.p - center(r2.time())).length(); 
    rec.normal = (rec.p - center(r2.time())) / radius;
    
    //Interaction information
    Float zRadius = std::sqrt(rec.p.x() * rec.p.x()  + rec.p.z()  * rec.p.z() );
    Float invZRadius = 1 / zRadius;
    Float cosPhi = rec.p.x() * invZRadius;
    Float sinPhi = rec.p.z() * invZRadius;
    Float theta = std::acos(clamp(rec.p.z() / radius, -1, 1));
    rec.dpdu = 2 * M_PI * vec3f(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = 2 * M_PI * vec3f(rec.p.z() * cosPhi, rec.p.z() * sinPhi, -radius * std::sin(theta));
    
    if(bump_tex) {
      point3f bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      vec3f o_u = cross(vec3f(rec.normal.x(),rec.normal.y(),rec.normal.z()), rec.dpdu);
      vec3f o_v = cross(vec3f(rec.normal.x(),rec.normal.y(),rec.normal.z()), rec.dpdv);
      rec.bump_normal = rec.normal + normal3f(bvbu.x() * o_v - bvbu.y() * o_u); 
      rec.bump_normal.make_unit_vector();
      rec.bump_normal = (*ObjectToWorld)(rec.bump_normal);
      
    }
    
    get_sphere_uv(rec.normal, rec.u, rec.v);
    if(!is_hit) {
      rec.normal = -rec.normal;
      rec.bump_normal = -rec.bump_normal;
    }
    rec.p = (*ObjectToWorld)(rec.p);
    rec.normal = !reverseOrientation ? (*ObjectToWorld)(rec.normal) : -(*ObjectToWorld)(rec.normal);
    rec.mat_ptr = mat_ptr.get();
    return(true);
  }
  return(false);
}


vec3f moving_sphere::random(const point3f& o, random_gen& rng, Float time) {
  vec3f direction = center(time) - o;
  Float distance_squared = direction.squared_length();
  onb uvw;
  uvw.build_from_w(direction);
  return(uvw.local_to_world(rng.random_to_sphere(radius,distance_squared)));
}

vec3f moving_sphere::random(const point3f& o, Sampler* sampler, Float time) {
  vec3f direction = center(time) - o;
  Float distance_squared = direction.squared_length();
  onb uvw;
  uvw.build_from_w(direction);
  return(uvw.local_to_world(rand_to_sphere(radius,distance_squared, sampler->Get2D())));
}


bool moving_sphere::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
  ray r2 = (*WorldToObject)(r);
  vec3f oc = r2.origin() - center(r2.time());
  Float a = dot(r2.direction(), r2.direction());
  Float b = 2 * dot(oc, r2.direction()); 
  Float c = dot(oc,oc) - radius * radius;
  Float temp1, temp2;
  if (!quadratic(a, b, c, &temp1, &temp2)) {
    return(false);
  }
  bool is_hit = true;
  bool second_is_hit = true;
  if(alpha_mask) {
    Float u;
    Float v;
    if(temp1 < t_max && temp1 > t_min) {
      point3f p1 = r2.point_at_parameter(temp1);
      p1 *= radius / p1.length(); 
      vec3f normal = (p1 - center(r2.time())) / radius;
      get_sphere_uv(normal, u, v);
      if(alpha_mask->value(u, v, rec.p).x() < sampler->Get1D()) {
        is_hit = false;
      }
    }
    if(temp2 < t_max && temp2 > t_min) {
      point3f p2 = r2.point_at_parameter(temp2);
      p2 *= radius / p2.length(); 
      vec3f normal = (p2 - center(r2.time())) / radius;
      get_sphere_uv(normal, u, v);
      if(alpha_mask->value(u, v, rec.p).x() < sampler->Get1D()) {
        if(!is_hit) {
          return(false);
        }
        second_is_hit = false;
      } 
    }
  }
  if(temp1 < t_max && temp1 > t_min && is_hit) {
    rec.t = temp1;
    rec.p = r2.point_at_parameter(rec.t);
    rec.p *= radius / (rec.p- center(r2.time())).length();
    rec.normal = (rec.p - center(r2.time())) / radius;
    
    //Interaction information
    Float zRadius = std::sqrt(rec.p.x() * rec.p.x()  + rec.p.z()  * rec.p.z() );
    Float invZRadius = 1 / zRadius;
    Float cosPhi = rec.p.x() * invZRadius;
    Float sinPhi = rec.p.z() * invZRadius;
    Float theta = std::acos(clamp(rec.p.z() / radius, -1, 1));
    rec.dpdu = 2 * M_PI * vec3f(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = 2 * M_PI * vec3f(rec.p.z() * cosPhi, rec.p.z() * sinPhi, -radius * std::sin(theta));
    get_sphere_uv(rec.normal, rec.u, rec.v);
    rec.has_bump = bump_tex ? true : false;
    
    if(bump_tex) {
      point3f bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      rec.bump_normal = rec.normal + normal3f(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv); 
      rec.bump_normal.make_unit_vector();
      rec.bump_normal = (*ObjectToWorld)(rec.bump_normal);
      
    }
    rec.p = (*ObjectToWorld)(rec.p);
    rec.normal = !reverseOrientation ? (*ObjectToWorld)(rec.normal) : -(*ObjectToWorld)(rec.normal);
    rec.mat_ptr = mat_ptr.get();
    return(true);
  }
  if(temp2 < t_max && temp2 > t_min && second_is_hit) {
    rec.t = temp2;
    rec.p = r2.point_at_parameter(rec.t);
    rec.p *= radius / (rec.p - center(r2.time())).length();
    rec.normal = (rec.p - center(r2.time())) / radius;
    
    //Interaction information
    Float zRadius = std::sqrt(rec.p.x() * rec.p.x()  + rec.p.z()  * rec.p.z() );
    Float invZRadius = 1 / zRadius;
    Float cosPhi = rec.p.x() * invZRadius;
    Float sinPhi = rec.p.z() * invZRadius;
    Float theta = std::acos(clamp(rec.p.z() / radius, -1, 1));
    rec.dpdu = 2 * M_PI * vec3f(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = 2 * M_PI * vec3f(rec.p.z() * cosPhi, rec.p.z() * sinPhi, -radius * std::sin(theta));
    get_sphere_uv(rec.normal, rec.u, rec.v);
    rec.has_bump = bump_tex ? true : false;
    
    if(bump_tex) {
      point3f bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      rec.bump_normal = rec.normal +  normal3f(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv); 
      rec.bump_normal.make_unit_vector();
      rec.bump_normal = (*ObjectToWorld)(rec.bump_normal);
      
    }
    
    if(alpha_mask) {
      rec.normal = -rec.normal;
      rec.bump_normal = -rec.bump_normal;
    }
    rec.p = (*ObjectToWorld)(rec.p);
    rec.normal = !reverseOrientation ? (*ObjectToWorld)(rec.normal) : -(*ObjectToWorld)(rec.normal);
    rec.mat_ptr = mat_ptr.get();
    return(true);
  }
  return(false);
}

Float moving_sphere::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) {
  hit_record rec;
  if(this->hit(ray(o,v,time), 0.001, FLT_MAX, rec, rng)) {
    point3f o2 = (*WorldToObject)(o);
    
    Float cos_theta_max = sqrt(1 - radius * radius/(center(time) - o2).squared_length());
    Float solid_angle = 2 * M_PI * (1-cos_theta_max);
    return(1/solid_angle);
  } else {
    return(0);
  }
}


Float moving_sphere::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) {
  hit_record rec;
  if(this->hit(ray(o,v,time), 0.001, FLT_MAX, rec, sampler)) {
    point3f o2 = (*WorldToObject)(o);
    
    Float cos_theta_max = sqrt(1 - radius * radius/(center(time) - o2).squared_length());
    Float solid_angle = 2 * M_PI * (1-cos_theta_max);
    return(1/solid_angle);
  } else {
    return(0);
  }
}

