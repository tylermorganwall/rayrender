#include "ellipsoid.h"


bool ellipsoid::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  ray r2 = (*WorldToObject)(r);
  ray scaled_ray(r2.origin() * point3f(inv_axes) + -center, r2.direction() * inv_axes);
  Float a = dot(scaled_ray.direction(), scaled_ray.direction());
  Float b = 2 * dot(scaled_ray.origin(), scaled_ray.direction()); 
  Float c = dot(scaled_ray.origin(),scaled_ray.origin()) - 1;
  Float temp1, temp2;
  if (!quadratic(a, b, c, &temp1, &temp2)) {
    return(false);
  }
  bool is_hit = true;
  bool second_is_hit = true;
  bool alpha_miss = false;
  if(alpha_mask) {
    Float u;
    Float v;
    if(temp1 < t_max && temp1 > t_min) {
      point3f p1 = scaled_ray.point_at_parameter(temp1) ;
      p1 *= 1/p1.length() * axes;
      vec3f normal = (p1 - center) * inv_axes;
      normal.make_unit_vector();
      get_sphere_uv(normal, u, v);
      if(alpha_mask->value(u, v, rec.p) < rng.unif_rand()) {
        is_hit = false;
      }
    }
    if(temp2 < t_max && temp2 > t_min) {
      point3f p2 = scaled_ray.point_at_parameter(temp2) ;
      p2 *= 1/p2.length() * axes;
      vec3f normal = (p2 - center) * inv_axes;
      normal.make_unit_vector();
      get_sphere_uv(normal, u, v);
      if(alpha_mask->value(u, v, rec.p) < rng.unif_rand()) {
        if(!is_hit) {
          alpha_miss = true;
        }
        second_is_hit = false;
      } 
    }
  }
  if(temp1 < t_max && temp1 > t_min && is_hit) {
    rec.t = temp1;
    rec.p = scaled_ray.point_at_parameter(rec.t) ;
    rec.normal = (rec.p - center);
    rec.mat_ptr = mat_ptr.get();
    normal3f trans_normal = rec.normal *  normal3f(inv_axes);
    trans_normal.make_unit_vector();
    
    get_sphere_uv(trans_normal, rec.u, rec.v);
    
    //Interaction information
    Float zRadius = std::sqrt(rec.p.x() * rec.p.x()  + rec.p.z()  * rec.p.z() );
    Float invZRadius = 1 / zRadius;
    Float cosPhi = rec.p.x() * invZRadius;
    Float sinPhi = rec.p.z() * invZRadius;
    Float theta = std::acos(clamp(rec.p.z(), -1, 1));
    rec.dpdu = 2 * M_PI * vec3f(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = 2 * M_PI * vec3f(rec.p.z() * cosPhi, rec.p.z() * sinPhi, -std::sin(theta));
    rec.has_bump = bump_tex ? true : false;
    
    if(bump_tex) {
      point3f bvbu = bump_tex->value(rec.u, rec.v, rec.p);
      rec.bump_normal = rec.normal + normal3f(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv); 
      rec.bump_normal *= inv_axes;
      rec.bump_normal.make_unit_vector();
    }
    rec.normal *= inv_axes;
    rec.normal.make_unit_vector();
    rec.p *= 1/rec.p.length() * axes;
    
    rec.pError = gamma(5) * Abs(rec.p);
    rec = (*ObjectToWorld)(rec);
    rec.normal *= reverseOrientation  ? -1 : 1;
    rec.bump_normal *= reverseOrientation  ? -1 : 1;
    rec.shape = this;
    rec.alpha_miss = alpha_miss;
    
    return(true);
  }
  if(temp2 < t_max && temp2 > t_min && second_is_hit) {
    rec.t = temp2;
    rec.p = scaled_ray.point_at_parameter(rec.t) ;
    rec.normal = (rec.p - center);
    rec.mat_ptr = mat_ptr.get();
    normal3f trans_normal = rec.normal *  normal3f(inv_axes);
    trans_normal.make_unit_vector();
    
    get_sphere_uv(trans_normal, rec.u, rec.v);
    //Interaction information
    Float zRadius = std::sqrt(rec.p.x() * rec.p.x()  + rec.p.z()  * rec.p.z() );
    Float invZRadius = 1 / zRadius;
    Float cosPhi = rec.p.x() * invZRadius;
    Float sinPhi = rec.p.z() * invZRadius;
    Float theta = std::acos(clamp(rec.p.z(), -1, 1));
    rec.dpdu = 2 * M_PI * vec3f(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = 2 * M_PI * vec3f(rec.p.z() * cosPhi, rec.p.z() * sinPhi, -std::sin(theta));
    rec.has_bump = bump_tex ? true : false;
    
    if(bump_tex) {
      point3f bvbu = bump_tex->value(rec.u, rec.v, rec.p);
      rec.bump_normal = rec.normal + normal3f(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv); 
      rec.bump_normal *= inv_axes;
      rec.bump_normal.make_unit_vector();
    }
    rec.normal *= inv_axes;
    rec.normal.make_unit_vector();
    rec.p *= 1/rec.p.length() * axes;
    
    if(alpha_mask) {
      rec.normal = -rec.normal;
      rec.bump_normal = -rec.bump_normal;
    }
    rec.pError = gamma(5) * Abs(rec.p);
    rec = (*ObjectToWorld)(rec);
    rec.normal *= reverseOrientation  ? -1 : 1;
    rec.bump_normal *= reverseOrientation  ? -1 : 1;
    rec.shape = this;
    rec.alpha_miss = alpha_miss;
    
    return(true);
  }
  return(false);
}


bool ellipsoid::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
  ray r2 = (*WorldToObject)(r);
  
  ray scaled_ray(r2.origin() * point3f(inv_axes) + -center, r2.direction() * inv_axes);
  Float a = dot(scaled_ray.direction(), scaled_ray.direction());
  Float b = 2 * dot(scaled_ray.origin(), scaled_ray.direction()); 
  Float c = dot(scaled_ray.origin(),scaled_ray.origin()) - 1;
  Float temp1, temp2;
  if (!quadratic(a, b, c, &temp1, &temp2)) {
    return(false);
  }
  bool is_hit = true;
  bool second_is_hit = true;
  bool alpha_miss = false;
  
  if(alpha_mask) {
    Float u;
    Float v;
    if(temp1 < t_max && temp1 > t_min) {
      point3f p1 = scaled_ray.point_at_parameter(temp1) ;
      p1 *= 1/p1.length() * axes;
      vec3f normal = (p1 - center) * inv_axes;
      normal.make_unit_vector();
      get_sphere_uv(normal, u, v);
      if(alpha_mask->value(u, v, rec.p) < sampler->Get1D()) {
        is_hit = false;
      }
    }
    if(temp2 < t_max && temp2 > t_min) {
      point3f p2 = scaled_ray.point_at_parameter(temp2) ;
      p2 *= 1/p2.length() * axes;
      vec3f normal = (p2 - center) * inv_axes;
      normal.make_unit_vector();
      get_sphere_uv(normal, u, v);
      if(alpha_mask->value(u, v, rec.p) < sampler->Get1D()) {
        if(!is_hit) {
          alpha_miss = true;
          return(false);
        }
        second_is_hit = false;
      } 
    }
  }
  if(temp1 < t_max && temp1 > t_min && is_hit) {
    rec.t = temp1;
    rec.p = scaled_ray.point_at_parameter(rec.t) ;
    rec.normal = (rec.p - center);
    rec.mat_ptr = mat_ptr.get();
    normal3f trans_normal = rec.normal *  normal3f(inv_axes);
    trans_normal.make_unit_vector();
    
    get_sphere_uv(trans_normal, rec.u, rec.v);
    
    //Interaction information
    Float zRadius = std::sqrt(rec.p.x() * rec.p.x()  + rec.p.z()  * rec.p.z() );
    Float invZRadius = 1 / zRadius;
    Float cosPhi = rec.p.x() * invZRadius;
    Float sinPhi = rec.p.z() * invZRadius;
    Float theta = std::acos(clamp(rec.p.z(), -1, 1));
    rec.dpdu = 2 * M_PI * vec3f(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = 2 * M_PI * vec3f(rec.p.z() * cosPhi, rec.p.z() * sinPhi, -std::sin(theta));
    rec.has_bump = bump_tex ? true : false;
    
    if(bump_tex) {
      point3f bvbu = bump_tex->value(rec.u, rec.v, rec.p);
      rec.bump_normal = rec.normal + normal3f(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv); 
      rec.bump_normal *= inv_axes;
      rec.bump_normal.make_unit_vector();
    }
    rec.normal *= inv_axes;
    rec.normal.make_unit_vector();
    rec.p *= 1/rec.p.length() * axes;
    
    rec.pError = gamma(5) * Abs(rec.p);
    rec = (*ObjectToWorld)(rec);
    rec.normal *= reverseOrientation  ? -1 : 1;
    rec.bump_normal *= reverseOrientation  ? -1 : 1;
    rec.shape = this;
    rec.alpha_miss = alpha_miss;
    
    return(true);
  }
  if(temp2 < t_max && temp2 > t_min && second_is_hit) {
    rec.t = temp2;
    rec.p = scaled_ray.point_at_parameter(rec.t) ;
    rec.normal = (rec.p - center);
    rec.mat_ptr = mat_ptr.get();
    normal3f trans_normal = rec.normal *  normal3f(inv_axes);
    trans_normal.make_unit_vector();
    
    get_sphere_uv(trans_normal, rec.u, rec.v);
    //Interaction information
    Float zRadius = std::sqrt(rec.p.x() * rec.p.x()  + rec.p.z()  * rec.p.z() );
    Float invZRadius = 1 / zRadius;
    Float cosPhi = rec.p.x() * invZRadius;
    Float sinPhi = rec.p.z() * invZRadius;
    Float theta = std::acos(clamp(rec.p.z(), -1, 1));
    rec.dpdu = 2 * M_PI * vec3f(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = 2 * M_PI * vec3f(rec.p.z() * cosPhi, rec.p.z() * sinPhi, -std::sin(theta));
    rec.has_bump = bump_tex ? true : false;
    
    if(bump_tex) {
      point3f bvbu = bump_tex->value(rec.u, rec.v, rec.p);
      rec.bump_normal = rec.normal + normal3f(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv); 
      rec.bump_normal *= inv_axes;
      rec.bump_normal.make_unit_vector();
    }
    rec.normal *= inv_axes;
    rec.normal.make_unit_vector();
    rec.p *= 1/rec.p.length() * axes;
    
    if(alpha_mask) {
      rec.normal = -rec.normal;
      rec.bump_normal = -rec.bump_normal;
    }
    rec.pError = gamma(5) * Abs(rec.p);
    rec = (*ObjectToWorld)(rec);
    rec.normal *= reverseOrientation  ? -1 : 1;
    rec.bump_normal *= reverseOrientation  ? -1 : 1;
    rec.shape = this;
    rec.alpha_miss = alpha_miss;
    
    
    return(true);
  }
  return(false);
}

//Not great
Float ellipsoid::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, rng)) {
    point3f o2 = (*WorldToObject)(o);
    Float cos_theta_max = sqrt(1 - 1/(center - o2).squared_length());
    Float solid_angle = 2 * M_PI * (1-cos_theta_max) * largest_proj_axis ;
    return(1/solid_angle);
  } else {
    return(0);
  }
}


//Not great
Float ellipsoid::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, sampler)) {
    point3f o2 = (*WorldToObject)(o);
    Float cos_theta_max = sqrt(1 - 1/(center - o2).squared_length());
    Float solid_angle = 2 * M_PI * (1-cos_theta_max) * largest_proj_axis ;
    return(1/solid_angle);
  } else {
    return(0);
  }
}

vec3f ellipsoid::random(const point3f& o, random_gen& rng, Float time) {
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
  Float phi = u.e[1] * 2 * M_PI;
  
  // Compute surface normal and sampled point on sphere
  vec3f nWorld = SphericalDirection(sinAlpha, cosAlpha, phi, -uvw.u(), -uvw.v(), -uvw.w()) * inv_axes;
  point3f pWorld = pCenter + radius * point3f(nWorld.x(), nWorld.y(), nWorld.z());
  return (pWorld-o);
}

vec3f ellipsoid::random(const point3f& o, Sampler* sampler, Float time) {
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
  Float phi = u.e[1] * 2 * M_PI;
  
  // Compute surface normal and sampled point on sphere
  vec3f nWorld = SphericalDirection(sinAlpha, cosAlpha, phi, -uvw.u(), -uvw.v(), -uvw.w()) * inv_axes;
  point3f pWorld = pCenter + radius * point3f(nWorld.x(), nWorld.y(), nWorld.z());
  return (pWorld-o);
}

bool ellipsoid::bounding_box(Float t0, Float t1, aabb& box) const {
  box = (*ObjectToWorld)(aabb(center + -point3f(radius,radius,radius) * axes, center + point3f(radius,radius,radius) * axes));
  return(true);
}
