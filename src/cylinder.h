#ifndef CYLINDERH
#define CYLINDERH

#include "hitable.h"
#include "material.h"

class cylinder: public hitable {
public:
  cylinder() {}
  cylinder(Float r, Float len, Float phi_min, Float phi_max, bool has_caps,
           material *mat, alpha_texture *alpha_mask, bump_texture* bump_tex) : 
  radius(r), length(len), phi_min(phi_min), phi_max(phi_max), has_caps(has_caps), mat_ptr(mat), 
  alpha_mask(alpha_mask), bump_tex(bump_tex) {};
  ~cylinder() {
    delete mat_ptr;
    delete alpha_mask;
    delete bump_tex;
  }
  virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng);
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng);
  virtual vec3 random(const vec3& o, random_gen& rng);
  virtual vec3 random(const vec3& o, Sampler* sampler);
  void get_cylinder_uv(const vec3& p, Float& u, Float& v) {
    Float phi = atan2(p.z(),p.x());
    // if (phi < 0) phi += 2 * M_PI;
    u = 1 - (phi + M_PI) / (2*M_PI);
    v = (p.y() + length/2)/length;
  };
  Float radius;
  Float length;
  Float phi_min;
  Float phi_max;
  bool has_caps;
  material *mat_ptr;
  alpha_texture *alpha_mask;
  bump_texture *bump_tex;
};

bool cylinder::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  vec3 oc = r.origin();
  vec3 dir = r.direction();
  dir.e[1] = 0;
  oc.e[1] = 0;
  Float a = dot(dir, dir);
  Float b = 2 * dot(oc, dir); 
  Float c = dot(oc, oc) - radius * radius;
  Float temp1, temp2;
  if (!quadratic(a, b, c, &temp1, &temp2)) {
    return(false);
  }
  bool is_hit = true;
  bool second_is_hit = true;
  if(alpha_mask) {
    vec3 temppoint = r.point_at_parameter(temp1);
    Float phi = atan2(temppoint.z(),temppoint.x());
    phi = phi < 0 ? phi + 2 * M_PI : phi;
    Float u;
    Float v;
    if(temp1 < t_max && temp1 > t_min && 
       temppoint.y() > -length/2 && temppoint.y() < length/2 && phi <= phi_max && phi >= phi_min) {
      Float hitRad = std::sqrt(temppoint.x() * temppoint.x() + temppoint.z() * temppoint.z());
      temppoint.e[0] *= radius / hitRad;
      temppoint.e[2] *= radius / hitRad;
      get_cylinder_uv(temppoint, u, v);
      if(alpha_mask->value(u, v, rec.p).x() < rng.unif_rand()) {
        is_hit = false;
      }
    }
    temppoint = r.point_at_parameter(temp2);
    phi = atan2(temppoint.z(),temppoint.x());
    phi = phi < 0 ? phi + 2 * M_PI : phi;
    if(temp2 < t_max && temp2 > t_min && 
       temppoint.y() > -length/2 && temppoint.y() < length/2 && phi <= phi_max && phi >= phi_min) {
      Float hitRad = std::sqrt(temppoint.x() * temppoint.x() + temppoint.z() * temppoint.z());
      temppoint.e[0] *= radius / hitRad;
      temppoint.e[2] *= radius / hitRad;
      get_cylinder_uv(temppoint, u, v);
      if(alpha_mask->value(u, v, rec.p).x() < rng.unif_rand()) {
        if(!is_hit) {
          return(false);
        }
        second_is_hit = false;
      } 
    }
  }
  vec3 temppoint = r.point_at_parameter(temp1);
  Float phi = atan2(temppoint.z(),temppoint.x());
  phi = phi < 0 ? phi + 2 * M_PI : phi;
  if(is_hit && temp1 < t_max && temp1 > t_min && 
     temppoint.y() > -length/2 && temppoint.y() < length/2 && phi <= phi_max && phi >= phi_min) {
    Float hitRad = std::sqrt(temppoint.x() * temppoint.x() + temppoint.z() * temppoint.z());
    temppoint.e[0] *= radius / hitRad;
    temppoint.e[2] *= radius / hitRad;
    rec.t = temp1;
    rec.p = temppoint;
    
    temppoint.e[1] = 0;
    rec.normal = dot(temppoint, dir) > 0 ? -temppoint / radius : temppoint / radius;
    get_cylinder_uv(rec.p, rec.u, rec.v);
    
    //Interaction information
    rec.dpdu = vec3(-temppoint.z(),0,  temppoint.x());
    rec.dpdv = vec3(0, length, 0);
    rec.has_bump = bump_tex ? true : false;
    
    if(bump_tex) {
      vec3 bvbu = bump_tex->value(rec.u, rec.v, rec.p);
      rec.bump_normal = rec.normal + bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv; 
      rec.bump_normal.make_unit_vector();
      rec.bump_normal *= dot(temppoint, dir) > 0 ? -1 : 1;
      
    }
    rec.mat_ptr = mat_ptr;
    return(true);
  }
  Float t_cyl = -(r.origin().y()-length/2) / r.direction().y();
  Float t_cyl2 = -(r.origin().y()+length/2) / r.direction().y();
  Float x = r.origin().x() + t_cyl*r.direction().x();
  Float z = r.origin().z() + t_cyl*r.direction().z();
  
  Float phi2 = atan2(z,x);
  phi2 = phi2 < 0 ? phi2 + 2 * M_PI : phi2;
  Float radHit2 = x*x + z*z;
  if(has_caps && t_cyl < temp2 && t_cyl > t_min && t_cyl < t_max && t_cyl < t_cyl2 && 
     radHit2 <= radius * radius && phi2 <= phi_max && phi2 >= phi_min) {
    vec3 p = r.point_at_parameter(t_cyl);
    p.e[1] = length/2;
    
    Float u = p.x() / (2.0 * radius) + 0.5;
    Float v = p.z() / (2.0 * radius) + 0.5;
    u = 1 - u;
    if(alpha_mask) {
      if(alpha_mask->value(u, v, rec.p).x() < 1) {
        return(false);
      }
    }
    rec.p = p;
    rec.normal = vec3(0,1,0);
    rec.t = t_cyl;
    rec.mat_ptr = mat_ptr;
    rec.u = u;
    rec.v = v;
    rec.dpdu = vec3(1, 0, 0);
    rec.dpdv = vec3(0, 0, 1);
    rec.has_bump = bump_tex ? true : false;
    
    if(bump_tex) {
      vec3 bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      rec.bump_normal = rec.normal + bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv;
      rec.bump_normal.make_unit_vector();
    }
    return(true);
  }
  Float x2 = r.origin().x() + t_cyl2*r.direction().x();
  Float z2 = r.origin().z() + t_cyl2*r.direction().z();
  
  Float phi3 = atan2(z2,x2);
  phi3 = phi3 < 0 ? phi3 + 2 * M_PI : phi3;
  Float radHit3 = x2*x2 + z2*z2;
  if(has_caps && t_cyl2 < temp2 && t_cyl2 > t_min && t_cyl2 < t_max && radHit3 <= radius * radius && phi3 <= phi_max && phi3 >= phi_min) {
    vec3 p = r.point_at_parameter(t_cyl2);
    p.e[1] = -length/2;
    
    Float u = p.x() / (2.0 * radius) + 0.5;
    Float v = p.z() / (2.0 * radius) + 0.5;
    u = 1 - u;
    if(alpha_mask) {
      if(alpha_mask->value(u, v, rec.p).x() < 1) {
        return(false);
      }
    }
    rec.p = p;
    rec.normal = vec3(0,-1,0);
    rec.t = t_cyl2;
    rec.mat_ptr = mat_ptr;
    rec.u = u;
    rec.v = v;
    rec.dpdu = vec3(1, 0, 0);
    rec.dpdv = vec3(0, 0, 1);
    rec.has_bump = bump_tex ? true : false;
    
    if(bump_tex) {
      vec3 bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      rec.bump_normal = rec.normal + bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv;
      rec.bump_normal.make_unit_vector();
    }
    return(true);
  }
  temppoint = r.point_at_parameter(temp2);
  phi = atan2(temppoint.z(),temppoint.x());
  phi = phi < 0 ? phi + 2 * M_PI : phi;
  if(second_is_hit && temp2 < t_max && temp2 > t_min && 
     temppoint.y() > -length/2 && temppoint.y() < length/2 && phi <= phi_max && phi >= phi_min) {
    Float hitRad = std::sqrt(temppoint.x() * temppoint.x() + temppoint.z() * temppoint.z());
    temppoint.e[0] *= radius / hitRad;
    temppoint.e[2] *= radius / hitRad;
    rec.t = temp2;
    rec.p = temppoint;
    
    temppoint.e[1] = 0;
    rec.normal = dot(temppoint,dir) > 0 ? -temppoint / radius : temppoint / radius;
    get_cylinder_uv(rec.p, rec.u, rec.v);
    
    //Interaction information
    rec.dpdu = vec3(-phi_max * temppoint.z(), 0,  phi_max * temppoint.x());
    rec.dpdv = vec3(0, length, 0);
    rec.has_bump = bump_tex ? true : false;
    
    if(bump_tex) {
      vec3 bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      rec.bump_normal = rec.normal + bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv; 
      rec.bump_normal.make_unit_vector();
    }
    
    rec.mat_ptr = mat_ptr;
    return(true);
  }
  return(false);
}

Float cylinder::pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, rng)) {
    Float area = length * radius * (phi_max - phi_min);
    Float distance_squared = rec.t * rec.t * v.squared_length();
    Float cosine = fabs(dot(v,rec.normal)/v.length());
    return(distance_squared / (cosine * area));
  } else {
    return(0);
  }
}

vec3 cylinder::random(const vec3& o, random_gen& rng) {
  Float r1 = rng.unif_rand();
  Float y1 = length*(rng.unif_rand()-0.5);
  Float phi = (phi_max - phi_min) * r1 + phi_min;
  Float x = radius * cos(phi);
  Float z = radius * sin(phi);
  return(vec3(x,y1,z)-o);
}

vec3 cylinder::random(const vec3& o, Sampler* sampler) {
  vec2 u = sampler->Get2D();
  Float r1 = u.x();
  Float y1 = length*(u.y()-0.5);
  Float phi = (phi_max - phi_min) * r1 + phi_min;
  Float x = radius * cos(phi);
  Float z = radius * sin(phi);
  return(vec3(x,y1,z)-o);
}


bool cylinder::bounding_box(Float t0, Float t1, aabb& box) const {
  box = aabb(-vec3(radius,length/2,radius), vec3(radius,length/2,radius));
  return(true);
}

#endif
