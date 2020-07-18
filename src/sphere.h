#ifndef SPHEREH
#define SPHEREH

#include "hitable.h"
#include "material.h"
#include "mathinline.h"

class sphere: public hitable {
  public:
    sphere() {}
    ~sphere() {
      delete mat_ptr;
      delete alpha_mask;
      delete bump_tex;
    }
    sphere(vec3 cen, Float r, material *mat, alpha_texture *alpha_mask, bump_texture* bump_tex) : center(cen), radius(r), 
           mat_ptr(mat), alpha_mask(alpha_mask), bump_tex(bump_tex) {};
    virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng);
    virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
    virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng);
    virtual vec3 random(const vec3& o, random_gen& rng);
    virtual vec3 random(const vec3& o, Sampler* sampler);
    vec3 center;
    Float radius;
    material *mat_ptr;
    alpha_texture *alpha_mask;
    bump_texture *bump_tex;
};

bool sphere::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  vec3 oc = r.origin() - center;
  Float a = dot(r.direction(), r.direction());
  Float b = 2 * dot(oc, r.direction()); 
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
      vec3 p1 = r.point_at_parameter(temp1);
      p1 *= radius / p1.length(); 
      vec3 normal = (p1 - center) / radius;
      get_sphere_uv(normal, u, v);
      if(alpha_mask->value(u, v, rec.p).x() < rng.unif_rand()) {
        is_hit = false;
      }
    }
    if(temp2 < t_max && temp2 > t_min) {
      vec3 p2 = r.point_at_parameter(temp2);
      p2 *= radius / p2.length(); 
      vec3 normal = (p2 - center) / radius;
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
    rec.p = r.point_at_parameter(rec.t);
    rec.p *= radius / rec.p.length(); 
    rec.normal = (rec.p - center) / radius;
    
    //Interaction information
    Float zRadius = std::sqrt(rec.p.x() * rec.p.x()  + rec.p.z()  * rec.p.z() );
    Float invZRadius = 1 / zRadius;
    Float cosPhi = rec.p.x() * invZRadius;
    Float sinPhi = rec.p.z() * invZRadius;
    Float theta = std::acos(clamp(rec.p.z() / radius, -1, 1));
    rec.dpdu = 2 * M_PI * vec3(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = 2 * M_PI * vec3(rec.p.z() * cosPhi, rec.p.z() * sinPhi, -radius * std::sin(theta));
    get_sphere_uv(rec.normal, rec.u, rec.v);
    rec.has_bump = bump_tex ? true : false;
    
    if(bump_tex) {
      vec3 bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      rec.bump_normal = rec.normal + bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv; 
      rec.bump_normal.make_unit_vector();
    }
    
    rec.mat_ptr = mat_ptr;
    return(true);
  }
  if(temp2 < t_max && temp2 > t_min && second_is_hit) {
    rec.t = temp2;
    rec.p = r.point_at_parameter(rec.t);
    rec.p *= radius / rec.p.length();
    rec.normal = (rec.p - center) / radius;
    
    //Interaction information
    Float zRadius = std::sqrt(rec.p.x() * rec.p.x()  + rec.p.z()  * rec.p.z() );
    Float invZRadius = 1 / zRadius;
    Float cosPhi = rec.p.x() * invZRadius;
    Float sinPhi = rec.p.z() * invZRadius;
    Float theta = std::acos(clamp(rec.p.z() / radius, -1, 1));
    rec.dpdu = 2 * M_PI * vec3(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = 2 * M_PI * vec3(rec.p.z() * cosPhi, rec.p.z() * sinPhi, -radius * std::sin(theta));
    get_sphere_uv(rec.normal, rec.u, rec.v);
    rec.has_bump = bump_tex ? true : false;
    
    if(bump_tex) {
      vec3 bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      rec.bump_normal = rec.normal +  bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv; 
      rec.bump_normal.make_unit_vector();
    }
    
    if(alpha_mask) {
      rec.normal = -rec.normal;
      rec.bump_normal = -rec.bump_normal;
    }
    rec.mat_ptr = mat_ptr;
    return(true);
  }
  return(false);
}

Float sphere::pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, rng)) {
    Float cos_theta_max = sqrt(1 - radius * radius/(center - o).squared_length());
    Float solid_angle = 2 * M_PI * (1-cos_theta_max);
    return(1/solid_angle);
  } else {
    return(0);
  }
}

vec3 sphere::random(const vec3& o, random_gen& rng) {
  vec3 direction = center - o;
  Float distance_squared = direction.squared_length();
  onb uvw;
  uvw.build_from_w(direction);
  return(uvw.local_to_world(rng.random_to_sphere(radius,distance_squared)));
}

vec3 sphere::random(const vec3& o, Sampler* sampler) {
  vec3 direction = center - o;
  Float distance_squared = direction.squared_length();
  onb uvw;
  uvw.build_from_w(direction);
  return(uvw.local_to_world(rand_to_sphere(radius,distance_squared, sampler->Get2D())));
}

bool sphere::bounding_box(Float t0, Float t1, aabb& box) const {
  box = aabb(center - vec3(radius,radius,radius), center + vec3(radius,radius,radius));
  return(true);
}

class moving_sphere: public hitable {
  public:
    moving_sphere() {}
    moving_sphere(vec3 cen0, vec3 cen1, Float t0, Float t1, Float r, 
                  material *mat,alpha_texture *alpha_mask, bump_texture* bump_tex) : 
                  center0(cen0), center1(cen1), time0(t0), time1(t1), radius(r), 
                  mat_ptr(mat), alpha_mask(alpha_mask), bump_tex(bump_tex) {};
    virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng);
    virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
    vec3 center(Float time) const;
    vec3 center0, center1;
    Float time0, time1;
    Float radius;
    material *mat_ptr;
    alpha_texture *alpha_mask;
    bump_texture *bump_tex;
};

vec3 moving_sphere::center(Float time) const {
  return(center0 + ((time - time0) / (time1 - time0)) * (center1 - center0));
}

bool moving_sphere::bounding_box(Float t0, Float t1, aabb& box) const {
  aabb box0(center(t0) - vec3(radius,radius,radius), center(t0) + vec3(radius,radius,radius));
  aabb box1(center(t1) - vec3(radius,radius,radius), center(t1) + vec3(radius,radius,radius));
  box = surrounding_box(box0, box1);
  return(true);
}

bool moving_sphere::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  vec3 oc = r.origin() - center(r.time());
  Float a = dot(r.direction(), r.direction());
  Float b = 2 * dot(oc, r.direction()); 
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
      vec3 p1 = r.point_at_parameter(temp1);
      p1 *= radius / p1.length(); 
      vec3 normal = (p1 - center(r.time())) / radius;
      get_sphere_uv(normal, u, v);
      if(alpha_mask->value(u, v, rec.p).x() < rng.unif_rand()) {
        is_hit = false;
      }
    }
    if(temp2 < t_max && temp2 > t_min) {
      vec3 p2 = r.point_at_parameter(temp2);
      p2 *= radius / p2.length(); 
      vec3 normal = (p2 - center(r.time())) / radius;
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
    rec.p = r.point_at_parameter(rec.t);
    rec.p *= radius / rec.p.length(); 
    rec.normal = (rec.p - center(r.time())) / radius;
    
    //Interaction information
    Float zRadius = std::sqrt(rec.p.x() * rec.p.x()  + rec.p.z()  * rec.p.z() );
    Float invZRadius = 1 / zRadius;
    Float cosPhi = rec.p.x() * invZRadius;
    Float sinPhi = rec.p.z() * invZRadius;
    Float theta = std::acos(clamp(rec.p.z() / radius, -1, 1));
    rec.dpdu = 2 * M_PI * vec3(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = 2 * M_PI * vec3(rec.p.z() * cosPhi, rec.p.z() * sinPhi, -radius * std::sin(theta));
    
    if(bump_tex) {
      vec3 bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      rec.bump_normal = rec.normal + bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv; 
      rec.bump_normal.make_unit_vector();
    }
    
    get_sphere_uv(rec.normal, rec.u, rec.v);
    rec.mat_ptr = mat_ptr;
    return(true);
  }
  if(temp2 < t_max && temp2 > t_min) {
    rec.t = temp2;
    rec.p = r.point_at_parameter(rec.t);
    rec.p *= radius / rec.p.length(); 
    rec.normal = (rec.p - center(r.time())) / radius;
    
    //Interaction information
    Float zRadius = std::sqrt(rec.p.x() * rec.p.x()  + rec.p.z()  * rec.p.z() );
    Float invZRadius = 1 / zRadius;
    Float cosPhi = rec.p.x() * invZRadius;
    Float sinPhi = rec.p.z() * invZRadius;
    Float theta = std::acos(clamp(rec.p.z() / radius, -1, 1));
    rec.dpdu = 2 * M_PI * vec3(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = 2 * M_PI * vec3(rec.p.z() * cosPhi, rec.p.z() * sinPhi, -radius * std::sin(theta));
    
    if(bump_tex) {
      vec3 bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      vec3 o_u = cross(rec.normal, rec.dpdu);
      vec3 o_v = cross(rec.normal, rec.dpdv);
      rec.bump_normal = rec.normal + bvbu.x() * o_v - bvbu.y() * o_u; 
      rec.bump_normal.make_unit_vector();
    }
    
    get_sphere_uv(rec.normal, rec.u, rec.v);
    if(!is_hit) {
      rec.normal = -rec.normal;
      rec.bump_normal = -rec.bump_normal;
    }
    rec.mat_ptr = mat_ptr;
    return(true);
  }
  return(false);
}

#endif
