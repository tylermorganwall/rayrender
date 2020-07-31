#ifndef CONEH
#define CONEH

#include "hitable.h"
#include "material.h"
#include "mathinline.h"

class cone: public hitable {
public:
  cone() {}
  ~cone() {
    delete mat_ptr;
    delete alpha_mask;
    delete bump_tex;
  }
  cone(Float r, Float h, material *mat, alpha_texture *alpha_mask, bump_texture* bump_tex) : 
     radius(r), height(h), mat_ptr(mat), alpha_mask(alpha_mask), bump_tex(bump_tex) {};
  virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng);
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng);
  virtual vec3 random(const vec3& o, random_gen& rng);
  virtual vec3 random(const vec3& o, Sampler* sampler);
  Float radius;
  Float height;
  //Float phi2; //unwrapped cone in a circle, angle around origin (radians)
  material *mat_ptr;
  alpha_texture *alpha_mask;
  bump_texture *bump_tex;
};

bool cone::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  vec3 oc = r.origin() - vec3(0,height, 0);
  Float k = radius / height;
  k = k*k;
  vec3 kvec = vec3(1,-k,1);
  Float a = dot(r.direction(), r.direction() * kvec);
  Float b = 2 * dot(oc, r.direction() * kvec); 
  Float c = dot(oc,oc * kvec);
  Float temp1, temp2;
  if (!quadratic(a, b, c, &temp1, &temp2)) {
    return(false);
  }
  Float t_cyl = -r.origin().y() / r.direction().y();
  Float phi;
  // Float phi;
  // bool is_hit = true;
  // bool second_is_hit = true;
  // bool base_is_hit = true;
  // if(alpha_mask) {
  //   Float u;
  //   Float v;
  //   vec3 temp_point = r.point_at_parameter(temp1);
  //   if(temp1 < t_max && temp1 > t_min && temp_point.y() > 0.0 && temp_point.y() < height) {
  //     phi = std::atan2(temp_point.z(), temp_point.x());
  //     phi += phi < 0 ? 2 * M_PI : 0;
  // 
  //     u = phi / (2*M_PI);
  //     v = temp_point.y() / height;
  //     u = 1-u;
  //     if(alpha_mask->value(u, v, temp_point).x() < rng.unif_rand()) {
  //       is_hit = false;
  //     }
  //   }
  //   Float t_cyl = -r.origin().y() / r.direction().y();
  //   Float x = r.origin().x() + t_cyl*r.direction().x();
  //   Float z = r.origin().z() + t_cyl*r.direction().z();
  //   Float radHit2 = x*x + z*z;
  //   if(t_cyl > t_min && t_cyl < t_max && radHit2 <= radius * radius) {
  //     vec3 p = r.point_at_parameter(t_cyl);
  //     u = p.x() / (2.0 * radius) + 0.5;
  //     v = p.z() / (2.0 * radius) + 0.5;
  //     u = 1 - u;
  //     if(alpha_mask->value(u, v, p).x() < rng.unif_rand()) {
  //       base_is_hit = false;
  //     }
  //   }
  //   temp_point = r.point_at_parameter(temp2);
  //   if(temp2 < t_max && temp2 > t_min && temp_point.y() > 0.0 && temp_point.y() < height) {
  //     phi = std::atan2(temp_point.z(), temp_point.x());
  //     phi += phi < 0 ? 2 * M_PI : 0;
  // 
  //     u = phi / (2*M_PI);
  //     v = temp_point.y() / height;
  //     u = 1-u;
  //     if(alpha_mask->value(u, v, temp_point).x() < rng.unif_rand()) {
  //       second_is_hit = false;
  //     }
  //   }
  // }
  vec3 temp_point = r.point_at_parameter(temp1);
  vec3 temp_point2 = r.point_at_parameter(temp2);
  
  Float x = r.origin().x() + t_cyl*r.direction().x();
  Float z = r.origin().z() + t_cyl*r.direction().z();
  Float radHit2 = x*x + z*z;
  bool hit_first  = temp1 < t_max && temp1 > t_min && 
                   temp_point.y() > 0.0 && temp_point.y() < height;
  bool hit_second = temp2 < t_max && temp2 > t_min && 
                   temp_point2.y() > 0.0 && temp_point2.y() < height;
  bool hit_base   = t_cyl < t_max && t_cyl > t_min && 
                    radHit2 <= radius * radius;
  // bool already_inside = (t_cyl < 0 || temp1 < 0) && r.origin().y() > 0.0 && r.origin().y() < height;
  
  // if(temp1 < t_max && temp1 > t_min && is_hit && (t_cyl > temp1 || !base_is_hit) && temp_point.y() > 0.0 && temp_point.y() < height) {
  if(hit_first && (t_cyl > temp1 || !hit_base)) {
    rec.t = temp1;
    rec.p = temp_point;
    phi = std::atan2(rec.p.z(), rec.p.x());
    phi += phi < 0 ? 2 * M_PI : 0;

    rec.u = phi / (2*M_PI);
    rec.v = rec.p.y() / height;
    rec.u = 1 - rec.u;
    //Interaction information
    rec.dpdu = 2 * M_PI * vec3(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = vec3(-rec.p.x() / (1 - rec.v), height, -rec.p.z()/ (1 - rec.v));
    rec.normal = unit_vector(-cross(rec.dpdu,rec.dpdv ));
    
    rec.has_bump = bump_tex ? true : false;

    if(bump_tex) {
      vec3 bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      rec.bump_normal = rec.normal + bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv;
      rec.bump_normal.make_unit_vector();
    }
    // if((!base_is_hit && t_cyl < temp1) || (already_inside && alpha_mask)) {
    //   rec.normal = -rec.normal;
    //   rec.bump_normal = -rec.bump_normal;
    // }

    rec.mat_ptr = mat_ptr;
    return(true);
  }
  // if((t_cyl < temp2 || !second_is_hit) && t_cyl > t_min && t_cyl < t_max && radHit2 <= radius * radius && base_is_hit) {
  if(hit_base && t_cyl < temp2) {
    vec3 p = r.point_at_parameter(t_cyl);
    p.e[1] = 0;

    Float u = p.x() / (2.0 * radius) + 0.5;
    Float v = p.z() / (2.0 * radius) + 0.5;
    u = 1 - u;
    rec.p = p;
    rec.normal = vec3(0,-1,0);
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
    // if((!second_is_hit && t_cyl > temp2) || (!is_hit && t_cyl > temp1) || (already_inside && alpha_mask)) {
    //   rec.normal = -rec.normal;
    //   rec.bump_normal = -rec.bump_normal;
    // }
    return(true);
  }
  // // if(temp2 < t_max && temp2 > t_min && second_is_hit && temp_point.y() >= 0.0 && temp_point.y() <= height) {
  if(hit_second) {
    rec.t = temp2;
    rec.p = temp_point2;
    phi = std::atan2(rec.p.z(), rec.p.x());
    phi += phi < 0 ? 2 * M_PI : 0;

    rec.u = phi / (2*M_PI);
    rec.v = rec.p.y() / height;
    rec.u = 1 - rec.u;

    //Interaction information
    rec.dpdu = 2 * M_PI * vec3(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = vec3(-rec.p.x() / (1 - rec.v), height, -rec.p.z()/ (1 - rec.v));
    rec.normal = unit_vector(-cross(rec.dpdu,rec.dpdv ));
    rec.has_bump = bump_tex ? true : false;

    if(bump_tex) {
      vec3 bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      rec.bump_normal = rec.normal + bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv;
      rec.bump_normal.make_unit_vector();
    }
    // if((!is_hit && temp2 < t_cyl) || (already_inside && alpha_mask)) {// || (!base_is_hit && temp2 > t_cyl)) {
    //   rec.normal = -rec.normal;
    //   rec.bump_normal = -rec.bump_normal;
    // }
    rec.mat_ptr = mat_ptr;
    return(true);
  }
  return(false);
}

Float cone::pdf_value(const vec3& o, const vec3& v, random_gen& rng) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, rng)) {
    Float maxval = ffmax(radius, 0.5f*height);
    Float cos_theta_max = sqrt(1 - maxval * maxval/o.squared_length());
    Float solid_angle = 2 * M_PI * (1-cos_theta_max);
    return(1/solid_angle);
  } else {
    return(0);
  }
}

vec3 cone::random(const vec3& o, random_gen& rng) {
  Float r1 = sqrt(1.0 - rng.unif_rand());
  Float phi_val = 2*M_PI*rng.unif_rand();
  Float height_val = r1 * height;
  Float radius_val = height_val * radius / height;
  Float x = radius_val * cos(phi_val);
  Float y = height_val;
  Float z = radius_val * sin(phi_val);
  return(vec3(x,y,z)-o);
}

vec3 cone::random(const vec3& o, Sampler* sampler) {
  vec2 u = sampler->Get2D();
  Float r1 = sqrt(1.0 - u.x());
  Float phi_val = 2*M_PI*u.y();
  Float height_val = r1 * height;
  Float radius_val = height_val * radius / height;
  Float x = radius_val * cos(phi_val);
  Float y = height_val;
  Float z = radius_val * sin(phi_val);
  return(vec3(x,y,z)-o);
}

bool cone::bounding_box(Float t0, Float t1, aabb& box) const {
  box = aabb(vec3(-radius,0,-radius), vec3(radius,height,radius));
  return(true);
}

#endif
