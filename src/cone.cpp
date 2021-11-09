#include "cone.h"
#include "efloat.h"

bool cone::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  ray r2 = (*WorldToObject)(r);
  vec3f oc = r2.origin() - point3f(0,height, 0);
  Float k = radius / height;
  k = k*k;
  vec3f kvec = vec3f(1,-k,1);
  Float a = dot(r2.direction(), r2.direction() * kvec);
  Float b = 2 * dot(oc, r2.direction() * kvec); 
  Float c = dot(oc,oc * kvec);
  Float temp1, temp2;
  if (!quadratic(a, b, c, &temp1, &temp2)) {
    return(false);
  }
  Float t_cyl = -r2.origin().y() / r2.direction().y();
  Float phi;
  // Float phi;
  // bool is_hit = true;
  // bool second_is_hit = true;
  // bool base_is_hit = true;
  // if(alpha_mask) {
  //   Float u;
  //   Float v;
  //   vec3f temp_point = r.point_at_parameter(temp1);
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
  //     vec3f p = r.point_at_parameter(t_cyl);
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
  point3f temp_point = r2.point_at_parameter(temp1);
  point3f temp_point2 = r2.point_at_parameter(temp2);
  
  Float x = r2.origin().x() + t_cyl*r2.direction().x();
  Float z = r2.origin().z() + t_cyl*r2.direction().z();
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
    rec.dpdu = 2 * M_PI * vec3f(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = vec3f(-rec.p.x() / (1 - rec.v), height, -rec.p.z()/ (1 - rec.v));
    rec.normal = unit_vector(-cross(rec.dpdu,rec.dpdv ));
    
    rec.has_bump = bump_tex ? true : false;
    
    if(bump_tex) {
      point3f bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      rec.bump_normal = rec.normal + normal3f(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv);
      rec.bump_normal.make_unit_vector();

    }
    // if((!base_is_hit && t_cyl < temp1) || (already_inside && alpha_mask)) {
    //   rec.normal = -rec.normal;
    //   rec.bump_normal = -rec.bump_normal;
    // }
    EFloat px = r2.origin().x() + rec.t * r2.direction().x();
    EFloat py = r2.origin().y() + rec.t * r2.direction().y();
    EFloat pz = r2.origin().z() + rec.t * r2.direction().z();
    rec.pError = vec3f(px.GetAbsoluteError(), py.GetAbsoluteError(),
                       pz.GetAbsoluteError());
    
    rec = (*ObjectToWorld)(rec);
    rec.normal *= reverseOrientation  ? -1 : 1;
    rec.bump_normal *= reverseOrientation  ? -1 : 1;
    rec.shape = this;
    rec.alpha_miss = false;
    
    rec.mat_ptr = mat_ptr.get();
    return(true);
  }
  // if((t_cyl < temp2 || !second_is_hit) && t_cyl > t_min && t_cyl < t_max && radHit2 <= radius * radius && base_is_hit) {
  if(hit_base && t_cyl < temp2) {
    point3f p = r2.point_at_parameter(t_cyl);
    p.e[1] = 0;
    
    Float u = p.x() / (2.0 * radius) + 0.5;
    Float v = p.z() / (2.0 * radius) + 0.5;
    u = 1 - u;
    rec.p = p;
    rec.normal = vec3f(0,-1,0);
    rec.t = t_cyl;
    rec.mat_ptr = mat_ptr.get();
    rec.u = u;
    rec.v = v;
    rec.dpdu = vec3f(1, 0, 0);
    rec.dpdv = vec3f(0, 0, 1);
    rec.has_bump = bump_tex ? true : false;
    
    if(bump_tex) {
      point3f bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      rec.bump_normal = rec.normal + normal3f(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv);
      rec.bump_normal.make_unit_vector();

    }
    // if((!second_is_hit && t_cyl > temp2) || (!is_hit && t_cyl > temp1) || (already_inside && alpha_mask)) {
    //   rec.normal = -rec.normal;
    //   rec.bump_normal = -rec.bump_normal;
    // }
    EFloat px = r2.origin().x() + rec.t * r2.direction().x();
    EFloat py = r2.origin().y() + rec.t * r2.direction().y();
    EFloat pz = r2.origin().z() + rec.t * r2.direction().z();
    rec.pError = vec3f(px.GetAbsoluteError(), py.GetAbsoluteError(),
                       pz.GetAbsoluteError());
    
    rec = (*ObjectToWorld)(rec);
    rec.normal *= reverseOrientation  ? -1 : 1;
    rec.bump_normal *= reverseOrientation  ? -1 : 1;
    rec.shape = this;
    rec.alpha_miss = false;
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
    rec.dpdu = 2 * M_PI * vec3f(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = vec3f(-rec.p.x() / (1 - rec.v), height, -rec.p.z()/ (1 - rec.v));
    rec.normal = unit_vector(-cross(rec.dpdu,rec.dpdv ));
    rec.has_bump = bump_tex ? true : false;
    
    if(bump_tex) {
      point3f bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      rec.bump_normal = rec.normal + normal3f(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv);
      rec.bump_normal.make_unit_vector();

    }
    // if((!is_hit && temp2 < t_cyl) || (already_inside && alpha_mask)) {// || (!base_is_hit && temp2 > t_cyl)) {
    //   rec.normal = -rec.normal;
    //   rec.bump_normal = -rec.bump_normal;
    // }
    EFloat px = r2.origin().x() + rec.t * r2.direction().x();
    EFloat py = r2.origin().y() + rec.t * r2.direction().y();
    EFloat pz = r2.origin().z() + rec.t * r2.direction().z();
    rec.pError = vec3f(px.GetAbsoluteError(), py.GetAbsoluteError(),
                       pz.GetAbsoluteError());
    
    rec = (*ObjectToWorld)(rec);
    rec.normal *= reverseOrientation  ? -1 : 1;
    rec.bump_normal *= reverseOrientation  ? -1 : 1;
    rec.shape = this;
    rec.alpha_miss = false;
    
    rec.mat_ptr = mat_ptr.get();
    return(true);
  }
  return(false);
}


bool cone::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
  ray r2 = (*WorldToObject)(r);
  
  vec3f oc = r2.origin() - point3f(0,height, 0);
  Float k = radius / height;
  k = k*k;
  vec3f kvec = vec3f(1,-k,1);
  Float a = dot(r2.direction(), r2.direction() * kvec);
  Float b = 2 * dot(oc, r2.direction() * kvec); 
  Float c = dot(oc,oc * kvec);
  Float temp1, temp2;
  if (!quadratic(a, b, c, &temp1, &temp2)) {
    return(false);
  }
  Float t_cyl = -r2.origin().y() / r2.direction().y();
  Float phi;
  // Float phi;
  // bool is_hit = true;
  // bool second_is_hit = true;
  // bool base_is_hit = true;
  // if(alpha_mask) {
  //   Float u;
  //   Float v;
  //   vec3f temp_point = r.point_at_parameter(temp1);
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
  //     vec3f p = r.point_at_parameter(t_cyl);
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
  point3f temp_point = r2.point_at_parameter(temp1);
  point3f temp_point2 = r2.point_at_parameter(temp2);
  
  Float x = r2.origin().x() + t_cyl*r2.direction().x();
  Float z = r2.origin().z() + t_cyl*r2.direction().z();
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
    rec.dpdu = 2 * M_PI * vec3f(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = vec3f(-rec.p.x() / (1 - rec.v), height, -rec.p.z()/ (1 - rec.v));
    rec.normal = unit_vector(-cross(rec.dpdu,rec.dpdv ));
    
    rec.has_bump = bump_tex ? true : false;
    
    if(bump_tex) {
      point3f bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      rec.bump_normal = rec.normal + normal3f(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv);
      rec.bump_normal.make_unit_vector();

    }
    // if((!base_is_hit && t_cyl < temp1) || (already_inside && alpha_mask)) {
    //   rec.normal = -rec.normal;
    //   rec.bump_normal = -rec.bump_normal;
    // }
    
    EFloat px = r2.origin().x() + rec.t * r2.direction().x();
    EFloat py = r2.origin().y() + rec.t * r2.direction().y();
    EFloat pz = r2.origin().z() + rec.t * r2.direction().z();
    rec.pError = vec3f(px.GetAbsoluteError(), py.GetAbsoluteError(),
                       pz.GetAbsoluteError());
    
    rec = (*ObjectToWorld)(rec);
    rec.normal *= reverseOrientation  ? -1 : 1;
    rec.bump_normal *= reverseOrientation  ? -1 : 1;
    rec.shape = this;
    rec.alpha_miss = false;
    
    rec.mat_ptr = mat_ptr.get();
    return(true);
  }
  // if((t_cyl < temp2 || !second_is_hit) && t_cyl > t_min && t_cyl < t_max && radHit2 <= radius * radius && base_is_hit) {
  if(hit_base && t_cyl < temp2) {
    point3f p = r2.point_at_parameter(t_cyl);
    p.e[1] = 0;
    
    Float u = p.x() / (2.0 * radius) + 0.5;
    Float v = p.z() / (2.0 * radius) + 0.5;
    u = 1 - u;
    rec.p = p;
    rec.normal = vec3f(0,-1,0);
    rec.t = t_cyl;
    rec.mat_ptr = mat_ptr.get();
    rec.u = u;
    rec.v = v;
    rec.dpdu = vec3f(1, 0, 0);
    rec.dpdv = vec3f(0, 0, 1);
    rec.has_bump = bump_tex ? true : false;
    
    if(bump_tex) {
      point3f bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      rec.bump_normal = rec.normal + normal3f(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv);
      rec.bump_normal.make_unit_vector();

    }
    // if((!second_is_hit && t_cyl > temp2) || (!is_hit && t_cyl > temp1) || (already_inside && alpha_mask)) {
    //   rec.normal = -rec.normal;
    //   rec.bump_normal = -rec.bump_normal;
    // }
    EFloat px = r2.origin().x() + rec.t * r2.direction().x();
    EFloat py = r2.origin().y() + rec.t * r2.direction().y();
    EFloat pz = r2.origin().z() + rec.t * r2.direction().z();
    rec.pError = vec3f(px.GetAbsoluteError(), py.GetAbsoluteError(),
                       pz.GetAbsoluteError());
    
    rec = (*ObjectToWorld)(rec);
    rec.normal *= reverseOrientation  ? -1 : 1;
    rec.bump_normal *= reverseOrientation  ? -1 : 1;
    rec.shape = this;
    rec.alpha_miss = false;
    
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
    rec.dpdu = 2 * M_PI * vec3f(-rec.p.z(), 0, rec.p.x());
    rec.dpdv = vec3f(-rec.p.x() / (1 - rec.v), height, -rec.p.z()/ (1 - rec.v));
    rec.normal = unit_vector(-cross(rec.dpdu,rec.dpdv ));
    rec.has_bump = bump_tex ? true : false;
    
    if(bump_tex) {
      point3f bvbu = bump_tex->value(rec.u,rec.v, rec.p);
      rec.bump_normal = rec.normal + normal3f(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv);
      rec.bump_normal.make_unit_vector();
      rec.bump_normal = (*ObjectToWorld)(rec.bump_normal);
      
    }
    // if((!is_hit && temp2 < t_cyl) || (already_inside && alpha_mask)) {// || (!base_is_hit && temp2 > t_cyl)) {
    //   rec.normal = -rec.normal;
    //   rec.bump_normal = -rec.bump_normal;
    // }
    EFloat px = r2.origin().x() + rec.t * r2.direction().x();
    EFloat py = r2.origin().y() + rec.t * r2.direction().y();
    EFloat pz = r2.origin().z() + rec.t * r2.direction().z();
    rec.pError = vec3f(px.GetAbsoluteError(), py.GetAbsoluteError(),
                       pz.GetAbsoluteError());
    
    rec = (*ObjectToWorld)(rec);
    rec.normal *= reverseOrientation  ? -1 : 1;
    rec.bump_normal *= reverseOrientation  ? -1 : 1;
    rec.shape = this;
    rec.alpha_miss = false;
    
    rec.mat_ptr = mat_ptr.get();
    return(true);
  }
  return(false);
}

Float cone::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, rng)) {
    point3f o2 = (*WorldToObject)(o);
    
    Float maxval = ffmax(radius, 0.5f*height);
    Float cos_theta_max = sqrt(1 - maxval * maxval/o2.squared_length());
    Float solid_angle = 2 * M_PI * (1-cos_theta_max);
    return(1/solid_angle);
  } else {
    return(0);
  }
}


Float cone::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, sampler)) {
    point3f o2 = (*WorldToObject)(o);
    Float maxval = ffmax(radius, 0.5f*height);
    Float cos_theta_max = sqrt(1 - maxval * maxval/o2.squared_length());
    Float solid_angle = 2 * M_PI * (1-cos_theta_max);
    return(1/solid_angle);
  } else {
    return(0);
  }
}

vec3f cone::random(const point3f& o, random_gen& rng, Float time) {
  Float r1 = sqrt(1.0 - rng.unif_rand());
  Float phi_val = 2*M_PI*rng.unif_rand();
  Float height_val = r1 * height;
  Float radius_val = height_val * radius / height;
  Float x = radius_val * cos(phi_val);
  Float y = height_val;
  Float z = radius_val * sin(phi_val);
  return((*ObjectToWorld)(point3f(x,y,z)-o));
}

vec3f cone::random(const point3f& o, Sampler* sampler, Float time) {
  vec2f u = sampler->Get2D();
  Float r1 = sqrt(1.0 - u.x());
  Float phi_val = 2*M_PI*u.y();
  Float height_val = r1 * height;
  Float radius_val = height_val * radius / height;
  Float x = radius_val * cos(phi_val);
  Float y = height_val;
  Float z = radius_val * sin(phi_val);
  return((*ObjectToWorld)(point3f(x,y,z)-o));
}

bool cone::bounding_box(Float t0, Float t1, aabb& box) const {
  box = (*ObjectToWorld)(aabb(vec3f(-radius,0,-radius), vec3f(radius,height,radius)));
  return(true);
}
