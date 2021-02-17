#include "rectangle.h"

bool xy_rect::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  Float t = (k-r.origin().z()) / r.direction().z();
  if(t < t_min || t > t_max) {
    return(false);
  }
  Float x = r.origin().x() + t*r.direction().x();
  Float y = r.origin().y() + t*r.direction().y();
  if(x < x0 || x > x1 || y < y0 || y > y1) {
    return(false);
  }
  Float u = (x-x0)/(x1-x0);
  Float v = (y-y0)/(y1-y0);
  if(flipped) {
    u = 1 - u;
  }
  if(alpha_mask) {
    if(alpha_mask->value(u, v, rec.p).x() < rng.unif_rand()) {
      return(false);
    }
    rec.normal = dot(r.direction(),vec3(0,0,1)) < 0 ? vec3(0,0,1) : vec3(0,0,-1);
  } else {
    rec.normal = flipped ? vec3(0,0,-1) : vec3(0,0,1);
  }
  rec.u = u;
  rec.v = v;
  rec.t = t;
  
  //Interaction information
  rec.dpdu = vec3(-1, 0, 0);
  rec.dpdv = vec3(0, 1, 0);
  rec.has_bump = bump_tex ? true : false;
  
  if(bump_tex) {
    vec3 bvbu = bump_tex->value(u,v, rec.p);
    bvbu.e[0] *= flipped ? -1 : 1;
    rec.bump_normal = rec.normal + bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv; 
    rec.bump_normal.make_unit_vector();
  }
  
  rec.mat_ptr = mp.get();
  rec.p = r.point_at_parameter(t);
  rec.p.e[2] = k;
  
  return(true);
}

bool xy_rect::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
  Float t = (k-r.origin().z()) / r.direction().z();
  if(t < t_min || t > t_max) {
    return(false);
  }
  Float x = r.origin().x() + t*r.direction().x();
  Float y = r.origin().y() + t*r.direction().y();
  if(x < x0 || x > x1 || y < y0 || y > y1) {
    return(false);
  }
  Float u = (x-x0)/(x1-x0);
  Float v = (y-y0)/(y1-y0);
  if(flipped) {
    u = 1 - u;
  }
  if(alpha_mask) {
    if(alpha_mask->value(u, v, rec.p).x() < sampler->Get1D()) {
      return(false);
    }
    rec.normal = dot(r.direction(),vec3(0,0,1)) < 0 ? vec3(0,0,1) : vec3(0,0,-1);
  } else {
    rec.normal = flipped ? vec3(0,0,-1) : vec3(0,0,1);
  }
  rec.u = u;
  rec.v = v;
  rec.t = t;
  //Interaction information
  rec.dpdu = vec3(-1, 0, 0);
  rec.dpdv = vec3(0, 1, 0);
  rec.has_bump = bump_tex ? true : false;
  
  if(bump_tex) {
    vec3 bvbu = bump_tex->value(u,v, rec.p);
    bvbu.e[0] *= flipped ? -1 : 1;
    rec.bump_normal = rec.normal + bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv; 
    rec.bump_normal.make_unit_vector();
  }
  
  rec.mat_ptr = mp.get();
  rec.p = r.point_at_parameter(t);
  rec.p.e[2] = k;
  
  return(true);
}

bool xy_rect::bounding_box(Float t0, Float t1, aabb& box) const {
  box = aabb(vec3(x0,y0,k-0.001), vec3(x1,y1,k+0.001));
  return(true);
}

Float xy_rect::pdf_value(const vec3& o, const vec3& v, random_gen& rng, Float time) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, rng)) {
    Float area = (x1-x0)*(y1-y0);
    Float distance_squared = rec.t * rec.t * v.squared_length();
    Float cosine = fabs(dot(v,rec.normal)/v.length());
    return(distance_squared / (cosine * area));
  } else {
    return(0);
  }
}

Float xy_rect::pdf_value(const vec3& o, const vec3& v, Sampler* sampler, Float time) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, sampler)) {
    Float area = (x1-x0)*(y1-y0);
    Float distance_squared = rec.t * rec.t * v.squared_length();
    Float cosine = fabs(dot(v,rec.normal)/v.length());
    return(distance_squared / (cosine * area));
  } else {
    return(0);
  }
}

bool xz_rect::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  Float t = (k-r.origin().y()) / r.direction().y();
  if(t < t_min || t > t_max) {
    return(false);
  }
  Float x = r.origin().x() + t*r.direction().x();
  Float z = r.origin().z() + t*r.direction().z();
  if(x < x0 || x > x1 || z < z0 || z > z1) {
    return(false);
  }
  Float u = 1-(x-x0)/(x1-x0);
  Float v = (z-z0)/(z1-z0);
  if(flipped) {
    u = 1 - u;
  }
  if(alpha_mask) {
    if(alpha_mask->value(u, v, rec.p).x() < rng.unif_rand()) {
      return(false);
    }
    rec.normal =  dot(r.direction(),vec3(0,1,0)) < 0 ? vec3(0,1,0) : vec3(0,-1,0);
  } else {
    rec.normal =  flipped ? vec3(0,-1,0) : vec3(0,1,0);
  }
  rec.u = u;
  rec.v = v;
  rec.t = t;
  
  //Interaction information
  rec.dpdu = vec3(1, 0, 0);
  rec.dpdv = vec3(0, 0, 1);
  rec.has_bump = bump_tex ? true : false;
  
  if(bump_tex) {
    vec3 bvbu = bump_tex->value(u,v, rec.p);
    bvbu.e[0] *= flipped ? -1 : 1;
    rec.bump_normal = rec.normal + bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv; 
    rec.bump_normal.make_unit_vector();
  }
  
  rec.mat_ptr = mp.get();
  rec.p = r.point_at_parameter(t);
  rec.p.e[1] = k;
  
  return(true);
}


bool xz_rect::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
  Float t = (k-r.origin().y()) / r.direction().y();
  if(t < t_min || t > t_max) {
    return(false);
  }
  Float x = r.origin().x() + t*r.direction().x();
  Float z = r.origin().z() + t*r.direction().z();
  if(x < x0 || x > x1 || z < z0 || z > z1) {
    return(false);
  }
  Float u = 1-(x-x0)/(x1-x0);
  Float v = (z-z0)/(z1-z0);
  if(flipped) {
    u = 1 - u;
  }
  if(alpha_mask) {
    if(alpha_mask->value(u, v, rec.p).x() < sampler->Get1D()) {
      return(false);
    }
    rec.normal =  dot(r.direction(),vec3(0,1,0)) < 0 ? vec3(0,1,0) : vec3(0,-1,0);
  } else {
    rec.normal =  flipped ? vec3(0,-1,0) : vec3(0,1,0);
  }
  rec.u = u;
  rec.v = v;
  rec.t = t;
  
  //Interaction information
  rec.dpdu = vec3(1, 0, 0);
  rec.dpdv = vec3(0, 0, 1);
  rec.has_bump = bump_tex ? true : false;
  
  if(bump_tex) {
    vec3 bvbu = bump_tex->value(u,v, rec.p);
    bvbu.e[0] *= flipped ? -1 : 1;
    rec.bump_normal = rec.normal + bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv; 
    rec.bump_normal.make_unit_vector();
  }
  
  rec.mat_ptr = mp.get();
  rec.p = r.point_at_parameter(t);
  rec.p.e[1] = k;
  
  return(true);
}

bool yz_rect::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  Float t = (k-r.origin().x()) / r.direction().x();
  if(t < t_min || t > t_max) {
    return(false);
  }
  Float z = r.origin().z() + t*r.direction().z();
  Float y = r.origin().y() + t*r.direction().y();
  if(z < z0 || z > z1 || y < y0 || y > y1) {
    return(false);
  }
  Float u = 1-(z-z0)/(z1-z0);
  Float v = (y-y0)/(y1-y0);
  if(flipped) {
    u = 1 - u;
  }
  if(alpha_mask) {
    if(alpha_mask->value(u, v, rec.p).x() < rng.unif_rand()) {
      return(false);
    }
    rec.normal =  dot(r.direction(),vec3(1,0,0)) < 0 ? vec3(1,0,0) : vec3(-1,0,0);
  } else {
    rec.normal =  flipped ? vec3(-1,0,0) : vec3(1,0,0);
  }
  rec.u = u;
  rec.v = v;
  rec.t = t;
  
  //Interaction information
  rec.dpdu = vec3(0, 0, 1);
  rec.dpdv = vec3(0, 1, 0);
  rec.has_bump = bump_tex ? true : false;
  if(bump_tex) {
    vec3 bvbu = bump_tex->value(u,v, rec.p);
    bvbu.e[0] *= flipped ? -1 : 1;
    rec.bump_normal = rec.normal + bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv; 
    rec.bump_normal.make_unit_vector();
  }
  
  rec.mat_ptr = mp.get();
  rec.p = r.point_at_parameter(t);
  rec.p.e[0] = k;
  
  return(true);
}


bool yz_rect::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
  Float t = (k-r.origin().x()) / r.direction().x();
  if(t < t_min || t > t_max) {
    return(false);
  }
  Float z = r.origin().z() + t*r.direction().z();
  Float y = r.origin().y() + t*r.direction().y();
  if(z < z0 || z > z1 || y < y0 || y > y1) {
    return(false);
  }
  Float u = 1-(z-z0)/(z1-z0);
  Float v = (y-y0)/(y1-y0);
  if(flipped) {
    u = 1 - u;
  }
  if(alpha_mask) {
    if(alpha_mask->value(u, v, rec.p).x() < sampler->Get1D()) {
      return(false);
    }
    rec.normal =  dot(r.direction(),vec3(1,0,0)) < 0 ? vec3(1,0,0) : vec3(-1,0,0);
  } else {
    rec.normal =  flipped ? vec3(-1,0,0) : vec3(1,0,0);
  }
  rec.u = u;
  rec.v = v;
  rec.t = t;
  
  //Interaction information
  rec.dpdu = vec3(0, 0, 1);
  rec.dpdv = vec3(0, 1, 0);
  rec.has_bump = bump_tex ? true : false;
  if(bump_tex) {
    vec3 bvbu = bump_tex->value(u,v, rec.p);
    bvbu.e[0] *= flipped ? -1 : 1;
    rec.bump_normal = rec.normal + bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv; 
    rec.bump_normal.make_unit_vector();
  }
  
  rec.mat_ptr = mp.get();
  rec.p = r.point_at_parameter(t);
  rec.p.e[0] = k;
  
  return(true);
}
