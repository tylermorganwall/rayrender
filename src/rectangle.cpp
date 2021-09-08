#include "rectangle.h"

bool xy_rect::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  ray r2 = (*WorldToObject)(r);
  
  Float t = (k-r2.origin().z()) / r2.direction().z();
  if(t < t_min || t > t_max) {
    return(false);
  }
  Float x = r2.origin().x() + t*r2.direction().x();
  Float y = r2.origin().y() + t*r2.direction().y();
  if(x < x0 || x > x1 || y < y0 || y > y1) {
    return(false);
  }
  Float u = (x-x0)/(x1-x0);
  Float v = (y-y0)/(y1-y0);
  if(reverseOrientation) {
    u = 1 - u;
  }
  if(alpha_mask) {
    if(alpha_mask->value(u, v, rec.p).x() < rng.unif_rand()) {
      return(false);
    }
    rec.normal = dot(r2.direction(),vec3f(0,0,1)) < 0 ? vec3f(0,0,1) : vec3f(0,0,-1);
  } else {
    rec.normal = vec3f(0,0,1);
  }
  rec.u = u;
  rec.v = v;
  rec.t = t;
  
  //Interaction information
  rec.dpdu = vec3f(-1, 0, 0);
  rec.dpdv = vec3f(0, 1, 0);
  rec.has_bump = bump_tex ? true : false;
  
  if(bump_tex) {
    point3f bvbu = bump_tex->value(u,v, rec.p);
    rec.bump_normal = rec.normal + normal3f(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv); 
    rec.bump_normal.make_unit_vector();
    rec.bump_normal = !reverseOrientation ? (*ObjectToWorld)(rec.bump_normal) : -(*ObjectToWorld)(rec.bump_normal);
  }
  
  rec.mat_ptr = mp.get();
  rec.p = r2.point_at_parameter(t);
  rec.p.e[2] = k;
  rec = (*ObjectToWorld)(rec);
  rec.normal = !reverseOrientation ? (rec.normal) : -(rec.normal);
  
  return(true);
}

bool xy_rect::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
  ray r2 = (*WorldToObject)(r);
  
  Float t = (k-r2.origin().z()) / r2.direction().z();
  if(t < t_min || t > t_max) {
    return(false);
  }
  Float x = r2.origin().x() + t*r2.direction().x();
  Float y = r2.origin().y() + t*r2.direction().y();
  if(x < x0 || x > x1 || y < y0 || y > y1) {
    return(false);
  }
  Float u = (x-x0)/(x1-x0);
  Float v = (y-y0)/(y1-y0);
  if(reverseOrientation) {
    u = 1 - u;
  }
  if(alpha_mask) {
    if(alpha_mask->value(u, v, rec.p).x() < sampler->Get1D()) {
      return(false);
    }
    rec.normal = dot(r2.direction(),vec3f(0,0,1)) < 0 ? vec3f(0,0,1) : vec3f(0,0,-1);
  } else {
    rec.normal = vec3f(0,0,1);
  }
  rec.u = u;
  rec.v = v;
  rec.t = t;
  //Interaction information
  rec.dpdu = vec3f(-1, 0, 0);
  rec.dpdv = vec3f(0, 1, 0);
  rec.has_bump = bump_tex ? true : false;
  
  if(bump_tex) {
    point3f bvbu = bump_tex->value(u,v, rec.p);
    rec.bump_normal = rec.normal + normal3f(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv); 
    rec.bump_normal.make_unit_vector();
    rec.bump_normal = !reverseOrientation ? (*ObjectToWorld)(rec.bump_normal) : -(*ObjectToWorld)(rec.bump_normal);
  }
  
  rec.mat_ptr = mp.get();
  rec.p = r2.point_at_parameter(t);
  rec.p.e[2] = k;
  rec = (*ObjectToWorld)(rec);
  rec.normal = !reverseOrientation ? (rec.normal) : -(rec.normal);
  
  return(true);
}

bool xy_rect::bounding_box(Float t0, Float t1, aabb& box) const {
  box = (*ObjectToWorld)(aabb(vec3f(x0,y0,k-0.001), vec3f(x1,y1,k+0.001)));
  return(true);
}

Float xy_rect::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) {
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

Float xy_rect::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) {
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
  ray r2 = (*WorldToObject)(r);
  
  Float t = (k-r2.origin().y()) / r2.direction().y();
  if(t < t_min || t > t_max) {
    return(false);
  }
  Float x = r2.origin().x() + t*r2.direction().x();
  Float z = r2.origin().z() + t*r2.direction().z();
  if(x < x0 || x > x1 || z < z0 || z > z1) {
    return(false);
  }
  Float u = 1-(x-x0)/(x1-x0);
  Float v = (z-z0)/(z1-z0);
  if(reverseOrientation) {
    u = 1 - u;
  }
  if(alpha_mask) {
    if(alpha_mask->value(u, v, rec.p).x() < rng.unif_rand()) {
      return(false);
    }
    rec.normal =  dot(r2.direction(),vec3f(0,1,0)) < 0 ? vec3f(0,1,0) : vec3f(0,-1,0);
  } else {
    rec.normal =  vec3f(0,1,0);
  }
           
  rec.u = u;
  rec.v = v;
  rec.t = t;
  
  //Interaction information
  rec.dpdu = vec3f(1, 0, 0);
  rec.dpdv = vec3f(0, 0, 1);
  rec.has_bump = bump_tex ? true : false;
  
  if(bump_tex) {
    point3f bvbu = bump_tex->value(u,v, rec.p);
    rec.bump_normal = rec.normal + normal3f(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv); 
    rec.bump_normal.make_unit_vector();
    rec.bump_normal = !reverseOrientation ? (*ObjectToWorld)(rec.bump_normal) : -(*ObjectToWorld)(rec.bump_normal);
  }
  
  rec.mat_ptr = mp.get();
  rec.p = r2.point_at_parameter(t);
  rec.p.e[1] = k;
  rec = (*ObjectToWorld)(rec);
  rec.normal = !reverseOrientation ? (rec.normal) : -(rec.normal);
  return(true);
}


bool xz_rect::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
  ray r2 = (*WorldToObject)(r);
  
  Float t = (k-r2.origin().y()) / r2.direction().y();
  if(t < t_min || t > t_max) {
    return(false);
  }
  Float x = r2.origin().x() + t*r2.direction().x();
  Float z = r2.origin().z() + t*r2.direction().z();
  if(x < x0 || x > x1 || z < z0 || z > z1) {
    return(false);
  }
  Float u = 1-(x-x0)/(x1-x0);
  Float v = (z-z0)/(z1-z0);
  if(reverseOrientation) {
    u = 1 - u;
  }
  if(alpha_mask) {
    if(alpha_mask->value(u, v, rec.p).x() < sampler->Get1D()) {
      return(false);
    }
    rec.normal =  dot(r2.direction(),vec3f(0,1,0)) < 0 ? vec3f(0,1,0) : vec3f(0,-1,0);
  } else {
    rec.normal =  vec3f(0,1,0);
  }
  rec.u = u;
  rec.v = v;
  rec.t = t;
  
  //Interaction information
  rec.dpdu = vec3f(1, 0, 0);
  rec.dpdv = vec3f(0, 0, 1);
  rec.has_bump = bump_tex ? true : false;
  
  if(bump_tex) {
    point3f bvbu = bump_tex->value(u,v, rec.p);
    rec.bump_normal = rec.normal + normal3f(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv); 
    rec.bump_normal.make_unit_vector();
    rec.bump_normal = !reverseOrientation ? (*ObjectToWorld)(rec.bump_normal) : -(*ObjectToWorld)(rec.bump_normal);
  }
  
  rec.mat_ptr = mp.get();
  rec.p = r2.point_at_parameter(t);
  rec.p.e[1] = k;
  rec = (*ObjectToWorld)(rec);
  rec.normal = !reverseOrientation ? (rec.normal) : -(rec.normal);
  
  return(true);
}

bool yz_rect::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  ray r2 = (*WorldToObject)(r);
  
  Float t = (k-r2.origin().x()) / r2.direction().x();
  if(t < t_min || t > t_max) {
    return(false);
  }
  Float z = r2.origin().z() + t*r2.direction().z();
  Float y = r2.origin().y() + t*r2.direction().y();
  if(z < z0 || z > z1 || y < y0 || y > y1) {
    return(false);
  }
  Float u = 1-(z-z0)/(z1-z0);
  Float v = (y-y0)/(y1-y0);
  if(reverseOrientation) {
    u = 1 - u;
  }
  if(alpha_mask) {
    if(alpha_mask->value(u, v, rec.p).x() < rng.unif_rand()) {
      return(false);
    }
    rec.normal =  dot(r2.direction(),vec3f(1,0,0)) < 0 ? vec3f(1,0,0) : vec3f(-1,0,0);
  } else {
    rec.normal =  vec3f(1,0,0);
  }
  rec.u = u;
  rec.v = v;
  rec.t = t;
  
  //Interaction information
  rec.dpdu = vec3f(0, 0, 1);
  rec.dpdv = vec3f(0, 1, 0);
  rec.has_bump = bump_tex ? true : false;
  if(bump_tex) {
    point3f bvbu = bump_tex->value(u,v, rec.p);
    rec.bump_normal = rec.normal + normal3f(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv); 
    rec.bump_normal.make_unit_vector();
    rec.bump_normal = !reverseOrientation ? (*ObjectToWorld)(rec.bump_normal) : -(*ObjectToWorld)(rec.bump_normal);
  }
  
  rec.mat_ptr = mp.get();
  rec.p = r2.point_at_parameter(t);
  rec.p.e[0] = k;
  rec = (*ObjectToWorld)(rec);
  rec.normal = !reverseOrientation ? (rec.normal) : -(rec.normal);
  
  return(true);
}


bool yz_rect::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
  ray r2 = (*WorldToObject)(r);
  
  Float t = (k-r2.origin().x()) / r2.direction().x();
  if(t < t_min || t > t_max) {
    return(false);
  }
  Float z = r2.origin().z() + t*r2.direction().z();
  Float y = r2.origin().y() + t*r2.direction().y();
  if(z < z0 || z > z1 || y < y0 || y > y1) {
    return(false);
  }
  Float u = 1-(z-z0)/(z1-z0);
  Float v = (y-y0)/(y1-y0);
  if(reverseOrientation) {
    u = 1 - u;
  }
  if(alpha_mask) {
    if(alpha_mask->value(u, v, rec.p).x() < sampler->Get1D()) {
      return(false);
    }
    rec.normal =  dot(r2.direction(),vec3f(1,0,0)) < 0 ? vec3f(1,0,0) : vec3f(-1,0,0);
  } else {
    rec.normal =  vec3f(1,0,0);
  }
  rec.u = u;
  rec.v = v;
  rec.t = t;
  
  //Interaction information
  rec.dpdu = vec3f(0, 0, 1);
  rec.dpdv = vec3f(0, 1, 0);
  rec.has_bump = bump_tex ? true : false;
  if(bump_tex) {
    point3f bvbu = bump_tex->value(u,v, rec.p);
    rec.bump_normal = rec.normal + normal3f(bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv); 
    rec.bump_normal.make_unit_vector();
    rec.bump_normal = !reverseOrientation ? (*ObjectToWorld)(rec.bump_normal) : -(*ObjectToWorld)(rec.bump_normal);
  }
  
  rec.mat_ptr = mp.get();
  rec.p = r2.point_at_parameter(t);
  rec.p.e[0] = k;
  rec = (*ObjectToWorld)(rec);
  rec.normal = !reverseOrientation ? (rec.normal) : -(rec.normal);
  
  return(true);
}
