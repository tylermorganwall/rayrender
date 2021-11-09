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
  bool alpha_miss = false;
  if(alpha_mask) {
    if(alpha_mask->value(u, v, rec.p).x() < rng.unif_rand()) {
      alpha_miss = true;
    }
    rec.normal = dot(r2.direction(),normal3f(0,0,1)) < 0 ? normal3f(0,0,1) : normal3f(0,0,-1);
  } else {
    rec.normal = normal3f(0,0,1);
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
  }
  
  rec.mat_ptr = mp.get();
  rec.p = r2.point_at_parameter(t);
  rec.p.e[2] = k;
  rec.pError = vec3f(0,0,0);
  
  rec = (*ObjectToWorld)(rec);
  if(!alpha_mask) {
    rec.normal *= reverseOrientation  ? -1 : 1;
    rec.bump_normal *= reverseOrientation  ? -1 : 1;
  }
  rec.shape = this;
  rec.alpha_miss = alpha_miss;
  
  
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
  bool alpha_miss = false;
  
  if(alpha_mask) {
    if(alpha_mask->value(u, v, rec.p).x() < sampler->Get1D()) {
      alpha_miss = true;
    }
    rec.normal = dot(r2.direction(),vec3f(0,0,1)) < 0 ? vec3f(0,0,1) : vec3f(0,0,-1);
  } else {
    rec.normal = normal3f(0,0,1);
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
  }
  
  rec.mat_ptr = mp.get();
  rec.p = r2.point_at_parameter(t);
  rec.p.e[2] = k;
  rec.pError = vec3f(0,0,0);
  
  rec = (*ObjectToWorld)(rec);
  if(!alpha_mask) {
    rec.normal *= reverseOrientation  ? -1 : 1;
    rec.bump_normal *= reverseOrientation  ? -1 : 1;
  }
  rec.shape = this;
  rec.alpha_miss = alpha_miss;
  
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

vec3f xy_rect::random(const point3f& o, random_gen& rng, Float time) {
  point3f random_point = (*ObjectToWorld)(point3f(x0 + rng.unif_rand() * (x1 - x0), y0 + rng.unif_rand() * (y1-y0),k));
  return(random_point - o);
}

vec3f xy_rect::random(const point3f& o, Sampler* sampler, Float time) {
  vec2f u = sampler->Get2D();
  point3f random_point = (*ObjectToWorld)(point3f(x0 + u.x() * (x1 - x0), y0 + u.y() * (y1-y0),k));
  return(random_point - o);
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
  bool alpha_miss = false;
  
  if(alpha_mask) {
    if(alpha_mask->value(u, v, rec.p).x() < rng.unif_rand()) {
      alpha_miss = true;
    }
    rec.normal =  dot(r2.direction(),normal3f(0,1,0)) < 0 ? normal3f(0,1,0) : normal3f(0,-1,0);
  } else {
    rec.normal =  normal3f(0,1,0);
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
  }
  
  rec.mat_ptr = mp.get();
  rec.p = r2.point_at_parameter(t);
  rec.p.e[1] = k;
  rec.pError = vec3f(0,0,0);
  rec = (*ObjectToWorld)(rec);
  if(!alpha_mask) {
    rec.normal *= reverseOrientation  ? -1 : 1;
    rec.bump_normal *= reverseOrientation  ? -1 : 1;
  }
  rec.shape = this;
  rec.alpha_miss = alpha_miss;
  
  
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
  bool alpha_miss = false;
  
  if(alpha_mask) {
    if(alpha_mask->value(u, v, rec.p).x() < sampler->Get1D()) {
      alpha_miss = true;
    }
    rec.normal =  dot(r2.direction(),normal3f(0,1,0)) < 0 ? normal3f(0,1,0) : normal3f(0,-1,0);
  } else {
    rec.normal =  normal3f(0,1,0);
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
  }
  
  rec.mat_ptr = mp.get();
  rec.p = r2.point_at_parameter(t);
  rec.p.e[1] = k;
  rec.pError = vec3f(0,0,0);
  
  rec = (*ObjectToWorld)(rec);
  if(!alpha_mask) {
    rec.normal *= reverseOrientation  ? -1 : 1;
    rec.bump_normal *= reverseOrientation  ? -1 : 1;
  }
  rec.shape = this;
  rec.alpha_miss = alpha_miss;
  
  return(true);
}

bool xz_rect::bounding_box(Float t0, Float t1, aabb& box) const {
  box = (*ObjectToWorld)(aabb(vec3f(x0,k-0.001,z0), vec3f(x1,k+0.001,z1)));
  return(true);
}
Float xz_rect::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, rng)) {
    vec3f v2 = (*WorldToObject)(v);
    Float area = (x1-x0)*(z1-z0);
    Float distance_squared = rec.t * rec.t * v2.squared_length();
    Float cosine = fabs(dot(v, rec.normal)/v2.length());
    return(distance_squared / (cosine * area));
  } else {
    return(0);
  }
}
Float xz_rect::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, sampler)) {
    vec3f v2 = (*WorldToObject)(v);
    
    Float area = (x1-x0)*(z1-z0);
    Float distance_squared = rec.t * rec.t * v2.squared_length();
    Float cosine = fabs(dot(v, rec.normal)/v2.length());
    return(distance_squared / (cosine * area));
  } else {
    return(0);
  }
}
vec3f xz_rect::random(const point3f& o, random_gen& rng, Float time) {
  point3f random_point = (*ObjectToWorld)(point3f(x0 + rng.unif_rand() * (x1 - x0), k, z0 + rng.unif_rand() * (z1-z0)));
  return(random_point - o);
}
vec3f xz_rect::random(const point3f& o, Sampler* sampler, Float time) {
  vec2f u = sampler->Get2D();
  point3f random_point = (*ObjectToWorld)(point3f(x0 + u.x() * (x1 - x0), k, z0 + u.y()  * (z1-z0)));
  return(random_point - o);
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
  bool alpha_miss = false;
  
  if(alpha_mask) {
    if(alpha_mask->value(u, v, rec.p).x() < rng.unif_rand()) {
      alpha_miss = true;
    }
    rec.normal =  dot(r2.direction(),normal3f(1,0,0)) < 0 ? normal3f(1,0,0) : normal3f(-1,0,0);
  } else {
    rec.normal =  normal3f(1,0,0);
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
  }
  
  rec.mat_ptr = mp.get();
  rec.p = r2.point_at_parameter(t);
  rec.p.e[0] = k;
  rec.pError = vec3f(0,0,0);
  
  rec = (*ObjectToWorld)(rec);
  if(!alpha_mask) {
    rec.normal *= reverseOrientation  ? -1 : 1;
    rec.bump_normal *= reverseOrientation  ? -1 : 1;
  }
  rec.shape = this;
  rec.alpha_miss = alpha_miss;
  
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
  bool alpha_miss = false;
  
  if(alpha_mask) {
    if(alpha_mask->value(u, v, rec.p).x() < sampler->Get1D()) {
      alpha_miss = true;
    }
    rec.normal =  dot(r2.direction(),normal3f(1,0,0)) < 0 ? normal3f(1,0,0) : normal3f(-1,0,0);
  } else {
    rec.normal =  normal3f(1,0,0);
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
  }
  
  rec.mat_ptr = mp.get();
  rec.p = r2.point_at_parameter(t);
  rec.p.e[0] = k;
  rec.pError = vec3f(0,0,0);
  
  rec = (*ObjectToWorld)(rec);
  if(!alpha_mask) {
    rec.normal *= reverseOrientation  ? -1 : 1;
    rec.bump_normal *= reverseOrientation  ? -1 : 1;
  }
  rec.shape = this;
  rec.alpha_miss = alpha_miss;
  
  return(true);
}

bool yz_rect::bounding_box(Float t0, Float t1, aabb& box) const {
  box = (*ObjectToWorld)(aabb(vec3f(k-0.001,y0,z0), vec3f(k+0.001,y1,z1)));
  return(true);
}

Float yz_rect::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, rng)) {
    Float area = (y1-y0)*(z1-z0);
    Float distance_squared = rec.t * rec.t * v.squared_length();
    Float cosine = fabs(dot(v,rec.normal)/v.length());
    return(distance_squared / (cosine * area));
  } else {
    return(0);
  }
}

Float yz_rect::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, sampler)) {
    Float area = (y1-y0)*(z1-z0);
    Float distance_squared = rec.t * rec.t * v.squared_length();
    Float cosine = fabs(dot(v,rec.normal)/v.length());
    return(distance_squared / (cosine * area));
  } else {
    return(0);
  }
}

vec3f yz_rect::random(const point3f& o, random_gen& rng, Float time) {
  point3f random_point = (*ObjectToWorld)(point3f(k, y0 + rng.unif_rand() * (y1 - y0), z0 + rng.unif_rand() * (z1-z0)));
  return(random_point-o);
}

vec3f yz_rect::random(const point3f& o, Sampler* sampler, Float time) {
  vec2f u = sampler->Get2D();
  point3f random_point = (*ObjectToWorld)(point3f(k, y0 + u.x() * (y1 - y0), z0 + u.y() * (z1-z0)));
  return(random_point - o);
}
