#include "hitable.h"

//Translate implementation

void get_sphere_uv(const vec3& p, Float& u, Float& v) {
  Float phi = atan2(p.z(),p.x());
  Float theta = asin(p.y());
  u = 1 - (phi + M_PI) / (2*M_PI);
  v = (theta + M_PI/2) / M_PI;
}


bool translate::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  ray moved_r(r.origin()-offset, r.direction(), r.time());
  if(ptr->hit(moved_r, t_min, t_max, rec, rng)) {
    rec.p += offset;
    return(true);
  } else {
    return(false);
  }
}

bool translate::bounding_box(Float t0, Float t1, aabb& box) const {
  if(ptr->bounding_box(t0,t1,box)) {
    box = aabb(box.min() + offset, box.max() + offset);
    return(true);
  } else {
    return(false);
  }
}

//Scale implementation

bool scale::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  ray scaled_r(r.origin() * inv_scale, r.direction() * inv_scale, r.time());
  if(ptr->hit(scaled_r, t_min, t_max, rec, rng)) {
    rec.p *= scale_factor;
    rec.normal *= scale_factor;
    rec.normal.make_unit_vector();
    if(rec.has_bump) {
      rec.bump_normal *= scale_factor;
      rec.bump_normal.make_unit_vector();
    }
    return(true);
  } else {
    return(false);
  }
}

bool scale::bounding_box(Float t0, Float t1, aabb& box) const {
  if(ptr->bounding_box(t0,t1,box)) {
    box = aabb(box.min() * scale_factor, box.max() * scale_factor);
    return(true);
  } else {
    return(false);
  }
}

//Rotate implementations


rotate_y::rotate_y(std::shared_ptr<hitable> p, Float angle) : ptr(p) {
  Float radians = (M_PI / 180.0) * angle;
  sin_theta = sin(radians);
  cos_theta = cos(radians);
  hasbox = ptr->bounding_box(0,1,bbox);
  vec3 min(FLT_MAX, FLT_MAX, FLT_MAX);
  vec3 max(-FLT_MAX, -FLT_MAX, -FLT_MAX);
  for(int i = 0; i < 2; i++) {
    for(int j = 0; j < 2; j++) {
      for(int k = 0; k < 2; k++) {
        Float x = i*bbox.max().x() + (1-i)*bbox.min().x();
        Float y = j*bbox.max().y() + (1-j)*bbox.min().y();
        Float z = k*bbox.max().z() + (1-k)*bbox.min().z();
        Float newx = cos_theta*x + sin_theta*z;
        Float newz = -sin_theta*x + cos_theta*z;
        vec3 tester(newx,y,newz);
        for(int c = 0; c < 3; c++) {
          if(tester[c] > max[c]) {
            max.e[c] = tester[c];
          }
          if(tester[c] < min[c]) {
            min.e[c] = tester[c];
          }
        }
      } 
    }
  }
  bbox = aabb(min, max);
}

bool rotate_y::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  vec3 origin = r.origin();
  vec3 direction = r.direction();
  origin.e[0] = cos_theta*r.origin()[0] - sin_theta*r.origin()[2];
  origin.e[2] = sin_theta*r.origin()[0] + cos_theta*r.origin()[2];
  direction.e[0] = cos_theta*r.direction()[0] - sin_theta*r.direction()[2];
  direction.e[2] = sin_theta*r.direction()[0] + cos_theta*r.direction()[2];
  ray rotated_r(origin, direction, r.time());
  if(ptr->hit(rotated_r, t_min, t_max, rec, rng)) {
    vec3 p = rec.p;
    vec3 normal = rec.normal;
    p.e[0] = cos_theta*rec.p.e[0] + sin_theta*rec.p.e[2];
    p.e[2] = -sin_theta*rec.p.e[0] + cos_theta*rec.p.e[2]; 
    normal.e[0] = cos_theta*rec.normal.e[0] + sin_theta*rec.normal.e[2];
    normal.e[2] = -sin_theta*rec.normal.e[0] + cos_theta*rec.normal.e[2]; 
    rec.p = p;
    rec.normal = normal;
    if(rec.has_bump) {
      normal = rec.bump_normal;
      normal.e[0] = cos_theta*rec.bump_normal.e[0] + sin_theta*rec.bump_normal.e[2];
      normal.e[2] = -sin_theta*rec.bump_normal.e[0] + cos_theta*rec.bump_normal.e[2]; 
      rec.bump_normal = normal;
    }
    return(true);
  } else {
    return(false);
  }
}


rotate_x::rotate_x(std::shared_ptr<hitable> p, Float angle) : ptr(p) {
  Float radians = (M_PI / 180.0) * angle;
  sin_theta = sin(radians);
  cos_theta = cos(radians);
  hasbox = ptr->bounding_box(0,1,bbox);
  vec3 min(FLT_MAX, FLT_MAX, FLT_MAX);
  vec3 max(-FLT_MAX, -FLT_MAX, -FLT_MAX);
  for(int i = 0; i < 2; i++) {
    for(int j = 0; j < 2; j++) {
      for(int k = 0; k < 2; k++) {
        Float x = i*bbox.max().x() + (1-i)*bbox.min().x();
        Float y = j*bbox.max().y() + (1-j)*bbox.min().y();
        Float z = k*bbox.max().z() + (1-k)*bbox.min().z();
        Float newy = cos_theta*y + sin_theta*z;
        Float newz = -sin_theta*y + cos_theta*z;
        vec3 tester(x,newy,newz);
        for(int c = 0; c < 3; c++) {
          if(tester[c] > max[c]) {
            max.e[c] = tester[c];
          }
          if(tester[c] < min[c]) {
            min.e[c] = tester[c];
          }
        }
      } 
    }
  }
  bbox = aabb(min, max);
}

bool rotate_x::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  vec3 origin = r.origin();
  vec3 direction = r.direction();
  origin.e[1] = cos_theta*r.origin()[1] - sin_theta*r.origin()[2];
  origin.e[2] = sin_theta*r.origin()[1] + cos_theta*r.origin()[2];
  direction.e[1] = cos_theta*r.direction()[1] - sin_theta*r.direction()[2];
  direction.e[2] = sin_theta*r.direction()[1] + cos_theta*r.direction()[2];
  ray rotated_r(origin, direction, r.time());
  if(ptr->hit(rotated_r, t_min, t_max, rec, rng)) {
    vec3 p = rec.p;
    vec3 normal = rec.normal;
    p.e[1] = cos_theta*rec.p.e[1] + sin_theta*rec.p.e[2];
    p.e[2] = -sin_theta*rec.p.e[1] + cos_theta*rec.p.e[2]; 
    normal.e[1] = cos_theta*rec.normal.e[1] + sin_theta*rec.normal.e[2];
    normal.e[2] = -sin_theta*rec.normal.e[1] + cos_theta*rec.normal.e[2]; 
    rec.p = p;
    rec.normal = normal;
    if(rec.has_bump) {
      normal = rec.bump_normal;
      normal.e[1] = cos_theta*rec.bump_normal.e[1] + sin_theta*rec.bump_normal.e[2];
      normal.e[2] = -sin_theta*rec.bump_normal.e[1] + cos_theta*rec.bump_normal.e[2]; 
      rec.bump_normal = normal;
    }
    return(true);
  } else {
    return(false);
  }
}


rotate_z::rotate_z(std::shared_ptr<hitable> p, Float angle) : ptr(p) {
  Float radians = (M_PI / 180.0) * angle;
  sin_theta = sin(radians);
  cos_theta = cos(radians);
  hasbox = ptr->bounding_box(0,1,bbox);
  vec3 min(FLT_MAX, FLT_MAX, FLT_MAX);
  vec3 max(-FLT_MAX, -FLT_MAX, -FLT_MAX);
  for(int i = 0; i < 2; i++) {
    for(int j = 0; j < 2; j++) {
      for(int k = 0; k < 2; k++) {
        Float x = i*bbox.max().x() + (1-i)*bbox.min().x();
        Float y = j*bbox.max().y() + (1-j)*bbox.min().y();
        Float z = k*bbox.max().z() + (1-k)*bbox.min().z();
        Float newx = cos_theta*x + sin_theta*y;
        Float newy = -sin_theta*x + cos_theta*y;
        vec3 tester(newx,newy,z);
        for(int c = 0; c < 3; c++) {
          if(tester[c] > max[c]) {
            max.e[c] = tester[c];
          }
          if(tester[c] < min[c]) {
            min.e[c] = tester[c];
          }
        }
      } 
    }
  }
  bbox = aabb(min, max);
}

bool rotate_z::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  vec3 origin = r.origin();
  vec3 direction = r.direction();
  origin.e[0] = cos_theta*r.origin()[0] - sin_theta*r.origin()[1];
  origin.e[1] = sin_theta*r.origin()[0] + cos_theta*r.origin()[1];
  direction.e[0] = cos_theta*r.direction()[0] - sin_theta*r.direction()[1];
  direction.e[1] = sin_theta*r.direction()[0] + cos_theta*r.direction()[1];
  ray rotated_r(origin, direction, r.time());
  if(ptr->hit(rotated_r, t_min, t_max, rec, rng)) {
    vec3 p = rec.p;
    vec3 normal = rec.normal;
    p.e[0] = cos_theta*rec.p.e[0] + sin_theta*rec.p.e[1];
    p.e[1] = -sin_theta*rec.p.e[0] + cos_theta*rec.p.e[1]; 
    normal.e[0] = cos_theta*rec.normal.e[0] + sin_theta*rec.normal.e[1];
    normal.e[1] = -sin_theta*rec.normal.e[0] + cos_theta*rec.normal.e[1]; 
    rec.p = p;
    rec.normal = normal;
    if(rec.has_bump) {
      normal = rec.bump_normal;
      normal.e[0] = cos_theta*rec.bump_normal.e[0] + sin_theta*rec.bump_normal.e[1];
      normal.e[1] = -sin_theta*rec.bump_normal.e[0] + cos_theta*rec.bump_normal.e[1]; 
      rec.bump_normal = normal;
    }
    return(true);
  } else {
    return(false);
  }
}