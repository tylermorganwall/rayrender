#include "hitable.h"

//Translate implementation

void get_sphere_uv(const vec3f& p, Float& u, Float& v) {
  Float phi = atan2(p.z(),p.x());
  Float theta = asin(p.y());
  u = 1 - (phi + M_PI) / (2*M_PI);
  v = (theta + M_PI/2) / M_PI;
}

void get_sphere_uv(const normal3f& p, Float& u, Float& v) {
  Float phi = atan2(p.z(),p.x());
  Float theta = asin(p.y());
  u = 1 - (phi + M_PI) / (2*M_PI);
  v = (theta + M_PI/2) / M_PI;
}



bool AnimatedHitable::bounding_box(Float t0, Float t1, aabb& box) const {
  primitive->bounding_box(t0, t1, box);
  box = PrimitiveToWorld.MotionBounds(box);
  return(true);
}

bool AnimatedHitable::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  Transform InterpolatedPrimToWorld;
  PrimitiveToWorld.Interpolate(r.time(), &InterpolatedPrimToWorld);
  ray ray_interp = Inverse(InterpolatedPrimToWorld)(r);
  
  if (!primitive->hit(ray_interp, t_min, t_max, rec, rng)) {
    return false;
  }
  r.tMax = ray_interp.tMax;
  if (!InterpolatedPrimToWorld.IsIdentity()) {
    rec = InterpolatedPrimToWorld(rec);
  }
  return true;
}

bool AnimatedHitable::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
  Transform InterpolatedPrimToWorld;
  PrimitiveToWorld.Interpolate(r.time(), &InterpolatedPrimToWorld);
  ray ray_interp = Inverse(InterpolatedPrimToWorld)(r);
  
  if (!primitive->hit(ray_interp, t_min, t_max, rec, sampler)) {
    return false;
  }
  r.tMax = ray_interp.tMax;
  if (!InterpolatedPrimToWorld.IsIdentity()) {
    rec = InterpolatedPrimToWorld(rec);
  }
  return true;
}