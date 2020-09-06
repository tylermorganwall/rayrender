#include "csg.h"

bool csg::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  constexpr Float maxDistance = 100; 
  constexpr Float threshold = 10e-6; 
  constexpr float delta = 10e-5; 
  Float t = 0; 
  uint32_t numSteps = 0; 
  csg_sphere temp(center, radius);
  // const ImplicitShape *isectShape = nullptr; 
  
  while (t < maxDistance) {
    Float minDistance = INFINITY; 
    vec3 dir = unit_vector(r.direction());
    vec3 from = r.origin() + t * dir; 
    // for (auto shape : scene) { 
    float d = temp.getDistance(from); 
    if (d < minDistance) {
      minDistance = d;
      // isectShape = shape;
    }
    // } 
    if (minDistance <= threshold * t) { 
      rec.p = r.origin() + t * dir;
      rec.normal = vec3( 
        temp.getDistance(rec.p + vec3(delta, 0, 0)) - temp.getDistance(rec.p + vec3(-delta, 0, 0)), 
        temp.getDistance(rec.p + vec3(0, delta, 0)) - temp.getDistance(rec.p + vec3(0, -delta, 0)), 
        temp.getDistance(rec.p + vec3(0, 0, delta)) - temp.getDistance(rec.p + vec3(0, 0, -delta))
      ); 
      rec.t = t;
      rec.normal.make_unit_vector(); 
      rec.u = 0.5;
      rec.v = 0.5;
      rec.dpdu = 0.5;
      rec.dpdv = 0.5;
      rec.mat_ptr = mat_ptr;
      rec.has_bump = false;
      rec.bump_normal =  rec.normal;
      return(true);
    } 
    t += minDistance; 
    numSteps++; 
  } 
  return(false);
}

bool csg::bounding_box(Float t0, Float t1, aabb& box) const {
  box = aabb(center - vec3(radius,radius,radius), center + vec3(radius,radius,radius));
  return(true);
}