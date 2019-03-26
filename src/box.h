#ifndef BOXH
#define BOXH

#include "hitablelist.h"
#include "xyrect.h"

class box : public hitable {
public:
  box() {}
  box(const vec3& p0, const vec3& p1, material *ptr);
  virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const;
  virtual bool bounding_box(float t0, float t1, aabb& box) const {
    box = aabb(pmin, pmax);
    return(true);
  }
  virtual float pdf_value(const vec3& o, const vec3& v) const {
    return(list_ptr.pdf_value(o,v));
  }
  virtual vec3 random(const vec3& o) const {
    return(list_ptr.random(o));
  }
  vec3 pmin, pmax;
  hitable_list list_ptr;
};

box::box(const vec3& p0, const vec3& p1, material *ptr) {
  pmin = p0;
  pmax = p1;
  hitable **list = new hitable*[6];
  list[0] = new xy_rect(p0.x(), p1.x(), p0.y(), p1.y(), p1.z(),ptr);
  list[1] = new flip_normals(new xy_rect(p0.x(), p1.x(), p0.y(), p1.y(), p0.z(), ptr));
  list[2] = new xz_rect(p0.x(), p1.x(), p0.z(), p1.z(), p1.y(), ptr);
  list[3] = new flip_normals(new xz_rect(p0.x(), p1.x(), p0.z(), p1.z(), p0.y(), ptr));
  list[4] = new yz_rect(p0.y(), p1.y(), p0.z(), p1.z(), p1.x(), ptr);
  list[5] = new flip_normals(new yz_rect(p0.y(), p1.y(), p0.z(), p1.z(), p0.x(), ptr));
  list_ptr = hitable_list(list,6);
}

bool box::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
  return(list_ptr.hit(r,t_min,t_max,rec));
}

#endif