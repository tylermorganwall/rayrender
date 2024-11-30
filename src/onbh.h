#ifndef ONBH
#define ONBH

#include "vectypes.h"

class onb {
public:
  onb() {}
  onb(vec3f u, vec3f v, vec3f w) {
    axis[0] = u;
    axis[1] = v;
    axis[2] = w;
  }
  onb(vec3f u, vec3f v, normal3f w) {
    axis[0] = u;
    axis[1] = v;
    axis[2] = convert_to_vec3(w);
  }
  inline vec3f operator[](int i) const {return(axis[i]);}
  vec3f u() const {return(axis[0]);}
  vec3f v() const {return(axis[1]);}
  vec3f w() const {return(axis[2]);}
  vec3f local(Float a, Float b, Float c) const {
    return(a*u() + b*v() + c*w());
  }
  vec3f local(const vec3f& a) const {
    return(a.x()*u() + a.y()*v() + a.z()*w());
  }
  vec3f local_to_world(const vec3f& a) const {
    return(a.x()*u() + a.y()*v() + a.z() * w());
  }
  vec3f world_to_local(const vec3f& a) const {
    return(vec3f(dot(a,u()),dot(a,v()),dot(a,w())));
  }
  void build_from_w(const vec3f& n) {
    axis[2] = unit_vector(n);
    vec3f a = std::fabs(w().x()) > 0.9999999 ? vec3f(0,1,0)  : vec3f(1,0,0);
    axis[1] = unit_vector(cross(w(),a));
    axis[0] = cross(w(), v());
  }
  void build_from_w(const normal3f& n) {
    axis[2] = unit_vector(convert_to_vec3(n));
    vec3f a = std::fabs(w().x()) > 0.9999999 ? vec3f(0,1,0)  : vec3f(1,0,0);
    axis[1] = unit_vector(cross(w(),a));
    axis[0] = cross(w(), v());
  }
  void build_from_w_normalized(const vec3f& n) {
    axis[2] = (n);
    vec3f a = std::abs(w().x()) > 0.9999999 ? vec3f(0,1,0)  : vec3f(1,0,0);
    axis[1] = unit_vector(cross(w(),a));
    axis[0] = cross(w(), v());
  }
  void build_from_w_normalized(const normal3f& n) {
    axis[2] = (convert_to_vec3(n));
    vec3f a = std::abs(w().x()) > 0.9999999 ? vec3f(0,1,0)  : vec3f(1,0,0);
    axis[1] = unit_vector(cross(w(),a));
    axis[0] = cross(w(), v());
  }

  void swap_yz() {
    vec3f tempaxis = axis[1];
    axis[1] = axis[2];
    axis[2] = tempaxis;
  }
  vec3f axis[3];
};


#endif
