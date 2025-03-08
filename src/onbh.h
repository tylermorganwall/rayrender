#ifndef ONBH
#define ONBH

#include "vectypes.h"
#include <string>
#include <cstdio>
#include <cmath>

class onb {
public:
  onb() {
    axis[0] = vec3f(1, 0, 0);
    axis[1] = vec3f(0, 1, 0);
    axis[2] = vec3f(0, 0, 1);
  }

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

  inline vec3f operator[](int i) const { return axis[i]; }
  vec3f u() const { return axis[0]; }
  vec3f v() const { return axis[1]; }
  vec3f w() const { return axis[2]; }

  vec3f local(Float a, Float b, Float c) const {
    return a*u() + b*v() + c*w();
  }

  vec3f local(const vec3f &a) const {
    return a.x*u() + a.y*v() + a.z*w();
  }

  vec3f local_to_world(const vec3f &a) const {
    return a.x*u() + a.y*v() + a.z*w();
  }
  
  vec3f world_to_local(const vec3f &a) const {
    return vec3f(dot(a, u()), dot(a, v()), dot(a, w()));
  }
  // Aliases for PBRT's naming.
  vec3f ToLocal(const vec3f &a) const { return world_to_local(a); }
  vec3f FromLocal(const vec3f &a) const { return local_to_world(a); }

  // Builds the basis from a given vector (w) and ensures right-handedness.
  void build_from_w(const vec3f &n) {
    axis[2] = unit_vector(n);
    vec3f a = (std::fabs(axis[2].x) > 0.9999999f) ? vec3f(0, 1, 0) : vec3f(1, 0, 0);
    axis[1] = unit_vector(cross(axis[2], a));
    axis[0] = cross(axis[1], axis[2]);
  }

  void build_from_w(const normal3f &n) {
    axis[2] = unit_vector(convert_to_vec3(n));
    vec3f a = (std::fabs(axis[2].x) > 0.9999999f) ? vec3f(0, 1, 0) : vec3f(1, 0, 0);
    axis[1] = unit_vector(cross(axis[2], a));
    axis[0] = cross(axis[1], axis[2]);
  }

  void build_from_w_normalized(const vec3f &n) {
    axis[2] = n;
    vec3f a = (std::fabs(axis[2].x) > 0.9999999f) ? vec3f(0, 1, 0) : vec3f(1, 0, 0);
    axis[1] = unit_vector(cross(axis[2], a));
    axis[0] = cross(axis[1], axis[2]);
  }

  void build_from_w_normalized(const normal3f &n) {
    axis[2] = unit_vector(convert_to_vec3(n));
    vec3f a = (std::fabs(axis[2].x) > 0.9999999f) ? vec3f(0, 1, 0) : vec3f(1, 0, 0);
    axis[1] = unit_vector(cross(axis[2], a));
    axis[0] = cross(axis[1], axis[2]);
  }

  void swap_yz() {
    vec3f tempaxis = axis[1];
    axis[1] = axis[2];
    axis[2] = tempaxis;
  }

  static onb FromXZ(const vec3f &x, const vec3f &z) {
    return onb(x, cross(z, x), z);
  }

  static onb FromXY(const vec3f &x, const vec3f &y) {
    return onb(x, y, cross(x, y));
  }

  static void coordinateSystem(const vec3f &v1, vec3f *v2, vec3f *v3) {
    Float sign = (v1.z >= 0) ? 1.0f : -1.0f;
    Float a = -1.0f / (sign + v1.z);
    Float b = v1.x * v1.y * a;
    *v2 = vec3f(1 + sign * (v1.x * v1.x) * a, sign * b, -sign * v1.x);
    *v3 = vec3f(b, sign + (v1.y * v1.y) * a, -v1.y);
  }

  static onb FromZ(const vec3f &z_) {
    vec3f x_, y_;
    coordinateSystem(z_, &x_, &y_);
    return onb(x_, y_, z_);
  }

  static onb FromX(const vec3f &x_) {
    vec3f y_, z_;
    coordinateSystem(x_, &y_, &z_);
    return onb(x_, y_, z_);
  }

  static onb FromY(const vec3f &y_) {
    vec3f z_, x_;
    coordinateSystem(y_, &z_, &x_);
    return onb(x_, y_, z_);
  }

  // Overloads for normal3f.
  static onb FromZ(const normal3f &nz) { return FromZ(convert_to_vec3(nz)); }
  static onb FromX(const normal3f &nx) { return FromX(convert_to_vec3(nx)); }
  static onb FromY(const normal3f &ny) { return FromY(convert_to_vec3(ny)); }

  union {
    vec3f axis[3];
    vec3f x, y, z;
  };
};

#endif
