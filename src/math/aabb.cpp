#include "../math/aabb.h"
#include "../utils/raylog.h"
#include "../math/simd.h"
#include <cmath>

const Float aabb::surface_area() const {
  point3f diag = Diag();
  return(min().xyz.x <= max().xyz.x ? 
         2*(diag.xyz.x * diag.xyz.y + diag.xyz.x * diag.xyz.z + diag.xyz.y * diag.xyz.z):
         10E20);
}

const Float aabb::Volume() const {
  point3f diag = Diag();
  return(min().xyz.x <= max().xyz.x ? 
         diag.xyz.x * diag.xyz.y * diag.xyz.z:
         10E20);
}

namespace {

inline bool hit_aabb_branchless(const aabb &box, const Ray &r, Float tmin,
                                Float tmax) {
  const int sx = r.inv_dir_is_neg[0] ? 1 : 0;
  const int sy = r.inv_dir_is_neg[1] ? 1 : 0;
  const int sz = r.inv_dir_is_neg[2] ? 1 : 0;

  const Float txNear =
      (box.bounds[sx].e[0] - r.o.e[0]) * r.inv_dir_pad.e[0];
  const Float txFar =
      (box.bounds[1 - sx].e[0] - r.o.e[0]) * r.inv_dir_pad.e[0];

  const Float tyNear =
      (box.bounds[sy].e[1] - r.o.e[1]) * r.inv_dir_pad.e[1];
  const Float tyFar =
      (box.bounds[1 - sy].e[1] - r.o.e[1]) * r.inv_dir_pad.e[1];

  const Float tzNear =
      (box.bounds[sz].e[2] - r.o.e[2]) * r.inv_dir_pad.e[2];
  const Float tzFar =
      (box.bounds[1 - sz].e[2] - r.o.e[2]) * r.inv_dir_pad.e[2];

  tmin = std::fmax(tmin, txNear);
  tmax = std::fmin(tmax, txFar);

  tmin = std::fmax(tmin, tyNear);
  tmax = std::fmin(tmax, tyFar);

  tmin = std::fmax(tmin, tzNear);
  tmax = std::fmin(tmax, tzFar);

  return tmin <= tmax;
}

} // namespace

const bool aabb::hit(const Ray &r, Float tmin, Float tmax,
                     random_gen &rng) const {
  return hit_aabb_branchless(*this, r, tmin, tmax);
}

const bool aabb::hit(const Ray &r, Float tmin, Float tmax,
                     Sampler *sampler) const {
  return hit_aabb_branchless(*this, r, tmin, tmax);
}

const point3f aabb::offset(const point3f p) const {
  point3f o = p + -min();
  if (max().xyz.x > min().xyz.x) {
    o.e[0] /= (max().xyz.x - min().xyz.x);
  }
  if (max().xyz.y > min().xyz.y) {
    o.e[1] /= (max().xyz.y - min().xyz.y);
  }
  if (max().xyz.z > min().xyz.z) {
    o.e[2] /= (max().xyz.z - min().xyz.z);
  }
  return(o);
}

const point3f aabb::offset(const vec3f p) const {
  point3f o = point3f(p.xyz.x,p.xyz.y,p.xyz.z) + -min();
  if (max().xyz.x > min().xyz.x) {
    o.e[0] /= (max().xyz.x - min().xyz.x);
  }
  if (max().xyz.y > min().xyz.y) {
    o.e[1] /= (max().xyz.y - min().xyz.y);
  }
  if (max().xyz.z > min().xyz.z) {
    o.e[2] /= (max().xyz.z - min().xyz.z);
  }
  return(o);
}

const point3f aabb::Corner(int corner) const {
  return point3f((*this).bounds[(corner & 1)].xyz.x,
                 (*this).bounds[(corner & 2) ? 1 : 0].xyz.y,
                 (*this).bounds[(corner & 4) ? 1 : 0].xyz.z);
}

const point3f aabb::Lerp(const point3f &t) const {
  return point3f(lerp(t.xyz.x, min().xyz.x, max().xyz.x),
                 lerp(t.xyz.y, min().xyz.y, max().xyz.y),
                 lerp(t.xyz.z, min().xyz.z, max().xyz.z));
}

const point2f aabb::Lerp(const point2f &t) const {
  return point2f(lerp(t.xy.x, min().xyz.x, max().xyz.x),
                 lerp(t.xy.y, min().xyz.y, max().xyz.y));
}

const point3f aabb::Centroid() const {
  return((bounds[0] + bounds[1]) / static_cast<Float>(2));
}

const point3f aabb::Diag() const {
  return((bounds[1] + (-bounds[0])));
}

int aabb::MaxDimension() const {
    point3f d = Diag();
    if (d.xyz.x > d.xyz.y && d.xyz.x > d.xyz.z)
      return 0;
    else if (d.xyz.y > d.xyz.z)
      return 1;
    else
      return 2;
}

void rayBBoxIntersect4Serial(const Ray& r,
                       const BBox4& bbox4,
                       Float tMin,
                       Float tMax,
                       IVec4& hits,
                       FVec4& tEnters) {
}
