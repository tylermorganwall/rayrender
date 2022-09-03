#include "aabb.h"

const Float aabb::surface_area() const {
  point3f diag = Diag();
  return(min().x() <= max().x() ? 
         2*(diag.x() * diag.y() + diag.x() * diag.z() + diag.y()*diag.z()):
         10E20);
}

const Float aabb::Volume() const {
  point3f diag = Diag();
  return(min().x() <= max().x() ? 
         diag.x() * diag.y() * diag.z():
         10E20);
}

bool aabb::hit(const ray &r, Float tmin, Float tmax, random_gen& rng) {
  Float txmin, txmax, tymin, tymax, tzmin, tzmax;
  txmin = (bounds[  r.sign[0]].x()-r.origin().x()) * r.inv_dir.x();
  txmax = (bounds[1-r.sign[0]].x()-r.origin().x()) * r.inv_dir_pad.x();
  tymin = (bounds[  r.sign[1]].y()-r.origin().y()) * r.inv_dir.y();
  tymax = (bounds[1-r.sign[1]].y()-r.origin().y()) * r.inv_dir_pad.y();
  tzmin = (bounds[  r.sign[2]].z()-r.origin().z()) * r.inv_dir.z();
  tzmax = (bounds[1-r.sign[2]].z()-r.origin().z()) * r.inv_dir_pad.z();
  tmin = ffmax(tzmin, ffmax(tymin, ffmax(txmin, tmin)));
  tmax = ffmin(tzmax, ffmin(tymax, ffmin(txmax, tmax)));
  return(tmin <= tmax);
}

bool aabb::hit(const ray &r, Float tmin, Float tmax, Sampler* sampler) {
  Float txmin, txmax, tymin, tymax, tzmin, tzmax;
  txmin = (bounds[  r.sign[0]].x()-r.origin().x()) * r.inv_dir.x();
  txmax = (bounds[1-r.sign[0]].x()-r.origin().x()) * r.inv_dir_pad.x();
  tymin = (bounds[  r.sign[1]].y()-r.origin().y()) * r.inv_dir.y();
  tymax = (bounds[1-r.sign[1]].y()-r.origin().y()) * r.inv_dir_pad.y();
  tzmin = (bounds[  r.sign[2]].z()-r.origin().z()) * r.inv_dir.z();
  tzmax = (bounds[1-r.sign[2]].z()-r.origin().z()) * r.inv_dir_pad.z();
  tmin = ffmax(tzmin, ffmax(tymin, ffmax(txmin, tmin)));
  tmax = ffmin(tzmax, ffmin(tymax, ffmin(txmax, tmax)));
  return(tmin <= tmax);
}

const point3f aabb::offset(const point3f p) {
  point3f o = p + -min();
  if (max().x() > min().x()) {
    o.e[0] /= (max().x() - min().x());
  }
  if (max().y() > min().y()) {
    o.e[1] /= (max().y() - min().y());
  }
  if (max().z() > min().z()) {
    o.e[2] /= (max().z() - min().z());
  }
  return(o);
}

const point3f aabb::offset(const vec3f p) {
  point3f o = point3f(p.x(),p.y(),p.z()) + -min();
  if (max().x() > min().x()) {
    o.e[0] /= (max().x() - min().x());
  }
  if (max().y() > min().y()) {
    o.e[1] /= (max().y() - min().y());
  }
  if (max().z() > min().z()) {
    o.e[2] /= (max().z() - min().z());
  }
  return(o);
}

const point3f aabb::Corner(int corner) const {
  return point3f((*this).bounds[(corner & 1)].x(),
                 (*this).bounds[(corner & 2) ? 1 : 0].y(),
                 (*this).bounds[(corner & 4) ? 1 : 0].z());
}

const point3f aabb::Lerp(const point3f &t) const {
  return point3f(lerp(t.x(), min().x(), max().x()),
                 lerp(t.y(), min().y(), max().y()),
                 lerp(t.z(), min().z(), max().z()));
}

const point2f aabb::Lerp(const point2f &t) const {
  return point2f(lerp(t.x(), min().x(), max().x()),
                 lerp(t.y(), min().y(), max().y()));
}

const point3f aabb::Centroid() const {
  return((bounds[0] + bounds[1])/2);
}

const point3f aabb::Diag() const {
  return((bounds[1] + bounds[0]));
}


