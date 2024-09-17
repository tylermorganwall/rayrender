#include "aabb.h"
#include "raylog.h"
#include "simd.h"

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

const bool aabb::hit(const ray &r, Float tmin, Float tmax, random_gen& rng) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("AABB");

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

const bool aabb::hit(const ray &r, Float tmin, Float tmax, Sampler* sampler) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("AABB");
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

const point3f aabb::offset(const point3f p) const {
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

const point3f aabb::offset(const vec3f p) const {
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
  return((bounds[1] - bounds[0]));
}

int aabb::MaxDimension() const {
    point3f d = Diag();
    if (d.x() > d.y() && d.x() > d.z())
      return 0;
    else if (d.y() > d.z())
      return 1;
    else
      return 2;
}



void rayBBoxIntersect4(const ray& ray,
                       const BBox4& bbox4,
                       Float tMin,
                       Float tMax,
                       IVec4& hits,
                       FVec4& tMins,
                       FVec4& tMaxs) {
    SCOPED_CONTEXT("Hit");
    SCOPED_TIMER_COUNTER("BBox4");
    // Prepare ray data
    FVec4 rayOriginX = simd_set1(ray.origin().x());
    FVec4 rayOriginY = simd_set1(ray.origin().y());
    FVec4 rayOriginZ = simd_set1(ray.origin().z());
    FVec4 rayInvDirX = simd_set1(ray.inv_dir_pad.x());
    FVec4 rayInvDirY = simd_set1(ray.inv_dir_pad.y());
    FVec4 rayInvDirZ = simd_set1(ray.inv_dir_pad.z());

    // Determine near and far indices based on ray direction signs
    int ixNear = ray.sign[0] ? 3 : 0;
    int iyNear = ray.sign[1] ? 4 : 1;
    int izNear = ray.sign[2] ? 5 : 2;
    int ixFar = ray.sign[0] ? 0 : 3;
    int iyFar = ray.sign[1] ? 1 : 4;
    int izFar = ray.sign[2] ? 2 : 5;

    // Compute t values for slabs
    FVec4 t0x = simd_mul(simd_sub(bbox4.corners[ixNear], rayOriginX), rayInvDirX);
    FVec4 t1x = simd_mul(simd_sub(bbox4.corners[ixFar], rayOriginX), rayInvDirX);

    FVec4 t0y = simd_mul(simd_sub(bbox4.corners[iyNear], rayOriginY), rayInvDirY);
    FVec4 t1y = simd_mul(simd_sub(bbox4.corners[iyFar], rayOriginY), rayInvDirY);

    FVec4 t0z = simd_mul(simd_sub(bbox4.corners[izNear], rayOriginZ), rayInvDirZ);
    FVec4 t1z = simd_mul(simd_sub(bbox4.corners[izFar], rayOriginZ), rayInvDirZ);

    // Compute tEnter and tExit
    FVec4 tEnter = simd_max(simd_max(t0x, t0y), t0z);
    FVec4 tExit = simd_min(simd_min(t1x, t1y), t1z);

    // Clamp tEnter and tExit with tMin and tMax
    tMins = simd_max(tEnter, simd_set1(tMin));
    tMaxs = simd_min(tExit, simd_set1(tMax));

    // Compute hit mask
    SimdMask hitMask = simd_less_equal(tMins, tMaxs);
    hits = simd_cast_to_int(hitMask);
}

void rayBBoxIntersect4Serial(const ray& ray,
                       const BBox4& bbox4,
                       Float tMin,
                       Float tMax,
                       IVec4& hits,
                       FVec4& tMins,
                       FVec4& tMaxs) {
    // Loop over each bounding box (up to 4)
    for (int i = 0; i < 4; ++i) {
        // Extract min and max coordinates for the i-th bounding box
        Float minX = bbox4.corners[0][i]; // bbox4.corners[0]: minX
        Float minY = bbox4.corners[1][i]; // bbox4.corners[1]: minY
        Float minZ = bbox4.corners[2][i]; // bbox4.corners[2]: minZ
        Float maxX = bbox4.corners[3][i]; // bbox4.corners[3]: maxX
        Float maxY = bbox4.corners[4][i]; // bbox4.corners[4]: maxY
        Float maxZ = bbox4.corners[5][i]; // bbox4.corners[5]: maxZ

        // Prepare ray data
        Float rayOriginX = ray.origin().x();
        Float rayOriginY = ray.origin().y();
        Float rayOriginZ = ray.origin().z();
        Float rayInvDirX = ray.inv_dir_pad.x();
        Float rayInvDirY = ray.inv_dir_pad.y();
        Float rayInvDirZ = ray.inv_dir_pad.z();

        // Determine near and far coordinates based on ray direction signs
        Float bboxNearX = ray.sign[0] ? maxX : minX;
        Float bboxFarX  = ray.sign[0] ? minX : maxX;
        Float bboxNearY = ray.sign[1] ? maxY : minY;
        Float bboxFarY  = ray.sign[1] ? minY : maxY;
        Float bboxNearZ = ray.sign[2] ? maxZ : minZ;
        Float bboxFarZ  = ray.sign[2] ? minZ : maxZ;

        // Compute t values for slabs
        Float t0x = (bboxNearX - rayOriginX) * rayInvDirX;
        Float t1x = (bboxFarX  - rayOriginX) * rayInvDirX;

        Float t0y = (bboxNearY - rayOriginY) * rayInvDirY;
        Float t1y = (bboxFarY  - rayOriginY) * rayInvDirY;

        Float t0z = (bboxNearZ - rayOriginZ) * rayInvDirZ;
        Float t1z = (bboxFarZ  - rayOriginZ) * rayInvDirZ;

        // Compute tEnter and tExit
        Float tEnter = std::max(std::max(t0x, t0y), t0z);
        Float tExit  = std::min(std::min(t1x, t1y), t1z);

        // Clamp tEnter and tExit with tMin and tMax
        Float tMin_i = std::max(tEnter, tMin);
        Float tMax_i = std::min(tExit, tMax);

        // Store tMins and tMaxs
        tMins[i] = tMin_i;
        tMaxs[i] = tMax_i;

        // Compute hit mask
        hits[i] = (tMin_i <= tMax_i) ? -1 : 0; // -1 indicates a hit (all bits set in integer)
    }
}
