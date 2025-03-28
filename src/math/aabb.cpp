#include "../math/aabb.h"
#include "../utils/raylog.h"
#include "../math/simd.h"

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

const bool aabb::hit(const Ray &r, Float tmin, Float tmax, random_gen& rng) const {
  // SCOPED_CONTEXT("Hit");
  // SCOPED_TIMER_COUNTER("AABB");

  // Float txmin, txmax, tymin, tymax, tzmin, tzmax;
  // txmin = (bounds[  r.sign[0]].xyz.x-r.origin().xyz.x) * r.inv_dir.xyz.x;
  // txmax = (bounds[1-r.sign[0]].xyz.x-r.origin().xyz.x) * r.inv_dir_pad.xyz.x;
  // tymin = (bounds[  r.sign[1]].xyz.y-r.origin().xyz.y) * r.inv_dir.xyz.y;
  // tymax = (bounds[1-r.sign[1]].xyz.y-r.origin().xyz.y) * r.inv_dir_pad.xyz.y;
  // tzmin = (bounds[  r.sign[2]].xyz.z-r.origin().xyz.z) * r.inv_dir.xyz.z;
  // tzmax = (bounds[1-r.sign[2]].xyz.z-r.origin().xyz.z) * r.inv_dir_pad.xyz.z;
  // tmin = ffmax(tzmin, ffmax(tymin, ffmax(txmin, tmin)));
  // tmax = ffmin(tzmax, ffmin(tymax, ffmin(txmax, tmax)));
  return(tmin <= tmax);
}

const bool aabb::hit(const Ray &r, Float tmin, Float tmax, Sampler* sampler) const {
  // SCOPED_CONTEXT("Hit");
  // SCOPED_TIMER_COUNTER("AABB");
  // Float txmin, txmax, tymin, tymax, tzmin, tzmax;
  // txmin = (bounds[  r.sign[0]].xyz.x-r.origin().xyz.x) * r.inv_dir.xyz.x;
  // txmax = (bounds[1-r.sign[0]].xyz.x-r.origin().xyz.x) * r.inv_dir_pad.xyz.x;
  // tymin = (bounds[  r.sign[1]].xyz.y-r.origin().xyz.y) * r.inv_dir.xyz.y;
  // tymax = (bounds[1-r.sign[1]].xyz.y-r.origin().xyz.y) * r.inv_dir_pad.xyz.y;
  // tzmin = (bounds[  r.sign[2]].xyz.z-r.origin().xyz.z) * r.inv_dir.xyz.z;
  // tzmax = (bounds[1-r.sign[2]].xyz.z-r.origin().xyz.z) * r.inv_dir_pad.xyz.z;
  // tmin = ffmax(tzmin, ffmax(tymin, ffmax(txmin, tmin)));
  // tmax = ffmin(tzmax, ffmin(tymax, ffmin(txmax, tmax)));
  return(tmin <= tmax);
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



// void rayBBoxIntersect4(const ray& ray,
//                        const BBox4& bbox4,
//                        Float tMin,
//                        Float tMax,
//                        IVec4& hits,
//                        FVec4& tMins,
//                        FVec4& tMaxs) {
//     // SCOPED_CONTEXT("Hit");
//     // SCOPED_TIMER_COUNTER("BBox4");
//     // Prepare ray data
//     FVec4 rayOriginX = simd_set1(ray.origin().xyz.x);
//     FVec4 rayOriginY = simd_set1(ray.origin().xyz.y);
//     FVec4 rayOriginZ = simd_set1(ray.origin().xyz.z);
//     FVec4 rayInvDirX = simd_set1(ray.inv_dir_pad.xyz.x);
//     FVec4 rayInvDirY = simd_set1(ray.inv_dir_pad.xyz.y);
//     FVec4 rayInvDirZ = simd_set1(ray.inv_dir_pad.xyz.z);

//     // Determine near and far indices based on ray direction signs
//     int ixNear = ray.sign[0] ? 3 : 0;
//     int iyNear = ray.sign[1] ? 4 : 1;
//     int izNear = ray.sign[2] ? 5 : 2;
//     int ixFar = ray.sign[0] ? 0 : 3;
//     int iyFar = ray.sign[1] ? 1 : 4;
//     int izFar = ray.sign[2] ? 2 : 5;

//     // Compute t values for slabs
//     FVec4 t0x = simd_mul(simd_sub(bbox4.corners[ixNear], rayOriginX), rayInvDirX);
//     FVec4 t1x = simd_mul(simd_sub(bbox4.corners[ixFar], rayOriginX), rayInvDirX);

//     FVec4 t0y = simd_mul(simd_sub(bbox4.corners[iyNear], rayOriginY), rayInvDirY);
//     FVec4 t1y = simd_mul(simd_sub(bbox4.corners[iyFar], rayOriginY), rayInvDirY);

//     FVec4 t0z = simd_mul(simd_sub(bbox4.corners[izNear], rayOriginZ), rayInvDirZ);
//     FVec4 t1z = simd_mul(simd_sub(bbox4.corners[izFar], rayOriginZ), rayInvDirZ);

//     // Compute tEnter and tExit
//     FVec4 tEnter = simd_max(simd_max(t0x, t0y), t0z);
//     FVec4 tExit = simd_min(simd_min(t1x, t1y), t1z);

//     // Clamp tEnter and tExit with tMin and tMax
//     tMins = simd_max(tEnter, simd_set1(tMin));
//     tMaxs = simd_min(tExit, simd_set1(tMax));

//     // Compute hit mask
//     SimdMask hitMask = simd_less_equal(tMins, tMaxs);
//     hits = simd_cast_to_int(hitMask);
// }

// inline void rayBBoxIntersect4(const ray& r,
//                        const BBox4& bbox4,
//                        Float tMin,
//                        Float tMax,
//                        IVec4& hits,
//                        FVec4& tEnters) {
//                        //FVec4& tExits) {
//     FVec4 bboxMinX = bbox4.getMinX();
//     FVec4 bboxMaxX = bbox4.getMaxX();
//     FVec4 bboxMinY = bbox4.getMinY();
//     FVec4 bboxMaxY = bbox4.getMaxY();
//     FVec4 bboxMinZ = bbox4.getMinZ();
//     FVec4 bboxMaxZ = bbox4.getMaxZ();
//     FVec4 origin0 = r.origin4[0];
//     FVec4 origin1 = r.origin4[1];
//     FVec4 origin2 = r.origin4[2];
//     FVec4 inv_pad0 = r.inv_dir_pad4[0];
//     FVec4 inv_pad1 = r.inv_dir_pad4[1];
//     FVec4 inv_pad2 = r.inv_dir_pad4[2];

//     // For X axis
//     FVec4 t0x = simd_mul(simd_sub(bboxMinX, origin0), inv_pad0);
//     FVec4 t1x = simd_mul(simd_sub(bboxMaxX, origin0), inv_pad0);

//     // For Y axis
//     FVec4 t0y = simd_mul(simd_sub(bboxMinY, origin1), inv_pad1);
//     FVec4 t1y = simd_mul(simd_sub(bboxMaxY, origin1), inv_pad1);

//     // For Z axis
//     FVec4 t0z = simd_mul(simd_sub(bboxMinZ, origin2), inv_pad2);
//     FVec4 t1z = simd_mul(simd_sub(bboxMaxZ, origin2), inv_pad2);

//     // Compute tEnter and tExit
//     tEnters = simd_max(simd_max(simd_min(t0x, t1x), simd_min(t0y, t1y)), simd_min(t0z, t1z));
//     FVec4 tExits  = simd_min(simd_min(simd_max(t0x, t1x), simd_max(t0y, t1y)), simd_max(t0z, t1z));

//     // Compute hit mask
//     FVec4 tmp_max = simd_max(tEnters, simd_set1(tMin));
//     FVec4 tmp_min = simd_min(tExits,  simd_set1(tMax));

//     SimdMask hitMask = simd_less_equal(tmp_max,tmp_min);

//     hits = simd_cast_to_int(hitMask);
//     // int hits = simd_extract_hitmask(simd_cast_to_int(hitMask));
//     // return(hits);
// }

void rayBBoxIntersect4Serial(const Ray& r,
                       const BBox4& bbox4,
                       Float tMin,
                       Float tMax,
                       IVec4& hits,
                       FVec4& tEnters) {
    // FVec4 bboxMinX = bbox4.getMinX();
    // FVec4 bboxMaxX = bbox4.getMaxX();
    // FVec4 bboxMinY = bbox4.getMinY();
    // FVec4 bboxMaxY = bbox4.getMaxY();
    // FVec4 bboxMinZ = bbox4.getMinZ();
    // FVec4 bboxMaxZ = bbox4.getMaxZ();
    
    // // Prepare ray data
    // Float rayOriginX = r.origin().xyz.x;
    // Float rayOriginY = r.origin().xyz.y;
    // Float rayOriginZ = r.origin().xyz.z;
    // Float rayInvDirX = r.inv_dir_pad.xyz.x;
    // Float rayInvDirY = r.inv_dir_pad.xyz.y;
    // Float rayInvDirZ = r.inv_dir_pad.xyz.z;
    // // Loop over each bounding box (up to 4)
    // for (int i = 0; i < 4; ++i) {
    //     // Extract min and max coordinates for the i-th bounding box
    //     Float minX = bboxMinX[i]; // bbox4.corners[0]: minX
    //     Float minY = bboxMinY[i]; // bbox4.corners[1]: minY
    //     Float minZ = bboxMinZ[i]; // bbox4.corners[2]: minZ
    //     Float maxX = bboxMaxX[i]; // bbox4.corners[3]: maxX
    //     Float maxY = bboxMaxY[i]; // bbox4.corners[4]: maxY
    //     Float maxZ = bboxMaxZ[i]; // bbox4.corners[5]: maxZ

    //     // Determine near and far coordinates based on ray direction signs
    //     Float bboxNearX = r.sign[0] ? maxX : minX;
    //     Float bboxFarX  = r.sign[0] ? minX : maxX;
    //     Float bboxNearY = r.sign[1] ? maxY : minY;
    //     Float bboxFarY  = r.sign[1] ? minY : maxY;
    //     Float bboxNearZ = r.sign[2] ? maxZ : minZ;
    //     Float bboxFarZ  = r.sign[2] ? minZ : maxZ;

    //     // Compute t values for slabs
    //     Float t0x = (bboxNearX - rayOriginX) * rayInvDirX;
    //     Float t1x = (bboxFarX  - rayOriginX) * rayInvDirX;

    //     Float t0y = (bboxNearY - rayOriginY) * rayInvDirY;
    //     Float t1y = (bboxFarY  - rayOriginY) * rayInvDirY;

    //     Float t0z = (bboxNearZ - rayOriginZ) * rayInvDirZ;
    //     Float t1z = (bboxFarZ  - rayOriginZ) * rayInvDirZ;

    //     // Compute tEnter and tExit
    //     Float tEnter = ffmax(ffmax(t0x, t0y), t0z);
    //     Float tExit  = ffmin(ffmin(t1x, t1y), t1z);

    //     // Clamp tEnter and tExit with tMin and tMax
    //     Float tMin_i = ffmax(tEnter, tMin);
    //     Float tMax_i = ffmin(tExit, tMax);

    //     // Store tMins and tMaxs
    //     tEnters[i] = tMin_i;

    //     // Compute hit mask
    //     hits[i] = (tMin_i <= tMax_i) ? -1 : 0; // -1 indicates a hit (all bits set in integer)
    // }
}
