#include "curve.h"

static point3f BlossomBezier(const point3f p[4], Float u0, Float u1, Float u2) {
  point3f a[3] = { lerp(u0, p[0], p[1]),
                lerp(u0, p[1], p[2]),
                lerp(u0, p[2], p[3]) };
  point3f b[2] = { lerp(u1, a[0], a[1]), lerp(u1, a[1], a[2]) };
  return(lerp(u2, b[0], b[1]));
}


static point3f EvalBezier(const point3f cp[4], Float u, vec3f *deriv = nullptr) {
  point3f cp1[3] = {lerp(u, cp[0], cp[1]), lerp(u, cp[1], cp[2]),
                    lerp(u, cp[2], cp[3])};
  point3f cp2[2] = {lerp(u, cp1[0], cp1[1]), lerp(u, cp1[1], cp1[2])};
  if (deriv) {
    if ((cp2[1] - cp2[0]).squared_length() > 0) {
      *deriv = 3 * (cp2[1] - cp2[0]);
    } else {
      // For a cubic Bezier, if the first three control points (say) are
      // coincident, then the derivative of the curve is legitimately (0,0,0)
      // at u=0.  This is problematic for us, though, since we'd like to be
      // able to compute a surface normal there.  In that case, just punt and
      // take the difference between the first and last control points, which
      // ain't great, but will hopefully do.
      *deriv = cp[3] - cp[0];
    }
  }
  return(lerp(u, cp2[0], cp2[1]));
}


inline void SubdivideBezier(const point3f cp[4], point3f cpSplit[7]) {
  cpSplit[0] = cp[0];
  cpSplit[1] = (cp[0] + cp[1]) / 2;
  cpSplit[2] = (cp[0] + 2 * cp[1] + cp[2]) / 4;
  cpSplit[3] = (cp[0] + 3 * cp[1] + 3 * cp[2] + cp[3]) / 8;
  cpSplit[4] = (cp[1] + 2 * cp[2] + cp[3]) / 4;
  cpSplit[5] = (cp[2] + cp[3]) / 2;
  cpSplit[6] = cp[3];
}


CurveCommon::CurveCommon(const vec3f c[4], Float width0, Float width1,
                         CurveType type, const vec3f *norm)
  : type(type), cpObj{c[0], c[1], c[2], c[3]}, width{width0, width1} {
  width[0] = width0;
  width[1] = width1;
  for (int i = 0; i < 4; ++i) {
    cpObj[i] = c[i];
  }
  if (norm) {
    n[0] = unit_vector(norm[0]);
    n[1] = unit_vector(norm[1]);
    normalAngle = std::acos(clamp(dot(n[0], n[1]), 0, 1));
    invSinNormalAngle = 1.0 / std::sin(normalAngle);
  }
}

bool curve::bounding_box(Float t0, Float t1, aabb& box) const {
  // Compute object-space control points for curve segment, cpObj
  vec3f cpObj[4];
  cpObj[0] = BlossomBezier(common->cpObj, uMin, uMin, uMin);
  cpObj[1] = BlossomBezier(common->cpObj, uMin, uMin, uMax);
  cpObj[2] = BlossomBezier(common->cpObj, uMin, uMax, uMax);
  cpObj[3] = BlossomBezier(common->cpObj, uMax, uMax, uMax);
  box = surrounding_box(aabb(cpObj[0], cpObj[1]), aabb(cpObj[2], cpObj[3]));
  Float width[2] = {lerp(uMin, common->width[0], common->width[1]),
                    lerp(uMax, common->width[0], common->width[1])};
  box = (*ObjectToWorld)(Expand(box, std::fmax(width[0], width[1]) * 0.5f));
  return(true);
}

bool curve::hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng) {
  ray r2 = (*WorldToObject)(r); 
  
  // Compute object-space control points for curve segment, cpObj
  point3f cpObj[4];
  cpObj[0] = BlossomBezier(common->cpObj, uMin, uMin, uMin);
  cpObj[1] = BlossomBezier(common->cpObj, uMin, uMin, uMax);
  cpObj[2] = BlossomBezier(common->cpObj, uMin, uMax, uMax);
  cpObj[3] = BlossomBezier(common->cpObj, uMax, uMax, uMax);
  
  // Project curve control points to plane perpendicular to ray
  vec3f unit_dir = unit_vector(r2.direction()); 
  vec3f dx = cross(unit_dir, cpObj[3] - cpObj[0]);
  if (dx.squared_length() == 0) {
    // If the ray and the vector between the first and last control
    // points are parallel, dx will be zero.  Generate an arbitrary xy
    // orientation for the ray coordinate system so that intersection
    // tests can proceed in this unusual case.
    onb uvw;
    uvw.build_from_w(unit_dir);
    dx = uvw.v();
  } else {
    dx.make_unit_vector();
  }
  Transform objectToRay = LookAt(r2.origin(), r2.origin() + unit_dir, dx);
  point3f cp[4] = {objectToRay(cpObj[0]), objectToRay(cpObj[1]),
                   objectToRay(cpObj[2]), objectToRay(cpObj[3])};
  // Before going any further, see if the ray's bounding box intersects
  // the curve's bounding box. We start with the y dimension, since the y
  // extent is generally the smallest (and is often tiny) due to our
  // careful orientation of the ray coordinate system above.
  Float maxWidth = std::fmax(lerp(uMin, common->width[0], common->width[1]),
                             lerp(uMax, common->width[0], common->width[1]));
  if (std::fmax(std::fmax(cp[0].y(), cp[1].y()), std::fmax(cp[2].y(), cp[3].y())) +
      0.5f * maxWidth < 0 ||
      std::fmin(std::fmin(cp[0].y(), cp[1].y()), std::fmin(cp[2].y(), cp[3].y())) -
      0.5f * maxWidth > 0) {
    return false;
  }
  
  // Check for non-overlap in x.
  if (std::fmax(std::fmax(cp[0].x(), cp[1].x()), std::fmax(cp[2].x(), cp[3].x())) +
      0.5f * maxWidth < 0 ||
      std::fmin(std::fmin(cp[0].x(), cp[1].x()), std::fmin(cp[2].x(), cp[3].x())) -
      0.5f * maxWidth > 0) {
    return false;
  }
  
  // Check for non-overlap in z.
  Float rayLength = r2.direction().length();
  Float zMax = rayLength * tmax;
  if (std::fmax(std::fmax(cp[0].z(), cp[1].z()), std::fmax(cp[2].z(), cp[3].z())) +
      0.5f * maxWidth < 0 ||
      std::fmin(std::fmin(cp[0].z(), cp[1].z()), std::fmin(cp[2].z(), cp[3].z())) -
      0.5f * maxWidth > zMax) {
    return false;
  }
  
  // Compute refinement depth for curve, maxDepth
  Float L0 = 0;
  for (int i = 0; i < 2; ++i) {
    L0 = std::fmax(
      L0, std::fmax(
          std::fmax(std::fabs(cp[i].x() - 2 * cp[i + 1].x() + cp[i + 2].x()),
                   std::fabs(cp[i].y() - 2 * cp[i + 1].y() + cp[i + 2].y())),
                   std::fabs(cp[i].z() - 2 * cp[i + 1].z() + cp[i + 2].z())));
  }
  
  Float eps = std::fmax(common->width[0], common->width[1]) * .05f;  // width / 20
  auto Log2 = [](float v) -> int {
    if (v < 1) {
      return 0;
    }
    uint32_t bits = FloatToBits(v);
    // https://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
    // (With an additional add so get round-to-nearest rather than
    // round down.)
    return (bits >> 23) - 127 + (bits & (1 << 22) ? 1 : 0);
  };
  // Compute log base 4 by dividing log2 in half.
  int r0 = Log2(1.41421356237f * 6.f * L0 / (8.f * eps)) / 2;
  int maxDepth = clamp(r0, 0, 10);
  return(recursiveIntersect(r2, tmin, tmax, rec, cp, uMin,
                            uMax, maxDepth, Inverse(objectToRay)));
}


bool curve::hit(const ray& r, Float tmin, Float tmax, hit_record& rec, Sampler* sampler) {
  ray r2 = (*WorldToObject)(r); 
  
  // Compute object-space control points for curve segment, cpObj
  point3f cpObj[4];
  cpObj[0] = BlossomBezier(common->cpObj, uMin, uMin, uMin);
  cpObj[1] = BlossomBezier(common->cpObj, uMin, uMin, uMax);
  cpObj[2] = BlossomBezier(common->cpObj, uMin, uMax, uMax);
  cpObj[3] = BlossomBezier(common->cpObj, uMax, uMax, uMax);
  
  // Project curve control points to plane perpendicular to ray
  vec3f unit_dir = unit_vector(r2.direction()); 
  vec3f dx = cross(unit_dir, cpObj[3] - cpObj[0]);
  if (dx.squared_length() == 0) {
    // If the ray and the vector between the first and last control
    // points are parallel, dx will be zero.  Generate an arbitrary xy
    // orientation for the ray coordinate system so that intersection
    // tests can proceed in this unusual case.
    onb uvw;
    uvw.build_from_w(unit_dir);
    dx = uvw.v();
  } else {
    dx.make_unit_vector();
  }

  Transform objectToRay = LookAt(r2.origin(), r2.origin() + unit_dir, dx);
  point3f cp[4] = {objectToRay(cpObj[0]), objectToRay(cpObj[1]),
                   objectToRay(cpObj[2]), objectToRay(cpObj[3])};
  // Before going any further, see if the ray's bounding box intersects
  // the curve's bounding box. We start with the y dimension, since the y
  // extent is generally the smallest (and is often tiny) due to our
  // careful orientation of the ray coordinate system above.
  Float maxWidth = std::fmax(lerp(uMin, common->width[0], common->width[1]),
                             lerp(uMax, common->width[0], common->width[1]));
  if (std::fmax(std::fmax(cp[0].y(), cp[1].y()), std::fmax(cp[2].y(), cp[3].y())) +
      0.5f * maxWidth < 0 ||
      std::fmin(std::fmin(cp[0].y(), cp[1].y()), std::fmin(cp[2].y(), cp[3].y())) -
      0.5f * maxWidth > 0) {
    return false;
  }
  
  // Check for non-overlap in x.
  if (std::fmax(std::fmax(cp[0].x(), cp[1].x()), std::fmax(cp[2].x(), cp[3].x())) +
      0.5f * maxWidth < 0 ||
      std::fmin(std::fmin(cp[0].x(), cp[1].x()), std::fmin(cp[2].x(), cp[3].x())) -
      0.5f * maxWidth > 0) {
    return false;
  }
  
  // Check for non-overlap in z.
  Float rayLength = r2.direction().length();
  Float zMax = rayLength * tmax;
  if (std::fmax(std::fmax(cp[0].z(), cp[1].z()), std::fmax(cp[2].z(), cp[3].z())) +
      0.5f * maxWidth < 0 ||
      std::fmin(std::fmin(cp[0].z(), cp[1].z()), std::fmin(cp[2].z(), cp[3].z())) -
      0.5f * maxWidth > zMax) {
    return false;
  }
  
  // Compute refinement depth for curve, maxDepth
  Float L0 = 0;
  for (int i = 0; i < 2; ++i) {
    L0 = std::fmax(
      L0, std::fmax(
          std::fmax(std::fabs(cp[i].x() - 2 * cp[i + 1].x() + cp[i + 2].x()),
                    std::fabs(cp[i].y() - 2 * cp[i + 1].y() + cp[i + 2].y())),
                    std::fabs(cp[i].z() - 2 * cp[i + 1].z() + cp[i + 2].z())));
  }
  
  Float eps = std::fmax(common->width[0], common->width[1]) * .05f;  // width / 20
  auto Log2 = [](float v) -> int {
    if (v < 1) {
      return 0;
    }
    uint32_t bits = FloatToBits(v);
    // https://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
    // (With an additional add so get round-to-nearest rather than
    // round down.)
    return (bits >> 23) - 127 + (bits & (1 << 22) ? 1 : 0);
  };
  // Compute log base 4 by dividing log2 in half.
  int r0 = Log2(1.41421356237f * 6.f * L0 / (8.f * eps)) / 2;
  int maxDepth = clamp(r0, 0, 10);
  return(recursiveIntersect(r2, tmin, tmax, rec, cp, uMin,
                            uMax, maxDepth, Inverse(objectToRay)));
}

bool curve::recursiveIntersect(const ray& r, Float tmin, Float tmax, hit_record& rec, 
                               const point3f cp[4], Float u0, Float u1, int depth,
                               const Transform &rayToObject) const {
  Float rayLength = r.direction().length();

  if (depth > 0) {
    // Split curve segment into sub-segments and test for intersection
    point3f cpSplit[7];
    SubdivideBezier(cp, cpSplit);

    // For each of the two segments, see if the ray's bounding box
    // overlaps the segment before recursively checking for
    // intersection with it.
    bool hit = false;
    Float u[3] = {u0, (u0 + u1) / 2.f, u1};
    // Pointer to the 4 control points for the current segment.
    const point3f *cps = cpSplit;
    for (int seg = 0; seg < 2; ++seg, cps += 3) {
      Float maxWidth =
        std::fmax(lerp(u[seg],     common->width[0], common->width[1]),
                  lerp(u[seg + 1], common->width[0], common->width[1]));

      // As above, check y first, since it most commonly lets us exit
      // out early.
      if (std::fmax(std::fmax(cps[0].y(), cps[1].y()),
                    std::fmax(cps[2].y(), cps[3].y())) +
                     0.5 * maxWidth < 0 ||
                     std::fmin(std::fmin(cps[0].y(), cps[1].y()),
                               std::fmin(cps[2].y(), cps[3].y())) -
                                0.5 * maxWidth > 0) {
        continue;
      }
      if (std::fmax(std::fmax(cps[0].x(), cps[1].x()),
                   std::fmax(cps[2].x(), cps[3].x())) +
                     0.5 * maxWidth < 0 ||
                     std::fmin(std::fmin(cps[0].x(), cps[1].x()),
                              std::fmin(cps[2].x(), cps[3].x())) -
                                0.5 * maxWidth > 0) {
        continue;
      }
      Float zMax = rayLength * tmax;
      if (std::fmax(std::fmax(cps[0].z(), cps[1].z()),
                    std::fmax(cps[2].z(), cps[3].z())) +
                     0.5 * maxWidth < 0 ||
                     std::fmin(std::fmin(cps[0].z(), cps[1].z()),
                               std::fmin(cps[2].z(), cps[3].z())) -
                                0.5 * maxWidth > zMax) {
        continue;
      }

      hit |= recursiveIntersect(r, r.time(), tmax, rec, cps,
                                u[seg], u[seg + 1], depth - 1, rayToObject);
      // If we found an intersection and this is a shadow ray,
      // we can exit out immediately.
      if (hit && rec.t > tmin && rec.t < tmax) {
        return(true);
      }
    }
    return(hit);
  } else {
    // Test sample point against tangent perpendicular at curve start
    Float edge = (cp[1].y() - cp[0].y()) * -cp[0].y() + cp[0].x() * (cp[0].x() - cp[1].x());
    if (edge < 0) {
      return(false);
    }

    // Test sample point against tangent perpendicular at curve end
    edge = (cp[2].y() - cp[3].y()) * -cp[3].y() + cp[3].x() * (cp[3].x() - cp[2].x());
    if (edge < 0) {
      return(false);
    }

    // Compute line w that gives minimum distance to sample point
    vec2f segmentDirection = vec2f(cp[3]) - vec2f(cp[0]);
    Float denom = segmentDirection.squared_length();
    if (denom == 0) {
      return(false);
    }
    Float w = dot(-vec2f(cp[0]), segmentDirection) / denom;

    // Compute u coordinate of curve intersection point and hitWidth
    Float u = clamp(lerp(w, u0, u1), u0, u1);
    Float hitWidth = lerp(u, common->width[0], common->width[1]);
    normal3f nHit; // specified in world space
    bool flipped_n = false;
    if (common->type == CurveType::Ribbon) {
      // Scale hitWidth based on ribbon orientation
      Float sin0 = std::sin((1 - u) * common->normalAngle) * common->invSinNormalAngle;
      Float sin1 = std::sin(u * common->normalAngle) * common->invSinNormalAngle;
      if(!std::isnan(sin0) && !std::isnan(sin1)) {
        nHit = sin0 * common->n[0] + sin1 * common->n[1];
      } else {
        nHit = common->n[0];
      }
      //FIX TO WORLD SPACE
      hitWidth *= AbsDot(nHit, r.direction()) / rayLength; //both specified in world space
      flipped_n = dot(nHit, r.direction()) > 0;
    }

    // Test intersection point against curve width
    vec3f dpcdw;
    point3f pc = EvalBezier(cp, clamp(w, 0, 1), &dpcdw);
    Float ptCurveDist2 = pc.x() * pc.x() + pc.y() * pc.y();
    if (ptCurveDist2 > hitWidth * hitWidth * .25) {
      return(false);
    }
    Float zMax = rayLength * tmax;
    if (pc.z() < 0 || pc.z() > zMax) {
      return(false);
    }

    // Compute v coordinate of curve intersection point
    Float ptCurveDist = std::sqrt(ptCurveDist2);
    Float edgeFunc = dpcdw.x() * -pc.y() + pc.x() * dpcdw.y();
    Float v = (edgeFunc > 0) ? 0.5f + ptCurveDist / hitWidth : 0.5f - ptCurveDist / hitWidth;
    Float tnew =  pc.z() / rayLength;

    // Compute hit t and partial derivatives for curve intersection
    rec.t = tnew;

    // Compute dpdu and dpdv for curve intersection
    vec3f dpdu, dpdv;
    EvalBezier(common->cpObj, u, &rec.dpdu);
    
    if (common->type == CurveType::Cylinder) {
      vec3f dpduPlane = (Inverse(rayToObject))(rec.dpdu);
      vec3f dpdvPlane = unit_vector(vec3f(-dpduPlane.y(), dpduPlane.x(), 0)) * hitWidth;
      Float theta = lerp(v, -90.0, 90.0);
      Transform rot = Rotate(theta, dpduPlane);
      dpdvPlane = rot(dpdvPlane);
      rec.dpdv = rayToObject(dpdvPlane);
      rec.normal = unit_vector(-cross(rec.dpdu ,rec.dpdv));
    } else if (common->type == CurveType::Ribbon) {
      rec.dpdv = unit_vector(cross(vec3f(nHit.x(),nHit.y(),nHit.z()), rec.dpdu)) * hitWidth;
      rec.normal = !flipped_n ? nHit : -nHit;
    } else if (common->type == CurveType::Flat) {
      //Compute curve dpdv for flat and cylinder curves
      //Assumes z-axis faces directly at viewer
      rec.dpdv = unit_vector(cross(rec.dpdu,-r.direction()));
      rec.normal = unit_vector(-r.direction());
    } 
    rec.dpdu.make_unit_vector();
    
    rec.u = u;
    rec.v = v;
    rec.mat_ptr = mat_ptr.get();
    rec.has_bump = false;
    rec.pError = vec3f(hitWidth, hitWidth, hitWidth);
    rec.p = r.point_at_parameter(rec.t);
    rec.shape = this;
    rec.alpha_miss = false;
    rec = (*ObjectToWorld)(rec);

    
    return(true);
  }
  return(false);
}

vec3f curve::random(const point3f& o, random_gen& rng, Float time) {
  return(vec3f(0,0,0));
};

vec3f curve::random(const point3f& o, Sampler* sampler, Float time) {
  return(vec3f(0,0,0));
};

Float curve::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) {
  return(0);
}

Float curve::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) {
  return(0);
}

