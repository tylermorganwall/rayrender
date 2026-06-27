#include "spectral_scene.h"

#include "../math/mathinline.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace rayrender {
namespace render {

namespace {

constexpr Float RayEpsilon = static_cast<Float>(1e-4);

Float DotNormal(const normal3f& lhs, const normal3f& rhs) {
  return lhs.xyz.x * rhs.xyz.x + lhs.xyz.y * rhs.xyz.y + lhs.xyz.z * rhs.xyz.z;
}

vec3f UnitOrZero(const vec3f& value) {
  Float length = value.length();
  if (!(length > 0) || !std::isfinite(length)) {
    return vec3f(0, 0, 0);
  }
  return value / length;
}

normal3f UnitOrZero(const normal3f& value) {
  Float length = value.length();
  if (!(length > 0) || !std::isfinite(length)) {
    return normal3f(0, 0, 0);
  }
  return value / length;
}

point3f TransformPoint(const std::array<Float, 16>& m, const point3f& p) {
  Float x = p.xyz.x;
  Float y = p.xyz.y;
  Float z = p.xyz.z;
  Float xp = m[0] * x + m[1] * y + m[2] * z + m[3];
  Float yp = m[4] * x + m[5] * y + m[6] * z + m[7];
  Float zp = m[8] * x + m[9] * y + m[10] * z + m[11];
  Float wp = m[12] * x + m[13] * y + m[14] * z + m[15];
  if (wp == 1 || wp == 0) {
    return point3f(xp, yp, zp);
  }
  return point3f(xp / wp, yp / wp, zp / wp);
}

vec3f TransformVector(const std::array<Float, 16>& m, const vec3f& v) {
  Float x = v.xyz.x;
  Float y = v.xyz.y;
  Float z = v.xyz.z;
  return vec3f(
    m[0] * x + m[1] * y + m[2] * z,
    m[4] * x + m[5] * y + m[6] * z,
    m[8] * x + m[9] * y + m[10] * z
  );
}

normal3f TransformNormal(const std::array<Float, 16>& mInv, const normal3f& n) {
  Float x = n.xyz.x;
  Float y = n.xyz.y;
  Float z = n.xyz.z;
  return UnitOrZero(normal3f(
    mInv[0] * x + mInv[4] * y + mInv[8] * z,
    mInv[1] * x + mInv[5] * y + mInv[9] * z,
    mInv[2] * x + mInv[6] * y + mInv[10] * z
  ));
}

vec3f TransformError(const Transform3f& transform, const vec3f& error) {
  return Abs(transform.ApplyVector(error));
}

Float ProjectTHit(const Ray& ray, const point3f& p) {
  vec3f delta = p - ray.origin();
  Float denom = ray.direction().squared_length();
  if (!(denom > 0)) {
    return 0;
  }
  return dot(delta, ray.direction()) / denom;
}

void GetSphereUV(const normal3f& n, Float& u, Float& v) {
  Float phi = std::atan2(n.xyz.z, n.xyz.x);
  Float theta = std::asin(std::clamp(n.xyz.y, static_cast<Float>(-1), static_cast<Float>(1)));
  u = 1 - (phi * static_cast<Float>(M_1_PI) + 1) * static_cast<Float>(0.5);
  v = theta * static_cast<Float>(M_1_PI) + static_cast<Float>(0.5);
}

bool BoundsIntersectRange(
  const Bounds3f& bounds,
  const Ray& ray,
  Float tMin,
  Float tMax,
  Float* entry,
  Float* exit
) {
  for (int axis = 0; axis < 3; ++axis) {
    if (ray.direction().e[axis] == 0) {
      if (ray.origin().e[axis] < bounds.pMin.e[axis] || ray.origin().e[axis] > bounds.pMax.e[axis]) {
        return false;
      }
      continue;
    }
    Float invD = static_cast<Float>(1) / ray.direction().e[axis];
    Float t0 = (bounds.pMin.e[axis] - ray.origin().e[axis]) * invD;
    Float t1 = (bounds.pMax.e[axis] - ray.origin().e[axis]) * invD;
    if (invD < 0) {
      std::swap(t0, t1);
    }
    tMin = t0 > tMin ? t0 : tMin;
    tMax = t1 < tMax ? t1 : tMax;
    if (tMax <= tMin) {
      return false;
    }
  }
  if (entry) {
    *entry = tMin;
  }
  if (exit) {
    *exit = tMax;
  }
  return true;
}

SurfaceInteraction TransformSurfaceInteraction(
  const SurfaceInteraction& local,
  const Transform3f& renderFromObject,
  const Ray& renderRay
) {
  SurfaceInteraction result = local;
  result.p = renderFromObject.ApplyPoint(local.p);
  result.pError = TransformError(renderFromObject, local.pError);
  result.n = renderFromObject.ApplyNormal(local.n);
  result.dpdu = renderFromObject.ApplyVector(local.dpdu);
  result.dpdv = renderFromObject.ApplyVector(local.dpdv);
  result.dndu = renderFromObject.ApplyNormal(local.dndu);
  result.dndv = renderFromObject.ApplyNormal(local.dndv);
  result.shadingNormal = renderFromObject.ApplyNormal(local.shadingNormal);
  result.shadingDpdu = renderFromObject.ApplyVector(local.shadingDpdu);
  result.shadingDpdv = renderFromObject.ApplyVector(local.shadingDpdv);
  result.shadingDndu = renderFromObject.ApplyNormal(local.shadingDndu);
  result.shadingDndv = renderFromObject.ApplyNormal(local.shadingDndv);
  result.bumpNormal = renderFromObject.ApplyNormal(local.bumpNormal);
  result.tHit = ProjectTHit(renderRay, result.p);
  return result;
}

class CallbackShape {
public:
  CallbackShape(std::string name, ShapeCapabilities capabilities, ShapeCallbacks callbacks)
    : name_(std::move(name)), capabilities_(capabilities), callbacks_(std::move(callbacks)) {}

  const std::string& Name() const {
    return name_;
  }

  const ShapeCapabilities& Capabilities() const {
    return capabilities_;
  }

  std::optional<ShapeIntersection> Intersect(const Ray& ray, Float tMin, Float tMax) const {
    if (!callbacks_.intersect) {
      return std::nullopt;
    }
    return callbacks_.intersect(ray, tMin, tMax);
  }

  bool HitP(const Ray& ray, Float tMin, Float tMax) const {
    if (callbacks_.hitP) {
      return callbacks_.hitP(ray, tMin, tMax);
    }
    return Intersect(ray, tMin, tMax).has_value();
  }

  bool Bounds(Bounds3f* bounds) const {
    if (!callbacks_.bounds) {
      return false;
    }
    return callbacks_.bounds(bounds);
  }

  Float Area() const {
    if (!callbacks_.area) {
      return 0;
    }
    return callbacks_.area();
  }

  std::optional<ShapeSample> SampleArea(point2f u) const {
    if (!callbacks_.sampleArea) {
      return std::nullopt;
    }
    return callbacks_.sampleArea(u);
  }

  Float PDFDirection(const Interaction& ref, const vec3f& wi) const {
    if (callbacks_.pdfDirection) {
      return callbacks_.pdfDirection(ref, wi);
    }
    Float area = Area();
    if (!(area > 0)) {
      return 0;
    }
    SpawnedRay ray = ref.SpawnRay(wi);
    std::optional<ShapeIntersection> intersection = Intersect(ray.ray, RayEpsilon, Infinity);
    if (!intersection) {
      return 0;
    }
    Float distanceSquared = (intersection->interaction.p - ref.p).squared_length();
    Float cosine = std::fabs(dot(intersection->interaction.n, -wi));
    if (!(cosine > 0)) {
      return 0;
    }
    return distanceSquared / (cosine * area);
  }

  std::optional<ShapeSample> SampleDirection(const Interaction& ref, point2f u) const {
    if (callbacks_.sampleDirection) {
      return callbacks_.sampleDirection(ref, u);
    }
    std::optional<ShapeSample> sample = SampleArea(u);
    if (!sample) {
      return std::nullopt;
    }
    vec3f wi = sample->interaction.p - ref.p;
    Float distanceSquared = wi.squared_length();
    wi = UnitOrZero(wi);
    Float cosine = std::fabs(dot(sample->interaction.n, -wi));
    if (!(cosine > 0) || !(sample->pdf > 0)) {
      return std::nullopt;
    }
    sample->pdf *= distanceSquared / cosine;
    return sample;
  }

  bool Contains(const point3f& p) const {
    if (!callbacks_.contains) {
      return false;
    }
    return callbacks_.contains(p);
  }

private:
  std::string name_;
  ShapeCapabilities capabilities_;
  ShapeCallbacks callbacks_;
};

} // namespace

struct Shape::Impl {
  explicit Impl(CallbackShape shape) : shape(std::move(shape)) {}

  CallbackShape shape;
};

Bounds3f::Bounds3f()
  : pMin(
      std::numeric_limits<Float>::max(),
      std::numeric_limits<Float>::max(),
      std::numeric_limits<Float>::max()
    ),
    pMax(
      std::numeric_limits<Float>::lowest(),
      std::numeric_limits<Float>::lowest(),
      std::numeric_limits<Float>::lowest()
    ) {}

Bounds3f::Bounds3f(point3f minPoint, point3f maxPoint)
  : pMin(
      std::min(minPoint.xyz.x, maxPoint.xyz.x),
      std::min(minPoint.xyz.y, maxPoint.xyz.y),
      std::min(minPoint.xyz.z, maxPoint.xyz.z)
    ),
    pMax(
      std::max(minPoint.xyz.x, maxPoint.xyz.x),
      std::max(minPoint.xyz.y, maxPoint.xyz.y),
      std::max(minPoint.xyz.z, maxPoint.xyz.z)
    ) {}

bool Bounds3f::IsValid() const {
  return pMin.xyz.x <= pMax.xyz.x &&
         pMin.xyz.y <= pMax.xyz.y &&
         pMin.xyz.z <= pMax.xyz.z;
}

point3f Bounds3f::Corner(int corner) const {
  return point3f(
    (corner & 1) ? pMax.xyz.x : pMin.xyz.x,
    (corner & 2) ? pMax.xyz.y : pMin.xyz.y,
    (corner & 4) ? pMax.xyz.z : pMin.xyz.z
  );
}

Bounds3f Union(const Bounds3f& lhs, const Bounds3f& rhs) {
  return Bounds3f(
    point3f(
      std::min(lhs.pMin.xyz.x, rhs.pMin.xyz.x),
      std::min(lhs.pMin.xyz.y, rhs.pMin.xyz.y),
      std::min(lhs.pMin.xyz.z, rhs.pMin.xyz.z)
    ),
    point3f(
      std::max(lhs.pMax.xyz.x, rhs.pMax.xyz.x),
      std::max(lhs.pMax.xyz.y, rhs.pMax.xyz.y),
      std::max(lhs.pMax.xyz.z, rhs.pMax.xyz.z)
    )
  );
}

bool Inside(const point3f& p, const Bounds3f& bounds) {
  return p.xyz.x >= bounds.pMin.xyz.x && p.xyz.x <= bounds.pMax.xyz.x &&
         p.xyz.y >= bounds.pMin.xyz.y && p.xyz.y <= bounds.pMax.xyz.y &&
         p.xyz.z >= bounds.pMin.xyz.z && p.xyz.z <= bounds.pMax.xyz.z;
}

bool MediumInterface::IsMediumTransition() const {
  return inside != outside;
}

bool MediumInterface::IsValid() const {
  return inside.IsValid() || outside.IsValid();
}

base::MediumHandle MediumInterface::Select(const vec3f& w, const normal3f& n) const {
  return dot(w, n) > 0 ? outside : inside;
}

bool Interaction::IsSurfaceInteraction() const {
  return n.squared_length() > 0;
}

point3f Interaction::OffsetRayOrigin(const vec3f& w) const {
  return ::OffsetRayOrigin(p, pError, n, w);
}

SpawnedRay Interaction::SpawnRay(const vec3f& d) const {
  base::MediumHandle selectedMedium = GetMedium(d);
  return SpawnedRay{Ray(OffsetRayOrigin(d), d, time), selectedMedium, selectedMedium.IsValid()};
}

SpawnedRay Interaction::SpawnRayTo(const point3f& p2) const {
  vec3f d = p2 - p;
  return SpawnRay(d);
}

base::MediumHandle Interaction::GetMedium(const vec3f& w) const {
  if (hasMediumInterface) {
    return mediumInterface.Select(w, n);
  }
  return hasMedium ? medium : base::MediumHandle::Invalid();
}

void SurfaceInteraction::SetShadingGeometry(
  normal3f ns,
  vec3f dpdus,
  vec3f dpdvs,
  normal3f dndus,
  normal3f dndvs,
  bool orientationIsAuthoritative
) {
  shadingNormal = UnitOrZero(ns);
  if (orientationIsAuthoritative) {
    if (DotNormal(n, shadingNormal) < 0) {
      n = -n;
    }
  } else if (DotNormal(shadingNormal, n) < 0) {
    shadingNormal = -shadingNormal;
  }
  shadingDpdu = dpdus;
  shadingDpdv = dpdvs;
  shadingDndu = dndus;
  shadingDndv = dndvs;
}

SurfaceInteraction SurfaceInteractionFromLegacyHit(
  const LegacySurfaceHit& hit,
  const Ray& ray,
  ShapeHandle shape,
  PrimitiveHandle primitive
) {
  SurfaceInteraction interaction;
  interaction.p = hit.p;
  interaction.tHit = hit.t;
  interaction.n = UnitOrZero(hit.normal);
  interaction.uv = point2f(hit.u, hit.v);
  interaction.wo = UnitOrZero(-ray.direction());
  interaction.time = ray.time();
  interaction.dpdu = hit.dpdu;
  interaction.dpdv = hit.dpdv;
  interaction.pError = hit.pError;
  interaction.hasBump = hit.hasBump;
  interaction.bumpNormal = UnitOrZero(hit.bumpNormal);
  interaction.alphaMiss = hit.alphaMiss;
  interaction.infiniteAreaHit = hit.infiniteAreaHit;
  interaction.faceIndex = hit.faceIndex;
  interaction.shape = shape;
  interaction.primitive = primitive;
  interaction.shadingNormal = interaction.hasBump ? interaction.bumpNormal : interaction.n;
  interaction.shadingDpdu = interaction.dpdu;
  interaction.shadingDpdv = interaction.dpdv;
  return interaction;
}

ShapeCapabilities ShapeCapabilities::LegacyHitable() {
  ShapeCapabilities capabilities;
  capabilities.canSampleArea = true;
  capabilities.canSampleDirection = true;
  capabilities.supportsAreaLight = true;
  capabilities.topologyClosed = false;
  capabilities.topologyManifold = false;
  return capabilities;
}

ShapeCapabilities ShapeCapabilities::Sphere() {
  ShapeCapabilities capabilities;
  capabilities.canSampleArea = true;
  capabilities.canSampleDirection = true;
  capabilities.supportsAreaLight = true;
  capabilities.topologyClosed = true;
  capabilities.topologyManifold = true;
  return capabilities;
}

ShapeCapabilities ShapeCapabilities::CSGImplicit() {
  ShapeCapabilities capabilities;
  capabilities.hasUV = false;
  capabilities.supportsNormalMap = false;
  capabilities.supportsDisplacement = false;
  capabilities.supportsAreaLight = false;
  capabilities.supportsDielectricRegion = false;
  capabilities.hasSignedDistance = true;
  capabilities.supportsContainment = true;
  capabilities.supportsFiniteDifferenceNormals = true;
  capabilities.topologyClosed = true;
  capabilities.topologyManifold = false;
  return capabilities;
}

std::vector<std::string> ValidateShapeCapabilities(
  const ShapeCapabilities& capabilities,
  const ShapeUsageRequirements& requirements
) {
  std::vector<std::string> errors;
  if (requirements.uvTexture && !capabilities.hasUV) {
    errors.push_back("shape does not provide UV coordinates required by UV texture");
  }
  if (requirements.normalMap && !capabilities.supportsNormalMap) {
    errors.push_back("shape does not support normal maps");
  }
  if (requirements.displacement && !capabilities.supportsDisplacement) {
    errors.push_back("shape does not support displacement");
  }
  if (requirements.areaLight && !capabilities.supportsAreaLight) {
    errors.push_back("shape does not support area lights");
  }
  if (requirements.dielectricRegionBoundary && !capabilities.supportsDielectricRegion) {
    errors.push_back("shape cannot be used as a dielectric region boundary");
  }
  return errors;
}

Transform3f::Transform3f()
  : m_({1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1}),
    mInv_(m_) {}

Transform3f Transform3f::Identity() {
  return Transform3f();
}

Transform3f Transform3f::Translate(const vec3f& delta) {
  std::array<Float, 16> matrix = {
    1, 0, 0, delta.xyz.x,
    0, 1, 0, delta.xyz.y,
    0, 0, 1, delta.xyz.z,
    0, 0, 0, 1
  };
  std::array<Float, 16> inverse = {
    1, 0, 0, -delta.xyz.x,
    0, 1, 0, -delta.xyz.y,
    0, 0, 1, -delta.xyz.z,
    0, 0, 0, 1
  };
  return FromMatrices(matrix, inverse);
}

Transform3f Transform3f::Scale(Float x, Float y, Float z) {
  if (x == 0 || y == 0 || z == 0) {
    throw std::invalid_argument("Transform3f::Scale requires non-zero scale factors");
  }
  std::array<Float, 16> matrix = {
    x, 0, 0, 0,
    0, y, 0, 0,
    0, 0, z, 0,
    0, 0, 0, 1
  };
  std::array<Float, 16> inverse = {
    static_cast<Float>(1) / x, 0, 0, 0,
    0, static_cast<Float>(1) / y, 0, 0,
    0, 0, static_cast<Float>(1) / z, 0,
    0, 0, 0, 1
  };
  return FromMatrices(matrix, inverse);
}

Transform3f Transform3f::FromMatrices(
  const std::array<Float, 16>& matrix,
  const std::array<Float, 16>& inverse
) {
  Transform3f transform;
  transform.m_ = matrix;
  transform.mInv_ = inverse;
  return transform;
}

point3f Transform3f::ApplyPoint(const point3f& p) const {
  return TransformPoint(m_, p);
}

vec3f Transform3f::ApplyVector(const vec3f& v) const {
  return TransformVector(m_, v);
}

normal3f Transform3f::ApplyNormal(const normal3f& n) const {
  return TransformNormal(mInv_, n);
}

Ray Transform3f::ApplyRay(const Ray& ray) const {
  return Ray(ApplyPoint(ray.origin()), ApplyVector(ray.direction()), ray.time(), ray.tMax);
}

point3f Transform3f::ApplyInversePoint(const point3f& p) const {
  return TransformPoint(mInv_, p);
}

vec3f Transform3f::ApplyInverseVector(const vec3f& v) const {
  return TransformVector(mInv_, v);
}

Ray Transform3f::ApplyInverseRay(const Ray& ray) const {
  return Ray(ApplyInversePoint(ray.origin()), ApplyInverseVector(ray.direction()), ray.time(), ray.tMax);
}

bool Transform3f::SwapsHandedness() const {
  Float det = m_[0] * (m_[5] * m_[10] - m_[6] * m_[9]) -
              m_[1] * (m_[4] * m_[10] - m_[6] * m_[8]) +
              m_[2] * (m_[4] * m_[9] - m_[5] * m_[8]);
  return det < 0;
}

Shape::Shape(std::shared_ptr<const Impl> impl) : impl_(std::move(impl)) {}

Shape Shape::Sphere(Float radius, point3f center) {
  if (!(radius > 0) || !std::isfinite(radius)) {
    throw std::invalid_argument("sphere shape radius must be positive finite");
  }

  ShapeCallbacks callbacks;
  callbacks.intersect = [radius, center](const Ray& ray, Float tMin, Float tMax) -> std::optional<ShapeIntersection> {
    vec3f oc = ray.origin() - center;
    Float a = dot(ray.direction(), ray.direction());
    Float halfB = dot(oc, ray.direction());
    Float c = dot(oc, oc) - radius * radius;
    Float discriminant = halfB * halfB - a * c;
    if (discriminant < 0) {
      return std::nullopt;
    }
    Float root = std::sqrt(discriminant);
    Float t = (-halfB - root) / a;
    if (t <= tMin || t >= tMax) {
      t = (-halfB + root) / a;
      if (t <= tMin || t >= tMax) {
        return std::nullopt;
      }
    }

    point3f p = ray(t);
    vec3f local = p - center;
    normal3f n = UnitOrZero(convert_to_normal3(local / radius));
    Float u = 0;
    Float v = 0;
    GetSphereUV(n, u, v);

    LegacySurfaceHit legacy;
    legacy.p = p;
    legacy.t = t;
    legacy.normal = n;
    legacy.pError = gamma(5) * Abs(local);
    legacy.u = u;
    legacy.v = v;

    Float zRadius = std::sqrt(local.xyz.x * local.xyz.x + local.xyz.z * local.xyz.z);
    if (zRadius > 0) {
      Float invZRadius = static_cast<Float>(1) / zRadius;
      Float cosPhi = local.xyz.x * invZRadius;
      Float sinPhi = local.xyz.z * invZRadius;
      Float theta = std::acos(std::clamp(local.xyz.z / radius, static_cast<Float>(-1), static_cast<Float>(1)));
      legacy.dpdu = static_cast<Float>(2) * static_cast<Float>(M_PI) *
                    vec3f(-local.xyz.z, 0, local.xyz.x);
      legacy.dpdv = static_cast<Float>(2) * static_cast<Float>(M_PI) *
                    vec3f(local.xyz.z * cosPhi, local.xyz.z * sinPhi, -radius * std::sin(theta));
    }

    ShapeIntersection intersection;
    intersection.tHit = t;
    intersection.interaction = SurfaceInteractionFromLegacyHit(legacy, ray);
    return intersection;
  };
  callbacks.hitP = [callbacks](const Ray& ray, Float tMin, Float tMax) {
    return callbacks.intersect(ray, tMin, tMax).has_value();
  };
  callbacks.bounds = [radius, center](Bounds3f* bounds) {
    *bounds = Bounds3f(
      center - vec3f(radius, radius, radius),
      center + vec3f(radius, radius, radius)
    );
    return true;
  };
  callbacks.area = [radius]() {
    return static_cast<Float>(4) * static_cast<Float>(M_PI) * radius * radius;
  };
  callbacks.sampleArea = [radius, center](point2f u) -> std::optional<ShapeSample> {
    Float z = static_cast<Float>(1) - static_cast<Float>(2) * u.xy.x;
    Float r = std::sqrt(std::max(static_cast<Float>(0), static_cast<Float>(1) - z * z));
    Float phi = static_cast<Float>(2) * static_cast<Float>(M_PI) * u.xy.y;
    vec3f local(radius * r * std::cos(phi), radius * r * std::sin(phi), radius * z);
    normal3f n = UnitOrZero(convert_to_normal3(local / radius));
    Float su = 0;
    Float sv = 0;
    GetSphereUV(n, su, sv);

    ShapeSample sample;
    sample.interaction.p = center + local;
    sample.interaction.pError = gamma(5) * Abs(local);
    sample.interaction.n = n;
    sample.interaction.uv = point2f(su, sv);
    sample.pdf = static_cast<Float>(1) /
                 (static_cast<Float>(4) * static_cast<Float>(M_PI) * radius * radius);
    return sample;
  };
  callbacks.contains = [radius, center](const point3f& p) {
    return (p - center).squared_length() <= radius * radius;
  };

  return FromCallbacks("Sphere", ShapeCapabilities::Sphere(), std::move(callbacks));
}

Shape Shape::Transformed(Shape shape, Transform3f renderFromObject) {
  if (!shape.IsValid()) {
    throw std::invalid_argument("cannot transform an invalid shape");
  }
  ShapeCapabilities capabilities = shape.Capabilities();
  ShapeCallbacks callbacks;
  callbacks.intersect = [shape, renderFromObject](const Ray& ray, Float tMin, Float tMax) -> std::optional<ShapeIntersection> {
    Ray objectRay = renderFromObject.ApplyInverseRay(ray);
    std::optional<ShapeIntersection> objectHit = shape.Intersect(objectRay, tMin, tMax);
    if (!objectHit) {
      return std::nullopt;
    }
    ShapeIntersection result = *objectHit;
    result.interaction = TransformSurfaceInteraction(objectHit->interaction, renderFromObject, ray);
    result.tHit = result.interaction.tHit;
    return result;
  };
  callbacks.hitP = [callbacks](const Ray& ray, Float tMin, Float tMax) {
    return callbacks.intersect(ray, tMin, tMax).has_value();
  };
  callbacks.bounds = [shape, renderFromObject](Bounds3f* bounds) {
    Bounds3f childBounds;
    if (!shape.Bounds(&childBounds)) {
      return false;
    }
    Bounds3f transformed(renderFromObject.ApplyPoint(childBounds.Corner(0)), renderFromObject.ApplyPoint(childBounds.Corner(0)));
    for (int i = 1; i < 8; ++i) {
      transformed = Union(
        transformed,
        Bounds3f(renderFromObject.ApplyPoint(childBounds.Corner(i)), renderFromObject.ApplyPoint(childBounds.Corner(i)))
      );
    }
    *bounds = transformed;
    return true;
  };
  callbacks.area = [shape]() {
    return shape.Area();
  };
  callbacks.contains = [shape, renderFromObject](const point3f& p) {
    return shape.Contains(renderFromObject.ApplyInversePoint(p));
  };
  return FromCallbacks("Transformed " + shape.Name(), capabilities, std::move(callbacks));
}

Shape Shape::CSGSignedDistance(
  SignedDistanceFunction signedDistance,
  CSGShapeOptions options,
  std::string name
) {
  if (!signedDistance) {
    throw std::invalid_argument("CSG shape requires a signed-distance function");
  }
  if (!options.bounds.IsValid()) {
    throw std::invalid_argument("CSG shape requires valid finite bounds");
  }
  if (!(options.hitThreshold > 0) || !(options.minStep > 0) || options.maxSteps <= 0) {
    throw std::invalid_argument("CSG shape requires positive hit threshold, min step, and max steps");
  }

  ShapeCallbacks callbacks;
  callbacks.intersect = [signedDistance, options](const Ray& ray, Float tMin, Float tMax) -> std::optional<ShapeIntersection> {
    Float entry = tMin;
    Float exit = tMax;
    if (!BoundsIntersectRange(options.bounds, ray, tMin, tMax, &entry, &exit)) {
      return std::nullopt;
    }
    Float t = std::max(entry, static_cast<Float>(0));
    Float end = std::min(exit, tMax);
    for (int i = 0; i < options.maxSteps && t < end; ++i) {
      point3f p = ray(t);
      if (!Inside(p, options.bounds)) {
        t += options.minStep;
        continue;
      }
      Float d = signedDistance(p);
      if (std::fabs(d) <= options.hitThreshold) {
        Float eps = std::max(options.hitThreshold * static_cast<Float>(4), static_cast<Float>(1e-4));
        normal3f n(
          signedDistance(point3f(p.xyz.x + eps, p.xyz.y, p.xyz.z)) -
            signedDistance(point3f(p.xyz.x - eps, p.xyz.y, p.xyz.z)),
          signedDistance(point3f(p.xyz.x, p.xyz.y + eps, p.xyz.z)) -
            signedDistance(point3f(p.xyz.x, p.xyz.y - eps, p.xyz.z)),
          signedDistance(point3f(p.xyz.x, p.xyz.y, p.xyz.z + eps)) -
            signedDistance(point3f(p.xyz.x, p.xyz.y, p.xyz.z - eps))
        );
        LegacySurfaceHit legacy;
        legacy.p = p;
        legacy.t = t;
        legacy.normal = UnitOrZero(n);
        legacy.pError = vec3f(options.hitThreshold, options.hitThreshold, options.hitThreshold);
        ShapeIntersection intersection;
        intersection.tHit = t;
        intersection.interaction = SurfaceInteractionFromLegacyHit(legacy, ray);
        return intersection;
      }
      t += std::max(std::fabs(d), options.minStep);
    }
    return std::nullopt;
  };
  callbacks.hitP = [callbacks](const Ray& ray, Float tMin, Float tMax) {
    return callbacks.intersect(ray, tMin, tMax).has_value();
  };
  callbacks.bounds = [options](Bounds3f* bounds) {
    *bounds = options.bounds;
    return true;
  };
  callbacks.contains = [signedDistance](const point3f& p) {
    return signedDistance(p) <= 0;
  };

  return FromCallbacks(std::move(name), ShapeCapabilities::CSGImplicit(), std::move(callbacks));
}

Shape Shape::FromCallbacks(
  std::string name,
  ShapeCapabilities capabilities,
  ShapeCallbacks callbacks
) {
  return Shape(std::make_shared<Impl>(CallbackShape(std::move(name), capabilities, std::move(callbacks))));
}

bool Shape::IsValid() const {
  return static_cast<bool>(impl_);
}

const std::string& Shape::Name() const {
  if (!impl_) {
    static const std::string invalid = "InvalidShape";
    return invalid;
  }
  return impl_->shape.Name();
}

const ShapeCapabilities& Shape::Capabilities() const {
  if (!impl_) {
    static const ShapeCapabilities invalid{};
    return invalid;
  }
  return impl_->shape.Capabilities();
}

std::optional<ShapeIntersection> Shape::Intersect(const Ray& ray, Float tMin, Float tMax) const {
  if (!impl_) {
    return std::nullopt;
  }
  return impl_->shape.Intersect(ray, tMin, tMax);
}

bool Shape::HitP(const Ray& ray, Float tMin, Float tMax) const {
  return impl_ && impl_->shape.HitP(ray, tMin, tMax);
}

bool Shape::Bounds(Bounds3f* bounds) const {
  return impl_ && impl_->shape.Bounds(bounds);
}

Float Shape::Area() const {
  return impl_ ? impl_->shape.Area() : 0;
}

std::optional<ShapeSample> Shape::SampleArea(point2f u) const {
  if (!impl_) {
    return std::nullopt;
  }
  return impl_->shape.SampleArea(u);
}

Float Shape::PDFDirection(const Interaction& ref, const vec3f& wi) const {
  return impl_ ? impl_->shape.PDFDirection(ref, wi) : 0;
}

std::optional<ShapeSample> Shape::SampleDirection(const Interaction& ref, point2f u) const {
  if (!impl_) {
    return std::nullopt;
  }
  return impl_->shape.SampleDirection(ref, u);
}

bool Shape::Contains(const point3f& p) const {
  return impl_ && impl_->shape.Contains(p);
}

GeometricPrimitive::GeometricPrimitive(
  Shape shape,
  PrimitiveBinding binding,
  ShapeHandle shapeHandle,
  PrimitiveHandle primitiveHandle
)
  : shape_(std::move(shape)),
    binding_(binding),
    shapeHandle_(shapeHandle),
    primitiveHandle_(primitiveHandle) {
  ShapeUsageRequirements requirements = binding_.shapeRequirements;
  requirements.areaLight = requirements.areaLight || binding_.areaLight.IsValid();
  std::vector<std::string> errors = ValidateShapeCapabilities(shape_.Capabilities(), requirements);
  if (!errors.empty()) {
    throw std::invalid_argument(errors.front());
  }
}

const Shape& GeometricPrimitive::GetShape() const {
  return shape_;
}

const PrimitiveBinding& GeometricPrimitive::Binding() const {
  return binding_;
}

ShapeHandle GeometricPrimitive::GetShapeHandle() const {
  return shapeHandle_;
}

PrimitiveHandle GeometricPrimitive::GetPrimitiveHandle() const {
  return primitiveHandle_;
}

std::optional<PrimitiveIntersection> GeometricPrimitive::Intersect(
  const Ray& ray,
  Float tMin,
  Float tMax
) const {
  std::optional<ShapeIntersection> shapeHit = shape_.Intersect(ray, tMin, tMax);
  if (!shapeHit) {
    return std::nullopt;
  }
  PrimitiveIntersection result;
  result.tHit = shapeHit->tHit;
  result.interaction = shapeHit->interaction;
  result.interaction.shape = shapeHandle_;
  result.interaction.primitive = primitiveHandle_;
  result.interaction.material = binding_.material;
  result.interaction.areaLight = binding_.areaLight;
  result.interaction.hasAreaLight = binding_.areaLight.IsValid();
  result.interaction.mediumInterface = binding_.mediumInterface;
  result.interaction.hasMediumInterface = binding_.hasMediumInterface;
  result.interaction.dielectricRegionId = binding_.dielectricRegionId;
  result.interaction.hasDielectricRegion = binding_.hasDielectricRegion;
  return result;
}

bool GeometricPrimitive::Bounds(Bounds3f* bounds) const {
  return shape_.Bounds(bounds);
}

TransformedPrimitive::TransformedPrimitive(
  GeometricPrimitive primitive,
  Transform3f renderFromPrimitive
)
  : primitive_(std::move(primitive)), renderFromPrimitive_(renderFromPrimitive) {}

std::optional<PrimitiveIntersection> TransformedPrimitive::Intersect(
  const Ray& ray,
  Float tMin,
  Float tMax
) const {
  Shape transformed = Shape::Transformed(primitive_.GetShape(), renderFromPrimitive_);
  GeometricPrimitive transformedPrimitive(
    transformed,
    primitive_.Binding(),
    primitive_.GetShapeHandle(),
    primitive_.GetPrimitiveHandle()
  );
  return transformedPrimitive.Intersect(ray, tMin, tMax);
}

bool TransformedPrimitive::Bounds(Bounds3f* bounds) const {
  Shape transformed = Shape::Transformed(primitive_.GetShape(), renderFromPrimitive_);
  return transformed.Bounds(bounds);
}

void Aggregate::Add(GeometricPrimitive primitive) {
  geometric_.push_back(std::move(primitive));
}

void Aggregate::Add(TransformedPrimitive primitive) {
  transformed_.push_back(std::move(primitive));
}

std::optional<PrimitiveIntersection> Aggregate::Intersect(
  const Ray& ray,
  Float tMin,
  Float tMax
) const {
  std::optional<PrimitiveIntersection> closest;
  Float closestT = tMax;
  for (const GeometricPrimitive& primitive : geometric_) {
    std::optional<PrimitiveIntersection> hit = primitive.Intersect(ray, tMin, closestT);
    if (hit) {
      closestT = hit->tHit;
      closest = hit;
    }
  }
  for (const TransformedPrimitive& primitive : transformed_) {
    std::optional<PrimitiveIntersection> hit = primitive.Intersect(ray, tMin, closestT);
    if (hit) {
      closestT = hit->tHit;
      closest = hit;
    }
  }
  return closest;
}

bool Aggregate::Bounds(Bounds3f* bounds) const {
  bool haveBounds = false;
  Bounds3f result;
  for (const GeometricPrimitive& primitive : geometric_) {
    Bounds3f primitiveBounds;
    if (!primitive.Bounds(&primitiveBounds)) {
      continue;
    }
    result = haveBounds ? Union(result, primitiveBounds) : primitiveBounds;
    haveBounds = true;
  }
  for (const TransformedPrimitive& primitive : transformed_) {
    Bounds3f primitiveBounds;
    if (!primitive.Bounds(&primitiveBounds)) {
      continue;
    }
    result = haveBounds ? Union(result, primitiveBounds) : primitiveBounds;
    haveBounds = true;
  }
  if (!haveBounds) {
    return false;
  }
  *bounds = result;
  return true;
}

std::size_t Aggregate::Size() const {
  return geometric_.size() + transformed_.size();
}

ShapeHandle Scene::AddShape(Shape shape) {
  if (!shape.IsValid()) {
    throw std::invalid_argument("cannot add invalid shape to spectral scene");
  }
  ShapeHandle handle = ShapeHandle::FromIndex(static_cast<ShapeHandle::IndexType>(shapes_.size()), 1);
  shapes_.push_back(std::move(shape));
  shapeGenerations_.push_back(handle.Generation());
  return handle;
}

PrimitiveHandle Scene::AddPrimitive(ShapeHandle shape, PrimitiveBinding binding) {
  const Shape& sceneShape = GetShape(shape);
  PrimitiveHandle handle = PrimitiveHandle::FromIndex(
    static_cast<PrimitiveHandle::IndexType>(primitiveBindings_.size()),
    1
  );
  GeometricPrimitive primitive(sceneShape, binding, shape, handle);
  primitiveBindings_.push_back(binding);
  primitiveGenerations_.push_back(handle.Generation());
  aggregate_.Add(std::move(primitive));
  return handle;
}

void Scene::AddInfiniteLight(base::LightHandle light) {
  if (light.IsValid()) {
    infiniteLights_.push_back(light);
  }
}

void Scene::AddLight(base::LightHandle light) {
  if (light.IsValid()) {
    lights_.push_back(light);
  }
}

const Shape& Scene::GetShape(ShapeHandle handle) const {
  if (!handle.IsValid() || handle.Index() >= shapes_.size() ||
      shapeGenerations_[handle.Index()] != handle.Generation()) {
    throw std::out_of_range("invalid spectral shape handle");
  }
  return shapes_[handle.Index()];
}

const Aggregate& Scene::GetAggregate() const {
  return aggregate_;
}

std::optional<PrimitiveIntersection> Scene::Intersect(const Ray& ray, Float tMin, Float tMax) const {
  return aggregate_.Intersect(ray, tMin, tMax);
}

bool Scene::Bounds(Bounds3f* bounds) const {
  return aggregate_.Bounds(bounds);
}

std::size_t Scene::ShapeCount() const {
  return shapes_.size();
}

std::size_t Scene::PrimitiveCount() const {
  return primitiveBindings_.size();
}

const std::vector<base::LightHandle>& Scene::Lights() const {
  return lights_;
}

const std::vector<base::LightHandle>& Scene::InfiniteLights() const {
  return infiniteLights_;
}

} // namespace render
} // namespace rayrender
