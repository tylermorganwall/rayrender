#ifndef RAYRENDER_RENDER_SPECTRAL_SCENE_H
#define RAYRENDER_RENDER_SPECTRAL_SCENE_H

#include "../base/base.h"
#include "../core/ray.h"
#include "../math/vectypes.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>

class hitable;

namespace rayrender {
namespace render {

template <typename Tag>
class SceneHandle {
public:
  using IndexType = std::uint32_t;
  using GenerationType = std::uint32_t;

  SceneHandle() = default;

  static constexpr SceneHandle Invalid() {
    return SceneHandle();
  }

  static constexpr SceneHandle FromIndex(IndexType index, GenerationType generation = 1) {
    return SceneHandle(index, generation);
  }

  static constexpr IndexType InvalidIndex() {
    return static_cast<IndexType>(-1);
  }

  constexpr bool IsValid() const {
    return index_ != InvalidIndex();
  }

  constexpr IndexType Index() const {
    return index_;
  }

  constexpr GenerationType Generation() const {
    return generation_;
  }

  constexpr bool operator==(const SceneHandle& other) const {
    return index_ == other.index_ && generation_ == other.generation_;
  }

  constexpr bool operator!=(const SceneHandle& other) const {
    return !(*this == other);
  }

private:
  constexpr SceneHandle(IndexType index, GenerationType generation) : index_(index), generation_(generation) {}

  IndexType index_ = InvalidIndex();
  GenerationType generation_ = 0;
};

struct ShapeHandleTag {};
struct PrimitiveHandleTag {};

using ShapeHandle = SceneHandle<ShapeHandleTag>;
using PrimitiveHandle = SceneHandle<PrimitiveHandleTag>;

struct Bounds3f {
  point3f pMin;
  point3f pMax;

  Bounds3f();
  Bounds3f(point3f minPoint, point3f maxPoint);

  bool IsValid() const;
  point3f Corner(int corner) const;
};

Bounds3f Union(const Bounds3f& lhs, const Bounds3f& rhs);
bool Inside(const point3f& p, const Bounds3f& bounds);

struct MediumInterface {
  base::MediumHandle inside = base::MediumHandle::Invalid();
  base::MediumHandle outside = base::MediumHandle::Invalid();

  bool IsMediumTransition() const;
  bool IsValid() const;
  base::MediumHandle Select(const vec3f& w, const normal3f& n) const;
};

struct SpawnedRay {
  Ray ray;
  base::MediumHandle medium = base::MediumHandle::Invalid();
  bool hasMedium = false;
};

struct Interaction {
  point3f p;
  vec3f pError;
  normal3f n;
  point2f uv;
  vec3f wo;
  Float time = 0;
  base::MediumHandle medium = base::MediumHandle::Invalid();
  bool hasMedium = false;
  MediumInterface mediumInterface;
  bool hasMediumInterface = false;

  bool IsSurfaceInteraction() const;
  point3f OffsetRayOrigin(const vec3f& w) const;
  SpawnedRay SpawnRay(const vec3f& d) const;
  SpawnedRay SpawnRayTo(const point3f& p2) const;
  base::MediumHandle GetMedium(const vec3f& w) const;
};

struct SurfaceInteraction : public Interaction {
  Float tHit = 0;
  vec3f dpdu;
  vec3f dpdv;
  normal3f dndu;
  normal3f dndv;
  normal3f shadingNormal;
  vec3f shadingDpdu;
  vec3f shadingDpdv;
  normal3f shadingDndu;
  normal3f shadingDndv;
  vec3f dpdx;
  vec3f dpdy;
  Float dudx = 0;
  Float dvdx = 0;
  Float dudy = 0;
  Float dvdy = 0;
  bool hasBump = false;
  normal3f bumpNormal;
  bool alphaMiss = false;
  bool infiniteAreaHit = false;
  int faceIndex = -1;
  ShapeHandle shape = ShapeHandle::Invalid();
  PrimitiveHandle primitive = PrimitiveHandle::Invalid();
  base::MaterialHandle material = base::MaterialHandle::Invalid();
  base::LightHandle areaLight = base::LightHandle::Invalid();
  bool hasAreaLight = false;
  int dielectricRegionId = -1;
  bool hasDielectricRegion = false;

  void SetShadingGeometry(
    normal3f ns,
    vec3f dpdus,
    vec3f dpdvs,
    normal3f dndus,
    normal3f dndvs,
    bool orientationIsAuthoritative
  );
};

struct LegacySurfaceHit {
  point3f p;
  Float t = 0;
  normal3f normal;
  vec3f dpdu;
  vec3f dpdv;
  vec3f pError;
  Float u = 0;
  Float v = 0;
  bool hasBump = false;
  normal3f bumpNormal;
  bool alphaMiss = false;
  bool infiniteAreaHit = false;
  int faceIndex = -1;
};

SurfaceInteraction SurfaceInteractionFromLegacyHit(
  const LegacySurfaceHit& hit,
  const Ray& ray,
  ShapeHandle shape = ShapeHandle::Invalid(),
  PrimitiveHandle primitive = PrimitiveHandle::Invalid()
);

struct ShapeIntersection {
  Float tHit = 0;
  SurfaceInteraction interaction;
};

struct ShapeSample {
  Interaction interaction;
  Float pdf = 0;
};

struct ShapeCapabilities {
  bool canIntersect = true;
  bool canBound = true;
  bool canSampleArea = false;
  bool canSampleDirection = false;
  bool hasUV = true;
  bool supportsNormalMap = true;
  bool supportsDisplacement = true;
  bool supportsAreaLight = false;
  bool supportsDielectricRegion = true;
  bool hasSignedDistance = false;
  bool supportsContainment = false;
  bool supportsFiniteDifferenceNormals = false;
  bool topologyClosed = false;
  bool topologyManifold = false;

  static ShapeCapabilities LegacyHitable();
  static ShapeCapabilities Sphere();
  static ShapeCapabilities CSGImplicit();
};

struct ShapeUsageRequirements {
  bool uvTexture = false;
  bool normalMap = false;
  bool displacement = false;
  bool areaLight = false;
  bool dielectricRegionBoundary = false;
};

std::vector<std::string> ValidateShapeCapabilities(
  const ShapeCapabilities& capabilities,
  const ShapeUsageRequirements& requirements
);

class Transform3f {
public:
  Transform3f();

  static Transform3f Identity();
  static Transform3f Translate(const vec3f& delta);
  static Transform3f Scale(Float x, Float y, Float z);
  static Transform3f FromMatrices(
    const std::array<Float, 16>& matrix,
    const std::array<Float, 16>& inverse
  );

  point3f ApplyPoint(const point3f& p) const;
  vec3f ApplyVector(const vec3f& v) const;
  normal3f ApplyNormal(const normal3f& n) const;
  Ray ApplyRay(const Ray& ray) const;
  point3f ApplyInversePoint(const point3f& p) const;
  vec3f ApplyInverseVector(const vec3f& v) const;
  Ray ApplyInverseRay(const Ray& ray) const;
  bool SwapsHandedness() const;

private:
  std::array<Float, 16> m_;
  std::array<Float, 16> mInv_;
};

struct CSGShapeOptions {
  Bounds3f bounds;
  Float hitThreshold = static_cast<Float>(1e-4);
  Float minStep = static_cast<Float>(1e-4);
  int maxSteps = 512;
};

using SignedDistanceFunction = std::function<Float(const point3f&)>;

struct ShapeCallbacks {
  std::function<std::optional<ShapeIntersection>(const Ray&, Float, Float)> intersect;
  std::function<bool(const Ray&, Float, Float)> hitP;
  std::function<bool(Bounds3f*)> bounds;
  std::function<Float()> area;
  std::function<std::optional<ShapeSample>(point2f)> sampleArea;
  std::function<Float(const Interaction&, const vec3f&)> pdfDirection;
  std::function<std::optional<ShapeSample>(const Interaction&, point2f)> sampleDirection;
  std::function<bool(const point3f&)> contains;
};

class Shape {
public:
  Shape() = default;

  static Shape Sphere(Float radius = 1, point3f center = point3f(0, 0, 0));
  static Shape Transformed(Shape shape, Transform3f renderFromObject);
  static Shape CSGSignedDistance(
    SignedDistanceFunction signedDistance,
    CSGShapeOptions options,
    std::string name = "CSG"
  );
  static Shape FromCallbacks(
    std::string name,
    ShapeCapabilities capabilities,
    ShapeCallbacks callbacks
  );
  static Shape LegacyHitable(
    std::shared_ptr<hitable> legacy,
    ShapeCapabilities capabilities = ShapeCapabilities::LegacyHitable()
  );

  bool IsValid() const;
  const std::string& Name() const;
  const ShapeCapabilities& Capabilities() const;
  std::optional<ShapeIntersection> Intersect(const Ray& ray, Float tMin, Float tMax) const;
  bool HitP(const Ray& ray, Float tMin, Float tMax) const;
  bool Bounds(Bounds3f* bounds) const;
  Float Area() const;
  std::optional<ShapeSample> SampleArea(point2f u) const;
  Float PDFDirection(const Interaction& ref, const vec3f& wi) const;
  std::optional<ShapeSample> SampleDirection(const Interaction& ref, point2f u) const;
  bool Contains(const point3f& p) const;

private:
  struct Impl;
  explicit Shape(std::shared_ptr<const Impl> impl);

  std::shared_ptr<const Impl> impl_;
};

struct PrimitiveBinding {
  base::MaterialHandle material = base::MaterialHandle::Invalid();
  base::LightHandle areaLight = base::LightHandle::Invalid();
  MediumInterface mediumInterface;
  bool hasMediumInterface = false;
  int dielectricRegionId = -1;
  bool hasDielectricRegion = false;
  ShapeUsageRequirements shapeRequirements;
};

struct PrimitiveIntersection {
  Float tHit = 0;
  SurfaceInteraction interaction;
};

class GeometricPrimitive {
public:
  GeometricPrimitive() = default;
  GeometricPrimitive(
    Shape shape,
    PrimitiveBinding binding,
    ShapeHandle shapeHandle = ShapeHandle::Invalid(),
    PrimitiveHandle primitiveHandle = PrimitiveHandle::Invalid()
  );

  const Shape& GetShape() const;
  const PrimitiveBinding& Binding() const;
  ShapeHandle GetShapeHandle() const;
  PrimitiveHandle GetPrimitiveHandle() const;
  std::optional<PrimitiveIntersection> Intersect(const Ray& ray, Float tMin, Float tMax) const;
  bool Bounds(Bounds3f* bounds) const;

private:
  Shape shape_;
  PrimitiveBinding binding_;
  ShapeHandle shapeHandle_ = ShapeHandle::Invalid();
  PrimitiveHandle primitiveHandle_ = PrimitiveHandle::Invalid();
};

class TransformedPrimitive {
public:
  TransformedPrimitive() = default;
  TransformedPrimitive(GeometricPrimitive primitive, Transform3f renderFromPrimitive);

  std::optional<PrimitiveIntersection> Intersect(const Ray& ray, Float tMin, Float tMax) const;
  bool Bounds(Bounds3f* bounds) const;

private:
  GeometricPrimitive primitive_;
  Transform3f renderFromPrimitive_;
};

class Aggregate {
public:
  void Add(GeometricPrimitive primitive);
  void Add(TransformedPrimitive primitive);
  std::optional<PrimitiveIntersection> Intersect(const Ray& ray, Float tMin, Float tMax) const;
  bool Bounds(Bounds3f* bounds) const;
  std::size_t Size() const;

private:
  std::vector<GeometricPrimitive> geometric_;
  std::vector<TransformedPrimitive> transformed_;
};

class Scene {
public:
  ShapeHandle AddShape(Shape shape);
  PrimitiveHandle AddPrimitive(ShapeHandle shape, PrimitiveBinding binding);
  void AddInfiniteLight(base::LightHandle light);
  void AddLight(base::LightHandle light);

  const Shape& GetShape(ShapeHandle handle) const;
  const Aggregate& GetAggregate() const;
  std::optional<PrimitiveIntersection> Intersect(const Ray& ray, Float tMin, Float tMax) const;
  bool Bounds(Bounds3f* bounds) const;
  std::size_t ShapeCount() const;
  std::size_t PrimitiveCount() const;
  const std::vector<base::LightHandle>& Lights() const;
  const std::vector<base::LightHandle>& InfiniteLights() const;

private:
  std::vector<Shape> shapes_;
  std::vector<ShapeHandle::GenerationType> shapeGenerations_;
  std::vector<PrimitiveBinding> primitiveBindings_;
  std::vector<PrimitiveHandle::GenerationType> primitiveGenerations_;
  Aggregate aggregate_;
  std::vector<base::LightHandle> lights_;
  std::vector<base::LightHandle> infiniteLights_;
};

} // namespace render
} // namespace rayrender

#endif
