#include "src/render/spectral_scene.h"

#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

using namespace rayrender::base;
using namespace rayrender::render;

namespace {

bool Approx(Float lhs, Float rhs, Float tolerance = static_cast<Float>(1e-5)) {
  return std::fabs(lhs - rhs) <= tolerance;
}

void Check(bool condition, const char* message) {
  if (!condition) {
    std::cerr << "PR10 Scene test failed: " << message << std::endl;
    std::exit(1);
  }
}

void CheckApprox(Float actual, Float expected, Float tolerance, const char* message) {
  if (!Approx(actual, expected, tolerance)) {
    std::cerr << "PR10 Scene test failed: " << message
              << " expected " << expected << " got " << actual << std::endl;
    std::exit(1);
  }
}

void CheckPointApprox(const point3f& actual, const point3f& expected, Float tolerance, const char* message) {
  CheckApprox(actual.xyz.x, expected.xyz.x, tolerance, message);
  CheckApprox(actual.xyz.y, expected.xyz.y, tolerance, message);
  CheckApprox(actual.xyz.z, expected.xyz.z, tolerance, message);
}

void CheckVectorApprox(const vec3f& actual, const vec3f& expected, Float tolerance, const char* message) {
  CheckApprox(actual.xyz.x, expected.xyz.x, tolerance, message);
  CheckApprox(actual.xyz.y, expected.xyz.y, tolerance, message);
  CheckApprox(actual.xyz.z, expected.xyz.z, tolerance, message);
}

void CheckNormalApprox(const normal3f& actual, const normal3f& expected, Float tolerance, const char* message) {
  CheckApprox(actual.xyz.x, expected.xyz.x, tolerance, message);
  CheckApprox(actual.xyz.y, expected.xyz.y, tolerance, message);
  CheckApprox(actual.xyz.z, expected.xyz.z, tolerance, message);
}

vec3f Normalize(const vec3f& value) {
  return value / value.length();
}

void TestLegacyHitConversion() {
  LegacySurfaceHit legacy;
  legacy.p = point3f(1, 2, 3);
  legacy.t = static_cast<Float>(4.5);
  legacy.normal = normal3f(0, 1, 0);
  legacy.dpdu = vec3f(1, 0, 0);
  legacy.dpdv = vec3f(0, 0, 1);
  legacy.pError = vec3f(static_cast<Float>(1e-4), static_cast<Float>(2e-4), static_cast<Float>(3e-4));
  legacy.u = static_cast<Float>(0.25);
  legacy.v = static_cast<Float>(0.75);
  legacy.hasBump = true;
  legacy.bumpNormal = normal3f(0, static_cast<Float>(0.8), static_cast<Float>(0.6));
  legacy.alphaMiss = true;
  legacy.infiniteAreaHit = true;
  legacy.faceIndex = 19;

  Ray ray(point3f(1, 2, 0), vec3f(0, 0, 1), static_cast<Float>(0.33));
  ShapeHandle shape = ShapeHandle::FromIndex(2, 1);
  PrimitiveHandle primitive = PrimitiveHandle::FromIndex(5, 1);
  SurfaceInteraction interaction = SurfaceInteractionFromLegacyHit(legacy, ray, shape, primitive);

  CheckPointApprox(interaction.p, legacy.p, static_cast<Float>(1e-7), "legacy hit p is preserved");
  CheckApprox(interaction.tHit, legacy.t, static_cast<Float>(1e-7), "legacy hit t is preserved");
  CheckNormalApprox(interaction.n, legacy.normal, static_cast<Float>(1e-7), "legacy hit normal is preserved");
  CheckVectorApprox(interaction.dpdu, legacy.dpdu, static_cast<Float>(1e-7), "legacy hit dpdu is preserved");
  CheckVectorApprox(interaction.dpdv, legacy.dpdv, static_cast<Float>(1e-7), "legacy hit dpdv is preserved");
  CheckVectorApprox(interaction.pError, legacy.pError, static_cast<Float>(1e-7), "legacy hit pError is preserved");
  CheckApprox(interaction.uv.xy.x, legacy.u, static_cast<Float>(1e-7), "legacy hit u is preserved");
  CheckApprox(interaction.uv.xy.y, legacy.v, static_cast<Float>(1e-7), "legacy hit v is preserved");
  Check(interaction.hasBump, "legacy bump flag is preserved");
  Check(interaction.alphaMiss, "legacy alpha flag is preserved");
  Check(interaction.infiniteAreaHit, "legacy infinite-area flag is preserved");
  Check(interaction.faceIndex == legacy.faceIndex, "legacy face index is preserved");
  Check(interaction.shape == shape, "shape handle is assigned");
  Check(interaction.primitive == primitive, "primitive handle is assigned");
  CheckVectorApprox(interaction.wo, vec3f(0, 0, -1), static_cast<Float>(1e-7), "wo is negative ray direction");
}

void TestSphereShapeIntersection() {
  Shape sphere = Shape::Sphere(1);
  Ray ray(point3f(0, 0, -3), vec3f(0, 0, 1), static_cast<Float>(0.2));
  std::optional<ShapeIntersection> hit = sphere.Intersect(ray, static_cast<Float>(0.001), Infinity);
  Check(hit.has_value(), "sphere intersects reference ray");
  CheckApprox(hit->tHit, 2, static_cast<Float>(1e-6), "sphere intersection t");
  CheckPointApprox(hit->interaction.p, point3f(0, 0, -1), static_cast<Float>(1e-6), "sphere intersection point");
  CheckNormalApprox(hit->interaction.n, normal3f(0, 0, -1), static_cast<Float>(1e-6), "sphere normal");
  CheckApprox(hit->interaction.uv.xy.x, static_cast<Float>(0.75), static_cast<Float>(1e-6), "sphere legacy u coordinate");
  CheckApprox(hit->interaction.uv.xy.y, static_cast<Float>(0.5), static_cast<Float>(1e-6), "sphere legacy v coordinate");
  Check(hit->interaction.dpdu.length() > 1, "sphere dpdu is populated");
  Check(hit->interaction.dpdv.length() > 1, "sphere dpdv is populated");
}

void TestShapeSampling() {
  Shape sphere = Shape::Sphere(2);
  CheckApprox(
    sphere.Area(),
    static_cast<Float>(16) * static_cast<Float>(M_PI),
    static_cast<Float>(1e-5),
    "sphere area"
  );
  std::optional<ShapeSample> sample = sphere.SampleArea(point2f(static_cast<Float>(0.25), static_cast<Float>(0.5)));
  Check(sample.has_value(), "sphere area sample exists");
  CheckApprox(
    sample->pdf,
    static_cast<Float>(1) / (static_cast<Float>(16) * static_cast<Float>(M_PI)),
    static_cast<Float>(1e-6),
    "sphere area sample pdf"
  );
  CheckApprox(sample->interaction.n.length(), 1, static_cast<Float>(1e-6), "sphere sample normal is unit length");

  Interaction ref;
  ref.p = point3f(0, 0, -6);
  ref.n = normal3f(0, 0, 1);
  ref.pError = vec3f(static_cast<Float>(1e-4));
  Float pdf = sphere.PDFDirection(ref, Normalize(point3f(0, 0, 0) - ref.p));
  Check(pdf > 0 && std::isfinite(pdf), "sphere directional pdf is finite positive");
}

void TestPrimitiveBindingsAndScene() {
  Scene scene;
  ShapeHandle sharedShape = scene.AddShape(Shape::Sphere(1));

  PrimitiveBinding materialOnly;
  materialOnly.material = MaterialHandle::FromIndex(11, 1);
  PrimitiveHandle first = scene.AddPrimitive(sharedShape, materialOnly);

  PrimitiveBinding emissiveBinding;
  emissiveBinding.material = MaterialHandle::FromIndex(12, 1);
  emissiveBinding.areaLight = LightHandle::FromIndex(21, 1);
  PrimitiveHandle second = scene.AddPrimitive(sharedShape, emissiveBinding);

  Check(scene.ShapeCount() == 1, "scene owns one shared shape");
  Check(scene.PrimitiveCount() == 2, "scene owns two primitives");
  Check(first != second, "shared shape primitives get distinct handles");

  GeometricPrimitive emissive(scene.GetShape(sharedShape), emissiveBinding, sharedShape, second);
  std::optional<PrimitiveIntersection> hit = emissive.Intersect(
    Ray(point3f(0, 0, -3), vec3f(0, 0, 1)),
    static_cast<Float>(0.001),
    Infinity
  );
  Check(hit.has_value(), "emissive primitive intersects");
  Check(hit->interaction.material == emissiveBinding.material, "primitive material binding is copied");
  Check(hit->interaction.hasAreaLight, "primitive area light binding is marked");
  Check(hit->interaction.areaLight == emissiveBinding.areaLight, "primitive area light handle is copied");
  Check(hit->interaction.material != MaterialHandle::Invalid(), "nonemissive material can coexist with area light");

  scene.AddLight(LightHandle::FromIndex(30, 1));
  scene.AddInfiniteLight(LightHandle::FromIndex(31, 1));
  Check(scene.Lights().size() == 1, "scene owns finite light handles");
  Check(scene.InfiniteLights().size() == 1, "scene owns infinite light handles");
}

void TestAggregateClosestHit() {
  Aggregate aggregate;
  PrimitiveBinding nearBinding;
  nearBinding.material = MaterialHandle::FromIndex(1, 1);
  PrimitiveBinding farBinding;
  farBinding.material = MaterialHandle::FromIndex(2, 1);

  aggregate.Add(GeometricPrimitive(
    Shape::Sphere(1, point3f(0, 0, 3)),
    farBinding,
    ShapeHandle::FromIndex(0, 1),
    PrimitiveHandle::FromIndex(0, 1)
  ));
  aggregate.Add(GeometricPrimitive(
    Shape::Sphere(1, point3f(0, 0, 0)),
    nearBinding,
    ShapeHandle::FromIndex(1, 1),
    PrimitiveHandle::FromIndex(1, 1)
  ));

  std::optional<PrimitiveIntersection> hit = aggregate.Intersect(
    Ray(point3f(0, 0, -5), vec3f(0, 0, 1)),
    static_cast<Float>(0.001),
    Infinity
  );
  Check(hit.has_value(), "aggregate finds closest hit");
  CheckApprox(hit->tHit, 4, static_cast<Float>(1e-6), "aggregate closest hit distance");
  Check(hit->interaction.material == nearBinding.material, "aggregate returns closest primitive binding");
}

void TestTransformedGeometry() {
  Shape translated = Shape::Transformed(Shape::Sphere(1), Transform3f::Translate(vec3f(2, 0, 0)));
  std::optional<ShapeIntersection> translatedHit = translated.Intersect(
    Ray(point3f(2, 0, -3), vec3f(0, 0, 1)),
    static_cast<Float>(0.001),
    Infinity
  );
  Check(translatedHit.has_value(), "translated shape intersects");
  CheckPointApprox(translatedHit->interaction.p, point3f(2, 0, -1), static_cast<Float>(1e-6), "translated hit point");
  CheckNormalApprox(translatedHit->interaction.n, normal3f(0, 0, -1), static_cast<Float>(1e-6), "translated normal");

  Transform3f mirror = Transform3f::Scale(static_cast<Float>(-1), 1, 1);
  Check(mirror.SwapsHandedness(), "negative scale reports handedness swap");
  Shape mirrored = Shape::Transformed(Shape::Sphere(1), mirror);
  std::optional<ShapeIntersection> mirroredHit = mirrored.Intersect(
    Ray(point3f(3, 0, 0), vec3f(-1, 0, 0)),
    static_cast<Float>(0.001),
    Infinity
  );
  Check(mirroredHit.has_value(), "mirrored shape intersects");
  CheckPointApprox(mirroredHit->interaction.p, point3f(1, 0, 0), static_cast<Float>(1e-6), "mirrored hit point");
  CheckNormalApprox(mirroredHit->interaction.n, normal3f(1, 0, 0), static_cast<Float>(1e-6), "mirrored final normal");
  Check(dot(mirroredHit->interaction.n, vec3f(1, 0, 0)) > 0, "mirrored normal keeps outward orientation");

  PrimitiveBinding binding;
  binding.material = MaterialHandle::FromIndex(7, 1);
  GeometricPrimitive primitive(Shape::Sphere(1), binding);
  TransformedPrimitive transformedPrimitive(primitive, Transform3f::Translate(vec3f(0, 0, 4)));
  std::optional<PrimitiveIntersection> primitiveHit = transformedPrimitive.Intersect(
    Ray(point3f(0, 0, 0), vec3f(0, 0, 1)),
    static_cast<Float>(0.001),
    Infinity
  );
  Check(primitiveHit.has_value(), "transformed primitive intersects");
  CheckApprox(primitiveHit->tHit, 3, static_cast<Float>(1e-6), "transformed primitive t");
}

void TestRayOriginOffsetsAndMedia() {
  Interaction interaction;
  interaction.p = point3f(1, 0, 0);
  interaction.pError = vec3f(static_cast<Float>(1e-3), static_cast<Float>(1e-3), static_cast<Float>(1e-3));
  interaction.n = normal3f(1, 0, 0);
  interaction.time = static_cast<Float>(0.4);
  interaction.mediumInterface.inside = MediumHandle::FromIndex(1, 1);
  interaction.mediumInterface.outside = MediumHandle::FromIndex(2, 1);
  interaction.hasMediumInterface = true;

  SpawnedRay outside = interaction.SpawnRay(vec3f(1, 0, 0));
  Check(outside.ray.origin().xyz.x > interaction.p.xyz.x, "spawned outside ray is offset along normal");
  Check(outside.medium == interaction.mediumInterface.outside, "outside medium selected");
  CheckApprox(outside.ray.time(), interaction.time, static_cast<Float>(1e-7), "spawned ray keeps time");

  SpawnedRay inside = interaction.SpawnRay(vec3f(-1, 0, 0));
  Check(inside.ray.origin().xyz.x < interaction.p.xyz.x, "spawned inside ray is offset opposite normal");
  Check(inside.medium == interaction.mediumInterface.inside, "inside medium selected");
}

void TestCSGCapabilitiesAndIntersection() {
  CSGShapeOptions options;
  options.bounds = Bounds3f(point3f(-2, -2, -2), point3f(2, 2, 2));
  Shape csg = Shape::CSGSignedDistance(
    [](const point3f& p) {
      return std::sqrt(p.xyz.x * p.xyz.x + p.xyz.y * p.xyz.y + p.xyz.z * p.xyz.z) - static_cast<Float>(1);
    },
    options,
    "SDF sphere"
  );

  Check(csg.Capabilities().hasSignedDistance, "CSG records signed-distance capability");
  Check(csg.Capabilities().supportsFiniteDifferenceNormals, "CSG records finite-difference normals");
  Check(csg.Contains(point3f(0, 0, 0)), "CSG containment uses signed distance");

  std::optional<ShapeIntersection> hit = csg.Intersect(
    Ray(point3f(0, 0, -3), vec3f(0, 0, 1)),
    static_cast<Float>(0.001),
    Infinity
  );
  Check(hit.has_value(), "CSG SDF intersects");
  CheckApprox(hit->tHit, 2, static_cast<Float>(2e-3), "CSG SDF t");
  CheckNormalApprox(hit->interaction.n, normal3f(0, 0, -1), static_cast<Float>(2e-3), "CSG finite-difference normal");

  ShapeUsageRequirements requirements;
  requirements.uvTexture = true;
  requirements.normalMap = true;
  requirements.displacement = true;
  requirements.areaLight = true;
  requirements.dielectricRegionBoundary = true;
  std::vector<std::string> errors = ValidateShapeCapabilities(csg.Capabilities(), requirements);
  Check(errors.size() == 5, "CSG rejects unsupported UV, normal-map, displacement, area-light, and region-boundary uses");

  PrimitiveBinding areaLightBinding;
  areaLightBinding.areaLight = LightHandle::FromIndex(8, 1);
  bool rejected = false;
  try {
    GeometricPrimitive invalid(csg, areaLightBinding);
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  Check(rejected, "primitive construction rejects area light on unsampleable CSG");
}

} // namespace

int main() {
  try {
    TestLegacyHitConversion();
    TestSphereShapeIntersection();
    TestShapeSampling();
    TestPrimitiveBindingsAndScene();
    TestAggregateClosestHit();
    TestTransformedGeometry();
    TestRayOriginOffsetsAndMedia();
    TestCSGCapabilitiesAndIntersection();
  } catch (const std::exception& error) {
    std::cerr << "PR10 Scene test failed with exception: " << error.what() << std::endl;
    return 1;
  }

  std::cout << "PR10 Scene tests passed" << std::endl;
  return 0;
}
