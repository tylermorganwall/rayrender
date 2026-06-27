#include "spectral_scene.h"

#include "../hitables/hitable.h"

#include <cfloat>
#include <stdexcept>
#include <utility>

namespace rayrender {
namespace render {

namespace {

LegacySurfaceHit LegacySurfaceHitFromRecord(const hit_record& record) {
  LegacySurfaceHit hit;
  hit.p = record.p;
  hit.t = record.t;
  hit.normal = record.normal;
  hit.dpdu = record.dpdu;
  hit.dpdv = record.dpdv;
  hit.pError = record.pError;
  hit.u = record.u;
  hit.v = record.v;
  hit.hasBump = record.has_bump;
  hit.bumpNormal = record.bump_normal;
  hit.alphaMiss = record.alpha_miss;
  hit.infiniteAreaHit = record.infinite_area_hit;
  return hit;
}

} // namespace

// We're capturing the shared ptr by value with [legacy] in the lambda to ensure it stays alive 
// as long as the callback exists

Shape Shape::LegacyHitable(
  std::shared_ptr<hitable> legacy,
  ShapeCapabilities capabilities
) {
  if (!legacy) {
    throw std::invalid_argument("LegacyHitable shape adapter requires a non-null hitable");
  }

  ShapeCallbacks callbacks;
  callbacks.intersect = [legacy](const Ray& ray, Float tMin, Float tMax) -> std::optional<ShapeIntersection> {
    random_gen rng(0);
    hit_record record;
    if (!legacy->hit(ray, tMin, tMax, record, rng)) {
      return std::nullopt;
    }
    ShapeIntersection intersection;
    intersection.tHit = record.t;
    intersection.interaction = SurfaceInteractionFromLegacyHit(
      LegacySurfaceHitFromRecord(record),
      ray
    );
    return intersection;
  };
  callbacks.hitP = [legacy](const Ray& ray, Float tMin, Float tMax) {
    random_gen rng(0);
    return legacy->HitP(ray, tMin, tMax, rng);
  };
  callbacks.bounds = [legacy](Bounds3f* bounds) {
    aabb legacyBounds;
    if (!legacy->bounding_box(0, 1, legacyBounds)) {
      return false;
    }
    *bounds = Bounds3f(legacyBounds.min(), legacyBounds.max());
    return true;
  };
  callbacks.area = []() {
    return static_cast<Float>(0);
  };
  callbacks.pdfDirection = [legacy](const Interaction& ref, const vec3f& wi) {
    random_gen rng(0);
    return legacy->pdf_value(ref.p, wi, rng, ref.time);
  };
  callbacks.sampleDirection = [legacy](const Interaction& ref, point2f) -> std::optional<ShapeSample> {
    random_gen rng(0);
    vec3f wi = legacy->random(ref.p, rng, ref.time);
    Float pdf = legacy->pdf_value(ref.p, wi, rng, ref.time);
    if (!(pdf > 0)) {
      return std::nullopt;
    }
    ShapeSample sample;
    sample.interaction = ref;
    sample.interaction.p = ref.p + wi;
    sample.pdf = pdf;
    return sample;
  };

  return Shape::FromCallbacks(legacy->GetName(), capabilities, std::move(callbacks));
}

} // namespace render
} // namespace rayrender
