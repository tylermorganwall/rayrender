#ifndef RAYRENDER_MEDIUM_BOUNDARY_H
#define RAYRENDER_MEDIUM_BOUNDARY_H
#include "../hitables/hitablelist.h"
#include "lights.h"
#include "medium.h"
#include <atomic>
#include <cstdlib>
#include <string>
#include <unordered_map>

class MediumBoundary final : public hitable {
public:
  MediumBoundary(std::shared_ptr<hitable> geometry, std::shared_ptr<const Medium> medium,
                 const Transform &object_to_world, bool keep_surface, uint64_t boundary_id);
  const bool hit(const Ray &, Float, Float, hit_record &, random_gen &) const override;
  const bool hit(const Ray &, Float, Float, hit_record &, Sampler *) const override;
  bool bounding_box(Float t0, Float t1, aabb &b) const override {
    return geometry->bounding_box(t0, t1, b);
  }
  Float pdf_value(const point3f &o, const vec3f &v, random_gen &r, Float t) override {
    return geometry->pdf_value(o, v, r, t);
  }
  Float pdf_value(const point3f &o, const vec3f &v, Sampler *r, Float t) override {
    return geometry->pdf_value(o, v, r, t);
  }
  vec3f random(const point3f &o, random_gen &r, Float t) override {
    return geometry->random(o, r, t);
  }
  vec3f random(const point3f &o, Sampler *r, Float t) override { return geometry->random(o, r, t); }
  std::string GetName() const override { return "MediumBoundary"; }
  size_t GetSize() override { return sizeof(*this) + geometry->GetSize(); }
  void hitable_info_bounds(Float t0, Float t1) const override {
    geometry->hitable_info_bounds(t0, t1);
  }
  std::shared_ptr<hitable> geometry;
  std::shared_ptr<const Medium> medium;
  Transform medium_to_world;
  bool keep_surface;
  Float orientation_sign = 1;

private:
  // ID in the owning scene; enclosing instances remap it into their parent scene.
  const uint64_t boundary_id;
  void Annotate(hit_record &) const;
};
struct MediumEntry {
  const MediumBoundary *boundary;
  uint64_t boundary_id;
  Transform medium_to_world;
};
struct VolumePathState {
  std::vector<MediumEntry> media;
  std::vector<dielectric *> glass;
  const MediumEntry *Active() const { return media.empty() ? nullptr : &media.back(); }
  void Cross(const hit_record &, const vec3f &direction);
  void SetRay(Ray &ray) {
    ray.segment_absorption = true;
    ray.pri_stack = &glass;
    ray.medium = Active() ? Active()->boundary->medium.get() : nullptr;
  }
};
struct VolumeStatistics {
  std::atomic<uint64_t> paths{0}, segments{0}, null_events{0}, scattering_events{0},
      shadow_candidates{0};
  void Reset() {
    paths = 0;
    segments = 0;
    null_events = 0;
    scattering_events = 0;
    shadow_candidates = 0;
  }
};
class VolumeScene {
public:
  // Local IDs are 1..BoundaryCount(), with zero reserved for ordinary surfaces.
  // Each instance reserves its child's entire range during scene construction.
  // Adding the returned offset to a child ID gives a unique parent-scene ID;
  // repeated additions identify the complete placement path without hashing.
  uint64_t ReserveBoundaryIds(uint64_t count);
  uint64_t NextBoundaryId() { return ReserveBoundaryIds(1) + 1; }
  uint64_t BoundaryCount() const { return boundary_count; }
  hitable_list boundaries;
  std::shared_ptr<VolumeLightSampler> light_sampler;
  using MediumCache = std::unordered_map<SEXP, std::shared_ptr<const Medium>>;
  std::shared_ptr<MediumCache> medium_cache = std::make_shared<MediumCache>();
  std::shared_ptr<const Medium> GetMedium(const Rcpp::List &description) {
    auto &result = (*medium_cache)[SEXP(description)];
    if (!result)
      result = LoadMedium(description);
    return result;
  }
  std::shared_ptr<hitable> boundary_bvh;
  std::vector<std::shared_ptr<VolumeScene>> children;
  bool has_media = false, has_emission = false, transparent_background = false;
  bool collect_statistics = std::getenv("RAYRENDER_VOLUME_STATS") &&
                            std::string(std::getenv("RAYRENDER_VOLUME_STATS")) == "true";
  mutable VolumeStatistics statistics;
  Rcpp::List Statistics() const;
  void Finish(Float t0, Float t1);
  VolumePathState InitialState(const Ray &, const std::atomic<bool> *cancel) const;

private:
  uint64_t boundary_count = 0;
};
// Returns false for ordinary geometry that cannot form a supported closed boundary.
bool ValidateMediumBoundary(hitable *geometry, bool required);
void ValidateMediumTransform(const Transform &);
point3f OffsetMediumOrigin(const hit_record &, const vec3f &);
#endif
