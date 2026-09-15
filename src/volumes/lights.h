#ifndef RAYRENDER_VOLUME_LIGHTS_H
#define RAYRENDER_VOLUME_LIGHTS_H
#include "../hitables/hitablelist.h"

// Shape adapters supply solid-angle direction proposals. The sampler can use
// their full mixture for atmospheric transport or retain a selected emitter
// and its joint probability for ordinary endpoint emission.
class VolumeLightAdapter {
public:
  virtual ~VolumeLightAdapter() = default;
  virtual vec3f Sample(const point3f &, Sampler *, Float time) const = 0;
  virtual Float Pdf(const point3f &, const vec3f &, random_gen &, Float time) const = 0;
};
class VolumeLightSampler {
public:
  struct Selection {
    vec3f wi{0};
    Float pdf = 0; // Selection PMF times the selected emitter's solid-angle PDF.
    size_t index = 0;
  };
  explicit VolumeLightSampler(const hitable_list &);
  ~VolumeLightSampler();
  // Marginal directional strategy, retained for clear-air in-scattering.
  vec3f Sample(const point3f &, Sampler *, Float time) const;
  Float Pdf(const point3f &, const vec3f &, random_gen &, Float time) const;

  Selection SampleEmitter(const point3f &, Sampler *, random_gen &, Float time) const;
  // Resolve the chosen surface from the offset shadow origin. This is only
  // needed for deterministic visibility; stochastic walks resolve it once
  // through the world, then use Matches() without resampling alpha masks.
  bool Endpoint(const Selection &, const Ray &, hit_record &, random_gen &) const;
  bool Matches(const Selection &, const hit_record &) const;
  Float EmitterPdf(const point3f &, const vec3f &, const hit_record &,
                   random_gen &, Float time) const;

private:
  struct Emitters;
  std::unique_ptr<Emitters> emitters;
  std::unique_ptr<VolumeLightAdapter> lights;
};
#endif
