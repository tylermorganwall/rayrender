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
  // Keep this context at the previous real scattering vertex for emitter-hit
  // MIS. Null collisions and visibility offsets must not replace it. A zero
  // normal denotes a medium; surfaces use an absolute cosine so transmission
  // and two-sided materials retain support on both sides.
  struct Context {
    point3f p{0};
    normal3f n{0};
    Float time = 0;
  };
  enum class SelectionMethod { BVH, Fixed };
  struct Selection {
    vec3f wi{0};
    Float pdf = 0; // Selection PMF times the selected emitter's solid-angle PDF.
    size_t index = 0;
  };
  explicit VolumeLightSampler(const hitable_list &, SelectionMethod = SelectionMethod::BVH);
  ~VolumeLightSampler();
  // Marginal directional strategy, retained for clear-air in-scattering.
  vec3f Sample(const point3f &, Sampler *, Float time) const;
  Float Pdf(const point3f &, const vec3f &, random_gen &, Float time) const;

  Selection SampleEmitter(const Context &, Sampler *, random_gen &) const;
  // The discrete probability is also useful when validating the selection
  // distribution independently of a shape's conditional direction sampler.
  double SelectionPmf(const Context &, size_t index) const;
  // Resolve the chosen surface from the offset shadow origin. This is only
  // needed for deterministic visibility; stochastic walks resolve it once
  // through the world, then use Matches() without resampling alpha masks.
  bool Endpoint(const Selection &, const Ray &, hit_record &, random_gen &) const;
  bool Matches(const Selection &, const hit_record &) const;
  Float EmitterPdf(const Context &, const vec3f &, const hit_record &, random_gen &) const;

private:
  struct Emitters;
  std::unique_ptr<Emitters> emitters;
  std::unique_ptr<VolumeLightAdapter> lights;
};
#endif
