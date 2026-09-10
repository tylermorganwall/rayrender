#ifndef RAYRENDER_ATMOSPHERE_H
#define RAYRENDER_ATMOSPHERE_H
#include "infinite.h"
#include <skymodelr/prague.hpp>
#include <array>

// Require the consolidated provider SDK, including in manual source builds.
// The wrapper validates the loaded provider's ABI and table size at runtime.
#if !defined(SKYMODELR_PRAGUE_ABI_VERSION) || SKYMODELR_PRAGUE_ABI_VERSION < 3
#error "rayrender requires the current Prague API: rebuild skymodelr 0.6.2 or later first."
#endif

struct AtmosphereSegment {
  point3f transmission{1}, radiance{0};
};
struct WeightedAtmosphereSegment {
  point3f transmission{1};
  std::array<double, 3> radiance{};
};
// Per-ray spectral endpoints shared by cumulative segment queries. The cache
// belongs to the path, so coefficient data stay immutable across workers.
struct AtmosphereSegmentCache {
  std::array<double, 16> radiance{}, zero_transmission{}, near_transmission{};
  bool ready = false, near_ready = false;
};
class Atmosphere {
public:
  virtual ~Atmosphere() = default;
  virtual AtmosphereSegment Segment(const point3f &, const vec3f &, double distance) const = 0;
  virtual AtmosphereSegment Segment(const point3f &p, const vec3f &w, double distance,
                                    AtmosphereSegmentCache &) const {
    return Segment(p, w, distance);
  }
  // Deferred callers supply independent uniforms for sampled source estimates.
  // Direct Segment queries always retain the complete deterministic evaluator.
  virtual bool SampleHazeCorrection() const { return false; }
  virtual AtmosphereSegment SampledSegment(const point3f &p, const vec3f &w, double d,
                                           AtmosphereSegmentCache &cache, double uniform) const {
    return Segment(p, w, d, cache);
  }
  virtual point3f Transmission(const point3f &, const vec3f &, double distance) const = 0;
  virtual point3f CelestialTransmission(const point3f &p, const vec3f &w,
                                        InfiniteLightSpectrum) const {
    return Transmission(p, w, INFINITY);
  }
  virtual point3f SkyRadiance(const point3f &, const vec3f &) const = 0;
  virtual bool MaySeeDisk(const point3f &, const vec3f &, double) const { return true; }
  virtual bool IntegrateInVolumes() const { return true; }
  virtual bool DeferredHaze() const { return false; }
};

// Track cumulative transport from a fixed ray origin. Taking differences and
// ratios of the same fitted field prevents null collisions or invisible
// boundaries from repeatedly applying the model's finite-distance fitting error.
// A direction change or a crossing into/out of a region with haze disabled
// starts a new atmospheric ray, so excluded distances never enter its integral.
class AtmosphereRay {
public:
  explicit AtmosphereRay(const Atmosphere *model) : model(model) {}
  void Start(const point3f &, const vec3f &, bool integrate);
  AtmosphereSegment Advance(const point3f &);
  // Deferred callers pass throughput/MIS before pending atmospheric extinction.
  // Flush applies that extinction before the next real interaction consumes it;
  // never reset or use eager Advance while a deferred span is pending.
  void Accumulate(const point3f &, const std::array<double, 3> &weight, double uniform);
  // Omitted variates retain exact evaluation for deterministic callers/tests.
  WeightedAtmosphereSegment Flush(double endpoint_uniform = -1, double correction_uniform = -1);
  point3f Origin(const point3f &fallback) const { return active ? origin : fallback; }
  point3f Remaining(const point3f &) const;
private:
  const Atmosphere *model;
  bool initialized = false, active = false;
  point3f origin{0};
  vec3f direction{0};
  double distance = 0;
  AtmosphereSegment cumulative;
  AtmosphereSegmentCache cache;
  // One reservoir stores weight changes along the pending span. The endpoint
  // term is always evaluated; a sampled correction estimates all interior terms.
  bool pending = false;
  double pending_distance = 0, selected_distance = 0, reservoir_mass = 0, selected_mass = 0;
  std::array<double, 3> pending_weight{}, selected_weight{};
};

// The fitted model owns immutable coefficients. No R API is called in a query.
class PragueInfiniteLight final : public InfiniteLight, public Atmosphere {
public:
  explicit PragueInfiniteLight(const Rcpp::List &, bool build_sampler = true);
  point3f Radiance(const point3f &, const vec3f &, Float) const override;
  vec3f Sample(const point3f &, vec2f, Float) const override;
  Float Pdf(const point3f &, const vec3f &, Float) const override;
  double SamplingWeight() const override { return sampling_weight; }
  size_t GetSize() const override;
  const Atmosphere *GetAtmosphere() const override { return this; }
  const Atmosphere *GetTransportAtmosphere() const override { return attenuation ? this : nullptr; }
  AtmosphereSegment Segment(const point3f &, const vec3f &, double) const override;
  AtmosphereSegment Segment(const point3f &, const vec3f &, double,
                            AtmosphereSegmentCache &) const override;
  point3f Transmission(const point3f &, const vec3f &, double) const override;
  point3f CelestialTransmission(const point3f &, const vec3f &, InfiniteLightSpectrum) const override;
  bool MaySeeDisk(const point3f &, const vec3f &, double) const override;
  point3f SkyRadiance(const point3f &, const vec3f &) const override;
  bool IntegrateInVolumes() const override { return haze_in_volumes; }
  bool DeferredHaze() const override { return deferred_haze; }
  bool SampleHazeCorrection() const override { return haze_correction_probability < 1; }
  AtmosphereSegment SampledSegment(const point3f &, const vec3f &, double,
                                   AtmosphereSegmentCache &, double) const override;

private:
  using Model = skymodelr::PragueSkyModel;
  using Vector = Model::Vector3;
  using SpectrumValues = std::array<double, 16>;
  std::shared_ptr<const Model> model;
  Transform light_to_environment, environment_to_light;
  std::array<double, 3> origin, gain;
  double meters_per_unit, altitude, elevation, azimuth, visibility, albedo, intensity;
  double sun_radius, sun_solid_angle, sun_fraction = 0, sampling_weight = 0;
  double haze_correction_probability = 1;
  bool include_sky, include_sun, attenuation = true, query_altitude = true,
       haze_in_volumes = true, deferred_haze = false;
  std::vector<double> wavelengths;
  std::vector<std::array<double, 3>> rgb_weights, transmission_weights;
  std::array<std::vector<std::array<double, 3>>, 2> celestial_weights;
  std::vector<double> proposal_altitudes;
  std::vector<std::unique_ptr<Distribution2D>> proposals;
  Vector Position(const point3f &) const;
  Vector LightingPosition(const point3f &) const;
  Vector Direction(const vec3f &) const;
  vec3f World(const Vector &) const;
  Model::Parameters Parameters(Vector, Vector) const;
  double Height(Vector) const;
  bool PlanetOccludes(Vector, Vector, double distance) const;
  SpectrumValues Spectrum(Vector, Vector, bool sun, bool sky, bool smooth = true) const;
  point3f RGB(const SpectrumValues &) const;
  AtmosphereSegment EvaluateSegment(const point3f &, const vec3f &, double,
                                    AtmosphereSegmentCache &, double) const;
  std::pair<size_t, double> Proposal(Vector) const;
  Float SkyPdf(Vector, Vector) const;
};
#endif
