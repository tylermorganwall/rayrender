#ifndef RAYRENDER_ATMOSPHERE_H
#define RAYRENDER_ATMOSPHERE_H
#include "infinite.h"
#include <skymodelr/prague.hpp>
#include <array>

struct AtmosphereSegment {
  point3f transmission{1}, radiance{0};
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
  virtual point3f Transmission(const point3f &, const vec3f &, double distance) const = 0;
  virtual point3f CelestialTransmission(const point3f &p, const vec3f &w,
                                        InfiniteLightSpectrum) const {
    return Transmission(p, w, INFINITY);
  }
  virtual point3f SkyRadiance(const point3f &, const vec3f &) const = 0;
  virtual bool MaySeeDisk(const point3f &, const vec3f &, double) const { return true; }
};

// Track cumulative transport from a fixed ray origin. Taking differences and
// ratios of the same fitted field prevents null collisions or invisible
// boundaries from repeatedly applying the model's finite-distance fitting error.
// A change of direction or entry/exit from glass starts a new atmospheric ray.
class AtmosphereRay {
public:
  explicit AtmosphereRay(const Atmosphere *model) : model(model) {}
  void Start(const point3f &, const vec3f &, bool outside_glass);
  AtmosphereSegment Advance(const point3f &);
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
  AtmosphereSegment Segment(const point3f &, const vec3f &, double) const override;
  AtmosphereSegment Segment(const point3f &, const vec3f &, double,
                            AtmosphereSegmentCache &) const override;
  point3f Transmission(const point3f &, const vec3f &, double) const override;
  point3f CelestialTransmission(const point3f &, const vec3f &, InfiniteLightSpectrum) const override;
  bool MaySeeDisk(const point3f &, const vec3f &, double) const override;
  point3f SkyRadiance(const point3f &, const vec3f &) const override;

private:
  using Model = skymodelr::PragueSkyModel;
  using Vector = Model::Vector3;
  using SpectrumValues = std::array<double, 16>;
  std::shared_ptr<const Model> model;
  Transform light_to_environment, environment_to_light;
  std::array<double, 3> origin, gain;
  double meters_per_unit, altitude, elevation, azimuth, visibility, albedo, intensity;
  double sun_radius, sun_solid_angle, sun_fraction = 0, sampling_weight = 0;
  bool include_sky, include_sun;
  std::vector<double> wavelengths;
  std::vector<std::array<double, 3>> rgb_weights, transmission_weights;
  std::array<std::vector<std::array<double, 3>>, 2> celestial_weights;
  std::vector<double> proposal_altitudes;
  std::vector<std::unique_ptr<Distribution2D>> proposals;
  Vector Position(const point3f &) const;
  Vector Direction(const vec3f &) const;
  vec3f World(const Vector &) const;
  Model::Parameters Parameters(Vector, Vector) const;
  double Height(Vector) const;
  bool PlanetOccludes(Vector, Vector, double distance) const;
  SpectrumValues Spectrum(Vector, Vector, bool sun, bool sky) const;
  point3f RGB(const SpectrumValues &) const;
  std::pair<size_t, double> Proposal(Vector) const;
  Float SkyPdf(Vector, Vector) const;
};
#endif
