#ifndef RAYRENDER_LIGHTS_INFINITE_H
#define RAYRENDER_LIGHTS_INFINITE_H

#include "../core/ray.h"
#include "../math/distributions.h"
#include "../math/transform.h"
#include "../materials/texture.h"
#include <Rcpp.h>
#include <memory>
#include <vector>

class TextureCache;
class Atmosphere;
enum class InfiniteLightSpectrum { RGB, Sun, Moon };

// Incident radiance and a directional proposal are separate operations. The
// position/time context permits future spatial atmosphere implementations.
class InfiniteLight {
public:
  virtual ~InfiniteLight() = default;
  virtual point3f Radiance(const point3f &p, const vec3f &wi, Float time) const = 0;
  virtual vec3f Sample(const point3f &p, vec2f u, Float time) const = 0;
  virtual Float Pdf(const point3f &p, const vec3f &wi, Float time) const = 0;
  virtual double SamplingWeight() const = 0;
  virtual size_t GetSize() const = 0;
  virtual const Atmosphere *GetAtmosphere() const { return nullptr; }
  virtual InfiniteLightSpectrum RadianceSpectrum() const { return InfiniteLightSpectrum::RGB; }
  virtual bool Available(const point3f &, const Atmosphere *) const { return true; }
  // Borrow the preview's mutable environment transform. Positions stay in world
  // space; only directions and atmosphere-relative offsets use this rotation.
  virtual void SetEnvironmentTransform(const Transform *to_world, const Transform *to_local) {
    environment_to_world = to_world;
    world_to_environment = to_local;
  }
protected:
  vec3f EnvironmentDirection(const vec3f &v) const {
    return world_to_environment ? (*world_to_environment)(v) : v;
  }
  vec3f WorldDirection(const vec3f &v) const {
    return environment_to_world ? (*environment_to_world)(v) : v;
  }
  const Transform *environment_to_world = nullptr, *world_to_environment = nullptr;
};

class ImageInfiniteLight final : public InfiniteLight {
public:
  ImageInfiniteLight(std::shared_ptr<texture> image, int width, int height,
                     Float rotation);
  point3f Radiance(const point3f &, const vec3f &, Float) const override;
  vec3f Sample(const point3f &, vec2f, Float) const override;
  Float Pdf(const point3f &, const vec3f &, Float) const override;
  double SamplingWeight() const override { return weight; }
  size_t GetSize() const override;

private:
  std::shared_ptr<texture> image;
  std::unique_ptr<Distribution2D> distribution;
  Transform light_to_world, world_to_light;
  double weight = 0;
};

class InfiniteLightMixture final : public InfiniteLight {
public:
  explicit InfiniteLightMixture(std::vector<std::shared_ptr<InfiniteLight>> lights);
  point3f Radiance(const point3f &, const vec3f &, Float) const override;
  vec3f Sample(const point3f &, vec2f, Float) const override;
  Float Pdf(const point3f &, const vec3f &, Float) const override;
  double SamplingWeight() const override { return weight; }
  size_t GetSize() const override;
  const Atmosphere *GetAtmosphere() const override;
  void SetEnvironmentTransform(const Transform *, const Transform *) override;

private:
  std::vector<std::shared_ptr<InfiniteLight>> lights;
  std::unique_ptr<Distribution1D> selection;
  double weight = 0;
  double ChoiceWeight(size_t index, const point3f &, const Atmosphere *) const;
  double ChoiceTotal(const point3f &, const Atmosphere *) const;
};

// A rectilinear image spanning a circular cone. Its texels are never resampled
// into the environment dome. Uniform solid-angle sampling remains efficient for
// small celestial disks and gives support to every phase/texture detail.
class DiskInfiniteLight final : public InfiniteLight {
public:
  DiskInfiniteLight(std::shared_ptr<texture> image, int width, int height,
                    vec3f direction, double angular_diameter, Float rotation,
                    bool clip_horizon = true,
                    InfiniteLightSpectrum spectrum = InfiniteLightSpectrum::RGB);
  point3f Radiance(const point3f &, const vec3f &, Float) const override;
  vec3f Sample(const point3f &, vec2f, Float) const override;
  Float Pdf(const point3f &, const vec3f &, Float) const override;
  double SamplingWeight() const override { return weight; }
  size_t GetSize() const override { return sizeof(*this); }
  InfiniteLightSpectrum RadianceSpectrum() const override { return spectrum; }
  bool Available(const point3f &, const Atmosphere *) const override;

private:
  bool Coordinates(const vec3f &, Float &u, Float &v) const;
  std::shared_ptr<texture> image;
  vec3<double> forward, right, up;
  double tan_radius, one_minus_cos_radius, solid_angle, weight = 0;
  bool clip_horizon;
  InfiniteLightSpectrum spectrum;
};

std::shared_ptr<InfiniteLight> BuildInfiniteLights(const Rcpp::List &descriptions,
                                                      TextureCache &textures);
#endif
