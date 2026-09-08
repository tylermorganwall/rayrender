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

private:
  std::vector<std::shared_ptr<InfiniteLight>> lights;
  std::unique_ptr<Distribution1D> selection;
  double weight = 0;
};

// A rectilinear image spanning a circular cone. Its texels are never resampled
// into the environment dome. Uniform solid-angle sampling remains efficient for
// small celestial disks and gives support to every phase/texture detail.
class DiskInfiniteLight final : public InfiniteLight {
public:
  DiskInfiniteLight(std::shared_ptr<texture> image, int width, int height,
                    vec3f direction, double angular_diameter, Float rotation,
                    bool clip_horizon = true);
  point3f Radiance(const point3f &, const vec3f &, Float) const override;
  vec3f Sample(const point3f &, vec2f, Float) const override;
  Float Pdf(const point3f &, const vec3f &, Float) const override;
  double SamplingWeight() const override { return weight; }
  size_t GetSize() const override { return sizeof(*this); }

private:
  bool Coordinates(const vec3f &, Float &u, Float &v) const;
  std::shared_ptr<texture> image;
  vec3<double> forward, right, up;
  double tan_radius, one_minus_cos_radius, solid_angle, weight = 0;
  bool clip_horizon;
};

std::shared_ptr<InfiniteLight> BuildInfiniteLights(const Rcpp::List &descriptions,
                                                      TextureCache &textures);
#endif
