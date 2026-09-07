#ifndef RAYRENDER_MEDIUM_H
#define RAYRENDER_MEDIUM_H
#include "../core/ray.h"
#include "../math/onbh.h"
#include "../math/transform.h"
#include <array>
#include <memory>
#include <optional>
#include <vector>

class texture;
struct PhaseFunctionSample {
  vec3f wi;
  Float p, pdf;
};
class HGPhaseFunction {
public:
  explicit HGPhaseFunction(Float g = 0) : g(g) {}
  Float p(const vec3f &wo, const vec3f &wi) const;
  PhaseFunctionSample Sample(const vec3f &wo, Float u, Float v) const;
  Float g;
};
struct MediumProperties {
  point3f sigma_a{0}, sigma_s{0}, Le{0};
  HGPhaseFunction phase;
};
struct MediumInteraction {
  point3f p;
  vec3f wo;
  Float time;
  MediumProperties properties;
};
struct RayMajorantSegment {
  double t_min, t_max;
  point3f sigma_maj;
};
struct MajorantGrid {
  point3f lo{-0.5f}, hi{0.5f};
  int resolution = 16;
  std::vector<Float> density;
  Float Get(int x, int y, int z) const { return density[x + resolution * (y + resolution * z)]; }
};
// Input direction is transformed, but t always measures world-space distance.
class RayMajorantIterator {
public:
  RayMajorantIterator(const Ray &medium_ray, double t_max, const point3f &sigma_t,
                      const MajorantGrid *grid);
  std::optional<RayMajorantSegment> Next();

private:
  const MajorantGrid *grid;
  point3f sigma_t;
  double t = 0, end = 0;
  std::array<int, 3> cell{}, step{};
  std::array<double, 3> next{}, delta{};
  bool finished = false;
};

struct SampledField {
  std::array<int, 3> dims{1, 1, 1};
  int channels = 1;
  std::vector<Float> values;
  Float Lookup(const point3f &normalized, int channel = 0) const;
};
class Medium {
public:
  explicit Medium(const Rcpp::List &description);
  virtual ~Medium() = default;
  virtual MediumProperties SamplePoint(const point3f &p) const;
  virtual RayMajorantIterator SampleRay(const Ray &r, double t_max) const;
  virtual bool IsEmissive() const;
  virtual size_t MemoryBytes() const { return sizeof(*this); }
  virtual bool IsHomogeneous() const { return true; }
  virtual Float Density(const point3f &) const { return 1; }
  virtual point3f Emission(const point3f &) const;
  point3f sigma_a, sigma_s, emission;
  Float g, emission_scale, temperature_scale, temperature_offset;
  bool has_temperature = false, has_rgb_emission = false;
  Float temperature = 0;
  Transform medium_to_object;
  std::shared_ptr<texture> legacy_albedo;

protected:
  MediumProperties Properties(Float density, const point3f &le, const point3f &p) const;
};
class GridMedium final : public Medium {
public:
  explicit GridMedium(const Rcpp::List &description);
  size_t MemoryBytes() const override;
  MediumProperties SamplePoint(const point3f &p) const override;
  RayMajorantIterator SampleRay(const Ray &r, double t_max) const override;
  bool IsHomogeneous() const override { return false; }
  Float Density(const point3f &p) const override;
  point3f Emission(const point3f &p) const override;

private:
  SampledField density, temperatures, emissions;
  MajorantGrid majorants;
  point3f Normalize(const point3f &p) const;
};
class NanoVDBMedium final : public Medium {
public:
  explicit NanoVDBMedium(const Rcpp::List &description);
  ~NanoVDBMedium();
  size_t MemoryBytes() const override;
  MediumProperties SamplePoint(const point3f &p) const override;
  RayMajorantIterator SampleRay(const Ray &r, double t_max) const override;
  bool IsHomogeneous() const override { return false; }
  Float Density(const point3f &p) const override;
  point3f Emission(const point3f &p) const override;

private:
  struct Data;
  std::unique_ptr<Data> data;
  MajorantGrid majorants;
};
std::shared_ptr<const Medium> LoadMedium(const Rcpp::List &description);
point3f BlackbodyRGB(Float kelvin);
#endif
