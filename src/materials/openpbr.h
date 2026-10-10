#ifndef RAYRENDER_OPENPBR_H
#define RAYRENDER_OPENPBR_H

#include "material.h"
#include <Rcpp.h>
#include <memory>
#include "../core/wavefront.h"

// The reference library and its math types stay in one translation unit. These
// types expose only rayrender's conventions: world directions pointing away
// from the surface, f times cosine, and camera-to-light transport weights.
struct OpenPBRSample {
  vec3f direction{0};
  point3f weight{0};
  Float pdf = 0;
  Float eta_squared = 1;
  bool specular = false, transmission = false;
};
struct OpenPBRWavefrontTextures {
  const texture *base;
  const roughness_texture *roughness;
  point2f repeat;
};

class OpenPBRInteraction {
public:
  struct Impl;
  explicit OpenPBRInteraction(std::unique_ptr<Impl> impl);
  ~OpenPBRInteraction();
  OpenPBRInteraction(OpenPBRInteraction&&) noexcept;
  point3f Evaluate(const vec3f& direction) const;
  Float Pdf(const vec3f& direction) const;
  OpenPBRSample Sample(Float branch, Float u, Float v) const;
  point3f Emission() const;
private:
  std::unique_ptr<Impl> impl;
};

// A separate uber shader. Inheriting the small dielectric identity lets the
// existing priority stack track its solid interior without substituting the
// legacy dielectric BSDF for its layered OpenPBR surface.
class OpenPBRMaterial final : public dielectric {
public:
  OpenPBRMaterial(const Rcpp::List& parameters, std::shared_ptr<texture> base,
                  std::shared_ptr<roughness_texture> roughness, bool has_roughness,
                  point2f texture_repeat = point2f(1, 1),
                  std::shared_ptr<const TextureNode> roughness_graph = nullptr);
  ~OpenPBRMaterial() override;
  OpenPBRInteraction Prepare(const Ray&, const hit_record&, Float exterior_ior = 1) const;
  bool is_dielectric() const override;
  bool is_delta_specular() const override { return false; }
  bool physical_normal_mapping() const override { return true; }
  point3f get_albedo(const hit_record&) const override;
  point3f emitted(const Ray&, const hit_record&, Float, Float, const point3f&, bool&) override;
  // OpenPBR uses the NEE integrator's joint sample/evaluate interface, including
  // mixed discrete and continuous lobes; legacy scatter must not silently apply.
  bool scatter(const Ray&, const hit_record&, scatter_record&, random_gen&) override;
  bool scatter(const Ray&, const hit_record&, scatter_record&, Sampler*) override;
  const std::string GetName() override { return "openpbr"; }
  size_t GetSize() override;
  Float EmissionEstimate() const;
  // Cold GPU export. The reference's GLM types and private storage stay here.
  OpenPBRWavefrontTextures ExportWavefront(std::vector<WFVector> &parameters) const;
private:
  struct Impl;
  std::unique_ptr<Impl> impl;
};

struct OpenPBRVolume {
  point3f absorption{0}, scattering{0};
  Float anisotropy = 0, ior = 1;
};
OpenPBRVolume OpenPBRInterior(const Rcpp::List& parameters);
#endif
