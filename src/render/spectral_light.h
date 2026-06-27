#ifndef RAYRENDER_RENDER_SPECTRAL_LIGHT_H
#define RAYRENDER_RENDER_SPECTRAL_LIGHT_H

#include "spectral_scene.h"

#include "../base/base.h"
#include "../core/ray.h"
#include "../math/vectypes.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace rayrender {
namespace render {

enum class LightType {
  Point,
  Spot,
  Distant,
  DiffuseArea,
  UniformInfinite,
  ImageInfinite
};

enum class LightFlags : std::uint32_t {
  None = 0,
  DeltaPosition = 1u << 0,
  DeltaDirection = 1u << 1,
  Area = 1u << 2,
  Infinite = 1u << 3
};

LightFlags operator|(LightFlags lhs, LightFlags rhs);
LightFlags operator&(LightFlags lhs, LightFlags rhs);
bool HasFlag(LightFlags flags, LightFlags flag);
const char* LightTypeName(LightType type);
bool IsDeltaLight(LightFlags flags);

enum class LightSampleMode {
  CompletePDF,
  AllowIncompletePDF
};

class LightSpectrum {
public:
  LightSpectrum();
  explicit LightSpectrum(base::Spectrum spectrum);

  static LightSpectrum Constant(Float value);
  static LightSpectrum FromRGBIlluminant(
    const base::RGBColorSpace& colorSpace,
    const base::RGB& rgb
  );

  bool IsValid() const;
  base::SampledSpectrum Sample(const base::SampledWavelengths& lambda) const;
  Float MaxValue() const;

private:
  struct NeutralIlluminant {
    std::shared_ptr<const base::DenselySampledSpectrum> illuminant;
    Float scale = 0;
  };

  std::variant<std::monostate, base::Spectrum, base::RGBIlluminantSpectrum, NeutralIlluminant> value_;
};

struct LightSampleContext {
  point3f p;
  normal3f n;
  normal3f ns;
  Float time = 0;
  base::MediumHandle medium = base::MediumHandle::Invalid();
  bool hasMedium = false;
  MediumInterface mediumInterface;
  bool hasMediumInterface = false;

  LightSampleContext() = default;
  explicit LightSampleContext(const Interaction& interaction);

  Interaction ToInteraction() const;
};

struct VisibilityTester {
  Interaction p0;
  Interaction p1;

  SpawnedRay SpawnRay() const;
  bool Unoccluded(const Aggregate& aggregate) const;
};

struct LightLiSample {
  base::SampledSpectrum L;
  vec3f wi;
  Float pdf = 0;
  Interaction pLight;
  VisibilityTester visibility;
  bool delta = false;
};

struct DiffuseAreaLightOptions {
  bool twoSided = false;
  bool allowUnsampledEmitter = false;
};

struct AreaLightBindingOptions {
  bool allowUnsampledEmitter = false;
};

struct LegacyEmissiveMaterialDescriptor {
  base::RGB color = base::RGB(1);
  Float intensity = 1;
  bool invisible = false;
};

struct LegacySpotLightMaterialDescriptor {
  point3f position;
  vec3f direction = vec3f(0, 0, -1);
  base::RGB color = base::RGB(1);
  Float intensity = 1;
  Float cosTotalWidth = static_cast<Float>(0.5);
  Float cosFalloffStart = static_cast<Float>(0.75);
  bool invisible = false;
};

class Light {
public:
  struct Concept;

  Light() = default;

  static Light Point(point3f position, LightSpectrum intensity);
  static Light Spot(
    point3f position,
    vec3f direction,
    LightSpectrum intensity,
    Float totalWidthDegrees,
    Float falloffStartDegrees
  );
  static Light SpotFromCosines(
    point3f position,
    vec3f direction,
    LightSpectrum intensity,
    Float cosTotalWidth,
    Float cosFalloffStart
  );
  static Light Distant(vec3f directionToLight, LightSpectrum radiance);
  static Light DiffuseArea(
    Shape shape,
    LightSpectrum radiance,
    DiffuseAreaLightOptions options = {}
  );
  static Light UniformInfinite(
    LightSpectrum radiance,
    Transform3f renderFromLight = Transform3f::Identity()
  );
  static Light ImageInfinite(
    int width,
    int height,
    std::vector<LightSpectrum> texels,
    Transform3f renderFromLight = Transform3f::Identity(),
    Float scale = 1
  );
  static Light ImageInfiniteFromRGB(
    int width,
    int height,
    const std::vector<base::RGB>& texels,
    const base::RGBColorSpace& colorSpace,
    Transform3f renderFromLight = Transform3f::Identity(),
    Float scale = 1
  );

  bool IsValid() const;
  LightType Type() const;
  LightFlags Flags() const;
  const std::string& Name() const;

  std::optional<LightLiSample> SampleLi(
    const LightSampleContext& ctx,
    point2f u,
    const base::SampledWavelengths& lambda,
    LightSampleMode mode = LightSampleMode::CompletePDF
  ) const;
  Float PDF_Li(
    const LightSampleContext& ctx,
    const vec3f& wi,
    LightSampleMode mode = LightSampleMode::CompletePDF
  ) const;
  base::SampledSpectrum Le(const Ray& ray, const base::SampledWavelengths& lambda) const;
  base::SampledSpectrum L(
    const Interaction& interaction,
    const vec3f& w,
    const base::SampledWavelengths& lambda
  ) const;
  base::SampledSpectrum Phi(const base::SampledWavelengths& lambda) const;
  void Preprocess(const Bounds3f& sceneBounds);
  Float PowerEstimate() const;

private:
  explicit Light(std::shared_ptr<Concept> impl);

  std::shared_ptr<Concept> impl_;
};

class SpectralLightTable {
public:
  base::LightHandle Add(Light light);
  const Light& Get(base::LightHandle handle) const;
  Light& Get(base::LightHandle handle);
  std::size_t Size() const;
  void Preprocess(const Bounds3f& sceneBounds);

private:
  void Validate(base::LightHandle handle) const;

  std::vector<Light> lights_;
};

struct SampledLight {
  base::LightHandle handle = base::LightHandle::Invalid();
  const Light* light = nullptr;
  Float pmf = 0;
};

class UniformLightSampler {
public:
  UniformLightSampler() = default;
  UniformLightSampler(const SpectralLightTable& table, std::vector<base::LightHandle> lights);

  std::optional<SampledLight> Sample(Float u) const;
  Float PMF(base::LightHandle light) const;
  Float PMFSum() const;
  std::size_t Size() const;

private:
  const SpectralLightTable* table_ = nullptr;
  std::vector<base::LightHandle> lights_;
};

class PowerLightSampler {
public:
  PowerLightSampler() = default;
  PowerLightSampler(const SpectralLightTable& table, std::vector<base::LightHandle> lights);

  std::optional<SampledLight> Sample(Float u) const;
  Float PMF(base::LightHandle light) const;
  Float PMFSum() const;
  std::size_t Size() const;

private:
  const SpectralLightTable* table_ = nullptr;
  std::vector<base::LightHandle> lights_;
  std::vector<Float> pmf_;
  std::vector<Float> cdf_;
};

PrimitiveBinding CompileAreaLightPrimitiveBinding(
  base::MaterialHandle material,
  base::LightHandle areaLight,
  const Shape& shape,
  AreaLightBindingOptions options = {}
);

Light ConvertLegacyEmissiveMaterialToAreaLight(
  const LegacyEmissiveMaterialDescriptor& legacy,
  const Shape& shape,
  const base::RGBColorSpace& colorSpace,
  DiffuseAreaLightOptions options = {}
);

Light ConvertLegacySpotLightMaterialToLight(
  const LegacySpotLightMaterialDescriptor& legacy,
  const base::RGBColorSpace& colorSpace
);

base::LightHandle RegisterLight(Scene& scene, SpectralLightTable& table, Light light);
void AttachLightToScene(Scene& scene, base::LightHandle handle, const Light& light);
bool ShouldCompileLegacyEnvironmentSphereToSpectralGeometry();
bool SpectralMaterialEmissionPathEnabled();

} // namespace render
} // namespace rayrender

#endif
