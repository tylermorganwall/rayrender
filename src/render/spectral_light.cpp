#include "spectral_light.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace rayrender {
namespace render {

namespace {

constexpr Float Pi = static_cast<Float>(3.14159265358979323846264338327950288);
constexpr Float Inv4Pi = static_cast<Float>(1) / (static_cast<Float>(4) * Pi);
constexpr Float ShadowEpsilon = static_cast<Float>(1e-4);

Float Clamp(Float value, Float low, Float high) {
  return std::max(low, std::min(value, high));
}

Float SafeSqrt(Float value) {
  return std::sqrt(std::max(static_cast<Float>(0), value));
}

Float Radians(Float degrees) {
  return degrees * Pi / static_cast<Float>(180);
}

Float DistanceSquared(const point3f& a, const point3f& b) {
  return (a - b).squared_length();
}

vec3f UnitOrZero(const vec3f& value) {
  Float length = value.length();
  if (!(length > 0) || !std::isfinite(length)) {
    return vec3f(0, 0, 0);
  }
  return value / length;
}

normal3f UnitOrZero(const normal3f& value) {
  Float length = value.length();
  if (!(length > 0) || !std::isfinite(length)) {
    return normal3f(0, 0, 0);
  }
  return value / length;
}

bool NearlyNeutral(const base::RGB& rgb) {
  constexpr Float tolerance = static_cast<Float>(1e-6);
  return std::fabs(rgb.r - rgb.g) <= tolerance && std::fabs(rgb.r - rgb.b) <= tolerance;
}

void ValidateFiniteNonNegative(Float value, const char* context) {
  if (!std::isfinite(value) || value < 0) {
    throw std::invalid_argument(std::string(context) + " must be finite and non-negative");
  }
}

void ValidateLightRGB(const base::RGB& rgb, const char* context) {
  if (!rgb.IsFinite() || rgb.r < 0 || rgb.g < 0 || rgb.b < 0) {
    throw std::invalid_argument(std::string(context) + " requires finite non-negative RGB values");
  }
}

Float SceneRadius(const Bounds3f& bounds) {
  if (!bounds.IsValid()) {
    return 1;
  }
  Float radius = static_cast<Float>(0.5) * (bounds.pMax - bounds.pMin).length();
  if (!(radius > 0) || !std::isfinite(radius)) {
    return 1;
  }
  return radius;
}

point3f SceneCenter(const Bounds3f& bounds) {
  if (!bounds.IsValid()) {
    return point3f(0, 0, 0);
  }
  return point3f(
    static_cast<Float>(0.5) * (bounds.pMin.xyz.x + bounds.pMax.xyz.x),
    static_cast<Float>(0.5) * (bounds.pMin.xyz.y + bounds.pMax.xyz.y),
    static_cast<Float>(0.5) * (bounds.pMin.xyz.z + bounds.pMax.xyz.z)
  );
}

vec3f SampleUniformSphere(point2f u) {
  Float z = static_cast<Float>(1) - static_cast<Float>(2) * Clamp(u[0], 0, 1);
  Float r = SafeSqrt(static_cast<Float>(1) - z * z);
  Float phi = static_cast<Float>(2) * Pi * Clamp(u[1], 0, 1);
  return vec3f(r * std::cos(phi), r * std::sin(phi), z);
}

Float UniformSpherePDF() {
  return Inv4Pi;
}

Float SmoothStep(Float x, Float edge0, Float edge1) {
  if (edge0 == edge1) {
    return x < edge0 ? 0 : 1;
  }
  Float t = Clamp((x - edge0) / (edge1 - edge0), 0, 1);
  return t * t * (static_cast<Float>(3) - static_cast<Float>(2) * t);
}

Float EvaluatePolynomial(Float x, Float c0, Float c1, Float c2, Float c3, Float c4, Float c5, Float c6) {
  return ((((((c6 * x + c5) * x + c4) * x + c3) * x + c2) * x + c1) * x + c0);
}

vec3f PbrtEqualAreaSquareToSphere(point2f p) {
  Float u = static_cast<Float>(2) * Clamp(p[0], 0, 1) - static_cast<Float>(1);
  Float v = static_cast<Float>(2) * Clamp(p[1], 0, 1) - static_cast<Float>(1);
  Float up = std::fabs(u);
  Float vp = std::fabs(v);
  Float signedDistance = static_cast<Float>(1) - (up + vp);
  Float d = std::fabs(signedDistance);
  Float r = static_cast<Float>(1) - d;
  Float phi = (r == 0 ? static_cast<Float>(1) : (vp - up) / r + static_cast<Float>(1)) * Pi /
              static_cast<Float>(4);
  Float z = std::copysign(static_cast<Float>(1) - r * r, signedDistance);
  Float scale = r * SafeSqrt(static_cast<Float>(2) - r * r);
  Float cosPhi = std::copysign(std::cos(phi), u);
  Float sinPhi = std::copysign(std::sin(phi), v);
  return vec3f(cosPhi * scale, sinPhi * scale, z);
}

point2f PbrtEqualAreaSphereToSquare(vec3f direction) {
  vec3f d = UnitOrZero(direction);
  Float x = std::fabs(d.xyz.x);
  Float y = std::fabs(d.xyz.y);
  Float z = std::fabs(d.xyz.z);

  Float r = SafeSqrt(static_cast<Float>(1) - z);
  Float a = std::max(x, y);
  Float b = std::min(x, y);
  b = a == 0 ? 0 : b / a;

  Float phi = EvaluatePolynomial(
    b,
    static_cast<Float>(0.406758566246788489601959989e-5),
    static_cast<Float>(0.636226545274016134946890922156),
    static_cast<Float>(0.61572017898280213493197203466e-2),
    static_cast<Float>(-0.247333733281268944196501420480),
    static_cast<Float>(0.881770664775316294736387951347e-1),
    static_cast<Float>(0.419038818029165735901852432784e-1),
    static_cast<Float>(-0.251390972343483509333252996350e-1)
  );
  if (x < y) {
    phi = static_cast<Float>(1) - phi;
  }

  Float v = phi * r;
  Float u = r - v;

  if (d.xyz.z < 0) {
    std::swap(u, v);
    u = static_cast<Float>(1) - u;
    v = static_cast<Float>(1) - v;
  }

  u = std::copysign(u, d.xyz.x);
  v = std::copysign(v, d.xyz.y);
  return point2f(
    static_cast<Float>(0.5) * (u + static_cast<Float>(1)),
    static_cast<Float>(0.5) * (v + static_cast<Float>(1))
  );
}

struct PiecewiseSample2D {
  point2f uv;
  Float pdf = 0;
};

class PiecewiseConstant2D {
public:
  PiecewiseConstant2D() = default;

  PiecewiseConstant2D(int width, int height, std::vector<Float> values)
    : width_(width), height_(height), values_(std::move(values)) {
    if (width_ <= 0 || height_ <= 0) {
      throw std::invalid_argument("PiecewiseConstant2D requires positive dimensions");
    }
    if (values_.size() != static_cast<std::size_t>(width_ * height_)) {
      throw std::invalid_argument("PiecewiseConstant2D value count does not match dimensions");
    }
    for (Float& value : values_) {
      if (!std::isfinite(value) || value < 0) {
        throw std::invalid_argument("PiecewiseConstant2D values must be finite and non-negative");
      }
    }
    sum_ = std::accumulate(values_.begin(), values_.end(), static_cast<Float>(0));
    if (!(sum_ > 0)) {
      std::fill(values_.begin(), values_.end(), static_cast<Float>(1));
      sum_ = static_cast<Float>(values_.size());
    }
    cdf_.reserve(values_.size());
    Float running = 0;
    for (Float value : values_) {
      running += value / sum_;
      cdf_.push_back(running);
    }
    if (!cdf_.empty()) {
      cdf_.back() = 1;
    }
  }

  std::optional<PiecewiseSample2D> Sample(point2f u) const {
    if (cdf_.empty()) {
      return std::nullopt;
    }
    Float selection = Clamp(u[0], 0, std::nextafter(static_cast<Float>(1), static_cast<Float>(0)));
    auto iter = std::lower_bound(cdf_.begin(), cdf_.end(), selection);
    std::size_t index = static_cast<std::size_t>(iter - cdf_.begin());
    if (index >= values_.size()) {
      index = values_.size() - 1;
    }
    int x = static_cast<int>(index % static_cast<std::size_t>(width_));
    int y = static_cast<int>(index / static_cast<std::size_t>(width_));
    Float previousCDF = index == 0 ? 0 : cdf_[index - 1];
    Float cellProbability = cdf_[index] - previousCDF;
    Float localX = cellProbability > 0 ? (selection - previousCDF) / cellProbability : static_cast<Float>(0.5);
    Float localY = Clamp(u[1], 0, std::nextafter(static_cast<Float>(1), static_cast<Float>(0)));
    point2f uv(
      (static_cast<Float>(x) + localX) / static_cast<Float>(width_),
      (static_cast<Float>(y) + localY) / static_cast<Float>(height_)
    );
    return PiecewiseSample2D{uv, PDF(uv)};
  }

  Float PDF(point2f uv) const {
    if (values_.empty() || !(sum_ > 0)) {
      return 0;
    }
    int x = std::clamp(static_cast<int>(std::floor(Clamp(uv[0], 0, std::nextafter(static_cast<Float>(1), static_cast<Float>(0))) * width_)), 0, width_ - 1);
    int y = std::clamp(static_cast<int>(std::floor(Clamp(uv[1], 0, std::nextafter(static_cast<Float>(1), static_cast<Float>(0))) * height_)), 0, height_ - 1);
    Float value = values_[static_cast<std::size_t>(y * width_ + x)];
    return value * static_cast<Float>(width_ * height_) / sum_;
  }

private:
  int width_ = 0;
  int height_ = 0;
  std::vector<Float> values_;
  std::vector<Float> cdf_;
  Float sum_ = 0;
};

Interaction MakeEndpoint(point3f p, normal3f n, Float time = 0) {
  Interaction interaction;
  interaction.p = p;
  interaction.n = UnitOrZero(n);
  interaction.pError = vec3f(ShadowEpsilon);
  interaction.time = time;
  return interaction;
}

} // namespace

LightFlags operator|(LightFlags lhs, LightFlags rhs) {
  return static_cast<LightFlags>(
    static_cast<std::uint32_t>(lhs) | static_cast<std::uint32_t>(rhs)
  );
}

LightFlags operator&(LightFlags lhs, LightFlags rhs) {
  return static_cast<LightFlags>(
    static_cast<std::uint32_t>(lhs) & static_cast<std::uint32_t>(rhs)
  );
}

bool HasFlag(LightFlags flags, LightFlags flag) {
  return (static_cast<std::uint32_t>(flags) & static_cast<std::uint32_t>(flag)) != 0;
}

const char* LightTypeName(LightType type) {
  switch (type) {
  case LightType::Point:
    return "Point";
  case LightType::Spot:
    return "Spot";
  case LightType::Distant:
    return "Distant";
  case LightType::DiffuseArea:
    return "DiffuseArea";
  case LightType::UniformInfinite:
    return "UniformInfinite";
  case LightType::ImageInfinite:
    return "ImageInfinite";
  }
  return "Unknown";
}

bool IsDeltaLight(LightFlags flags) {
  return HasFlag(flags, LightFlags::DeltaPosition) || HasFlag(flags, LightFlags::DeltaDirection);
}

LightSpectrum::LightSpectrum() = default;

LightSpectrum::LightSpectrum(base::Spectrum spectrum) : value_(std::move(spectrum)) {}

LightSpectrum LightSpectrum::Constant(Float value) {
  ValidateFiniteNonNegative(value, "LightSpectrum constant");
  return LightSpectrum(base::Spectrum(base::ConstantSpectrum(value)));
}

LightSpectrum LightSpectrum::FromRGBIlluminant(
  const base::RGBColorSpace& colorSpace,
  const base::RGB& rgb
) {
  ValidateLightRGB(rgb, "LightSpectrum RGB illuminant");
  if (NearlyNeutral(rgb) && !colorSpace.HasSpectralReconstruction()) {
    if (!colorSpace.illuminant) {
      throw std::invalid_argument("RGB illuminant light reconstruction requires a color-space illuminant");
    }
    LightSpectrum result;
    result.value_ = NeutralIlluminant{colorSpace.illuminant, rgb.r};
    return result;
  }
  LightSpectrum result;
  result.value_ = base::RGBIlluminantSpectrum(colorSpace, rgb);
  return result;
}

bool LightSpectrum::IsValid() const {
  return !std::holds_alternative<std::monostate>(value_);
}

base::SampledSpectrum LightSpectrum::Sample(const base::SampledWavelengths& lambda) const {
  return std::visit(
    [&lambda](const auto& value) -> base::SampledSpectrum {
      using T = std::decay_t<decltype(value)>;
      if constexpr (std::is_same<T, std::monostate>::value) {
        return base::SampledSpectrum(0);
      } else if constexpr (std::is_same<T, base::Spectrum>::value) {
        return value.Sample(lambda);
      } else if constexpr (std::is_same<T, base::RGBIlluminantSpectrum>::value) {
        return value.Sample(lambda);
      } else {
        if (!value.illuminant) {
          return base::SampledSpectrum(0);
        }
        return value.illuminant->Sample(lambda) * value.scale;
      }
    },
    value_
  );
}

Float LightSpectrum::MaxValue() const {
  return std::visit(
    [](const auto& value) -> Float {
      using T = std::decay_t<decltype(value)>;
      if constexpr (std::is_same<T, std::monostate>::value) {
        return 0;
      } else if constexpr (std::is_same<T, base::Spectrum>::value) {
        return std::max(static_cast<Float>(0), value.MaxValue());
      } else if constexpr (std::is_same<T, base::RGBIlluminantSpectrum>::value) {
        return std::max(static_cast<Float>(0), value.MaxValue());
      } else {
        if (!value.illuminant) {
          return 0;
        }
        return std::max(static_cast<Float>(0), value.illuminant->MaxValue() * value.scale);
      }
    },
    value_
  );
}

LightSampleContext::LightSampleContext(const Interaction& interaction)
  : p(interaction.p),
    n(interaction.n),
    ns(interaction.n),
    time(interaction.time),
    medium(interaction.medium),
    hasMedium(interaction.hasMedium),
    mediumInterface(interaction.mediumInterface),
    hasMediumInterface(interaction.hasMediumInterface) {}

Interaction LightSampleContext::ToInteraction() const {
  Interaction interaction;
  interaction.p = p;
  interaction.n = n;
  interaction.pError = vec3f(ShadowEpsilon);
  interaction.time = time;
  interaction.medium = medium;
  interaction.hasMedium = hasMedium;
  interaction.mediumInterface = mediumInterface;
  interaction.hasMediumInterface = hasMediumInterface;
  return interaction;
}

SpawnedRay VisibilityTester::SpawnRay() const {
  SpawnedRay ray = p0.SpawnRayTo(p1.p);
  ray.ray.tMax = static_cast<Float>(1) - ShadowEpsilon;
  return ray;
}

bool VisibilityTester::Unoccluded(const Aggregate& aggregate) const {
  SpawnedRay ray = SpawnRay();
  return !aggregate.Intersect(ray.ray, ShadowEpsilon, ray.ray.tMax).has_value();
}

struct Light::Concept {
  virtual ~Concept() = default;
  virtual LightType Type() const = 0;
  virtual LightFlags Flags() const = 0;
  virtual const std::string& Name() const = 0;
  virtual std::optional<LightLiSample> SampleLi(
    const LightSampleContext& ctx,
    point2f u,
    const base::SampledWavelengths& lambda,
    LightSampleMode mode
  ) const = 0;
  virtual Float PDF_Li(const LightSampleContext& ctx, const vec3f& wi, LightSampleMode mode) const = 0;
  virtual base::SampledSpectrum Le(const Ray& ray, const base::SampledWavelengths& lambda) const = 0;
  virtual base::SampledSpectrum L(
    const Interaction& interaction,
    const vec3f& w,
    const base::SampledWavelengths& lambda
  ) const = 0;
  virtual base::SampledSpectrum Phi(const base::SampledWavelengths& lambda) const = 0;
  virtual void Preprocess(const Bounds3f& sceneBounds) = 0;
  virtual Float PowerEstimate() const = 0;
};

class PointLightImpl : public Light::Concept {
public:
  PointLightImpl(point3f position, LightSpectrum intensity)
    : position_(position), intensity_(std::move(intensity)) {}

  LightType Type() const override { return LightType::Point; }
  LightFlags Flags() const override { return LightFlags::DeltaPosition; }
  const std::string& Name() const override { return name_; }

  std::optional<LightLiSample> SampleLi(
    const LightSampleContext& ctx,
    point2f,
    const base::SampledWavelengths& lambda,
    LightSampleMode
  ) const override {
    vec3f toLight = position_ - ctx.p;
    Float distance2 = toLight.squared_length();
    if (!(distance2 > 0)) {
      return std::nullopt;
    }
    vec3f wi = toLight / std::sqrt(distance2);
    Interaction pLight = MakeEndpoint(position_, normal3f(-wi.xyz.x, -wi.xyz.y, -wi.xyz.z), ctx.time);
    LightLiSample sample;
    sample.L = intensity_.Sample(lambda) / distance2;
    sample.wi = wi;
    sample.pdf = 1;
    sample.pLight = pLight;
    sample.visibility = VisibilityTester{ctx.ToInteraction(), pLight};
    sample.delta = true;
    return sample;
  }

  Float PDF_Li(const LightSampleContext&, const vec3f&, LightSampleMode) const override {
    return 0;
  }

  base::SampledSpectrum Le(const Ray&, const base::SampledWavelengths&) const override {
    return base::SampledSpectrum(0);
  }

  base::SampledSpectrum L(
    const Interaction&,
    const vec3f&,
    const base::SampledWavelengths&
  ) const override {
    return base::SampledSpectrum(0);
  }

  base::SampledSpectrum Phi(const base::SampledWavelengths& lambda) const override {
    return intensity_.Sample(lambda) * (static_cast<Float>(4) * Pi);
  }

  void Preprocess(const Bounds3f&) override {}

  Float PowerEstimate() const override {
    return static_cast<Float>(4) * Pi * intensity_.MaxValue();
  }

private:
  point3f position_;
  LightSpectrum intensity_;
  std::string name_ = "PointLight";
};

class SpotLightImpl : public Light::Concept {
public:
  SpotLightImpl(
    point3f position,
    vec3f direction,
    LightSpectrum intensity,
    Float cosTotalWidth,
    Float cosFalloffStart
  )
    : position_(position),
      direction_(UnitOrZero(direction)),
      intensity_(std::move(intensity)),
      cosTotalWidth_(cosTotalWidth),
      cosFalloffStart_(cosFalloffStart) {
    if (!(direction_.squared_length() > 0)) {
      throw std::invalid_argument("SpotLight direction must be non-zero");
    }
    if (cosFalloffStart_ < cosTotalWidth_) {
      throw std::invalid_argument("SpotLight falloff start angle must be inside total width");
    }
  }

  LightType Type() const override { return LightType::Spot; }
  LightFlags Flags() const override { return LightFlags::DeltaPosition; }
  const std::string& Name() const override { return name_; }

  std::optional<LightLiSample> SampleLi(
    const LightSampleContext& ctx,
    point2f,
    const base::SampledWavelengths& lambda,
    LightSampleMode
  ) const override {
    vec3f toLight = position_ - ctx.p;
    Float distance2 = toLight.squared_length();
    if (!(distance2 > 0)) {
      return std::nullopt;
    }
    vec3f wi = toLight / std::sqrt(distance2);
    vec3f fromLight = -wi;
    base::SampledSpectrum radiance = intensity_.Sample(lambda) * Falloff(fromLight) / distance2;
    if (!radiance) {
      return std::nullopt;
    }
    Interaction pLight = MakeEndpoint(position_, normal3f(direction_.xyz.x, direction_.xyz.y, direction_.xyz.z), ctx.time);
    LightLiSample sample;
    sample.L = radiance;
    sample.wi = wi;
    sample.pdf = 1;
    sample.pLight = pLight;
    sample.visibility = VisibilityTester{ctx.ToInteraction(), pLight};
    sample.delta = true;
    return sample;
  }

  Float PDF_Li(const LightSampleContext&, const vec3f&, LightSampleMode) const override {
    return 0;
  }

  base::SampledSpectrum Le(const Ray&, const base::SampledWavelengths&) const override {
    return base::SampledSpectrum(0);
  }

  base::SampledSpectrum L(
    const Interaction&,
    const vec3f&,
    const base::SampledWavelengths&
  ) const override {
    return base::SampledSpectrum(0);
  }

  base::SampledSpectrum Phi(const base::SampledWavelengths& lambda) const override {
    Float coneFactor = static_cast<Float>(2) * Pi *
                       ((static_cast<Float>(1) - cosFalloffStart_) +
                        (cosFalloffStart_ - cosTotalWidth_) / static_cast<Float>(2));
    return intensity_.Sample(lambda) * coneFactor;
  }

  void Preprocess(const Bounds3f&) override {}

  Float PowerEstimate() const override {
    Float coneFactor = static_cast<Float>(2) * Pi *
                       ((static_cast<Float>(1) - cosFalloffStart_) +
                        (cosFalloffStart_ - cosTotalWidth_) / static_cast<Float>(2));
    return intensity_.MaxValue() * coneFactor;
  }

private:
  Float Falloff(const vec3f& w) const {
    Float cosTheta = dot(direction_, UnitOrZero(w));
    if (cosTheta < cosTotalWidth_) {
      return 0;
    }
    if (cosTheta > cosFalloffStart_) {
      return 1;
    }
    return SmoothStep(cosTheta, cosTotalWidth_, cosFalloffStart_);
  }

  point3f position_;
  vec3f direction_;
  LightSpectrum intensity_;
  Float cosTotalWidth_ = 0;
  Float cosFalloffStart_ = 0;
  std::string name_ = "SpotLight";
};

class DistantLightImpl : public Light::Concept {
public:
  DistantLightImpl(vec3f directionToLight, LightSpectrum radiance)
    : directionToLight_(UnitOrZero(directionToLight)), radiance_(std::move(radiance)) {
    if (!(directionToLight_.squared_length() > 0)) {
      throw std::invalid_argument("DistantLight direction must be non-zero");
    }
  }

  LightType Type() const override { return LightType::Distant; }
  LightFlags Flags() const override { return LightFlags::DeltaDirection; }
  const std::string& Name() const override { return name_; }

  std::optional<LightLiSample> SampleLi(
    const LightSampleContext& ctx,
    point2f,
    const base::SampledWavelengths& lambda,
    LightSampleMode
  ) const override {
    Interaction pLight = MakeEndpoint(
      ctx.p + directionToLight_ * (static_cast<Float>(2) * sceneRadius_),
      normal3f(-directionToLight_.xyz.x, -directionToLight_.xyz.y, -directionToLight_.xyz.z),
      ctx.time
    );
    LightLiSample sample;
    sample.L = radiance_.Sample(lambda);
    sample.wi = directionToLight_;
    sample.pdf = 1;
    sample.pLight = pLight;
    sample.visibility = VisibilityTester{ctx.ToInteraction(), pLight};
    sample.delta = true;
    return sample;
  }

  Float PDF_Li(const LightSampleContext&, const vec3f&, LightSampleMode) const override {
    return 0;
  }

  base::SampledSpectrum Le(const Ray&, const base::SampledWavelengths&) const override {
    return base::SampledSpectrum(0);
  }

  base::SampledSpectrum L(
    const Interaction&,
    const vec3f&,
    const base::SampledWavelengths&
  ) const override {
    return base::SampledSpectrum(0);
  }

  base::SampledSpectrum Phi(const base::SampledWavelengths& lambda) const override {
    return radiance_.Sample(lambda) * (Pi * sceneRadius_ * sceneRadius_);
  }

  void Preprocess(const Bounds3f& sceneBounds) override {
    sceneCenter_ = SceneCenter(sceneBounds);
    sceneRadius_ = SceneRadius(sceneBounds);
  }

  Float PowerEstimate() const override {
    return radiance_.MaxValue() * Pi * sceneRadius_ * sceneRadius_;
  }

private:
  vec3f directionToLight_;
  LightSpectrum radiance_;
  point3f sceneCenter_ = point3f(0, 0, 0);
  Float sceneRadius_ = 1;
  std::string name_ = "DistantLight";
};

class DiffuseAreaLightImpl : public Light::Concept {
public:
  DiffuseAreaLightImpl(Shape shape, LightSpectrum radiance, DiffuseAreaLightOptions options)
    : shape_(std::move(shape)), radiance_(std::move(radiance)), options_(options) {
    ValidateShape();
  }

  LightType Type() const override { return LightType::DiffuseArea; }
  LightFlags Flags() const override { return LightFlags::Area; }
  const std::string& Name() const override { return name_; }

  std::optional<LightLiSample> SampleLi(
    const LightSampleContext& ctx,
    point2f u,
    const base::SampledWavelengths& lambda,
    LightSampleMode
  ) const override {
    std::optional<ShapeSample> shapeSample = shape_.SampleDirection(ctx.ToInteraction(), u);
    if (!shapeSample || !(shapeSample->pdf > 0) || DistanceSquared(shapeSample->interaction.p, ctx.p) == 0) {
      return std::nullopt;
    }
    vec3f wi = UnitOrZero(shapeSample->interaction.p - ctx.p);
    base::SampledSpectrum emitted = L(shapeSample->interaction, -wi, lambda);
    if (!emitted) {
      return std::nullopt;
    }
    LightLiSample sample;
    sample.L = emitted;
    sample.wi = wi;
    sample.pdf = shapeSample->pdf;
    sample.pLight = shapeSample->interaction;
    sample.visibility = VisibilityTester{ctx.ToInteraction(), shapeSample->interaction};
    sample.delta = false;
    return sample;
  }

  Float PDF_Li(const LightSampleContext& ctx, const vec3f& wi, LightSampleMode) const override {
    return shape_.PDFDirection(ctx.ToInteraction(), wi);
  }

  base::SampledSpectrum Le(const Ray&, const base::SampledWavelengths&) const override {
    return base::SampledSpectrum(0);
  }

  base::SampledSpectrum L(
    const Interaction& interaction,
    const vec3f& w,
    const base::SampledWavelengths& lambda
  ) const override {
    if (!options_.twoSided && dot(interaction.n, w) <= 0) {
      return base::SampledSpectrum(0);
    }
    return radiance_.Sample(lambda);
  }

  base::SampledSpectrum Phi(const base::SampledWavelengths& lambda) const override {
    Float sideScale = options_.twoSided ? static_cast<Float>(2) : static_cast<Float>(1);
    return radiance_.Sample(lambda) * (sideScale * Pi * shape_.Area());
  }

  void Preprocess(const Bounds3f&) override {}

  Float PowerEstimate() const override {
    Float sideScale = options_.twoSided ? static_cast<Float>(2) : static_cast<Float>(1);
    return radiance_.MaxValue() * sideScale * Pi * shape_.Area();
  }

private:
  void ValidateShape() const {
    if (!shape_.IsValid()) {
      throw std::invalid_argument("DiffuseAreaLight requires a valid shape");
    }
    const ShapeCapabilities& capabilities = shape_.Capabilities();
    bool sampleable = capabilities.canSampleDirection || capabilities.canSampleArea;
    if (!capabilities.supportsAreaLight || !sampleable) {
      if (!options_.allowUnsampledEmitter) {
        throw std::invalid_argument(
          "Area lights require sampleable shapes; CSG area lights require csg_mesh_sampler(mode = \"render_mesh\") or advanced unsampled-emitter mode"
        );
      }
    }
  }

  Shape shape_;
  LightSpectrum radiance_;
  DiffuseAreaLightOptions options_;
  std::string name_ = "DiffuseAreaLight";
};

class UniformInfiniteLightImpl : public Light::Concept {
public:
  UniformInfiniteLightImpl(LightSpectrum radiance, Transform3f renderFromLight)
    : radiance_(std::move(radiance)), renderFromLight_(renderFromLight) {}

  LightType Type() const override { return LightType::UniformInfinite; }
  LightFlags Flags() const override { return LightFlags::Infinite; }
  const std::string& Name() const override { return name_; }

  std::optional<LightLiSample> SampleLi(
    const LightSampleContext& ctx,
    point2f u,
    const base::SampledWavelengths& lambda,
    LightSampleMode mode
  ) const override {
    if (mode == LightSampleMode::AllowIncompletePDF) {
      return std::nullopt;
    }
    vec3f wLight = SampleUniformSphere(u);
    vec3f wi = UnitOrZero(renderFromLight_.ApplyVector(wLight));
    Interaction pLight = MakeEndpoint(
      ctx.p + wi * (static_cast<Float>(2) * sceneRadius_),
      normal3f(-wi.xyz.x, -wi.xyz.y, -wi.xyz.z),
      ctx.time
    );
    LightLiSample sample;
    sample.L = radiance_.Sample(lambda);
    sample.wi = wi;
    sample.pdf = UniformSpherePDF();
    sample.pLight = pLight;
    sample.visibility = VisibilityTester{ctx.ToInteraction(), pLight};
    sample.delta = false;
    return sample;
  }

  Float PDF_Li(const LightSampleContext&, const vec3f&, LightSampleMode mode) const override {
    if (mode == LightSampleMode::AllowIncompletePDF) {
      return 0;
    }
    return UniformSpherePDF();
  }

  base::SampledSpectrum Le(const Ray&, const base::SampledWavelengths& lambda) const override {
    return radiance_.Sample(lambda);
  }

  base::SampledSpectrum L(
    const Interaction&,
    const vec3f&,
    const base::SampledWavelengths&
  ) const override {
    return base::SampledSpectrum(0);
  }

  base::SampledSpectrum Phi(const base::SampledWavelengths& lambda) const override {
    return radiance_.Sample(lambda) *
           (static_cast<Float>(4) * Pi * Pi * sceneRadius_ * sceneRadius_);
  }

  void Preprocess(const Bounds3f& sceneBounds) override {
    sceneCenter_ = SceneCenter(sceneBounds);
    sceneRadius_ = SceneRadius(sceneBounds);
  }

  Float PowerEstimate() const override {
    return radiance_.MaxValue() * static_cast<Float>(4) * Pi * Pi * sceneRadius_ * sceneRadius_;
  }

private:
  LightSpectrum radiance_;
  Transform3f renderFromLight_;
  point3f sceneCenter_ = point3f(0, 0, 0);
  Float sceneRadius_ = 1;
  std::string name_ = "UniformInfiniteLight";
};

class ImageInfiniteLightImpl : public Light::Concept {
public:
  ImageInfiniteLightImpl(
    int width,
    int height,
    std::vector<LightSpectrum> texels,
    Transform3f renderFromLight,
    Float scale
  )
    : width_(width),
      height_(height),
      texels_(std::move(texels)),
      renderFromLight_(renderFromLight),
      scale_(scale) {
    if (width_ <= 0 || height_ <= 0) {
      throw std::invalid_argument("ImageInfiniteLight requires positive dimensions");
    }
    if (texels_.size() != static_cast<std::size_t>(width_ * height_)) {
      throw std::invalid_argument("ImageInfiniteLight texel count does not match dimensions");
    }
    ValidateFiniteNonNegative(scale_, "ImageInfiniteLight scale");
    std::vector<Float> weights;
    weights.reserve(texels_.size());
    for (const LightSpectrum& texel : texels_) {
      weights.push_back(texel.MaxValue());
    }
    distribution_ = PiecewiseConstant2D(width_, height_, std::move(weights));
  }

  LightType Type() const override { return LightType::ImageInfinite; }
  LightFlags Flags() const override { return LightFlags::Infinite; }
  const std::string& Name() const override { return name_; }

  std::optional<LightLiSample> SampleLi(
    const LightSampleContext& ctx,
    point2f u,
    const base::SampledWavelengths& lambda,
    LightSampleMode
  ) const override {
    std::optional<PiecewiseSample2D> uv = distribution_.Sample(u);
    if (!uv || !(uv->pdf > 0)) {
      return std::nullopt;
    }
    vec3f wLight = PbrtEqualAreaSquareToSphere(uv->uv);
    vec3f wi = UnitOrZero(renderFromLight_.ApplyVector(wLight));
    Interaction pLight = MakeEndpoint(
      ctx.p + wi * (static_cast<Float>(2) * sceneRadius_),
      normal3f(-wi.xyz.x, -wi.xyz.y, -wi.xyz.z),
      ctx.time
    );
    LightLiSample sample;
    sample.L = ImageLe(uv->uv, lambda);
    sample.wi = wi;
    sample.pdf = uv->pdf * Inv4Pi;
    sample.pLight = pLight;
    sample.visibility = VisibilityTester{ctx.ToInteraction(), pLight};
    sample.delta = false;
    return sample;
  }

  Float PDF_Li(const LightSampleContext&, const vec3f& wi, LightSampleMode) const override {
    vec3f wLight = UnitOrZero(renderFromLight_.ApplyInverseVector(wi));
    point2f uv = PbrtEqualAreaSphereToSquare(wLight);
    return distribution_.PDF(uv) * Inv4Pi;
  }

  base::SampledSpectrum Le(const Ray& ray, const base::SampledWavelengths& lambda) const override {
    vec3f wLight = UnitOrZero(renderFromLight_.ApplyInverseVector(ray.direction()));
    return ImageLe(PbrtEqualAreaSphereToSquare(wLight), lambda);
  }

  base::SampledSpectrum L(
    const Interaction&,
    const vec3f&,
    const base::SampledWavelengths&
  ) const override {
    return base::SampledSpectrum(0);
  }

  base::SampledSpectrum Phi(const base::SampledWavelengths& lambda) const override {
    base::SampledSpectrum sum(0);
    for (const LightSpectrum& texel : texels_) {
      sum += texel.Sample(lambda);
    }
    Float invCount = static_cast<Float>(1) / static_cast<Float>(texels_.size());
    return sum * invCount * scale_ *
           (static_cast<Float>(4) * Pi * Pi * sceneRadius_ * sceneRadius_);
  }

  void Preprocess(const Bounds3f& sceneBounds) override {
    sceneCenter_ = SceneCenter(sceneBounds);
    sceneRadius_ = SceneRadius(sceneBounds);
  }

  Float PowerEstimate() const override {
    Float sum = 0;
    for (const LightSpectrum& texel : texels_) {
      sum += texel.MaxValue();
    }
    Float average = texels_.empty() ? 0 : sum / static_cast<Float>(texels_.size());
    return average * scale_ * static_cast<Float>(4) * Pi * Pi * sceneRadius_ * sceneRadius_;
  }

private:
  const LightSpectrum& Texel(point2f uv) const {
    int x = std::clamp(static_cast<int>(std::floor(Clamp(uv[0], 0, std::nextafter(static_cast<Float>(1), static_cast<Float>(0))) * width_)), 0, width_ - 1);
    int y = std::clamp(static_cast<int>(std::floor(Clamp(uv[1], 0, std::nextafter(static_cast<Float>(1), static_cast<Float>(0))) * height_)), 0, height_ - 1);
    return texels_[static_cast<std::size_t>(y * width_ + x)];
  }

  base::SampledSpectrum ImageLe(point2f uv, const base::SampledWavelengths& lambda) const {
    return Texel(uv).Sample(lambda) * scale_;
  }

  int width_ = 0;
  int height_ = 0;
  std::vector<LightSpectrum> texels_;
  Transform3f renderFromLight_;
  Float scale_ = 1;
  PiecewiseConstant2D distribution_;
  point3f sceneCenter_ = point3f(0, 0, 0);
  Float sceneRadius_ = 1;
  std::string name_ = "ImageInfiniteLight";
};

Light::Light(std::shared_ptr<Concept> impl) : impl_(std::move(impl)) {}

Light Light::Point(point3f position, LightSpectrum intensity) {
  return Light(std::make_shared<PointLightImpl>(position, std::move(intensity)));
}

Light Light::Spot(
  point3f position,
  vec3f direction,
  LightSpectrum intensity,
  Float totalWidthDegrees,
  Float falloffStartDegrees
) {
  if (falloffStartDegrees > totalWidthDegrees) {
    throw std::invalid_argument("Spot falloff start angle must be less than or equal to total width");
  }
  return SpotFromCosines(
    position,
    direction,
    std::move(intensity),
    std::cos(Radians(totalWidthDegrees)),
    std::cos(Radians(falloffStartDegrees))
  );
}

Light Light::SpotFromCosines(
  point3f position,
  vec3f direction,
  LightSpectrum intensity,
  Float cosTotalWidth,
  Float cosFalloffStart
) {
  return Light(std::make_shared<SpotLightImpl>(
    position,
    direction,
    std::move(intensity),
    cosTotalWidth,
    cosFalloffStart
  ));
}

Light Light::Distant(vec3f directionToLight, LightSpectrum radiance) {
  return Light(std::make_shared<DistantLightImpl>(directionToLight, std::move(radiance)));
}

Light Light::DiffuseArea(Shape shape, LightSpectrum radiance, DiffuseAreaLightOptions options) {
  return Light(std::make_shared<DiffuseAreaLightImpl>(std::move(shape), std::move(radiance), options));
}

Light Light::UniformInfinite(LightSpectrum radiance, Transform3f renderFromLight) {
  return Light(std::make_shared<UniformInfiniteLightImpl>(std::move(radiance), renderFromLight));
}

Light Light::ImageInfinite(
  int width,
  int height,
  std::vector<LightSpectrum> texels,
  Transform3f renderFromLight,
  Float scale
) {
  return Light(std::make_shared<ImageInfiniteLightImpl>(
    width,
    height,
    std::move(texels),
    renderFromLight,
    scale
  ));
}

Light Light::ImageInfiniteFromRGB(
  int width,
  int height,
  const std::vector<base::RGB>& texels,
  const base::RGBColorSpace& colorSpace,
  Transform3f renderFromLight,
  Float scale
) {
  std::vector<LightSpectrum> spectra;
  spectra.reserve(texels.size());
  for (const base::RGB& texel : texels) {
    spectra.push_back(LightSpectrum::FromRGBIlluminant(colorSpace, texel));
  }
  return ImageInfinite(width, height, std::move(spectra), renderFromLight, scale);
}

bool Light::IsValid() const {
  return static_cast<bool>(impl_);
}

LightType Light::Type() const {
  if (!impl_) {
    throw std::runtime_error("Invalid Light has no type");
  }
  return impl_->Type();
}

LightFlags Light::Flags() const {
  return impl_ ? impl_->Flags() : LightFlags::None;
}

const std::string& Light::Name() const {
  if (!impl_) {
    static const std::string invalid = "InvalidLight";
    return invalid;
  }
  return impl_->Name();
}

std::optional<LightLiSample> Light::SampleLi(
  const LightSampleContext& ctx,
  point2f u,
  const base::SampledWavelengths& lambda,
  LightSampleMode mode
) const {
  if (!impl_) {
    return std::nullopt;
  }
  return impl_->SampleLi(ctx, u, lambda, mode);
}

Float Light::PDF_Li(const LightSampleContext& ctx, const vec3f& wi, LightSampleMode mode) const {
  return impl_ ? impl_->PDF_Li(ctx, wi, mode) : 0;
}

base::SampledSpectrum Light::Le(const Ray& ray, const base::SampledWavelengths& lambda) const {
  return impl_ ? impl_->Le(ray, lambda) : base::SampledSpectrum(0);
}

base::SampledSpectrum Light::L(
  const Interaction& interaction,
  const vec3f& w,
  const base::SampledWavelengths& lambda
) const {
  return impl_ ? impl_->L(interaction, w, lambda) : base::SampledSpectrum(0);
}

base::SampledSpectrum Light::Phi(const base::SampledWavelengths& lambda) const {
  return impl_ ? impl_->Phi(lambda) : base::SampledSpectrum(0);
}

void Light::Preprocess(const Bounds3f& sceneBounds) {
  if (impl_) {
    impl_->Preprocess(sceneBounds);
  }
}

Float Light::PowerEstimate() const {
  return impl_ ? impl_->PowerEstimate() : 0;
}

base::LightHandle SpectralLightTable::Add(Light light) {
  if (!light.IsValid()) {
    throw std::invalid_argument("Cannot add an invalid Light to SpectralLightTable");
  }
  base::LightHandle handle = base::LightHandle::FromIndex(static_cast<base::LightHandle::IndexType>(lights_.size()), 1);
  lights_.push_back(std::move(light));
  return handle;
}

const Light& SpectralLightTable::Get(base::LightHandle handle) const {
  Validate(handle);
  return lights_[handle.Index()];
}

Light& SpectralLightTable::Get(base::LightHandle handle) {
  Validate(handle);
  return lights_[handle.Index()];
}

std::size_t SpectralLightTable::Size() const {
  return lights_.size();
}

void SpectralLightTable::Preprocess(const Bounds3f& sceneBounds) {
  for (Light& light : lights_) {
    light.Preprocess(sceneBounds);
  }
}

void SpectralLightTable::Validate(base::LightHandle handle) const {
  if (!handle.IsValid() || handle.Generation() != 1 || handle.Index() >= lights_.size()) {
    throw std::out_of_range("Invalid LightHandle");
  }
}

UniformLightSampler::UniformLightSampler(
  const SpectralLightTable& table,
  std::vector<base::LightHandle> lights
)
  : table_(&table), lights_(std::move(lights)) {
  for (base::LightHandle light : lights_) {
    (void)table_->Get(light);
  }
}

std::optional<SampledLight> UniformLightSampler::Sample(Float u) const {
  if (!table_ || lights_.empty()) {
    return std::nullopt;
  }
  Float clamped = Clamp(u, 0, std::nextafter(static_cast<Float>(1), static_cast<Float>(0)));
  std::size_t index = std::min(
    static_cast<std::size_t>(clamped * static_cast<Float>(lights_.size())),
    lights_.size() - 1
  );
  base::LightHandle handle = lights_[index];
  return SampledLight{handle, &table_->Get(handle), static_cast<Float>(1) / static_cast<Float>(lights_.size())};
}

Float UniformLightSampler::PMF(base::LightHandle light) const {
  if (lights_.empty()) {
    return 0;
  }
  return std::find(lights_.begin(), lights_.end(), light) != lights_.end()
           ? static_cast<Float>(1) / static_cast<Float>(lights_.size())
           : 0;
}

Float UniformLightSampler::PMFSum() const {
  Float sum = 0;
  for (base::LightHandle light : lights_) {
    sum += PMF(light);
  }
  return sum;
}

std::size_t UniformLightSampler::Size() const {
  return lights_.size();
}

PowerLightSampler::PowerLightSampler(
  const SpectralLightTable& table,
  std::vector<base::LightHandle> lights
)
  : table_(&table), lights_(std::move(lights)) {
  if (lights_.empty()) {
    return;
  }
  std::vector<Float> weights;
  weights.reserve(lights_.size());
  Float sum = 0;
  for (base::LightHandle light : lights_) {
    Float weight = std::max(static_cast<Float>(0), table_->Get(light).PowerEstimate());
    weights.push_back(weight);
    sum += weight;
  }
  if (!(sum > 0)) {
    weights.assign(lights_.size(), static_cast<Float>(1));
    sum = static_cast<Float>(lights_.size());
  }
  Float running = 0;
  pmf_.reserve(weights.size());
  cdf_.reserve(weights.size());
  for (Float weight : weights) {
    Float p = weight / sum;
    pmf_.push_back(p);
    running += p;
    cdf_.push_back(running);
  }
  cdf_.back() = 1;
}

std::optional<SampledLight> PowerLightSampler::Sample(Float u) const {
  if (!table_ || lights_.empty()) {
    return std::nullopt;
  }
  Float clamped = Clamp(u, 0, std::nextafter(static_cast<Float>(1), static_cast<Float>(0)));
  auto iter = std::lower_bound(cdf_.begin(), cdf_.end(), clamped);
  std::size_t index = static_cast<std::size_t>(iter - cdf_.begin());
  if (index >= lights_.size()) {
    index = lights_.size() - 1;
  }
  base::LightHandle handle = lights_[index];
  return SampledLight{handle, &table_->Get(handle), pmf_[index]};
}

Float PowerLightSampler::PMF(base::LightHandle light) const {
  for (std::size_t i = 0; i < lights_.size(); ++i) {
    if (lights_[i] == light) {
      return pmf_[i];
    }
  }
  return 0;
}

Float PowerLightSampler::PMFSum() const {
  Float sum = 0;
  for (Float p : pmf_) {
    sum += p;
  }
  return sum;
}

std::size_t PowerLightSampler::Size() const {
  return lights_.size();
}

PrimitiveBinding CompileAreaLightPrimitiveBinding(
  base::MaterialHandle material,
  base::LightHandle areaLight,
  const Shape& shape,
  AreaLightBindingOptions options
) {
  if (!areaLight.IsValid()) {
    throw std::invalid_argument("Area light primitive binding requires a valid LightHandle");
  }
  const ShapeCapabilities& capabilities = shape.Capabilities();
  bool sampleable = capabilities.canSampleDirection || capabilities.canSampleArea;
  if (!capabilities.supportsAreaLight || !sampleable) {
    if (!options.allowUnsampledEmitter) {
      throw std::invalid_argument(
        "Area lights require sampleable shapes; CSG area lights require csg_mesh_sampler(mode = \"render_mesh\") or advanced unsampled-emitter mode"
      );
    }
  }
  PrimitiveBinding binding;
  binding.material = material;
  binding.areaLight = areaLight;
  binding.shapeRequirements.areaLight = true;
  return binding;
}

Light ConvertLegacyEmissiveMaterialToAreaLight(
  const LegacyEmissiveMaterialDescriptor& legacy,
  const Shape& shape,
  const base::RGBColorSpace& colorSpace,
  DiffuseAreaLightOptions options
) {
  base::RGB scaled = legacy.color * legacy.intensity;
  return Light::DiffuseArea(shape, LightSpectrum::FromRGBIlluminant(colorSpace, scaled), options);
}

Light ConvertLegacySpotLightMaterialToLight(
  const LegacySpotLightMaterialDescriptor& legacy,
  const base::RGBColorSpace& colorSpace
) {
  base::RGB scaled = legacy.color * legacy.intensity;
  return Light::SpotFromCosines(
    legacy.position,
    legacy.direction,
    LightSpectrum::FromRGBIlluminant(colorSpace, scaled),
    legacy.cosTotalWidth,
    legacy.cosFalloffStart
  );
}

base::LightHandle RegisterLight(Scene& scene, SpectralLightTable& table, Light light) {
  base::LightHandle handle = table.Add(std::move(light));
  AttachLightToScene(scene, handle, table.Get(handle));
  return handle;
}

void AttachLightToScene(Scene& scene, base::LightHandle handle, const Light& light) {
  if (HasFlag(light.Flags(), LightFlags::Infinite)) {
    scene.AddInfiniteLight(handle);
  } else {
    scene.AddLight(handle);
  }
}

bool ShouldCompileLegacyEnvironmentSphereToSpectralGeometry() {
  return false;
}

bool SpectralMaterialEmissionPathEnabled() {
  return false;
}

} // namespace render
} // namespace rayrender
