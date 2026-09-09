#include "infinite.h"
#include "atmosphere.h"
#include "../materials/texturecache.h"
#include "../math/mathinline.h"
#include "../hitables/hitable.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>


// Infinite lights provide both incident radiance and a directional sampling
// distribution. Images cover the whole sphere, celestial textures occupy small
// cones, and mixtures combine their emission while sharing atmospheric filtering.
namespace {
// Bilinear reconstruction can put radiance in a texel whose tabulated weight
// is zero. A small uniform component guarantees support there and at seams.
constexpr Float uniform_fraction = .001f;
constexpr Float sphere_pdf = 1 / (4 * M_PI);

// Uniform solid-angle sampling in the renderer's Y-up environment frame.
vec3f uniform_direction(vec2f u) {
  Float y = 1 - 2 * u[0], r = std::sqrt(std::max(Float(0), 1 - y * y));
  Float phi = 2 * M_PI * u[1];
  return vec3f(r * std::cos(phi), y, r * std::sin(phi));
}


// Invert the latitude-longitude convention used by get_sphere_uv_z(), including
// its azimuth sign, so sampling and texture lookup agree at the seam.
vec3f image_direction(vec2f uv) {
  Float theta = (1 - uv[1]) * M_PI, phi = -uv[0] * 2 * M_PI;
  return vec3f(std::sin(theta) * std::sin(phi), std::cos(theta),
               std::sin(theta) * std::cos(phi));
}
}


// Build an importance map from emitted luminance. sin(theta) accounts for the
// smaller solid angle of latitude-longitude texels near the poles.
ImageInfiniteLight::ImageInfiniteLight(std::shared_ptr<texture> source, int width,
                                     int height, Float rotation)
    : image(std::move(source)), light_to_world(RotateY(rotation)),
      world_to_light(Inverse(light_to_world)) {
  if (width <= 0 || height <= 0)
    throw std::runtime_error("Infinite light image dimensions must be positive.");
  std::vector<Float> values(size_t(width) * height);
  double max_value = 0;
  // Normalize in double precision before building the Float CDF so bright HDR
  // images cannot overflow its integral.
  for (int y = 0; y < height; ++y) {
    if (y % 64 == 0) Rcpp::checkUserInterrupt();
    Float v = Float(y + .5) / height;
    double sin_theta = std::sin(M_PI * v);
    for (int x = 0; x < width; ++x) {
      point3f rgb = image->value(Float(x + .5) / width, v, point3f(0));
      for (int c = 0; c < 3; ++c)
        if (!std::isfinite(rgb[c]))
          throw std::runtime_error("Infinite light radiance must be finite; check image values and intensity.");
      double lum = std::max(0.0, .212671 * rgb[0] + .715160 * rgb[1] + .072169 * rgb[2]);
      double value = lum * sin_theta;
      values[size_t(y) * width + x] = Float(value);
      max_value = std::max(max_value, value);
      weight += value / double(values.size());
    }
  }


  // A common scale preserves sampling probabilities while keeping the CDF
  // manageable. With no positive luminance, leave sampling to the uniform path.
  if (max_value > 0) {
    for (size_t i = 0; i < values.size(); ++i)
      values[i] = double(values[i]) / max_value;
    distribution = std::make_unique<Distribution2D>(values.data(), width, height);
  }


  // Integrate over solid angle, so image and small-disk selection weights have
  // the same units. This common factor leaves image-only mixtures unchanged.
  weight *= 2 * M_PI * M_PI;
}


// Undo the enclosing environment transform and this image's own rotation before
// looking up radiance. The importance map does not replace texture reconstruction.
point3f ImageInfiniteLight::Radiance(const point3f &p, const vec3f &wi, Float) const {
  Float u, v;
  get_sphere_uv_z(unit_vector(world_to_light(EnvironmentDirection(wi))), u, v);
  return image->value(u, v, p);
}


// Mix a uniform sphere proposal with the image proposal. Remapping u.x within
// each branch leaves a uniform variate for that branch's directional sampler.
vec3f ImageInfiniteLight::Sample(const point3f &, vec2f u, Float) const {
  if (!distribution) return WorldDirection(uniform_direction(u));
  if (u[0] < uniform_fraction) {
    u.xy.x /= uniform_fraction;
    return WorldDirection(uniform_direction(u));
  }


  u.xy.x = (u[0] - uniform_fraction) / (1 - uniform_fraction);
  Float pdf;
  return WorldDirection(light_to_world(image_direction(distribution->SampleContinuous(u, &pdf))));
}


// Convert the image's UV density to solid angle with 2*pi^2*sin(theta), then
// include the uniform branch. Pole queries retain the uniform density alone.
Float ImageInfiniteLight::Pdf(const point3f &, const vec3f &wi, Float) const {
  if (!distribution) return sphere_pdf;
  Float u, v;
  get_sphere_uv_z(unit_vector(world_to_light(EnvironmentDirection(wi))), u, v);
  Float sin_theta = std::sin(v * M_PI);
  Float map_pdf = v > 0 && v < 1 && sin_theta > 0 ? distribution->Pdf(vec2f(u, v)) /
                                    (2 * M_PI * M_PI * sin_theta) : 0;
  return (1 - uniform_fraction) * map_pdf + uniform_fraction * sphere_pdf;
}


size_t ImageInfiniteLight::GetSize() const {
  return sizeof(*this) + (distribution ? distribution->GetSize() : 0);
}


// Choose sources using their estimated integrated radiance. Scaling by the
// brightest source keeps relative weights representable in the Float CDF.
InfiniteLightMixture::InfiniteLightMixture(std::vector<std::shared_ptr<InfiniteLight>> sources)
    : lights(std::move(sources)) {
  if (lights.empty()) throw std::runtime_error("An infinite light mixture must contain lights.");
  double maximum = 0;
  for (const auto &light : lights) {
    weight += light->SamplingWeight();
    maximum = std::max(maximum, light->SamplingWeight());
  }
  std::vector<Float> weights;
  for (const auto &light : lights)
    weights.push_back(maximum > 0 ? light->SamplingWeight() / maximum : 1);


  // Retain a useful sky proposal even when the intrinsic solar disk is much
  // brighter than a dusk sky. A fully occluded disk is excluded per position.
  bool celestial = false;
  for (const auto &light : lights)
    celestial |= light->RadianceSpectrum() != InfiniteLightSpectrum::RGB;
  if (celestial) for (size_t i = 0; i < lights.size(); ++i) {
    if (!lights[i]->GetAtmosphere() || !(weights[i] > 0)) continue;
    Float other = 0;
    for (size_t j = 0; j < weights.size(); ++j) if (i != j) other += weights[j];
    weights[i] = std::max(weights[i], other);
  }
  selection = std::make_unique<Distribution1D>(weights.data(), weights.size());
}


// Emission adds across all sources, independently of which source was sampled.
// The atmosphere's own radiance is already filtered; other lights receive the
// transmission appropriate to RGB, solar, or lunar reference spectra.
point3f InfiniteLightMixture::Radiance(const point3f &p, const vec3f &wi, Float time) const {
  point3f result(0);
  const Atmosphere *atmosphere = GetAtmosphere();

  // All sources with the same reference spectrum share one transmission query
  // for this position and direction, even when several textures contribute.
  point3f transmission[3];
  bool evaluated[3] = {false, false, false};
  for (const auto &light : lights) {
    point3f value = light->Radiance(p, wi, time);
    // Most directions miss a disk. Do not evaluate atmospheric transmission
    // until a source contributes in this direction.
    if (atmosphere && !light->GetAtmosphere() &&
        (value[0] != 0 || value[1] != 0 || value[2] != 0)) {
      auto spectrum = light->RadianceSpectrum();
      int index = int(spectrum);
      if (!evaluated[index]) {
        transmission[index] = spectrum == InfiniteLightSpectrum::RGB
          ? atmosphere->Transmission(p, wi, INFINITY)
          : atmosphere->CelestialTransmission(p, wi, spectrum);
        evaluated[index] = true;
      }
      value *= transmission[index];
    }
    result += value;
  }
  return result;
}


// The light factory permits at most one atmospheric sky in the mixture. Other lights
// use that same object for visibility and transmission at their shading points.
const Atmosphere *InfiniteLightMixture::GetAtmosphere() const {
  for (const auto &source : lights)
    if (const auto *atmosphere = source->GetAtmosphere()) return atmosphere;
  return nullptr;
}


// Propagate the enclosing environment rotation to every source so preview
// changes affect sampling, radiance lookup, and atmospheric filtering together.
void InfiniteLightMixture::SetEnvironmentTransform(const Transform *to_world,
                                                   const Transform *to_local) {
  InfiniteLight::SetEnvironmentTransform(to_world, to_local);
  for (auto &source : lights) source->SetEnvironmentTransform(to_world, to_local);
}


// Reweight source selection at the shading point when the atmosphere can hide
// entire celestial disks. Sample and Pdf must use the same available-source mass.
vec3f InfiniteLightMixture::Sample(const point3f &p, vec2f u, Float time) const {
  const Atmosphere *atmosphere = GetAtmosphere();
  double total = atmosphere ? ChoiceTotal(p, atmosphere) : 0;
  if (total > 0) {
    // Walk the available weights and reuse the residual distance within the
    // chosen interval as that source's first sampling variate.
    double target = std::clamp(double(u[0]), 0.0, std::nextafter(1.0, 0.0)) * total;
    size_t last = 0;
    for (size_t i = 0; i < lights.size(); ++i) {
      double weight = ChoiceWeight(i, p, atmosphere);
      if (!(weight > 0)) continue;
      last = i;
      if (target < weight) {
        u.xy.x = std::min(Float(target / weight), std::nextafter(Float(1), Float(0)));
        return lights[i]->Sample(p, u, time);
      }
      target -= weight;
    }

    // Rounding can leave the target at the end of the last positive interval.
    u.xy.x = std::nextafter(Float(1), Float(0));
    return lights[last]->Sample(p, u, time);
  }


  // Use the static distribution without position-dependent weights, including
  // the fallback when those weights have no positive total.
  Float remapped;
  int index = selection->SampleDiscrete(u[0], nullptr, &remapped);
  u.xy.x = remapped;
  return lights[index]->Sample(p, u, time);
}


// A direction can be proposed by several lights. Its density is the weighted
// sum of all proposals, not just the density of the source selected by Sample().
Float InfiniteLightMixture::Pdf(const point3f &p, const vec3f &wi, Float time) const {
  double result = 0;
  const Atmosphere *atmosphere = GetAtmosphere();
  double total = atmosphere ? ChoiceTotal(p, atmosphere) : 0;
  for (size_t i = 0; i < lights.size(); ++i)
    result += (total > 0 ? ChoiceWeight(i, p, atmosphere) / total : selection->DiscretePDF(i)) *
              lights[i]->Pdf(p, wi, time);
  return Float(result);
}


// Availability removes whole hidden sources without changing the relative
// weights of the remaining ones. ChoiceTotal supplies their normalization.
double InfiniteLightMixture::ChoiceWeight(size_t i, const point3f &p, const Atmosphere *atmosphere) const {
  return lights[i]->Available(p, atmosphere) ? selection->DiscretePDF(i) : 0;
}


double InfiniteLightMixture::ChoiceTotal(const point3f &p, const Atmosphere *atmosphere) const {
  double total = 0;
  for (size_t i = 0; i < lights.size(); ++i) total += ChoiceWeight(i, p, atmosphere);
  return total;
}


size_t InfiniteLightMixture::GetSize() const {
  size_t result = sizeof(*this) + selection->GetSize();
  for (const auto &light : lights) result += light->GetSize();
  return result;
}


// Keep a celestial image in its own rectilinear projection instead of baking it
// into a sky map. A local orthonormal frame maps between its pixels and directions.
DiskInfiniteLight::DiskInfiniteLight(std::shared_ptr<texture> source, int width,
                                   int height, vec3f direction, double diameter,
                                   Float rotation, bool horizon, InfiniteLightSpectrum source_spectrum)
    : image(std::move(source)), clip_horizon(horizon), spectrum(source_spectrum) {
  if (width <= 0 || height <= 0 || !std::isfinite(diameter) ||
      diameter <= 0 || diameter >= 180)
    throw std::runtime_error("Celestial disk dimensions and angular diameter must be positive; diameter must be below 180 degrees.");
  vec3<double> d(direction[0], direction[1], direction[2]);
  if (!std::isfinite(d.length()) || d.length() == 0)
    throw std::runtime_error("Celestial disk direction must be finite and nonzero.");


  // Use world up to orient the disk, with a second reference axis at the poles
  // where the usual cross product would become degenerate.
  forward = unit_vector(d);
  vec3<double> vertical = std::abs(forward[1]) < .999999 ?
      vec3<double>(0, 1, 0) : vec3<double>(0, 0, 1);
  right = unit_vector(cross(forward, vertical));
  up = cross(right, forward);


  // Rotate the entire frame rather than reconstructing it after rotation.
  double a = rotation * M_PI / 180, c = std::cos(a), s = std::sin(a);
  auto rotate = [&](vec3<double> v) {
    return vec3<double>(c * v[0] + s * v[2], v[1], -s * v[0] + c * v[2]);
  };
  forward = rotate(forward); right = rotate(right); up = rotate(up);


  // Store both the tangent-plane radius for texture projection and the spherical
  // cap area for sampling. They describe the same apparent angular diameter.
  double radius = diameter * M_PI / 360;
  tan_radius = std::tan(radius);
  // Avoid cancellation for the small angles used by the Sun and Moon.
  one_minus_cos_radius = 2 * std::pow(std::sin(radius / 2), 2);
  solid_angle = 2 * M_PI * one_minus_cos_radius;
  if (!std::isfinite(Float(1 / solid_angle)))
    throw std::runtime_error("Celestial disk angular diameter is too small for the renderer's PDF precision.");


  // Estimate this source's integrated brightness for mixture selection. Mask
  // texels outside the disk and, when requested, below the fixed local horizon.
  for (int y = 0; y < height; ++y) {
    if (y % 64 == 0) Rcpp::checkUserInterrupt();
    double v = (y + .5) / height, yy = 2 * v - 1;
    for (int x = 0; x < width; ++x) {
      double u = (x + .5) / width, xx = 2 * u - 1;
      if (xx * xx + yy * yy > 1) continue;
      if (clip_horizon && forward[1] + tan_radius * (xx * right[1] + yy * up[1]) < 0)
        continue;
      point3f rgb = image->value(Float(u), Float(v), point3f(0));
      for (int c = 0; c < 3; ++c)
        if (!std::isfinite(rgb[c]))
          throw std::runtime_error("Celestial disk radiance must be finite.");

      // Rectilinear pixels subtend different solid angles across the disk;
      // this projection Jacobian converts their UV areas to spherical areas.
      double jacobian = 4 * tan_radius * tan_radius /
          std::pow(1 + tan_radius * tan_radius * (xx * xx + yy * yy), 1.5);
      weight += std::max(0.0, .212671 * rgb[0] + .715160 * rgb[1] + .072169 * rgb[2]) *
          jacobian / (double(width) * height);
    }
  }
}


// Project a direction onto the disk's tangent plane, reject points outside its
// circular support, then map the surviving [-1,1] coordinates into texture UVs.
bool DiskInfiniteLight::Coordinates(const vec3f &wi, Float &u, Float &v) const {
  vec3<double> w(wi[0], wi[1], wi[2]);
  double z = dot(w, forward);
  if (z <= 0) return false;
  double x = dot(w, right) / (z * tan_radius);
  double y = dot(w, up) / (z * tan_radius);
  // Testing transverse coordinates avoids acos/dot cancellation at small angles.
  if (x * x + y * y > 1) return false;
  u = Float(.5 + .5 * x); v = Float(.5 + .5 * y);
  return true;
}


// The optional local horizon mask affects emission here. When a native atmosphere
// is present, its spherical-Earth filtering is applied by the enclosing mixture.
point3f DiskInfiniteLight::Radiance(const point3f &p, const vec3f &wi, Float) const {
  vec3f local = EnvironmentDirection(wi);
  Float u, v;
  if ((clip_horizon && local[1] < 0) || !Coordinates(local, u, v)) return point3f(0);
  return image->value(u, v, p);
}


// Sample the entire cone uniformly in solid angle, including dark phase texels
// and clipped portions. Keeping that support fixed makes the PDF simply 1/area.
// delta = 1-cos(theta) preserves precision for the Sun and Moon's small radii.
vec3f DiskInfiniteLight::Sample(const point3f &, vec2f u, Float) const {
  double delta = double(u[0]) * one_minus_cos_radius;
  double sin_theta = std::sqrt(delta * (2 - delta)), phi = 2 * M_PI * u[1];
  vec3<double> w = (1 - delta) * forward +
      sin_theta * (std::cos(phi) * right + std::sin(phi) * up);
  return WorldDirection(vec3f(w[0], w[1], w[2]));
}


Float DiskInfiniteLight::Pdf(const point3f &, const vec3f &wi, Float) const {
  Float u, v;
  return Coordinates(EnvironmentDirection(wi), u, v) ? Float(1 / solid_angle) : 0;
}


// Cull an entire celestial proposal only when the atmosphere says even its upper
// limb is hidden. Individual directions still receive atmospheric filtering afterward.
bool DiskInfiniteLight::Available(const point3f &p, const Atmosphere *atmosphere) const {
  return !atmosphere || spectrum == InfiniteLightSpectrum::RGB ||
    atmosphere->MaySeeDisk(p, WorldDirection(vec3f(forward[0], forward[1], forward[2])),
                          std::atan(tan_radius));
}


// Turn prepared R descriptions into native lights before rendering starts.
// Image storage comes from the shared texture cache; atmosphere construction
// instead loads model coefficients and its directional sampling proposals.
std::shared_ptr<InfiniteLight> BuildInfiniteLights(const Rcpp::List &descriptions,
                                                      TextureCache &textures) {
  std::vector<std::shared_ptr<InfiniteLight>> lights;
  for (R_xlen_t i = 0; i < descriptions.size(); ++i) {
    Rcpp::checkUserInterrupt();
    Rcpp::List item = descriptions[i];
    std::string filename = Rcpp::as<std::string>(item["filename"]);
    std::string type = Rcpp::as<std::string>(item["type"]);

    // One sky supplies the scene's atmospheric transport for all other sources.
    if (type == "prague") {
      for (const auto &source : lights)
        if (source->GetAtmosphere())
          throw std::runtime_error("Only one atmospheric sky can be used in a scene.");
      lights.push_back(std::make_shared<PragueInfiniteLight>(item));
      continue;
    }


    // Both image-based types share linear HDR loading and intensity validation.
    // Reject invalid texels here rather than allowing worker queries to see them.
    if (type != "image" && type != "disk")
      throw std::runtime_error("Unsupported infinite light type.");
    Float intensity = Rcpp::as<Float>(item["intensity"]);
    double rotation = Rcpp::as<double>(item["rotation"]);
    if (!std::isfinite(intensity) || intensity < 0 || !std::isfinite(rotation))
      throw std::runtime_error("Infinite light intensity and rotation must be finite and representable; intensity must be nonnegative.");
    int width, height, channels;
    Float *data = textures.LookupFloat(filename, width, height, channels, 3);
    for (size_t j = 0; j < size_t(width) * height * channels; ++j) {
      if (j % 1048576 == 0) Rcpp::checkUserInterrupt();
      if (!std::isfinite(data[j]))
        throw std::runtime_error("Infinite light image contains non-finite pixels: " + filename);
    }
    auto image = std::make_shared<image_texture_float>(
        data, width, height, channels, 1, 1, intensity);


    // A disk also carries its geometric support and reference spectrum. The
    // spectrum tag selects the appropriate atmospheric transmission calculation.
    // R's public rotation convention is the negative of the renderer's RotateY.
    if (type == "disk") {
      Rcpp::NumericVector d = item["direction"];
      if (d.size() != 3) throw std::runtime_error("Celestial disk direction must have three components.");
      std::string spectrum = item.containsElementNamed("radiance_spectrum")
                               ? Rcpp::as<std::string>(item["radiance_spectrum"]) : "rgb";
      if (spectrum != "rgb" && spectrum != "sun" && spectrum != "moon")
        throw std::runtime_error("Unknown celestial radiance spectrum.");
      auto source_spectrum = spectrum == "sun" ? InfiniteLightSpectrum::Sun
                            : spectrum == "moon" ? InfiniteLightSpectrum::Moon : InfiniteLightSpectrum::RGB;
      lights.push_back(std::make_shared<DiskInfiniteLight>(
          image, width, height, vec3f(d[0], d[1], d[2]),
          Rcpp::as<double>(item["angular_diameter"]), -Float(std::fmod(rotation, 360.0)),
          Rcpp::as<bool>(item["clip_horizon"]), source_spectrum));
    } else {
      lights.push_back(std::make_shared<ImageInfiniteLight>(
          image, width, height, -Float(std::fmod(rotation, 360.0))));
    }
  }


  // Avoid an extra dispatch layer for a single source; otherwise expose one
  // additive light with a consistent sampling mixture to the integrators.
  if (lights.size() == 1) return lights.front();
  return std::make_shared<InfiniteLightMixture>(std::move(lights));
}
