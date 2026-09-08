#include "infinite.h"
#include "../materials/texturecache.h"
#include "../math/mathinline.h"
#include "../hitables/hitable.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {
// Bilinear reconstruction can put radiance in a texel whose tabulated weight
// is zero. A small uniform component guarantees support there and at seams.
constexpr Float uniform_fraction = .001f;
constexpr Float sphere_pdf = 1 / (4 * M_PI);
vec3f uniform_direction(vec2f u) {
  Float y = 1 - 2 * u[0], r = std::sqrt(std::max(Float(0), 1 - y * y));
  Float phi = 2 * M_PI * u[1];
  return vec3f(r * std::cos(phi), y, r * std::sin(phi));
}
vec3f image_direction(vec2f uv) {
  Float theta = (1 - uv[1]) * M_PI, phi = -uv[0] * 2 * M_PI;
  return vec3f(std::sin(theta) * std::sin(phi), std::cos(theta),
               std::sin(theta) * std::cos(phi));
}
}

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
  if (max_value > 0) {
    for (size_t i = 0; i < values.size(); ++i)
      values[i] = double(values[i]) / max_value;
    distribution = std::make_unique<Distribution2D>(values.data(), width, height);
  }
  // Integrate over solid angle, so image and small-disk selection weights have
  // the same units. This common factor leaves image-only mixtures unchanged.
  weight *= 2 * M_PI * M_PI;
}

point3f ImageInfiniteLight::Radiance(const point3f &p, const vec3f &wi, Float) const {
  Float u, v;
  get_sphere_uv_z(unit_vector(world_to_light(wi)), u, v);
  return image->value(u, v, p);
}

vec3f ImageInfiniteLight::Sample(const point3f &, vec2f u, Float) const {
  if (!distribution) return uniform_direction(u);
  if (u[0] < uniform_fraction) {
    u.xy.x /= uniform_fraction;
    return uniform_direction(u);
  }
  u.xy.x = (u[0] - uniform_fraction) / (1 - uniform_fraction);
  Float pdf;
  return light_to_world(image_direction(distribution->SampleContinuous(u, &pdf)));
}

Float ImageInfiniteLight::Pdf(const point3f &, const vec3f &wi, Float) const {
  if (!distribution) return sphere_pdf;
  Float u, v;
  get_sphere_uv_z(unit_vector(world_to_light(wi)), u, v);
  Float sin_theta = std::sin(v * M_PI);
  Float map_pdf = v > 0 && v < 1 && sin_theta > 0 ? distribution->Pdf(vec2f(u, v)) /
                                    (2 * M_PI * M_PI * sin_theta) : 0;
  return (1 - uniform_fraction) * map_pdf + uniform_fraction * sphere_pdf;
}

size_t ImageInfiniteLight::GetSize() const {
  return sizeof(*this) + (distribution ? distribution->GetSize() : 0);
}

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
  selection = std::make_unique<Distribution1D>(weights.data(), weights.size());
}

point3f InfiniteLightMixture::Radiance(const point3f &p, const vec3f &wi, Float time) const {
  point3f result(0);
  for (const auto &light : lights) result += light->Radiance(p, wi, time);
  return result;
}

vec3f InfiniteLightMixture::Sample(const point3f &p, vec2f u, Float time) const {
  Float remapped;
  int index = selection->SampleDiscrete(u[0], nullptr, &remapped);
  u.xy.x = remapped;
  return lights[index]->Sample(p, u, time);
}

Float InfiniteLightMixture::Pdf(const point3f &p, const vec3f &wi, Float time) const {
  double result = 0;
  for (size_t i = 0; i < lights.size(); ++i)
    result += selection->DiscretePDF(i) * lights[i]->Pdf(p, wi, time);
  return Float(result);
}

size_t InfiniteLightMixture::GetSize() const {
  size_t result = sizeof(*this) + selection->GetSize();
  for (const auto &light : lights) result += light->GetSize();
  return result;
}

DiskInfiniteLight::DiskInfiniteLight(std::shared_ptr<texture> source, int width,
                                   int height, vec3f direction, double diameter,
                                   Float rotation, bool horizon)
    : image(std::move(source)), clip_horizon(horizon) {
  if (width <= 0 || height <= 0 || !std::isfinite(diameter) ||
      diameter <= 0 || diameter >= 180)
    throw std::runtime_error("Celestial disk dimensions and angular diameter must be positive; diameter must be below 180 degrees.");
  vec3<double> d(direction[0], direction[1], direction[2]);
  if (!std::isfinite(d.length()) || d.length() == 0)
    throw std::runtime_error("Celestial disk direction must be finite and nonzero.");
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
  double radius = diameter * M_PI / 360;
  tan_radius = std::tan(radius);
  // Avoid cancellation for the small angles used by the Sun and Moon.
  one_minus_cos_radius = 2 * std::pow(std::sin(radius / 2), 2);
  solid_angle = 2 * M_PI * one_minus_cos_radius;
  if (!std::isfinite(Float(1 / solid_angle)))
    throw std::runtime_error("Celestial disk angular diameter is too small for the renderer's PDF precision.");
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
      double jacobian = 4 * tan_radius * tan_radius /
          std::pow(1 + tan_radius * tan_radius * (xx * xx + yy * yy), 1.5);
      weight += std::max(0.0, .212671 * rgb[0] + .715160 * rgb[1] + .072169 * rgb[2]) *
          jacobian / (double(width) * height);
    }
  }
}

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

point3f DiskInfiniteLight::Radiance(const point3f &p, const vec3f &wi, Float) const {
  Float u, v;
  if ((clip_horizon && wi[1] < 0) || !Coordinates(wi, u, v)) return point3f(0);
  return image->value(u, v, p);
}

vec3f DiskInfiniteLight::Sample(const point3f &, vec2f u, Float) const {
  double delta = double(u[0]) * one_minus_cos_radius;
  double sin_theta = std::sqrt(delta * (2 - delta)), phi = 2 * M_PI * u[1];
  vec3<double> w = (1 - delta) * forward +
      sin_theta * (std::cos(phi) * right + std::sin(phi) * up);
  return vec3f(w[0], w[1], w[2]);
}

Float DiskInfiniteLight::Pdf(const point3f &, const vec3f &wi, Float) const {
  Float u, v;
  return Coordinates(wi, u, v) ? Float(1 / solid_angle) : 0;
}

std::shared_ptr<InfiniteLight> BuildInfiniteLights(const Rcpp::List &descriptions,
                                                      TextureCache &textures) {
  std::vector<std::shared_ptr<InfiniteLight>> lights;
  for (R_xlen_t i = 0; i < descriptions.size(); ++i) {
    Rcpp::checkUserInterrupt();
    Rcpp::List item = descriptions[i];
    std::string filename = Rcpp::as<std::string>(item["filename"]);
    std::string type = Rcpp::as<std::string>(item["type"]);
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
    // R's public rotation convention is the negative of the renderer's RotateY.
    if (type == "disk") {
      Rcpp::NumericVector d = item["direction"];
      if (d.size() != 3) throw std::runtime_error("Celestial disk direction must have three components.");
      lights.push_back(std::make_shared<DiskInfiniteLight>(
          image, width, height, vec3f(d[0], d[1], d[2]),
          Rcpp::as<double>(item["angular_diameter"]), -Float(std::fmod(rotation, 360.0)),
          Rcpp::as<bool>(item["clip_horizon"])));
    } else {
      lights.push_back(std::make_shared<ImageInfiniteLight>(
          image, width, height, -Float(std::fmod(rotation, 360.0))));
    }
  }
  if (lights.size() == 1) return lights.front();
  return std::make_shared<InfiniteLightMixture>(std::move(lights));
}
