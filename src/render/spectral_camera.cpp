#include "spectral_camera.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace rayrender {
namespace render {

namespace {

void ValidateFinite(Float value, const char* name) {
  if (!std::isfinite(value)) {
    throw std::invalid_argument(std::string(name) + " must be finite");
  }
}

void ValidatePoint(const point2f& value, const char* name) {
  ValidateFinite(value.xy.x, name);
  ValidateFinite(value.xy.y, name);
}

void ValidatePoint(const point3f& value, const char* name) {
  ValidateFinite(value.xyz.x, name);
  ValidateFinite(value.xyz.y, name);
  ValidateFinite(value.xyz.z, name);
}

void ValidateVector(const vec3f& value, const char* name) {
  ValidateFinite(value.xyz.x, name);
  ValidateFinite(value.xyz.y, name);
  ValidateFinite(value.xyz.z, name);
}

Float DegreesToRadians(Float degrees) {
  return degrees * static_cast<Float>(M_PI) / static_cast<Float>(180);
}

Float SafeLength(const vec3f& value) {
  return std::sqrt(value.squared_length());
}

vec3f SafeNormalize(const vec3f& value, const char* name) {
  Float length = SafeLength(value);
  if (!(length > 0) || !std::isfinite(length)) {
    throw std::invalid_argument(std::string(name) + " must have non-zero finite length");
  }
  return value / length;
}

void BuildCameraFrame(
  const point3f& lookfrom,
  const point3f& lookat,
  const vec3f& upHint,
  vec3f& forward,
  vec3f& right,
  vec3f& up
) {
  ValidatePoint(lookfrom, "camera lookfrom");
  ValidatePoint(lookat, "camera lookat");
  ValidateVector(upHint, "camera up");

  forward = SafeNormalize(lookat - lookfrom, "camera view direction");
  vec3f upCandidate = upHint.length() > 0 ? SafeNormalize(upHint, "camera up") : vec3f(0, 1, 0);
  right = cross(upCandidate, forward);
  if (right.squared_length() < static_cast<Float>(1e-12)) {
    upCandidate = std::abs(forward.xyz.y) < static_cast<Float>(0.999) ? vec3f(0, 1, 0) : vec3f(1, 0, 0);
    right = cross(upCandidate, forward);
    if (right.squared_length() < static_cast<Float>(1e-12)) {
      upCandidate = vec3f(0, 0, 1);
      right = cross(upCandidate, forward);
    }
  }
  right = SafeNormalize(right, "camera right");
  up = SafeNormalize(cross(forward, right), "camera resolved up");
}

SpectralCameraOptions NormalizeOptions(SpectralCameraOptions options) {
  if (options.filmWidth <= 0 || options.filmHeight <= 0) {
    throw std::invalid_argument("camera film dimensions must be positive");
  }
  ValidateFinite(options.shutterOpen, "camera shutterOpen");
  ValidateFinite(options.shutterClose, "camera shutterClose");
  if (options.shutterClose < options.shutterOpen) {
    std::swap(options.shutterOpen, options.shutterClose);
  }
  return options;
}

point2f SampleUniformDiskConcentric(const point2f& u) {
  Float uOffsetX = static_cast<Float>(2) * u.xy.x - static_cast<Float>(1);
  Float uOffsetY = static_cast<Float>(2) * u.xy.y - static_cast<Float>(1);
  if (uOffsetX == 0 && uOffsetY == 0) {
    return point2f(0, 0);
  }

  Float theta;
  Float radius;
  if (std::abs(uOffsetX) > std::abs(uOffsetY)) {
    radius = uOffsetX;
    theta = static_cast<Float>(M_PI_4) * (uOffsetY / uOffsetX);
  } else {
    radius = uOffsetY;
    theta = static_cast<Float>(M_PI_2) - static_cast<Float>(M_PI_4) * (uOffsetX / uOffsetY);
  }
  return point2f(radius * std::cos(theta), radius * std::sin(theta));
}

bool WavelengthsUnchanged(
  const base::SampledWavelengths& before,
  const base::SampledWavelengths& after
) {
  for (int i = 0; i < base::NSpectrumSamples; ++i) {
    if (before[i] != after[i] || before.PDF(i) != after.PDF(i)) {
      return false;
    }
  }
  return true;
}

} // namespace

SpectralCamera SpectralCamera::Perspective(const PerspectiveCameraParameters& params) {
  if (!(params.vfov > 0 && params.vfov < static_cast<Float>(180))) {
    throw std::invalid_argument("perspective camera vfov must be in (0, 180)");
  }
  if (!(params.aspect > 0) || !std::isfinite(params.aspect)) {
    throw std::invalid_argument("perspective camera aspect must be positive finite");
  }
  if (!(params.focusDistance > 0) || !std::isfinite(params.focusDistance)) {
    throw std::invalid_argument("perspective camera focusDistance must be positive finite");
  }
  if (params.aperture < 0 || !std::isfinite(params.aperture)) {
    throw std::invalid_argument("perspective camera aperture must be finite non-negative");
  }

  Geometry geometry;
  geometry.origin = params.lookfrom;
  BuildCameraFrame(params.lookfrom, params.lookat, params.up, geometry.forward, geometry.right, geometry.up);
  Float halfHeight = std::tan(DegreesToRadians(params.vfov) / static_cast<Float>(2));
  Float halfWidth = params.aspect * halfHeight;
  geometry.focusDistance = params.focusDistance;
  geometry.lensRadius = params.aperture / static_cast<Float>(2);
  geometry.lowerLeft = geometry.origin -
                       halfWidth * geometry.focusDistance * geometry.right -
                       halfHeight * geometry.focusDistance * geometry.up +
                       geometry.focusDistance * geometry.forward;
  geometry.horizontal = static_cast<Float>(2) * halfWidth * geometry.focusDistance * geometry.right;
  geometry.vertical = static_cast<Float>(2) * halfHeight * geometry.focusDistance * geometry.up;
  return Create(CameraType::Perspective, geometry, NormalizeOptions(params.options));
}

SpectralCamera SpectralCamera::Orthographic(const OrthographicCameraParameters& params) {
  if (!(params.width > 0) || !std::isfinite(params.width)) {
    throw std::invalid_argument("orthographic camera width must be positive finite");
  }
  if (!(params.height > 0) || !std::isfinite(params.height)) {
    throw std::invalid_argument("orthographic camera height must be positive finite");
  }
  if (!(params.focusDistance > 0) || !std::isfinite(params.focusDistance)) {
    throw std::invalid_argument("orthographic camera focusDistance must be positive finite");
  }
  if (params.aperture < 0 || !std::isfinite(params.aperture)) {
    throw std::invalid_argument("orthographic camera aperture must be finite non-negative");
  }

  Geometry geometry;
  geometry.origin = params.lookfrom;
  BuildCameraFrame(params.lookfrom, params.lookat, params.up, geometry.forward, geometry.right, geometry.up);
  geometry.focusDistance = params.focusDistance;
  geometry.lensRadius = params.aperture / static_cast<Float>(2);
  geometry.lowerLeft = geometry.origin -
                       params.width / static_cast<Float>(2) * geometry.right -
                       params.height / static_cast<Float>(2) * geometry.up;
  geometry.horizontal = params.width * geometry.right;
  geometry.vertical = params.height * geometry.up;
  return Create(CameraType::Orthographic, geometry, NormalizeOptions(params.options));
}

SpectralCamera SpectralCamera::Create(
  CameraType type,
  Geometry geometry,
  SpectralCameraOptions options
) {
  SpectralCamera camera;
  camera.type_ = type;
  camera.geometry_ = geometry;
  camera.options_ = std::move(options);
  return camera;
}

std::optional<CameraRay> SpectralCamera::GenerateRay(
  const CameraSample& sample,
  base::SampledWavelengths& lambda
) const {
  base::SampledWavelengths before = lambda;
  std::optional<CameraRay> ray = GenerateBaseRay(sample);
  if (!WavelengthsUnchanged(before, lambda)) {
    throw std::logic_error("nondispersive spectral camera modified sampled wavelengths");
  }
  return ray;
}

std::optional<CameraRay> SpectralCamera::GenerateRayDifferential(
  const CameraSample& sample,
  base::SampledWavelengths& lambda
) const {
  base::SampledWavelengths before = lambda;
  std::optional<CameraRay> ray = GenerateBaseRay(sample);
  if (!ray || !options_.enableDifferentials) {
    return ray;
  }

  CameraSample xSample = sample;
  xSample.pFilm.xy.x += 1;
  CameraSample ySample = sample;
  ySample.pFilm.xy.y += 1;
  std::optional<CameraRay> rx = GenerateBaseRay(xSample);
  std::optional<CameraRay> ry = GenerateBaseRay(ySample);
  ray->hasDifferentials = rx && ry;
  if (ray->hasDifferentials) {
    ray->rx = rx->ray;
    ray->ry = ry->ray;
  }
  if (!WavelengthsUnchanged(before, lambda)) {
    throw std::logic_error("nondispersive spectral camera modified sampled wavelengths");
  }
  return ray;
}

CameraType SpectralCamera::Type() const {
  return type_;
}

const SpectralCameraOptions& SpectralCamera::Options() const {
  return options_;
}

point3f SpectralCamera::Origin() const {
  return geometry_.origin;
}

vec3f SpectralCamera::Forward() const {
  return geometry_.forward;
}

vec3f SpectralCamera::Right() const {
  return geometry_.right;
}

vec3f SpectralCamera::Up() const {
  return geometry_.up;
}

std::optional<CameraRay> SpectralCamera::GenerateBaseRay(const CameraSample& sample) const {
  ValidatePoint(sample.pFilm, "camera sample pFilm");
  ValidatePoint(sample.pLens, "camera sample pLens");
  ValidateFinite(sample.time, "camera sample time");

  Float time = SampleTime(sample.time);
  point3f filmPoint = FilmPoint(sample.pFilm);
  point2f disk = geometry_.lensRadius > 0 ? SampleUniformDiskConcentric(sample.pLens) : point2f(0, 0);
  point3f lensPoint = geometry_.origin +
                      geometry_.lensRadius * disk.xy.x * geometry_.right +
                      geometry_.lensRadius * disk.xy.y * geometry_.up;

  if (type_ == CameraType::Orthographic) {
    if (geometry_.lensRadius > 0) {
      point3f focusPoint = filmPoint + geometry_.focusDistance * geometry_.forward;
      vec3f direction = SafeNormalize(focusPoint - lensPoint, "orthographic camera lens ray direction");
      return AttachCameraState(Ray(lensPoint, direction, time));
    }
    return AttachCameraState(Ray(filmPoint, geometry_.forward, time));
  }

  point3f rayOrigin = geometry_.lensRadius > 0 ? lensPoint : geometry_.origin;
  vec3f direction = SafeNormalize(filmPoint - rayOrigin, "perspective camera ray direction");
  return AttachCameraState(Ray(rayOrigin, direction, time));
}

point3f SpectralCamera::FilmPoint(const point2f& pFilm) const {
  Float s = pFilm.xy.x / static_cast<Float>(options_.filmWidth);
  Float t = pFilm.xy.y / static_cast<Float>(options_.filmHeight);
  return geometry_.lowerLeft + s * geometry_.horizontal + t * geometry_.vertical;
}

Float SpectralCamera::SampleTime(Float u) const {
  return (static_cast<Float>(1) - u) * options_.shutterOpen + u * options_.shutterClose;
}

CameraRay SpectralCamera::AttachCameraState(Ray ray) const {
  CameraRay cameraRay;
  cameraRay.ray = ray;
  cameraRay.weight = base::SampledSpectrum(1);
  cameraRay.medium = options_.medium;
  cameraRay.hasInitialMedium = options_.medium.IsValid();
  cameraRay.regionInitialization = options_.regionInitialization;
  return cameraRay;
}

CameraHandle SpectralCameraTable::Add(SpectralCamera camera) {
  CameraHandle handle = CameraHandle::FromIndex(static_cast<CameraHandle::IndexType>(cameras_.size()), 1);
  cameras_.push_back(std::move(camera));
  generations_.push_back(handle.Generation());
  return handle;
}

const SpectralCamera& SpectralCameraTable::Get(CameraHandle handle) const {
  if (!handle.IsValid() || handle.Index() >= cameras_.size() ||
      generations_[handle.Index()] != handle.Generation()) {
    throw std::out_of_range("invalid spectral camera handle");
  }
  return cameras_[handle.Index()];
}

std::optional<CameraRay> SpectralCameraTable::GenerateRay(
  CameraHandle handle,
  const CameraSample& sample,
  base::SampledWavelengths& lambda
) const {
  return Get(handle).GenerateRay(sample, lambda);
}

std::optional<CameraRay> SpectralCameraTable::GenerateRayDifferential(
  CameraHandle handle,
  const CameraSample& sample,
  base::SampledWavelengths& lambda
) const {
  return Get(handle).GenerateRayDifferential(sample, lambda);
}

std::size_t SpectralCameraTable::Size() const {
  return cameras_.size();
}

CameraFilmSample GenerateCameraRayFromFilm(
  const Film& film,
  const SpectralCamera& camera,
  const CameraSample& sample,
  Float wavelengthSample,
  bool generateDifferentials
) {
  CameraFilmSample result;
  result.wavelengths = film.SampleWavelengths(wavelengthSample);
  if (generateDifferentials) {
    result.cameraRay = camera.GenerateRayDifferential(sample, result.wavelengths);
  } else {
    result.cameraRay = camera.GenerateRay(sample, result.wavelengths);
  }
  return result;
}

} // namespace render
} // namespace rayrender
