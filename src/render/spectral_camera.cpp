#include "spectral_camera.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
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

Float Sqr(Float value) {
  return value * value;
}

Float RadicalInverse(int base, std::uint64_t index) {
  Float invBase = static_cast<Float>(1) / static_cast<Float>(base);
  Float invBaseN = 1;
  Float reversed = 0;
  while (index > 0) {
    std::uint64_t digit = index % static_cast<std::uint64_t>(base);
    reversed += static_cast<Float>(digit) * (invBaseN *= invBase);
    index /= static_cast<std::uint64_t>(base);
  }
  return reversed;
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

bool Quadratic(Float a, Float b, Float c, Float* t0, Float* t1) {
  if (a == 0) {
    if (b == 0) {
      return false;
    }
    *t0 = *t1 = -c / b;
    return true;
  }

  Float discriminant = b * b - static_cast<Float>(4) * a * c;
  if (discriminant < 0) {
    return false;
  }
  Float rootDiscriminant = std::sqrt(discriminant);

  Float q = b < 0
    ? static_cast<Float>(-0.5) * (b - rootDiscriminant)
    : static_cast<Float>(-0.5) * (b + rootDiscriminant);
  if (q == 0) {
    *t0 = *t1 = 0;
    return true;
  }

  *t0 = q / a;
  *t1 = c / q;
  if (*t0 > *t1) {
    std::swap(*t0, *t1);
  }
  return true;
}

bool RefractCamera(const vec3f& wi, normal3f n, Float eta, vec3f* wt) {
  Float cosThetaI = dot(n, wi);
  if (cosThetaI < 0) {
    eta = static_cast<Float>(1) / eta;
    cosThetaI = -cosThetaI;
    n = -n;
  }

  Float sin2ThetaI = std::max(static_cast<Float>(0), static_cast<Float>(1) - Sqr(cosThetaI));
  Float sin2ThetaT = sin2ThetaI / Sqr(eta);
  if (sin2ThetaT >= 1) {
    return false;
  }

  Float cosThetaT = std::sqrt(std::max(static_cast<Float>(0), static_cast<Float>(1) - sin2ThetaT));
  *wt = -wi / eta + (cosThetaI / eta - cosThetaT) * convert_to_vec3(n);
  return true;
}

point3f LensPointFromCamera(const point3f& p) {
  return point3f(p.xyz.x, p.xyz.y, -p.xyz.z);
}

vec3f LensVectorFromCamera(const vec3f& v) {
  return vec3f(v.xyz.x, v.xyz.y, -v.xyz.z);
}

Ray LensRayFromCamera(const Ray& ray) {
  return Ray(LensPointFromCamera(ray.origin()), LensVectorFromCamera(ray.direction()), ray.time());
}

Ray CameraRayFromLens(const Ray& ray) {
  return Ray(LensPointFromCamera(ray.origin()), LensVectorFromCamera(ray.direction()), ray.time());
}

bool IntersectSphericalElement(Float radius, Float zCenter, const Ray& ray, Float* t, normal3f* n) {
  vec3f direction = ray.direction();
  point3f origin = ray.origin();
  point3f o = origin - vec3f(0, 0, zCenter);

  Float a = dot(direction, direction);
  Float b = static_cast<Float>(2) * dot(direction, o);
  Float c = dot(o, o) - radius * radius;
  Float t0 = 0;
  Float t1 = 0;
  if (!Quadratic(a, b, c, &t0, &t1)) {
    return false;
  }

  bool useCloserT = (direction.xyz.z > 0) ^ (radius < 0);
  *t = useCloserT ? std::min(t0, t1) : std::max(t0, t1);
  if (*t < 0) {
    return false;
  }

  *n = Faceforward(unit_vector(convert_to_normal3(o + *t * direction)), -direction);
  return true;
}

void ValidateLensElement(const RealisticCameraLensElement& element, std::size_t index) {
  ValidateFinite(element.curvatureRadius, "realistic camera lens curvatureRadius");
  ValidateFinite(element.thickness, "realistic camera lens thickness");
  ValidateFinite(element.apertureRadius, "realistic camera lens apertureRadius");
  if (element.thickness < 0) {
    throw std::invalid_argument("realistic camera lens thickness must be non-negative");
  }
  if (!(element.apertureRadius > 0)) {
    throw std::invalid_argument("realistic camera lens apertureRadius must be positive");
  }
  if (element.IsStop()) {
    return;
  }

  const Float samples[] = {base::LambdaMin, static_cast<Float>(550), base::LambdaMax};
  for (Float lambda : samples) {
    try {
      (void)element.Eta(lambda);
    } catch (const std::exception& error) {
      throw std::invalid_argument(
        "realistic camera lens element " + std::to_string(index) +
        " has invalid eta spectrum: " + error.what()
      );
    }
  }
}

} // namespace

RealisticCameraLensElement RealisticCameraLensElement::Spherical(
  Float curvatureRadius,
  Float thickness,
  Float etaValue,
  Float apertureRadius
) {
  return Spherical(curvatureRadius, thickness, ConstantEtaSpectrum(etaValue), apertureRadius);
}

RealisticCameraLensElement RealisticCameraLensElement::Spherical(
  Float curvatureRadius,
  Float thickness,
  EtaSpectrumHandle etaValue,
  Float apertureRadius
) {
  RealisticCameraLensElement element;
  element.curvatureRadius = curvatureRadius;
  element.thickness = thickness;
  element.eta = std::move(etaValue);
  element.apertureRadius = apertureRadius;
  return element;
}

RealisticCameraLensElement RealisticCameraLensElement::ApertureStop(
  Float thickness,
  Float apertureRadius
) {
  RealisticCameraLensElement element;
  element.curvatureRadius = 0;
  element.thickness = thickness;
  element.eta = ConstantEtaSpectrum(1);
  element.apertureRadius = apertureRadius;
  return element;
}

bool RealisticCameraLensElement::IsStop() const {
  return curvatureRadius == 0;
}

bool RealisticCameraLensElement::IsDispersive() const {
  return !IsStop() && !EtaSpectrumIsConstant(eta);
}

Float RealisticCameraLensElement::Eta(Float lambdaNm) const {
  return IsStop() ? static_cast<Float>(1) : EvaluateEtaSpectrum(eta, lambdaNm);
}

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

SpectralCamera SpectralCamera::Realistic(const RealisticCameraParameters& params) {
  if (!(params.filmDiagonal > 0) || !std::isfinite(params.filmDiagonal)) {
    throw std::invalid_argument("realistic camera filmDiagonal must be positive finite");
  }
  if (params.lensElements.empty()) {
    throw std::invalid_argument("realistic camera requires at least one lens element");
  }
  if (!(params.lensElements.back().thickness > 0)) {
    throw std::invalid_argument("realistic camera rear element thickness must place the lens in front of the film");
  }

  Geometry geometry;
  geometry.origin = params.lookfrom;
  BuildCameraFrame(params.lookfrom, params.lookat, params.up, geometry.forward, geometry.right, geometry.up);
  geometry.filmDiagonal = params.filmDiagonal;
  geometry.lensElements = params.lensElements;
  for (std::size_t i = 0; i < geometry.lensElements.size(); ++i) {
    ValidateLensElement(geometry.lensElements[i], i);
    geometry.hasDispersiveLens = geometry.hasDispersiveLens || geometry.lensElements[i].IsDispersive();
  }
  SpectralCamera camera = Create(CameraType::Realistic, geometry, NormalizeOptions(params.options));
  camera.PrecomputeExitPupilBounds();
  return camera;
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

void SpectralCamera::PrecomputeExitPupilBounds() {
  if (type_ != CameraType::Realistic || geometry_.lensElements.empty()) {
    return;
  }

  constexpr int nBounds = 64;
  geometry_.exitPupilBounds.clear();
  geometry_.exitPupilBounds.reserve(nBounds);
  Float halfFilmDiagonal = geometry_.filmDiagonal / static_cast<Float>(2);
  for (int i = 0; i < nBounds; ++i) {
    Float r0 = static_cast<Float>(i) / static_cast<Float>(nBounds) * halfFilmDiagonal;
    Float r1 = static_cast<Float>(i + 1) / static_cast<Float>(nBounds) * halfFilmDiagonal;
    geometry_.exitPupilBounds.push_back(BoundExitPupil(r0, r1));
  }
}

std::optional<CameraRay> SpectralCamera::GenerateRay(
  const CameraSample& sample,
  base::SampledWavelengths& lambda
) const {
  base::SampledWavelengths before = lambda;
  std::optional<CameraRay> ray = GenerateBaseRay(sample, lambda);
  if (!geometry_.hasDispersiveLens && !WavelengthsUnchanged(before, lambda)) {
    throw std::logic_error("nondispersive spectral camera modified sampled wavelengths");
  }
  return ray;
}

std::optional<CameraRay> SpectralCamera::GenerateRayDifferential(
  const CameraSample& sample,
  base::SampledWavelengths& lambda
) const {
  base::SampledWavelengths before = lambda;
  std::optional<CameraRay> ray = GenerateBaseRay(sample, lambda);
  if (!ray || !options_.enableDifferentials) {
    return ray;
  }

  CameraSample xSample = sample;
  xSample.pFilm.xy.x += 1;
  CameraSample ySample = sample;
  ySample.pFilm.xy.y += 1;
  std::optional<CameraRay> rx = GenerateBaseRay(xSample, lambda);
  std::optional<CameraRay> ry = GenerateBaseRay(ySample, lambda);
  ray->hasDifferentials = rx && ry;
  if (ray->hasDifferentials) {
    ray->rx = rx->ray;
    ray->ry = ry->ray;
  }
  if (!geometry_.hasDispersiveLens && !WavelengthsUnchanged(before, lambda)) {
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

std::optional<CameraRay> SpectralCamera::GenerateBaseRay(
  const CameraSample& sample,
  base::SampledWavelengths& lambda
) const {
  if (type_ == CameraType::Realistic) {
    return GenerateRealisticRay(sample, lambda);
  }
  return GenerateProjectiveRay(sample);
}

std::optional<CameraRay> SpectralCamera::GenerateProjectiveRay(const CameraSample& sample) const {
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

std::optional<CameraRay> SpectralCamera::GenerateRealisticRay(
  const CameraSample& sample,
  base::SampledWavelengths& lambda
) const {
  ValidatePoint(sample.pFilm, "camera sample pFilm");
  ValidatePoint(sample.pLens, "camera sample pLens");
  ValidateFinite(sample.time, "camera sample time");
  if (geometry_.lensElements.empty()) {
    return std::nullopt;
  }
  if (geometry_.hasDispersiveLens) {
    lambda.TerminateSecondary();
  }

  Float s = sample.pFilm.xy.x / static_cast<Float>(options_.filmWidth);
  Float t = sample.pFilm.xy.y / static_cast<Float>(options_.filmHeight);
  point2f pFilm2 = PhysicalExtent().Lerp(point2f(s, t));
  point3f pFilm(-pFilm2.xy.x, pFilm2.xy.y, 0);

  Float rearZ = LensRearZ();
  if (!(RearElementRadius() > 0) || !(rearZ > 0)) {
    return std::nullopt;
  }
  Float exitPupilBoundsArea = 0;
  point3f pRear = SampleExitPupil(point2f(pFilm.xyz.x, pFilm.xyz.y), sample.pLens, &exitPupilBoundsArea);
  if (!(exitPupilBoundsArea > 0) || !std::isfinite(exitPupilBoundsArea)) {
    return std::nullopt;
  }
  vec3f filmDirection = SafeNormalize(pRear - pFilm, "realistic camera film ray direction");

  Float time = SampleTime(sample.time);
  Ray rFilm(pFilm, filmDirection, time);
  Ray cameraSpaceRay;
  if (!TraceLensesFromFilm(rFilm, lambda[0], &cameraSpaceRay)) {
    return std::nullopt;
  }

  vec3f renderDirection = SafeNormalize(
    CameraToRenderVector(cameraSpaceRay.direction()),
    "realistic camera render ray direction"
  );
  CameraRay result = AttachCameraState(
    Ray(CameraToRenderPoint(cameraSpaceRay.origin()), renderDirection, time)
  );

  Float cosTheta = filmDirection.xyz.z;
  Float cos4Theta = Sqr(Sqr(cosTheta));
  Float weight = cos4Theta * exitPupilBoundsArea / (rearZ * rearZ);
  if (!(weight >= 0) || !std::isfinite(weight)) {
    return std::nullopt;
  }
  result.weight = base::SampledSpectrum(weight);
  return result;
}

bool SpectralCamera::TraceLensesFromFilm(const Ray& ray, Float lambdaNm, Ray* out) const {
  Float elementZ = 0;
  Ray lensRay = LensRayFromCamera(ray);

  for (int i = static_cast<int>(geometry_.lensElements.size()) - 1; i >= 0; --i) {
    const RealisticCameraLensElement& element = geometry_.lensElements[static_cast<std::size_t>(i)];
    elementZ -= element.thickness;

    Float t = 0;
    normal3f n;
    if (element.IsStop()) {
      if (lensRay.direction().xyz.z >= 0) {
        return false;
      }
      t = (elementZ - lensRay.origin().xyz.z) / lensRay.direction().xyz.z;
      if (t < 0 || !std::isfinite(t)) {
        return false;
      }
    } else {
      Float zCenter = elementZ + element.curvatureRadius;
      if (!IntersectSphericalElement(element.curvatureRadius, zCenter, lensRay, &t, &n)) {
        return false;
      }
    }

    point3f pHit = lensRay(t);
    Float hitRadius2 = pHit.xyz.x * pHit.xyz.x + pHit.xyz.y * pHit.xyz.y;
    if (hitRadius2 > element.apertureRadius * element.apertureRadius) {
      return false;
    }
    lensRay.o = pHit;

    if (!element.IsStop()) {
      Float etaI = element.Eta(lambdaNm);
      Float etaT = 1;
      if (i > 0) {
        const RealisticCameraLensElement& previous =
          geometry_.lensElements[static_cast<std::size_t>(i - 1)];
        etaT = previous.IsStop() ? static_cast<Float>(1) : previous.Eta(lambdaNm);
      }

      vec3f refracted;
      vec3f incident = SafeNormalize(-lensRay.direction(), "realistic camera incident lens direction");
      if (!RefractCamera(incident, n, etaT / etaI, &refracted)) {
        return false;
      }
      lensRay.d = refracted;
    }
  }

  if (out != nullptr) {
    *out = CameraRayFromLens(Ray(lensRay.origin(), lensRay.direction(), ray.time()));
  }
  return true;
}

Bounds2f SpectralCamera::BoundExitPupil(Float pFilmX0, Float pFilmX1) const {
  Bounds2f pupilBounds;
  constexpr int nSamples = 4096;
  int nExitingRays = 0;

  Float rearRadius = RearElementRadius();
  Bounds2f projectedRearBounds(
    point2f(static_cast<Float>(-1.5) * rearRadius, static_cast<Float>(-1.5) * rearRadius),
    point2f(static_cast<Float>(1.5) * rearRadius, static_cast<Float>(1.5) * rearRadius)
  );

  for (int i = 0; i < nSamples; ++i) {
    Float uFilm = (static_cast<Float>(i) + static_cast<Float>(0.5)) / static_cast<Float>(nSamples);
    point3f pFilm(
      base::Lerp(uFilm, pFilmX0, pFilmX1),
      0,
      0
    );
    point2f u(RadicalInverse(2, static_cast<std::uint64_t>(i)), RadicalInverse(3, static_cast<std::uint64_t>(i)));
    point2f pRear2 = projectedRearBounds.Lerp(u);
    point3f pRear(pRear2.xy.x, pRear2.xy.y, LensRearZ());
    point2f pRearPoint(pRear.xyz.x, pRear.xyz.y);

    if (Inside(pRearPoint, pupilBounds) ||
        TraceLensesFromFilm(Ray(pFilm, SafeNormalize(pRear - pFilm, "exit pupil probe ray")), static_cast<Float>(550), nullptr)) {
      pupilBounds = UnionB(pupilBounds, pRearPoint);
      ++nExitingRays;
    }
  }

  if (nExitingRays == 0) {
    return projectedRearBounds;
  }
  return Expand(pupilBounds, static_cast<Float>(2) * projectedRearBounds.Diagonal().length() / std::sqrt(static_cast<Float>(nSamples)));
}

point3f SpectralCamera::SampleExitPupil(
  const point2f& pFilm,
  const point2f& lensSample,
  Float* sampleBoundsArea
) const {
  Float rFilm = std::sqrt(pFilm.xy.x * pFilm.xy.x + pFilm.xy.y * pFilm.xy.y);
  int rIndex = 0;
  if (!geometry_.exitPupilBounds.empty() && geometry_.filmDiagonal > 0) {
    rIndex = static_cast<int>(rFilm / (geometry_.filmDiagonal / static_cast<Float>(2)) *
                              static_cast<Float>(geometry_.exitPupilBounds.size()));
    rIndex = std::max(0, std::min(static_cast<int>(geometry_.exitPupilBounds.size()) - 1, rIndex));
  }

  Bounds2f pupilBounds = geometry_.exitPupilBounds.empty()
    ? Bounds2f(
        point2f(-RearElementRadius(), -RearElementRadius()),
        point2f(RearElementRadius(), RearElementRadius())
      )
    : geometry_.exitPupilBounds[static_cast<std::size_t>(rIndex)];

  if (sampleBoundsArea != nullptr) {
    *sampleBoundsArea = pupilBounds.Area();
  }

  point2f pLens = pupilBounds.Lerp(lensSample);
  Float sinTheta = rFilm != 0 ? pFilm.xy.y / rFilm : 0;
  Float cosTheta = rFilm != 0 ? pFilm.xy.x / rFilm : 1;
  return point3f(
    cosTheta * pLens.xy.x - sinTheta * pLens.xy.y,
    sinTheta * pLens.xy.x + cosTheta * pLens.xy.y,
    LensRearZ()
  );
}

Bounds2f SpectralCamera::PhysicalExtent() const {
  Float aspect = static_cast<Float>(options_.filmHeight) / static_cast<Float>(options_.filmWidth);
  Float x = std::sqrt(geometry_.filmDiagonal * geometry_.filmDiagonal /
                      (static_cast<Float>(1) + aspect * aspect));
  Float y = aspect * x;
  return Bounds2f(point2f(-x / static_cast<Float>(2), -y / static_cast<Float>(2)),
                  point2f(x / static_cast<Float>(2), y / static_cast<Float>(2)));
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

Float SpectralCamera::LensRearZ() const {
  return geometry_.lensElements.empty() ? static_cast<Float>(0) : geometry_.lensElements.back().thickness;
}

Float SpectralCamera::RearElementRadius() const {
  return geometry_.lensElements.empty() ? static_cast<Float>(0) : geometry_.lensElements.back().apertureRadius;
}

point3f SpectralCamera::CameraToRenderPoint(const point3f& p) const {
  return geometry_.origin +
         p.xyz.x * geometry_.right +
         p.xyz.y * geometry_.up +
         p.xyz.z * geometry_.forward;
}

vec3f SpectralCamera::CameraToRenderVector(const vec3f& v) const {
  return v.xyz.x * geometry_.right +
         v.xyz.y * geometry_.up +
         v.xyz.z * geometry_.forward;
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
