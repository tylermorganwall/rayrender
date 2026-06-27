#ifndef RAYRENDER_RENDER_SPECTRAL_CAMERA_H
#define RAYRENDER_RENDER_SPECTRAL_CAMERA_H

#include "spectral_film.h"

#include "../base/base.h"
#include "../core/ray.h"
#include "../math/vectypes.h"

#include <cstdint>
#include <optional>
#include <vector>

namespace rayrender {
namespace render {

struct CameraSample {
  point2f pFilm;
  point2f pLens;
  Float time = 0;
  Float filterWeight = 1;
};

class CameraHandle {
public:
  using IndexType = std::uint32_t;
  using GenerationType = std::uint32_t;

  CameraHandle() = default;

  static constexpr CameraHandle Invalid() {
    return CameraHandle();
  }

  static constexpr CameraHandle FromIndex(IndexType index, GenerationType generation = 1) {
    return CameraHandle(index, generation);
  }

  static constexpr IndexType InvalidIndex() {
    return static_cast<IndexType>(-1);
  }

  constexpr bool IsValid() const {
    return index_ != InvalidIndex();
  }

  constexpr IndexType Index() const {
    return index_;
  }

  constexpr GenerationType Generation() const {
    return generation_;
  }

  constexpr bool operator==(const CameraHandle& other) const {
    return index_ == other.index_ && generation_ == other.generation_;
  }

  constexpr bool operator!=(const CameraHandle& other) const {
    return !(*this == other);
  }

private:
  constexpr CameraHandle(IndexType index, GenerationType generation) : index_(index), generation_(generation) {}

  IndexType index_ = InvalidIndex();
  GenerationType generation_ = 0;
};

enum class CameraType {
  Perspective,
  Orthographic
};

enum class CameraRegionInitializationMode {
  None,
  ExplicitRegions,
  DeferredContainment
};

struct CameraRegionInitialization {
  CameraRegionInitializationMode mode = CameraRegionInitializationMode::None;
  std::vector<int> explicitRegionIds;
  bool requiresContainmentQuery = false;
};

struct CameraRay {
  Ray ray;
  base::SampledSpectrum weight = base::SampledSpectrum(1);
  base::MediumHandle medium = base::MediumHandle::Invalid();
  bool hasInitialMedium = false;
  CameraRegionInitialization regionInitialization;
  bool hasDifferentials = false;
  Ray rx;
  Ray ry;
};

struct CameraFilmSample {
  base::SampledWavelengths wavelengths;
  std::optional<CameraRay> cameraRay;
};

struct SpectralCameraOptions {
  int filmWidth = 1;
  int filmHeight = 1;
  Float shutterOpen = 0;
  Float shutterClose = 1;
  base::MediumHandle medium = base::MediumHandle::Invalid();
  CameraRegionInitialization regionInitialization;
  bool enableDifferentials = true;
};

struct PerspectiveCameraParameters {
  point3f lookfrom;
  point3f lookat;
  vec3f up = vec3f(0, 1, 0);
  Float vfov = 40;
  Float aspect = 1;
  Float aperture = 0;
  Float focusDistance = 1;
  SpectralCameraOptions options;
};

struct OrthographicCameraParameters {
  point3f lookfrom;
  point3f lookat;
  vec3f up = vec3f(0, 1, 0);
  Float width = 2;
  Float height = 2;
  Float aperture = 0;
  Float focusDistance = static_cast<Float>(1e6);
  SpectralCameraOptions options;
};

class SpectralCamera {
public:
  static SpectralCamera Perspective(const PerspectiveCameraParameters& params);
  static SpectralCamera Orthographic(const OrthographicCameraParameters& params);

  std::optional<CameraRay> GenerateRay(
    const CameraSample& sample,
    base::SampledWavelengths& lambda
  ) const;

  std::optional<CameraRay> GenerateRayDifferential(
    const CameraSample& sample,
    base::SampledWavelengths& lambda
  ) const;

  CameraType Type() const;
  const SpectralCameraOptions& Options() const;
  point3f Origin() const;
  vec3f Forward() const;
  vec3f Right() const;
  vec3f Up() const;

private:
  struct Geometry {
    point3f origin;
    vec3f forward;
    vec3f right;
    vec3f up;
    point3f lowerLeft;
    vec3f horizontal;
    vec3f vertical;
    Float lensRadius = 0;
    Float focusDistance = 1;
  };

  static SpectralCamera Create(CameraType type, Geometry geometry, SpectralCameraOptions options);

  std::optional<CameraRay> GenerateBaseRay(const CameraSample& sample) const;
  point3f FilmPoint(const point2f& pFilm) const;
  Float SampleTime(Float u) const;
  CameraRay AttachCameraState(Ray ray) const;

  CameraType type_ = CameraType::Perspective;
  Geometry geometry_;
  SpectralCameraOptions options_;
};

class SpectralCameraTable {
public:
  CameraHandle Add(SpectralCamera camera);
  const SpectralCamera& Get(CameraHandle handle) const;

  std::optional<CameraRay> GenerateRay(
    CameraHandle handle,
    const CameraSample& sample,
    base::SampledWavelengths& lambda
  ) const;

  std::optional<CameraRay> GenerateRayDifferential(
    CameraHandle handle,
    const CameraSample& sample,
    base::SampledWavelengths& lambda
  ) const;

  std::size_t Size() const;

private:
  std::vector<SpectralCamera> cameras_;
  std::vector<CameraHandle::GenerationType> generations_;
};

CameraFilmSample GenerateCameraRayFromFilm(
  const Film& film,
  const SpectralCamera& camera,
  const CameraSample& sample,
  Float wavelengthSample,
  bool generateDifferentials = false
);

} // namespace render
} // namespace rayrender

#endif
