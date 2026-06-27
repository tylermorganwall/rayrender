#ifndef RAYRENDER_RENDER_SPECTRAL_FILM_H
#define RAYRENDER_RENDER_SPECTRAL_FILM_H

#include "../base/base.h"

#include <array>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace rayrender {
namespace render {

struct FilmPoint2i {
  int x = 0;
  int y = 0;
};

struct FilmPoint2f {
  Float x = 0;
  Float y = 0;
};

struct FilmBounds2i {
  FilmPoint2i min;
  FilmPoint2i max;

  static FilmBounds2i Full(int width, int height);

  int Width() const;
  int Height() const;
  int Area() const;
  bool IsEmpty() const;
  bool Contains(FilmPoint2i p) const;
};

FilmBounds2i Intersect(const FilmBounds2i& lhs, const FilmBounds2i& rhs);

enum class PixelSensorMode {
  CIE1931,
  RGBColorSpace,
  Measured
};

enum class WavelengthSamplingMode {
  Visible,
  Uniform
};

enum class FilmFilterType {
  Box,
  Triangle
};

enum class FilmOutputEncoding {
  Linear,
  SRGB
};

class SensorResponseCurve {
public:
  SensorResponseCurve() = default;
  SensorResponseCurve(int lambdaMin, int lambdaMax, std::vector<Float> values);

  static SensorResponseCurve FromSpectrum(
    const base::Spectrum& spectrum,
    int lambdaMin = static_cast<int>(base::LambdaMin),
    int lambdaMax = static_cast<int>(base::LambdaMax)
  );

  template <typename F>
  static SensorResponseCurve SampleFunction(
    F&& function,
    int lambdaMin = static_cast<int>(base::LambdaMin),
    int lambdaMax = static_cast<int>(base::LambdaMax)
  ) {
    std::vector<Float> values;
    values.reserve(static_cast<std::size_t>(lambdaMax - lambdaMin + 1));
    for (int lambda = lambdaMin; lambda <= lambdaMax; ++lambda) {
      values.push_back(static_cast<Float>(function(static_cast<Float>(lambda))));
    }
    return SensorResponseCurve(lambdaMin, lambdaMax, std::move(values));
  }

  bool IsValid() const;
  Float operator()(Float lambda) const;
  base::SampledSpectrum Sample(const base::SampledWavelengths& lambda) const;

  int LambdaMinValue() const;
  int LambdaMaxValue() const;
  const std::vector<Float>& Values() const;

private:
  int lambdaMin_ = static_cast<int>(base::LambdaMin);
  int lambdaMax_ = static_cast<int>(base::LambdaMax);
  std::vector<Float> values_;
};

struct SensorCalibrationSwatch {
  base::Spectrum reflectance;
};

base::ColorSpaceMatrix3x3 WhiteBalanceBradford(
  const base::Chromaticity& sourceWhite,
  const base::Chromaticity& targetWhite
);

base::RGB ApplyMatrixToRGB(const base::ColorSpaceMatrix3x3& matrix, const base::RGB& rgb);
base::XYZ ApplyMatrixToXYZ(const base::ColorSpaceMatrix3x3& matrix, const base::RGB& rgb);

class PixelSensor {
public:
  PixelSensor() = default;

  static PixelSensor CreateCIE1931(
    const base::NamedSpectrumRegistry& registry,
    const base::RGBColorSpace& outputColorSpace,
    const base::Spectrum* whiteBalanceIlluminant = nullptr,
    Float exposureTime = 1,
    Float iso = 100
  );

  static PixelSensor CreateRGBColorSpace(
    const base::NamedSpectrumRegistry& registry,
    const base::RGBColorSpace& sensorColorSpace,
    const base::RGBColorSpace& outputColorSpace,
    const base::Spectrum* whiteBalanceIlluminant = nullptr,
    Float exposureTime = 1,
    Float iso = 100
  );

  static PixelSensor CreateMeasured(
    SensorResponseCurve red,
    SensorResponseCurve green,
    SensorResponseCurve blue,
    base::ColorSpaceMatrix3x3 xyzFromSensorRGB,
    const base::NamedSpectrumRegistry& registry,
    const base::RGBColorSpace& outputColorSpace,
    const base::Spectrum* whiteBalanceIlluminant = nullptr,
    Float exposureTime = 1,
    Float iso = 100
  );

  static PixelSensor CreateCalibratedRGB(
    SensorResponseCurve red,
    SensorResponseCurve green,
    SensorResponseCurve blue,
    const std::vector<SensorCalibrationSwatch>& swatches,
    const base::Spectrum& sensorIlluminant,
    const base::RGBColorSpace& outputColorSpace,
    const base::NamedSpectrumRegistry& registry,
    Float exposureTime = 1,
    Float iso = 100
  );

  static base::ColorSpaceMatrix3x3 CalibrateXYZFromSensorRGB(
    const std::array<SensorResponseCurve, 3>& response,
    const std::vector<SensorCalibrationSwatch>& swatches,
    const base::Spectrum& sensorIlluminant,
    const base::Spectrum& outputIlluminant,
    const base::NamedSpectrumRegistry& registry
  );

  base::RGB ToSensorRGB(base::SampledSpectrum L, const base::SampledWavelengths& lambda) const;
  base::RGB DeterministicSensorRGB(const base::Spectrum& radiance) const;
  base::XYZ DeterministicXYZ(const base::Spectrum& radiance) const;

  bool IsValid() const;
  PixelSensorMode Mode() const;
  Float ImagingRatio() const;
  const base::ColorSpaceMatrix3x3& XYZFromSensorRGB() const;
  const std::array<SensorResponseCurve, 3>& Responses() const;
  const std::string& Name() const;

private:
  PixelSensor(
    std::array<SensorResponseCurve, 3> response,
    base::ColorSpaceMatrix3x3 xyzFromSensorRGB,
    PixelSensorMode mode,
    Float imagingRatio,
    std::string name
  );

  static Float ComputeImagingRatio(Float exposureTime, Float iso);

  std::array<SensorResponseCurve, 3> response_;
  base::ColorSpaceMatrix3x3 xyzFromSensorRGB_ = base::ColorSpaceMatrix3x3::Identity();
  PixelSensorMode mode_ = PixelSensorMode::CIE1931;
  Float imagingRatio_ = 1;
  std::string name_ = "cie1931";
};

struct FilmFilter {
  FilmFilterType type = FilmFilterType::Box;
  Float radiusX = static_cast<Float>(0.5);
  Float radiusY = static_cast<Float>(0.5);

  static FilmFilter Box(Float radiusX = static_cast<Float>(0.5), Float radiusY = static_cast<Float>(0.5));
  static FilmFilter Triangle(Float radiusX = 1, Float radiusY = 1);

  Float Evaluate(Float dx, Float dy) const;
};

struct SpectralVisibleSurface {
  bool valid = false;
  Float depth = 0;
  int primitiveId = -1;
  base::RGB albedo;
  bool hasAlbedo = false;
};

struct FilmOptions {
  int width = 0;
  int height = 0;
  PixelSensor sensor;
  base::RGBColorSpace outputColorSpace = base::RGBColorSpace::SRGB();
  FilmFilter filter = FilmFilter::Box();
  Float maxComponentValue = std::numeric_limits<Float>::infinity();
  WavelengthSamplingMode wavelengthSampling = WavelengthSamplingMode::Visible;
  bool deterministicSingleThread = true;
};

class FilmTile;

class Film {
public:
  explicit Film(FilmOptions options);

  base::SampledWavelengths SampleWavelengths(Float u) const;
  void AddSample(
    FilmPoint2i p,
    base::SampledSpectrum L,
    const base::SampledWavelengths& lambda,
    const SpectralVisibleSurface* visibleSurface,
    Float sampleWeight = 1
  );
  void AddFilteredSample(
    FilmPoint2f p,
    base::SampledSpectrum L,
    const base::SampledWavelengths& lambda,
    const SpectralVisibleSurface* visibleSurface,
    Float sampleWeight = 1
  );

  FilmTile CreateTile(FilmBounds2i bounds) const;
  void MergeTile(const FilmTile& tile);

  base::RGB ToOutputRGB(base::SampledSpectrum L, const base::SampledWavelengths& lambda) const;
  base::RGB GetPixelSensorRGB(FilmPoint2i p) const;
  base::RGB GetPixelLinearRGB(FilmPoint2i p) const;
  base::RGB GetPixelEncodedRGB(FilmPoint2i p, FilmOutputEncoding encoding) const;
  base::RGB GetPreviewLinearRGB(FilmPoint2i p) const;
  Float PixelWeightSum(FilmPoint2i p) const;
  const SpectralVisibleSurface* GetVisibleSurface(FilmPoint2i p) const;

  int Width() const;
  int Height() const;
  FilmBounds2i Bounds() const;
  const PixelSensor& Sensor() const;
  const base::RGBColorSpace& OutputColorSpace() const;
  const FilmFilter& Filter() const;
  bool DeterministicSingleThread() const;

private:
  friend class FilmTile;

  struct Pixel {
    double sensorSum[3] = {0, 0, 0};
    double weightSum = 0;
    SpectralVisibleSurface visibleSurface;
    bool hasVisibleSurface = false;
  };

  std::size_t PixelIndex(FilmPoint2i p) const;
  void AddSensorRGB(FilmPoint2i p, base::RGB sensorRGB, const SpectralVisibleSurface* visibleSurface, Float weight);
  base::RGB ClampSensorRGB(base::RGB sensorRGB) const;

  FilmOptions options_;
  base::ColorSpaceMatrix3x3 outputRGBFromSensorRGB_ = base::ColorSpaceMatrix3x3::Identity();
  std::vector<Pixel> pixels_;
};

class FilmTile {
public:
  FilmTile(
    FilmBounds2i bounds,
    const PixelSensor* sensor,
    FilmFilter filter,
    Float maxComponentValue
  );

  void AddSample(
    FilmPoint2i p,
    base::SampledSpectrum L,
    const base::SampledWavelengths& lambda,
    const SpectralVisibleSurface* visibleSurface,
    Float sampleWeight = 1
  );
  void AddFilteredSample(
    FilmPoint2f p,
    base::SampledSpectrum L,
    const base::SampledWavelengths& lambda,
    const SpectralVisibleSurface* visibleSurface,
    Float sampleWeight = 1
  );

  FilmBounds2i Bounds() const;

private:
  friend class Film;

  struct Pixel {
    double sensorSum[3] = {0, 0, 0};
    double weightSum = 0;
    SpectralVisibleSurface visibleSurface;
    bool hasVisibleSurface = false;
  };

  std::size_t PixelIndex(FilmPoint2i p) const;
  void AddSensorRGB(FilmPoint2i p, base::RGB sensorRGB, const SpectralVisibleSurface* visibleSurface, Float weight);
  base::RGB ClampSensorRGB(base::RGB sensorRGB) const;

  FilmBounds2i bounds_;
  const PixelSensor* sensor_ = nullptr;
  FilmFilter filter_;
  Float maxComponentValue_ = std::numeric_limits<Float>::infinity();
  std::vector<Pixel> pixels_;
};

} // namespace render
} // namespace rayrender

#endif
