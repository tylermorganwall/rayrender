#include "spectral_film.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace rayrender {
namespace render {

namespace {

base::RGB RGBFromArray(const std::array<Float, 3>& value) {
  return base::RGB(value[0], value[1], value[2]);
}

base::XYZ XYZFromArray(const std::array<Float, 3>& value) {
  return base::XYZ(value[0], value[1], value[2]);
}

void ValidateFinite(Float value, const char* context) {
  if (!std::isfinite(value)) {
    throw std::invalid_argument(std::string(context) + " must be finite");
  }
}

const base::Spectrum& CIESpectrum(const base::NamedSpectrumRegistry& registry, const char* name) {
  return registry.GetOrThrow(name);
}

std::array<SensorResponseCurve, 3> CIEResponses(const base::NamedSpectrumRegistry& registry) {
  return {
    SensorResponseCurve::FromSpectrum(CIESpectrum(registry, "cie-x")),
    SensorResponseCurve::FromSpectrum(CIESpectrum(registry, "cie-y")),
    SensorResponseCurve::FromSpectrum(CIESpectrum(registry, "cie-z"))
  };
}

base::XYZ SpectrumXYZ(
  const base::Spectrum& spectrum,
  const base::NamedSpectrumRegistry& registry
) {
  return base::SpectrumToXYZ(
    spectrum,
    CIESpectrum(registry, "cie-x"),
    CIESpectrum(registry, "cie-y"),
    CIESpectrum(registry, "cie-z")
  );
}

base::ColorSpaceMatrix3x3 WhiteBalanceFromIlluminant(
  const base::Spectrum* whiteBalanceIlluminant,
  const base::RGBColorSpace& outputColorSpace,
  const base::NamedSpectrumRegistry& registry
) {
  if (whiteBalanceIlluminant == nullptr) {
    return base::ColorSpaceMatrix3x3::Identity();
  }
  return WhiteBalanceBradford(SpectrumXYZ(*whiteBalanceIlluminant, registry).xy(), outputColorSpace.w);
}

base::Spectrum OutputIlluminantSpectrum(
  const base::RGBColorSpace& outputColorSpace,
  const base::Spectrum& fallback
) {
  if (outputColorSpace.illuminant) {
    return base::Spectrum(*outputColorSpace.illuminant);
  }
  return fallback;
}

base::RGB ProjectReflectanceSensor(
  const base::Spectrum& reflectance,
  const base::Spectrum& illuminant,
  const std::array<SensorResponseCurve, 3>& response
) {
  Float gIntegral = 0;
  base::RGB result;
  for (int lambda = static_cast<int>(base::LambdaMin); lambda <= static_cast<int>(base::LambdaMax); ++lambda) {
    Float lambdaValue = static_cast<Float>(lambda);
    Float illum = illuminant(lambdaValue);
    Float refl = reflectance(lambdaValue);
    gIntegral += response[1](lambdaValue) * illum;
    result.r += response[0](lambdaValue) * refl * illum;
    result.g += response[1](lambdaValue) * refl * illum;
    result.b += response[2](lambdaValue) * refl * illum;
  }
  if (!(gIntegral != 0) || !std::isfinite(gIntegral)) {
    throw std::invalid_argument("sensor calibration illuminant has invalid green-channel integral");
  }
  return result / gIntegral;
}

base::XYZ ProjectReflectanceXYZ(
  const base::Spectrum& reflectance,
  const base::Spectrum& illuminant,
  const base::NamedSpectrumRegistry& registry
) {
  const base::Spectrum& x = CIESpectrum(registry, "cie-x");
  const base::Spectrum& y = CIESpectrum(registry, "cie-y");
  const base::Spectrum& z = CIESpectrum(registry, "cie-z");

  Float yIntegral = 0;
  base::XYZ result;
  for (int lambda = static_cast<int>(base::LambdaMin); lambda <= static_cast<int>(base::LambdaMax); ++lambda) {
    Float lambdaValue = static_cast<Float>(lambda);
    Float illum = illuminant(lambdaValue);
    Float refl = reflectance(lambdaValue);
    yIntegral += y(lambdaValue) * illum;
    result.x += x(lambdaValue) * refl * illum;
    result.y += y(lambdaValue) * refl * illum;
    result.z += z(lambdaValue) * refl * illum;
  }
  if (!(yIntegral != 0) || !std::isfinite(yIntegral)) {
    throw std::invalid_argument("sensor calibration illuminant has invalid CIE Y integral");
  }
  return result / yIntegral;
}

std::array<double, 3> Solve3x3(
  std::array<std::array<double, 3>, 3> matrix,
  std::array<double, 3> rhs
) {
  for (int column = 0; column < 3; ++column) {
    int pivot = column;
    for (int row = column + 1; row < 3; ++row) {
      if (std::fabs(matrix[row][column]) > std::fabs(matrix[pivot][column])) {
        pivot = row;
      }
    }
    if (std::fabs(matrix[pivot][column]) < 1e-14) {
      throw std::invalid_argument("sensor calibration matrix is singular");
    }
    if (pivot != column) {
      std::swap(matrix[pivot], matrix[column]);
      std::swap(rhs[pivot], rhs[column]);
    }
    double invPivot = 1.0 / matrix[column][column];
    for (int col = column; col < 3; ++col) {
      matrix[column][col] *= invPivot;
    }
    rhs[column] *= invPivot;
    for (int row = 0; row < 3; ++row) {
      if (row == column) {
        continue;
      }
      double factor = matrix[row][column];
      for (int col = column; col < 3; ++col) {
        matrix[row][col] -= factor * matrix[column][col];
      }
      rhs[row] -= factor * rhs[column];
    }
  }
  return rhs;
}

base::ColorSpaceMatrix3x3 LinearLeastSquaresSensorToXYZ(
  const std::vector<base::RGB>& sensorRGB,
  const std::vector<base::XYZ>& xyz
) {
  if (sensorRGB.size() != xyz.size() || sensorRGB.size() < 3) {
    throw std::invalid_argument("sensor calibration requires at least three matched swatches");
  }

  std::array<std::array<double, 3>, 3> ata{};
  std::array<std::array<double, 3>, 3> atb{};

  for (std::size_t i = 0; i < sensorRGB.size(); ++i) {
    std::array<double, 3> a = {sensorRGB[i].r, sensorRGB[i].g, sensorRGB[i].b};
    std::array<double, 3> b = {xyz[i].x, xyz[i].y, xyz[i].z};
    for (int row = 0; row < 3; ++row) {
      for (int column = 0; column < 3; ++column) {
        ata[row][column] += a[row] * a[column];
      }
      for (int output = 0; output < 3; ++output) {
        atb[row][output] += a[row] * b[output];
      }
    }
  }

  base::ColorSpaceMatrix3x3 result;
  for (int output = 0; output < 3; ++output) {
    std::array<double, 3> rhs = {atb[0][output], atb[1][output], atb[2][output]};
    std::array<double, 3> row = Solve3x3(ata, rhs);
    for (int input = 0; input < 3; ++input) {
      result(output, input) = static_cast<Float>(row[input]);
    }
  }
  return result;
}

} // namespace

FilmBounds2i FilmBounds2i::Full(int width, int height) {
  return FilmBounds2i{{0, 0}, {width, height}};
}

int FilmBounds2i::Width() const {
  return std::max(0, max.x - min.x);
}

int FilmBounds2i::Height() const {
  return std::max(0, max.y - min.y);
}

int FilmBounds2i::Area() const {
  return Width() * Height();
}

bool FilmBounds2i::IsEmpty() const {
  return Width() == 0 || Height() == 0;
}

bool FilmBounds2i::Contains(FilmPoint2i p) const {
  return p.x >= min.x && p.x < max.x && p.y >= min.y && p.y < max.y;
}

FilmBounds2i Intersect(const FilmBounds2i& lhs, const FilmBounds2i& rhs) {
  FilmBounds2i result{
    {std::max(lhs.min.x, rhs.min.x), std::max(lhs.min.y, rhs.min.y)},
    {std::min(lhs.max.x, rhs.max.x), std::min(lhs.max.y, rhs.max.y)}
  };
  if (result.max.x < result.min.x) {
    result.max.x = result.min.x;
  }
  if (result.max.y < result.min.y) {
    result.max.y = result.min.y;
  }
  return result;
}

SensorResponseCurve::SensorResponseCurve(int lambdaMin, int lambdaMax, std::vector<Float> values)
  : lambdaMin_(lambdaMin),
    lambdaMax_(lambdaMax),
    values_(std::move(values)) {
  if (lambdaMax_ < lambdaMin_) {
    throw std::invalid_argument("sensor response lambdaMax must be >= lambdaMin");
  }
  if (values_.size() != static_cast<std::size_t>(lambdaMax_ - lambdaMin_ + 1)) {
    throw std::invalid_argument("sensor response value count must match wavelength interval");
  }
  for (Float value : values_) {
    ValidateFinite(value, "sensor response value");
  }
}

SensorResponseCurve SensorResponseCurve::FromSpectrum(
  const base::Spectrum& spectrum,
  int lambdaMin,
  int lambdaMax
) {
  return SampleFunction([&](Float lambda) { return spectrum(lambda); }, lambdaMin, lambdaMax);
}

bool SensorResponseCurve::IsValid() const {
  return !values_.empty();
}

Float SensorResponseCurve::operator()(Float lambda) const {
  if (!std::isfinite(lambda)) {
    throw std::invalid_argument("sensor response wavelength must be finite");
  }
  if (values_.empty()) {
    return 0;
  }
  int offset = static_cast<int>(std::lround(lambda)) - lambdaMin_;
  if (offset < 0 || offset >= static_cast<int>(values_.size())) {
    return 0;
  }
  return values_[static_cast<std::size_t>(offset)];
}

base::SampledSpectrum SensorResponseCurve::Sample(const base::SampledWavelengths& lambda) const {
  base::SampledSpectrum result;
  for (int i = 0; i < base::NSpectrumSamples; ++i) {
    result[i] = (*this)(lambda[i]);
  }
  return result;
}

int SensorResponseCurve::LambdaMinValue() const {
  return lambdaMin_;
}

int SensorResponseCurve::LambdaMaxValue() const {
  return lambdaMax_;
}

const std::vector<Float>& SensorResponseCurve::Values() const {
  return values_;
}

base::ColorSpaceMatrix3x3 WhiteBalanceBradford(
  const base::Chromaticity& sourceWhite,
  const base::Chromaticity& targetWhite
) {
  // pbrt uses the Bradford LMS transform for PixelSensor white balance.
  const base::ColorSpaceMatrix3x3 lmsFromXYZ = base::ColorSpaceMatrix3x3::FromRows(
    {static_cast<Float>(0.8951), static_cast<Float>(0.2664), static_cast<Float>(-0.1614)},
    {static_cast<Float>(-0.7502), static_cast<Float>(1.7135), static_cast<Float>(0.0367)},
    {static_cast<Float>(0.0389), static_cast<Float>(-0.0685), static_cast<Float>(1.0296)}
  );
  const base::ColorSpaceMatrix3x3 xyzFromLMS = base::ColorSpaceMatrix3x3::FromRows(
    {static_cast<Float>(0.986993), static_cast<Float>(-0.147054), static_cast<Float>(0.159963)},
    {static_cast<Float>(0.432305), static_cast<Float>(0.51836), static_cast<Float>(0.0492912)},
    {static_cast<Float>(-0.00852866), static_cast<Float>(0.0400428), static_cast<Float>(0.968487)}
  );
  base::XYZ sourceXYZ = base::XYZ::FromxyY(sourceWhite);
  base::XYZ targetXYZ = base::XYZ::FromxyY(targetWhite);
  std::array<Float, 3> sourceLMS = lmsFromXYZ.Apply({sourceXYZ.x, sourceXYZ.y, sourceXYZ.z});
  std::array<Float, 3> targetLMS = lmsFromXYZ.Apply({targetXYZ.x, targetXYZ.y, targetXYZ.z});
  for (int i = 0; i < 3; ++i) {
    if (sourceLMS[i] == 0 || !std::isfinite(sourceLMS[i]) || !std::isfinite(targetLMS[i])) {
      throw std::invalid_argument("white balance source and target whites must be finite and nonzero");
    }
  }
  base::ColorSpaceMatrix3x3 correction = base::ColorSpaceMatrix3x3::Diagonal(
    targetLMS[0] / sourceLMS[0],
    targetLMS[1] / sourceLMS[1],
    targetLMS[2] / sourceLMS[2]
  );
  return xyzFromLMS * correction * lmsFromXYZ;
}

base::RGB ApplyMatrixToRGB(const base::ColorSpaceMatrix3x3& matrix, const base::RGB& rgb) {
  return RGBFromArray(matrix.Apply({rgb.r, rgb.g, rgb.b}));
}

base::XYZ ApplyMatrixToXYZ(const base::ColorSpaceMatrix3x3& matrix, const base::RGB& rgb) {
  return XYZFromArray(matrix.Apply({rgb.r, rgb.g, rgb.b}));
}

PixelSensor PixelSensor::CreateCIE1931(
  const base::NamedSpectrumRegistry& registry,
  const base::RGBColorSpace& outputColorSpace,
  const base::Spectrum* whiteBalanceIlluminant,
  Float exposureTime,
  Float iso
) {
  base::ColorSpaceMatrix3x3 xyzFromSensorRGB =
    WhiteBalanceFromIlluminant(whiteBalanceIlluminant, outputColorSpace, registry);
  return PixelSensor(
    CIEResponses(registry),
    xyzFromSensorRGB,
    PixelSensorMode::CIE1931,
    ComputeImagingRatio(exposureTime, iso),
    "cie1931"
  );
}

PixelSensor PixelSensor::CreateRGBColorSpace(
  const base::NamedSpectrumRegistry& registry,
  const base::RGBColorSpace& sensorColorSpace,
  const base::RGBColorSpace& outputColorSpace,
  const base::Spectrum* whiteBalanceIlluminant,
  Float exposureTime,
  Float iso
) {
  const base::Spectrum& x = CIESpectrum(registry, "cie-x");
  const base::Spectrum& y = CIESpectrum(registry, "cie-y");
  const base::Spectrum& z = CIESpectrum(registry, "cie-z");
  const base::ColorSpaceMatrix3x3& rgbFromXYZ = sensorColorSpace.xyzToRGB;

  std::array<SensorResponseCurve, 3> response = {
    SensorResponseCurve::SampleFunction([&](Float lambda) {
      return rgbFromXYZ(0, 0) * x(lambda) + rgbFromXYZ(0, 1) * y(lambda) + rgbFromXYZ(0, 2) * z(lambda);
    }),
    SensorResponseCurve::SampleFunction([&](Float lambda) {
      return rgbFromXYZ(1, 0) * x(lambda) + rgbFromXYZ(1, 1) * y(lambda) + rgbFromXYZ(1, 2) * z(lambda);
    }),
    SensorResponseCurve::SampleFunction([&](Float lambda) {
      return rgbFromXYZ(2, 0) * x(lambda) + rgbFromXYZ(2, 1) * y(lambda) + rgbFromXYZ(2, 2) * z(lambda);
    })
  };

  base::ColorSpaceMatrix3x3 xyzFromSensorRGB = sensorColorSpace.rgbToXYZ;
  xyzFromSensorRGB =
    WhiteBalanceFromIlluminant(whiteBalanceIlluminant, outputColorSpace, registry) * xyzFromSensorRGB;

  return PixelSensor(
    std::move(response),
    xyzFromSensorRGB,
    PixelSensorMode::RGBColorSpace,
    ComputeImagingRatio(exposureTime, iso),
    std::string("rgb:") + sensorColorSpace.name
  );
}

PixelSensor PixelSensor::CreateMeasured(
  SensorResponseCurve red,
  SensorResponseCurve green,
  SensorResponseCurve blue,
  base::ColorSpaceMatrix3x3 xyzFromSensorRGB,
  const base::NamedSpectrumRegistry& registry,
  const base::RGBColorSpace& outputColorSpace,
  const base::Spectrum* whiteBalanceIlluminant,
  Float exposureTime,
  Float iso
) {
  xyzFromSensorRGB =
    WhiteBalanceFromIlluminant(whiteBalanceIlluminant, outputColorSpace, registry) * xyzFromSensorRGB;
  return PixelSensor(
    {std::move(red), std::move(green), std::move(blue)},
    xyzFromSensorRGB,
    PixelSensorMode::Measured,
    ComputeImagingRatio(exposureTime, iso),
    "measured"
  );
}

PixelSensor PixelSensor::CreateCalibratedRGB(
  SensorResponseCurve red,
  SensorResponseCurve green,
  SensorResponseCurve blue,
  const std::vector<SensorCalibrationSwatch>& swatches,
  const base::Spectrum& sensorIlluminant,
  const base::RGBColorSpace& outputColorSpace,
  const base::NamedSpectrumRegistry& registry,
  Float exposureTime,
  Float iso
) {
  std::array<SensorResponseCurve, 3> response = {std::move(red), std::move(green), std::move(blue)};
  base::Spectrum outputIlluminant = OutputIlluminantSpectrum(outputColorSpace, sensorIlluminant);
  base::ColorSpaceMatrix3x3 xyzFromSensorRGB =
    CalibrateXYZFromSensorRGB(response, swatches, sensorIlluminant, outputIlluminant, registry);
  return PixelSensor(
    std::move(response),
    xyzFromSensorRGB,
    PixelSensorMode::Measured,
    ComputeImagingRatio(exposureTime, iso),
    "calibrated"
  );
}

base::ColorSpaceMatrix3x3 PixelSensor::CalibrateXYZFromSensorRGB(
  const std::array<SensorResponseCurve, 3>& response,
  const std::vector<SensorCalibrationSwatch>& swatches,
  const base::Spectrum& sensorIlluminant,
  const base::Spectrum& outputIlluminant,
  const base::NamedSpectrumRegistry& registry
) {
  if (swatches.size() < 3) {
    throw std::invalid_argument("sensor calibration requires at least three reflectance swatches");
  }
  std::vector<base::RGB> sensorRGB;
  std::vector<base::XYZ> xyz;
  sensorRGB.reserve(swatches.size());
  xyz.reserve(swatches.size());

  Float sensorWhiteG = base::InnerProduct(sensorIlluminant, CIESpectrum(registry, "cie-y"));
  Float outputWhiteY = base::InnerProduct(outputIlluminant, CIESpectrum(registry, "cie-y"));
  if (!(sensorWhiteG > 0) || !(outputWhiteY > 0)) {
    throw std::invalid_argument("sensor calibration illuminants must have positive CIE Y integrals");
  }
  Float normalization = outputWhiteY / sensorWhiteG;

  for (const SensorCalibrationSwatch& swatch : swatches) {
    // Mirrors pbrt PixelSensor calibration: project ColorChecker reflectances
    // under the sensor illuminant and solve sensor RGB -> output XYZ.
    sensorRGB.push_back(ProjectReflectanceSensor(swatch.reflectance, sensorIlluminant, response));
    xyz.push_back(ProjectReflectanceXYZ(swatch.reflectance, outputIlluminant, registry) * normalization);
  }
  return LinearLeastSquaresSensorToXYZ(sensorRGB, xyz);
}

base::RGB PixelSensor::ToSensorRGB(base::SampledSpectrum L, const base::SampledWavelengths& lambda) const {
  if (!IsValid()) {
    return base::RGB(0, 0, 0);
  }
  // pbrt PixelSensor::ToSensorRGB centralizes wavelength-PDF division here.
  L = base::SafeDiv(L, lambda.PDF());
  return imagingRatio_ * base::RGB(
    (response_[0].Sample(lambda) * L).Average(),
    (response_[1].Sample(lambda) * L).Average(),
    (response_[2].Sample(lambda) * L).Average()
  );
}

base::RGB PixelSensor::DeterministicSensorRGB(const base::Spectrum& radiance) const {
  base::RGB result;
  for (int lambda = static_cast<int>(base::LambdaMin); lambda <= static_cast<int>(base::LambdaMax); ++lambda) {
    Float lambdaValue = static_cast<Float>(lambda);
    Float value = radiance(lambdaValue);
    result.r += response_[0](lambdaValue) * value;
    result.g += response_[1](lambdaValue) * value;
    result.b += response_[2](lambdaValue) * value;
  }
  return imagingRatio_ * result;
}

base::XYZ PixelSensor::DeterministicXYZ(const base::Spectrum& radiance) const {
  base::RGB sensorRGB = DeterministicSensorRGB(radiance);
  base::XYZ xyz = ApplyMatrixToXYZ(xyzFromSensorRGB_, sensorRGB);
  return xyz / base::CIEYIntegral;
}

bool PixelSensor::IsValid() const {
  return response_[0].IsValid() && response_[1].IsValid() && response_[2].IsValid();
}

PixelSensorMode PixelSensor::Mode() const {
  return mode_;
}

Float PixelSensor::ImagingRatio() const {
  return imagingRatio_;
}

const base::ColorSpaceMatrix3x3& PixelSensor::XYZFromSensorRGB() const {
  return xyzFromSensorRGB_;
}

const std::array<SensorResponseCurve, 3>& PixelSensor::Responses() const {
  return response_;
}

const std::string& PixelSensor::Name() const {
  return name_;
}

PixelSensor::PixelSensor(
  std::array<SensorResponseCurve, 3> response,
  base::ColorSpaceMatrix3x3 xyzFromSensorRGB,
  PixelSensorMode mode,
  Float imagingRatio,
  std::string name
)
  : response_(std::move(response)),
    xyzFromSensorRGB_(xyzFromSensorRGB),
    mode_(mode),
    imagingRatio_(imagingRatio),
    name_(std::move(name)) {
  if (!IsValid()) {
    throw std::invalid_argument("PixelSensor requires three valid response curves");
  }
  if (!xyzFromSensorRGB_.IsFinite()) {
    throw std::invalid_argument("PixelSensor sensor-to-XYZ matrix must be finite");
  }
  ValidateFinite(imagingRatio_, "PixelSensor imaging ratio");
}

Float PixelSensor::ComputeImagingRatio(Float exposureTime, Float iso) {
  if (!(exposureTime >= 0) || !(iso >= 0) || !std::isfinite(exposureTime) || !std::isfinite(iso)) {
    throw std::invalid_argument("PixelSensor exposure time and ISO must be finite non-negative values");
  }
  return exposureTime * iso / static_cast<Float>(100);
}

FilmFilter FilmFilter::Box(Float radiusXValue, Float radiusYValue) {
  return FilmFilter{FilmFilterType::Box, radiusXValue, radiusYValue};
}

FilmFilter FilmFilter::Triangle(Float radiusXValue, Float radiusYValue) {
  return FilmFilter{FilmFilterType::Triangle, radiusXValue, radiusYValue};
}

Float FilmFilter::Evaluate(Float dx, Float dy) const {
  if (radiusX < 0 || radiusY < 0 || !std::isfinite(radiusX) || !std::isfinite(radiusY)) {
    throw std::invalid_argument("FilmFilter radii must be finite non-negative values");
  }
  Float ax = std::fabs(dx);
  Float ay = std::fabs(dy);
  if (ax > radiusX || ay > radiusY) {
    return 0;
  }
  switch (type) {
  case FilmFilterType::Box:
    return 1;
  case FilmFilterType::Triangle:
    if (radiusX == 0 || radiusY == 0) {
      return 0;
    }
    return std::max(static_cast<Float>(0), radiusX - ax) / radiusX *
           std::max(static_cast<Float>(0), radiusY - ay) / radiusY;
  }
  return 0;
}

Film::Film(FilmOptions options)
  : options_(std::move(options)) {
  if (options_.width <= 0 || options_.height <= 0) {
    throw std::invalid_argument("Film dimensions must be positive");
  }
  if (!options_.sensor.IsValid()) {
    throw std::invalid_argument("Film requires a valid PixelSensor");
  }
  outputRGBFromSensorRGB_ = options_.outputColorSpace.xyzToRGB * options_.sensor.XYZFromSensorRGB();
  pixels_.resize(static_cast<std::size_t>(options_.width) * static_cast<std::size_t>(options_.height));
}

base::SampledWavelengths Film::SampleWavelengths(Float u) const {
  switch (options_.wavelengthSampling) {
  case WavelengthSamplingMode::Visible:
    return base::SampledWavelengths::SampleVisible(u);
  case WavelengthSamplingMode::Uniform:
    return base::SampledWavelengths::SampleUniform(u);
  }
  return base::SampledWavelengths::SampleVisible(u);
}

void Film::AddSample(
  FilmPoint2i p,
  base::SampledSpectrum L,
  const base::SampledWavelengths& lambda,
  const SpectralVisibleSurface* visibleSurface,
  Float sampleWeight
) {
  if (!Bounds().Contains(p) || sampleWeight == 0) {
    return;
  }
  AddSensorRGB(p, ClampSensorRGB(options_.sensor.ToSensorRGB(L, lambda)), visibleSurface, sampleWeight);
}

void Film::AddFilteredSample(
  FilmPoint2f p,
  base::SampledSpectrum L,
  const base::SampledWavelengths& lambda,
  const SpectralVisibleSurface* visibleSurface,
  Float sampleWeight
) {
  FilmPoint2f pDiscrete{p.x + static_cast<Float>(0.5), p.y + static_cast<Float>(0.5)};
  FilmBounds2i sampleBounds{
    {
      static_cast<int>(std::floor(pDiscrete.x - options_.filter.radiusX)),
      static_cast<int>(std::floor(pDiscrete.y - options_.filter.radiusY))
    },
    {
      static_cast<int>(std::floor(pDiscrete.x + options_.filter.radiusX)) + 1,
      static_cast<int>(std::floor(pDiscrete.y + options_.filter.radiusY)) + 1
    }
  };
  sampleBounds = Intersect(sampleBounds, Bounds());
  base::RGB sensorRGB = ClampSensorRGB(options_.sensor.ToSensorRGB(L, lambda));
  for (int y = sampleBounds.min.y; y < sampleBounds.max.y; ++y) {
    for (int x = sampleBounds.min.x; x < sampleBounds.max.x; ++x) {
      Float filterWeight = options_.filter.Evaluate(
        p.x - static_cast<Float>(x) - static_cast<Float>(0.5),
        p.y - static_cast<Float>(y) - static_cast<Float>(0.5)
      );
      if (filterWeight != 0) {
        AddSensorRGB({x, y}, sensorRGB, visibleSurface, sampleWeight * filterWeight);
      }
    }
  }
}

FilmTile Film::CreateTile(FilmBounds2i bounds) const {
  return FilmTile(Intersect(bounds, Bounds()), &options_.sensor, options_.filter, options_.maxComponentValue);
}

void Film::MergeTile(const FilmTile& tile) {
  FilmBounds2i tileBounds = Intersect(tile.Bounds(), Bounds());
  for (int y = tileBounds.min.y; y < tileBounds.max.y; ++y) {
    for (int x = tileBounds.min.x; x < tileBounds.max.x; ++x) {
      const FilmTile::Pixel& source = tile.pixels_[tile.PixelIndex({x, y})];
      Pixel& target = pixels_[PixelIndex({x, y})];
      for (int c = 0; c < 3; ++c) {
        target.sensorSum[c] += source.sensorSum[c];
      }
      target.weightSum += source.weightSum;
      if (source.hasVisibleSurface) {
        target.visibleSurface = source.visibleSurface;
        target.hasVisibleSurface = true;
      }
    }
  }
}

base::RGB Film::ToOutputRGB(base::SampledSpectrum L, const base::SampledWavelengths& lambda) const {
  return ApplyMatrixToRGB(outputRGBFromSensorRGB_, options_.sensor.ToSensorRGB(L, lambda));
}

base::RGB Film::GetPixelSensorRGB(FilmPoint2i p) const {
  const Pixel& pixel = pixels_.at(PixelIndex(p));
  base::RGB rgb(
    static_cast<Float>(pixel.sensorSum[0]),
    static_cast<Float>(pixel.sensorSum[1]),
    static_cast<Float>(pixel.sensorSum[2])
  );
  if (pixel.weightSum != 0) {
    rgb /= static_cast<Float>(pixel.weightSum);
  }
  return rgb;
}

base::RGB Film::GetPixelLinearRGB(FilmPoint2i p) const {
  return ApplyMatrixToRGB(outputRGBFromSensorRGB_, GetPixelSensorRGB(p));
}

base::RGB Film::GetPixelEncodedRGB(FilmPoint2i p, FilmOutputEncoding encoding) const {
  base::RGB linear = GetPixelLinearRGB(p);
  switch (encoding) {
  case FilmOutputEncoding::Linear:
    return linear;
  case FilmOutputEncoding::SRGB:
    return options_.outputColorSpace.Encode(linear);
  }
  return linear;
}

base::RGB Film::GetPreviewLinearRGB(FilmPoint2i p) const {
  return GetPixelLinearRGB(p);
}

Float Film::PixelWeightSum(FilmPoint2i p) const {
  return static_cast<Float>(pixels_.at(PixelIndex(p)).weightSum);
}

const SpectralVisibleSurface* Film::GetVisibleSurface(FilmPoint2i p) const {
  const Pixel& pixel = pixels_.at(PixelIndex(p));
  return pixel.hasVisibleSurface ? &pixel.visibleSurface : nullptr;
}

int Film::Width() const {
  return options_.width;
}

int Film::Height() const {
  return options_.height;
}

FilmBounds2i Film::Bounds() const {
  return FilmBounds2i::Full(options_.width, options_.height);
}

const PixelSensor& Film::Sensor() const {
  return options_.sensor;
}

const base::RGBColorSpace& Film::OutputColorSpace() const {
  return options_.outputColorSpace;
}

const FilmFilter& Film::Filter() const {
  return options_.filter;
}

bool Film::DeterministicSingleThread() const {
  return options_.deterministicSingleThread;
}

std::size_t Film::PixelIndex(FilmPoint2i p) const {
  if (!Bounds().Contains(p)) {
    throw std::out_of_range("Film pixel is outside bounds");
  }
  return static_cast<std::size_t>(p.y) * static_cast<std::size_t>(options_.width) + static_cast<std::size_t>(p.x);
}

void Film::AddSensorRGB(
  FilmPoint2i p,
  base::RGB sensorRGB,
  const SpectralVisibleSurface* visibleSurface,
  Float weight
) {
  Pixel& pixel = pixels_[PixelIndex(p)];
  pixel.sensorSum[0] += static_cast<double>(weight) * sensorRGB.r;
  pixel.sensorSum[1] += static_cast<double>(weight) * sensorRGB.g;
  pixel.sensorSum[2] += static_cast<double>(weight) * sensorRGB.b;
  pixel.weightSum += weight;
  if (visibleSurface != nullptr && visibleSurface->valid) {
    pixel.visibleSurface = *visibleSurface;
    pixel.hasVisibleSurface = true;
  }
}

base::RGB Film::ClampSensorRGB(base::RGB sensorRGB) const {
  Float maxComponent = std::max(sensorRGB.r, std::max(sensorRGB.g, sensorRGB.b));
  if (maxComponent > options_.maxComponentValue) {
    sensorRGB *= options_.maxComponentValue / maxComponent;
  }
  return sensorRGB;
}

FilmTile::FilmTile(
  FilmBounds2i bounds,
  const PixelSensor* sensor,
  FilmFilter filter,
  Float maxComponentValue
)
  : bounds_(bounds),
    sensor_(sensor),
    filter_(filter),
    maxComponentValue_(maxComponentValue),
    pixels_(static_cast<std::size_t>(bounds.Area())) {
  if (bounds_.IsEmpty()) {
    throw std::invalid_argument("FilmTile bounds must not be empty");
  }
  if (sensor_ == nullptr || !sensor_->IsValid()) {
    throw std::invalid_argument("FilmTile requires a valid PixelSensor");
  }
}

void FilmTile::AddSample(
  FilmPoint2i p,
  base::SampledSpectrum L,
  const base::SampledWavelengths& lambda,
  const SpectralVisibleSurface* visibleSurface,
  Float sampleWeight
) {
  if (!bounds_.Contains(p) || sampleWeight == 0) {
    return;
  }
  AddSensorRGB(p, ClampSensorRGB(sensor_->ToSensorRGB(L, lambda)), visibleSurface, sampleWeight);
}

void FilmTile::AddFilteredSample(
  FilmPoint2f p,
  base::SampledSpectrum L,
  const base::SampledWavelengths& lambda,
  const SpectralVisibleSurface* visibleSurface,
  Float sampleWeight
) {
  FilmPoint2f pDiscrete{p.x + static_cast<Float>(0.5), p.y + static_cast<Float>(0.5)};
  FilmBounds2i sampleBounds{
    {
      static_cast<int>(std::floor(pDiscrete.x - filter_.radiusX)),
      static_cast<int>(std::floor(pDiscrete.y - filter_.radiusY))
    },
    {
      static_cast<int>(std::floor(pDiscrete.x + filter_.radiusX)) + 1,
      static_cast<int>(std::floor(pDiscrete.y + filter_.radiusY)) + 1
    }
  };
  sampleBounds = Intersect(sampleBounds, bounds_);
  base::RGB sensorRGB = ClampSensorRGB(sensor_->ToSensorRGB(L, lambda));
  for (int y = sampleBounds.min.y; y < sampleBounds.max.y; ++y) {
    for (int x = sampleBounds.min.x; x < sampleBounds.max.x; ++x) {
      Float filterWeight = filter_.Evaluate(
        p.x - static_cast<Float>(x) - static_cast<Float>(0.5),
        p.y - static_cast<Float>(y) - static_cast<Float>(0.5)
      );
      if (filterWeight != 0) {
        AddSensorRGB({x, y}, sensorRGB, visibleSurface, sampleWeight * filterWeight);
      }
    }
  }
}

FilmBounds2i FilmTile::Bounds() const {
  return bounds_;
}

std::size_t FilmTile::PixelIndex(FilmPoint2i p) const {
  if (!bounds_.Contains(p)) {
    throw std::out_of_range("FilmTile pixel is outside bounds");
  }
  int x = p.x - bounds_.min.x;
  int y = p.y - bounds_.min.y;
  return static_cast<std::size_t>(y) * static_cast<std::size_t>(bounds_.Width()) + static_cast<std::size_t>(x);
}

void FilmTile::AddSensorRGB(
  FilmPoint2i p,
  base::RGB sensorRGB,
  const SpectralVisibleSurface* visibleSurface,
  Float weight
) {
  Pixel& pixel = pixels_[PixelIndex(p)];
  pixel.sensorSum[0] += static_cast<double>(weight) * sensorRGB.r;
  pixel.sensorSum[1] += static_cast<double>(weight) * sensorRGB.g;
  pixel.sensorSum[2] += static_cast<double>(weight) * sensorRGB.b;
  pixel.weightSum += weight;
  if (visibleSurface != nullptr && visibleSurface->valid) {
    pixel.visibleSurface = *visibleSurface;
    pixel.hasVisibleSurface = true;
  }
}

base::RGB FilmTile::ClampSensorRGB(base::RGB sensorRGB) const {
  Float maxComponent = std::max(sensorRGB.r, std::max(sensorRGB.g, sensorRGB.b));
  if (maxComponent > maxComponentValue_) {
    sensorRGB *= maxComponentValue_ / maxComponent;
  }
  return sensorRGB;
}

} // namespace render
} // namespace rayrender
