#include "src/render/spectral_film.h"

#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

using namespace rayrender::base;
using namespace rayrender::render;

namespace {

bool Approx(Float lhs, Float rhs, Float tolerance = static_cast<Float>(1e-5)) {
  return std::fabs(lhs - rhs) <= tolerance;
}

bool ApproxRel(
  Float lhs,
  Float rhs,
  Float relativeTolerance = static_cast<Float>(1e-4),
  Float absoluteTolerance = static_cast<Float>(1e-4)
) {
  Float scale = std::max(static_cast<Float>(1), std::max(std::fabs(lhs), std::fabs(rhs)));
  return std::fabs(lhs - rhs) <= std::max(absoluteTolerance, relativeTolerance * scale);
}

void Check(bool condition, const char* message) {
  if (!condition) {
    std::cerr << "PR8 Film/PixelSensor test failed: " << message << std::endl;
    std::exit(1);
  }
}

void CheckRGBApprox(const RGB& actual, const RGB& expected, Float tolerance, const char* message) {
  if (!Approx(actual.r, expected.r, tolerance) ||
      !Approx(actual.g, expected.g, tolerance) ||
      !Approx(actual.b, expected.b, tolerance)) {
    std::cerr << "PR8 Film/PixelSensor test failed: " << message
              << " expected [" << expected.r << ", " << expected.g << ", " << expected.b << "]"
              << " got [" << actual.r << ", " << actual.g << ", " << actual.b << "]"
              << std::endl;
    std::exit(1);
  }
}

void CheckRGBRel(const RGB& actual, const RGB& expected, Float tolerance, const char* message) {
  if (!ApproxRel(actual.r, expected.r, tolerance) ||
      !ApproxRel(actual.g, expected.g, tolerance) ||
      !ApproxRel(actual.b, expected.b, tolerance)) {
    std::cerr << "PR8 Film/PixelSensor test failed: " << message
              << " expected [" << expected.r << ", " << expected.g << ", " << expected.b << "]"
              << " got [" << actual.r << ", " << actual.g << ", " << actual.b << "]"
              << std::endl;
    std::exit(1);
  }
}

Spectrum MakeSmoothRadiance() {
  SpectrumDataPolicy policy;
  policy.validation = SpectrumValueValidation::NonNegative;
  return Spectrum(PiecewiseLinearSpectrum(
    {LambdaMin, 430, 510, 610, LambdaMax},
    {static_cast<Float>(0.15), static_cast<Float>(0.7), static_cast<Float>(1.1),
     static_cast<Float>(0.45), static_cast<Float>(0.2)},
    policy
  ));
}

Spectrum MakeScaledDenseSpectrum(const Spectrum& source, Float scale) {
  return Spectrum(DenselySampledSpectrum::SampleFunction([&](Float lambda) {
    return scale * source(lambda);
  }));
}

RGB EstimateSensor(
  const PixelSensor& sensor,
  const Spectrum& radiance,
  WavelengthSamplingMode mode,
  int sampleCount
) {
  RGB sum;
  for (int i = 0; i < sampleCount; ++i) {
    Float u = (static_cast<Float>(i) + static_cast<Float>(0.5)) / static_cast<Float>(sampleCount);
    SampledWavelengths lambda = mode == WavelengthSamplingMode::Visible
      ? SampledWavelengths::SampleVisible(u)
      : SampledWavelengths::SampleUniform(u);
    sum += sensor.ToSensorRGB(radiance.Sample(lambda), lambda);
  }
  return sum / static_cast<Float>(sampleCount);
}

void AddStratifiedSpectrum(Film& film, FilmPoint2i pixel, const Spectrum& radiance, int sampleCount) {
  for (int i = 0; i < sampleCount; ++i) {
    Float u = (static_cast<Float>(i) + static_cast<Float>(0.5)) / static_cast<Float>(sampleCount);
    SampledWavelengths lambda = film.SampleWavelengths(u);
    film.AddSample(pixel, radiance.Sample(lambda), lambda, nullptr, 1);
  }
}

SampledSpectrum ConstantPacket(Float value) {
  return SampledSpectrum(value);
}

FilmOptions MakeFilmOptions(
  int width,
  int height,
  const PixelSensor& sensor,
  const RGBColorSpace& outputColorSpace
) {
  FilmOptions options;
  options.width = width;
  options.height = height;
  options.sensor = sensor;
  options.outputColorSpace = outputColorSpace;
  options.filter = FilmFilter::Triangle(1, 1);
  options.wavelengthSampling = WavelengthSamplingMode::Visible;
  options.deterministicSingleThread = true;
  return options;
}

void TestSensorEstimator(const NamedSpectrumRegistry& registry, const RGBColorSpace& colorSpace) {
  PixelSensor sensor = PixelSensor::CreateCIE1931(registry, colorSpace);
  Spectrum radiance = MakeSmoothRadiance();
  RGB expected = sensor.DeterministicSensorRGB(radiance);

  RGB visible = EstimateSensor(sensor, radiance, WavelengthSamplingMode::Visible, 32768);
  RGB uniform = EstimateSensor(sensor, radiance, WavelengthSamplingMode::Uniform, 32768);

  CheckRGBRel(visible, expected, static_cast<Float>(4e-3), "visible wavelength estimator matches deterministic integral");
  CheckRGBRel(uniform, expected, static_cast<Float>(4e-3), "uniform wavelength estimator matches deterministic integral");
  CheckRGBRel(visible, uniform, static_cast<Float>(4e-3), "visible and uniform wavelength sampling agree");
}

void TestPDFDivisionExactlyOnce(const NamedSpectrumRegistry& registry, const RGBColorSpace& colorSpace) {
  PixelSensor sensor = PixelSensor::CreateCIE1931(registry, colorSpace);
  Spectrum radiance = MakeSmoothRadiance();
  SampledWavelengths lambda = SampledWavelengths::SampleUniform(static_cast<Float>(0.37));
  SampledSpectrum L = radiance.Sample(lambda);
  RGB correct = sensor.ToSensorRGB(L, lambda);
  RGB doubleDivided = sensor.ToSensorRGB(SafeDiv(L, lambda.PDF()), lambda);

  Check(doubleDivided.g > correct.g * 100, "deliberate double PDF division is detectably wrong");

  SampledWavelengths terminated = SampledWavelengths::SampleVisible(static_cast<Float>(0.23));
  Float originalPDF0 = terminated.PDF(0);
  terminated.TerminateSecondary();
  RGB terminatedRGB = sensor.ToSensorRGB(radiance.Sample(terminated), terminated);
  Float expectedY = registry.GetOrThrow("cie-y")(terminated[0]) *
                    radiance(terminated[0]) / originalPDF0;
  Check(
    ApproxRel(terminatedRGB.g, expectedY, static_cast<Float>(1e-5)),
    "terminated secondary wavelengths keep the single surviving PDF contribution"
  );
}

void TestRGBSensorMode(const NamedSpectrumRegistry& registry, const RGBColorSpace& colorSpace) {
  PixelSensor cie = PixelSensor::CreateCIE1931(registry, colorSpace);
  PixelSensor rgb = PixelSensor::CreateRGBColorSpace(registry, colorSpace, colorSpace);
  Spectrum radiance = MakeSmoothRadiance();
  SampledWavelengths lambda = SampledWavelengths::SampleUniform(static_cast<Float>(0.41));

  Film cieFilm(MakeFilmOptions(1, 1, cie, colorSpace));
  Film rgbFilm(MakeFilmOptions(1, 1, rgb, colorSpace));
  RGB cieOut = cieFilm.ToOutputRGB(radiance.Sample(lambda), lambda);
  RGB rgbOut = rgbFilm.ToOutputRGB(radiance.Sample(lambda), lambda);
  CheckRGBApprox(rgbOut, cieOut, static_cast<Float>(5e-4), "RGB color-space sensor converts back to output RGB");
  Check(rgb.Mode() == PixelSensorMode::RGBColorSpace, "RGB color-space sensor reports its mode");
}

void TestCalibration(const NamedSpectrumRegistry& registry, const RGBColorSpace& colorSpace) {
  std::array<SensorResponseCurve, 3> cie = {
    SensorResponseCurve::FromSpectrum(registry.GetOrThrow("cie-x")),
    SensorResponseCurve::FromSpectrum(registry.GetOrThrow("cie-y")),
    SensorResponseCurve::FromSpectrum(registry.GetOrThrow("cie-z"))
  };
  Spectrum equalEnergy(ConstantSpectrum(1));
  std::vector<SensorCalibrationSwatch> swatches = {
    {Spectrum(PiecewiseLinearSpectrum({LambdaMin, 450, LambdaMax}, {1, static_cast<Float>(0.2), static_cast<Float>(0.1)}))},
    {Spectrum(PiecewiseLinearSpectrum({LambdaMin, 550, LambdaMax}, {static_cast<Float>(0.1), 1, static_cast<Float>(0.15)}))},
    {Spectrum(PiecewiseLinearSpectrum({LambdaMin, 650, LambdaMax}, {static_cast<Float>(0.05), static_cast<Float>(0.25), 1}))},
    {Spectrum(ConstantSpectrum(static_cast<Float>(0.5)))}
  };

  ColorSpaceMatrix3x3 calibrated = PixelSensor::CalibrateXYZFromSensorRGB(
    cie,
    swatches,
    equalEnergy,
    equalEnergy,
    registry
  );
  Check(Approx(calibrated(0, 0), 1, static_cast<Float>(2e-4)), "calibration matrix identity xx");
  Check(Approx(calibrated(1, 1), 1, static_cast<Float>(2e-4)), "calibration matrix identity yy");
  Check(Approx(calibrated(2, 2), 1, static_cast<Float>(2e-4)), "calibration matrix identity zz");
  Check(std::fabs(calibrated(0, 1)) < static_cast<Float>(2e-4), "calibration matrix off diagonal xy");
  Check(std::fabs(calibrated(1, 2)) < static_cast<Float>(2e-4), "calibration matrix off diagonal yz");

  PixelSensor calibratedSensor = PixelSensor::CreateCalibratedRGB(
    cie[0],
    cie[1],
    cie[2],
    swatches,
    equalEnergy,
    colorSpace,
    registry
  );
  Check(calibratedSensor.Mode() == PixelSensorMode::Measured, "calibrated RGB sensor uses measured mode");
}

void TestFilmAccumulationAndTiles(const NamedSpectrumRegistry& registry, const RGBColorSpace& colorSpace) {
  PixelSensor sensor = PixelSensor::CreateCIE1931(registry, colorSpace);
  Film serial(MakeFilmOptions(3, 2, sensor, colorSpace));
  Film tiled(MakeFilmOptions(3, 2, sensor, colorSpace));

  SpectralVisibleSurface surface;
  surface.valid = true;
  surface.depth = static_cast<Float>(4.5);
  surface.primitiveId = 17;
  surface.albedo = RGB(static_cast<Float>(0.2), static_cast<Float>(0.3), static_cast<Float>(0.4));
  surface.hasAlbedo = true;

  std::vector<FilmPoint2i> pixels = {{0, 0}, {1, 0}, {2, 1}, {1, 1}};
  std::vector<Float> weights = {1, static_cast<Float>(2), static_cast<Float>(0.5), static_cast<Float>(3)};
  std::vector<Float> wavelengths = {static_cast<Float>(0.11), static_cast<Float>(0.31), static_cast<Float>(0.61), static_cast<Float>(0.89)};

  FilmTile tile = tiled.CreateTile(tiled.Bounds());
  for (std::size_t i = 0; i < pixels.size(); ++i) {
    SampledWavelengths lambda = SampledWavelengths::SampleUniform(wavelengths[i]);
    SampledSpectrum L(static_cast<Float>(1) + static_cast<Float>(i));
    const SpectralVisibleSurface* maybeSurface = i == 1 ? &surface : nullptr;
    serial.AddSample(pixels[i], L, lambda, maybeSurface, weights[i]);
    tile.AddSample(pixels[i], L, lambda, maybeSurface, weights[i]);
  }
  tiled.MergeTile(tile);

  for (int y = 0; y < 2; ++y) {
    for (int x = 0; x < 3; ++x) {
      FilmPoint2i p{x, y};
      CheckRGBApprox(tiled.GetPixelLinearRGB(p), serial.GetPixelLinearRGB(p), static_cast<Float>(1e-5), "tile merge equals serial film");
      Check(Approx(tiled.PixelWeightSum(p), serial.PixelWeightSum(p), static_cast<Float>(1e-6)), "tile merge preserves weights");
    }
  }

  const SpectralVisibleSurface* stored = tiled.GetVisibleSurface({1, 0});
  Check(stored != nullptr && stored->primitiveId == 17, "visible-surface hook stores metadata");
  CheckRGBApprox(
    tiled.GetPreviewLinearRGB({1, 0}),
    tiled.GetPixelLinearRGB({1, 0}),
    static_cast<Float>(1e-7),
    "preview reads accumulated Film channels"
  );
}

void TestFilteringAndWeights(const NamedSpectrumRegistry& registry, const RGBColorSpace& colorSpace) {
  PixelSensor sensor = PixelSensor::CreateCIE1931(registry, colorSpace);
  FilmOptions options = MakeFilmOptions(2, 2, sensor, colorSpace);
  options.filter = FilmFilter::Triangle(1, 1);
  Film film(options);

  SampledWavelengths lambda = SampledWavelengths::SampleUniform(static_cast<Float>(0.5));
  film.AddFilteredSample({static_cast<Float>(0.5), static_cast<Float>(0.5)}, ConstantPacket(1), lambda, nullptr, 2);

  Check(Approx(film.PixelWeightSum({0, 0}), 2), "center filtered sample stores expected weight");
  Check(Approx(film.PixelWeightSum({1, 0}), 0), "triangle edge sample has zero neighboring weight at radius");

  RGB first = film.GetPixelLinearRGB({0, 0});
  film.AddSample({0, 0}, ConstantPacket(1), lambda, nullptr, 6);
  RGB second = film.GetPixelLinearRGB({0, 0});
  CheckRGBApprox(second, first, static_cast<Float>(1e-5), "sample weighting normalizes repeated equal samples");
}

void TestWhitePointAndOutputEncoding(const NamedSpectrumRegistry& registry, const RGBColorSpace& colorSpace) {
  PixelSensor sensor = PixelSensor::CreateCIE1931(registry, colorSpace);
  Film film(MakeFilmOptions(1, 1, sensor, colorSpace));
  Spectrum d65 = MakeScaledDenseSpectrum(registry.GetOrThrow("stdillum-D65"), static_cast<Float>(0.01));
  AddStratifiedSpectrum(film, {0, 0}, d65, 32768);
  RGB linear = film.GetPixelLinearRGB({0, 0});
  Check(
    ApproxRel(linear.r, linear.g, static_cast<Float>(4e-3)) &&
      ApproxRel(linear.g, linear.b, static_cast<Float>(4e-3)),
    "D65 white point is neutral in sRGB output"
  );

  RGB encoded = film.GetPixelEncodedRGB({0, 0}, FilmOutputEncoding::SRGB);
  RGB expectedEncoded = colorSpace.Encode(linear);
  CheckRGBApprox(encoded, expectedEncoded, static_cast<Float>(1e-6), "output encoding is applied after Film accumulation");

  Spectrum equalEnergy = MakeScaledDenseSpectrum(Spectrum(ConstantSpectrum(1)), static_cast<Float>(0.01));
  PixelSensor balanced = PixelSensor::CreateCIE1931(registry, colorSpace, &equalEnergy);
  Film balancedFilm(MakeFilmOptions(1, 1, balanced, colorSpace));
  AddStratifiedSpectrum(balancedFilm, {0, 0}, equalEnergy, 32768);
  RGB balancedRGB = balancedFilm.GetPixelLinearRGB({0, 0});
  Check(
    ApproxRel(balancedRGB.r, balancedRGB.g, static_cast<Float>(4e-3)) &&
      ApproxRel(balancedRGB.g, balancedRGB.b, static_cast<Float>(4e-3)),
    "Bradford white balance maps source white to output white"
  );
}

void TestNegativeOutputIsNotClipped(const NamedSpectrumRegistry& registry, const RGBColorSpace& colorSpace) {
  PixelSensor rgb = PixelSensor::CreateRGBColorSpace(registry, colorSpace, colorSpace);
  Film film(MakeFilmOptions(1, 1, rgb, colorSpace));
  Float blueU = (static_cast<Float>(430) - LambdaMin) / (LambdaMax - LambdaMin);
  SampledWavelengths lambda = SampledWavelengths::SampleUniform(blueU);
  Spectrum narrowBlue(PiecewiseLinearSpectrum(
    {LambdaMin, 430, 440, LambdaMax},
    {0, 1, 0, 0}
  ));
  film.AddSample({0, 0}, narrowBlue.Sample(lambda), lambda, nullptr, 1);
  RGB linear = film.GetPixelLinearRGB({0, 0});
  Check(
    linear.r < 0 || linear.g < 0 || linear.b < 0,
    "negative sensor/output transform components are preserved before documented output policy"
  );
}

void TestPreviewDoesNotUsePacketComponents(const NamedSpectrumRegistry& registry, const RGBColorSpace& colorSpace) {
  PixelSensor sensor = PixelSensor::CreateCIE1931(registry, colorSpace);
  Film film(MakeFilmOptions(1, 1, sensor, colorSpace));
  SampledWavelengths lambda = SampledWavelengths::SampleUniform(static_cast<Float>(0.17));
  SampledSpectrum packet{1, 2, 3, 4};
  film.AddSample({0, 0}, packet, lambda, nullptr, 1);
  RGB preview = film.GetPreviewLinearRGB({0, 0});
  CheckRGBApprox(preview, film.GetPixelLinearRGB({0, 0}), static_cast<Float>(1e-7), "preview comes from Film output channels");
  Check(
    !Approx(preview.r, packet[0], static_cast<Float>(1e-3)) ||
      !Approx(preview.g, packet[1], static_cast<Float>(1e-3)) ||
      !Approx(preview.b, packet[2], static_cast<Float>(1e-3)),
    "preview does not reinterpret wavelength packet lanes as RGB"
  );
}

} // namespace

int main(int argc, char** argv) {
  try {
    std::string assetDirectory = argc > 1 ? argv[1] : FindSpectralAssetDirectory();
    NamedSpectrumRegistry registry = NamedSpectrumRegistry::LoadFromDirectory(assetDirectory);
    RGBColorSpace colorSpace = LoadSRGBColorSpace(assetDirectory);

    TestSensorEstimator(registry, colorSpace);
    TestPDFDivisionExactlyOnce(registry, colorSpace);
    TestRGBSensorMode(registry, colorSpace);
    TestCalibration(registry, colorSpace);
    TestFilmAccumulationAndTiles(registry, colorSpace);
    TestFilteringAndWeights(registry, colorSpace);
    TestWhitePointAndOutputEncoding(registry, colorSpace);
    TestNegativeOutputIsNotClipped(registry, colorSpace);
    TestPreviewDoesNotUsePacketComponents(registry, colorSpace);
  } catch (const std::exception& e) {
    std::cerr << "PR8 Film/PixelSensor test failed with exception: " << e.what() << std::endl;
    return 1;
  }

  std::cout << "PR8 Film/PixelSensor tests passed" << std::endl;
  return 0;
}
