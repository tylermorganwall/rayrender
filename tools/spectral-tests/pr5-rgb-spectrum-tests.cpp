#include "src/base/base.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace rayrender::base;

namespace {

bool Approx(Float lhs, Float rhs, Float tolerance = static_cast<Float>(1e-5)) {
  return std::fabs(lhs - rhs) <= tolerance;
}

void Check(bool condition, const char* message) {
  if (!condition) {
    std::cerr << "PR5 RGB spectrum test failed: " << message << std::endl;
    std::exit(1);
  }
}

template <typename F>
void CheckThrowsContaining(F&& f, const std::string& fragment, const char* message) {
  try {
    f();
  } catch (const std::exception& e) {
    if (std::string(e.what()).find(fragment) != std::string::npos) {
      return;
    }
    std::cerr << "PR5 RGB spectrum test failed: " << message
              << " (unexpected diagnostic: " << e.what() << ")" << std::endl;
    std::exit(1);
  }
  Check(false, message);
}

std::vector<unsigned char> ReadBinaryFile(const std::string& path) {
  std::ifstream input(path.c_str(), std::ios::binary);
  if (!input) {
    throw std::runtime_error("Unable to open test binary file: " + path);
  }
  return std::vector<unsigned char>(
    std::istreambuf_iterator<char>(input),
    std::istreambuf_iterator<char>()
  );
}

void WriteBinaryFile(const std::string& path, const std::vector<unsigned char>& bytes) {
  std::ofstream output(path.c_str(), std::ios::binary);
  if (!output) {
    throw std::runtime_error("Unable to write test binary file: " + path);
  }
  output.write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
}

std::string TempPath(const std::string& label) {
  static int counter = 0;
  const char* tmp = std::getenv("TMPDIR");
  std::ostringstream path;
  path << (tmp && tmp[0] != '\0' ? tmp : "/tmp")
       << "/rayrender-pr5-" << label << "-" << static_cast<unsigned long>(std::rand())
       << "-" << counter++ << ".bin";
  return path.str();
}

XYZ ReflectanceToXYZ(
  const RGBAlbedoSpectrum& spectrum,
  const DenselySampledSpectrum& illuminant,
  const Spectrum& x,
  const Spectrum& y,
  const Spectrum& z
) {
  Float xIntegral = 0;
  Float yIntegral = 0;
  Float zIntegral = 0;
  for (int lambda = static_cast<int>(LambdaMin); lambda <= static_cast<int>(LambdaMax); ++lambda) {
    Float lambdaValue = static_cast<Float>(lambda);
    Float value = spectrum(lambdaValue) * illuminant(lambdaValue);
    xIntegral += value * x(lambdaValue);
    yIntegral += value * y(lambdaValue);
    zIntegral += value * z(lambdaValue);
  }
  return XYZ(xIntegral / CIEYIntegral, yIntegral / CIEYIntegral, zIntegral / CIEYIntegral);
}

void CheckRGBApprox(const RGB& actual, const RGB& expected, Float tolerance, const char* message) {
  if (!Approx(actual.r, expected.r, tolerance) ||
      !Approx(actual.g, expected.g, tolerance) ||
      !Approx(actual.b, expected.b, tolerance)) {
    std::cerr << "PR5 RGB spectrum test failed: " << message
              << " expected [" << expected.r << ", " << expected.g << ", " << expected.b << "]"
              << " got [" << actual.r << ", " << actual.g << ", " << actual.b << "]"
              << std::endl;
    std::exit(1);
  }
}

void TestTransferFunctions() {
  Check(Approx(RGBColorEncoding::LinearToSRGB(0), 0), "linear transfer zero");
  Check(
    Approx(RGBColorEncoding::LinearToSRGB(0.0031308f), 0.040449936f, static_cast<Float>(1e-8)),
    "pbrt linear to sRGB toe"
  );
  Check(
    Approx(RGBColorEncoding::LinearToSRGB(0.18f), 0.46135613f, static_cast<Float>(1e-6)),
    "pbrt linear to sRGB mid gray"
  );
  Check(
    Approx(RGBColorEncoding::SRGBToLinear(0.25f), 0.05087606f, static_cast<Float>(1e-6)),
    "pbrt sRGB to linear quarter"
  );
  Check(
    Approx(RGBColorEncoding::SRGBToLinear(0.5f), 0.21404105f, static_cast<Float>(1e-6)),
    "pbrt sRGB to linear half"
  );

  RGBColorSpace srgb = RGBColorSpace::SRGB();
  RGB encoded(0.25f, 0.5f, 0.75f);
  CheckRGBApprox(srgb.Encode(srgb.Decode(encoded)), encoded, static_cast<Float>(2e-5), "sRGB transfer round trip");
}

void TestMatrixConstruction() {
  RGBColorSpace srgb = RGBColorSpace::SRGB();
  Check(Approx(srgb.r.x, 0.64f), "sRGB red chromaticity x");
  Check(Approx(srgb.g.y, 0.60f), "sRGB green chromaticity y");
  Check(Approx(srgb.b.x, 0.15f), "sRGB blue chromaticity x");
  Check(Approx(srgb.w.x, 0.3127f, static_cast<Float>(1e-4)), "sRGB D65 white x");

  Check(Approx(srgb.rgbToXYZ(0, 0), 0.4123908f, static_cast<Float>(5e-6)), "sRGB matrix row 0 col 0");
  Check(Approx(srgb.rgbToXYZ(1, 1), 0.7151687f, static_cast<Float>(5e-6)), "sRGB matrix row 1 col 1");
  Check(Approx(srgb.rgbToXYZ(2, 2), 0.9505321f, static_cast<Float>(5e-6)), "sRGB matrix row 2 col 2");

  RGB roundTrip = srgb.ToLinearRGB(srgb.ToXYZ(RGB(0.2f, 0.4f, 0.8f)));
  CheckRGBApprox(roundTrip, RGB(0.2f, 0.4f, 0.8f), static_cast<Float>(2e-5), "sRGB matrix round trip");
}

void TestTableLayoutAndInterpolation(const std::string& assetDirectory) {
  std::string tablePath = detail::JoinPath(assetDirectory, "rgb-to-spectrum-srgb-v1.bin");
  RGBToSpectrumTable table = RGBToSpectrumTable::LoadFromFile(tablePath, "sRGB");

  Check(table.ColorSpaceId() == "sRGB", "table color-space id");
  Check(Approx(table.ScaleNode(0), 0), "first table scale node");
  Check(Approx(table.ScaleNode(1), 1.67704457e-06f, static_cast<Float>(1e-12)), "second table scale node");
  Check(Approx(table.ScaleNode(31), 0.482147723f, static_cast<Float>(1e-7)), "middle table scale node");
  Check(Approx(table.ScaleNode(63), 1), "last table scale node");

  Check(
    Approx(table.Coefficient(0, 0, 0, 0, 0), -0.00058463082f, static_cast<Float>(1e-10)),
    "pbrt coefficient [0][0][0][0][0]"
  );
  Check(
    Approx(table.Coefficient(0, 0, 0, 0, 1), 0.846466124f, static_cast<Float>(1e-7)),
    "pbrt coefficient [0][0][0][0][1]"
  );
  Check(
    Approx(table.Coefficient(0, 0, 0, 0, 2), -310.922729f, static_cast<Float>(1e-4)),
    "pbrt coefficient [0][0][0][0][2]"
  );
  Check(
    Approx(table.Coefficient(1, 7, 11, 13, 0), -0.001190462f, static_cast<Float>(1e-9)),
    "pbrt interior coefficient c0"
  );
  Check(
    Approx(table.Coefficient(2, 63, 63, 63, 2), 309.9353f, static_cast<Float>(1e-4)),
    "pbrt final coefficient"
  );

  RGBSigmoidPolynomial coeffs = table(RGB(0.25f, 0.5f, 0.75f));
  Check(Approx(coeffs.C0(), 1.168388e-06f, static_cast<Float>(1e-10)), "pbrt interpolation c0");
  Check(Approx(coeffs.C1(), -0.006714183f, static_cast<Float>(1e-8)), "pbrt interpolation c1");
  Check(Approx(coeffs.C2(), 3.298782f, static_cast<Float>(1e-5)), "pbrt interpolation c2");

  RGBSigmoidPolynomial gray = table(RGB(0.5f, 0.5f, 0.5f));
  Check(Approx(gray(400), 0.5f), "neutral gray spectrum is constant");
  RGBSigmoidPolynomial black = table(RGB(0, 0, 0));
  RGBSigmoidPolynomial white = table(RGB(1, 1, 1));
  Check(black(500) == 0, "black neutral spectrum is zero");
  Check(white(500) == 1, "white neutral spectrum is one");
  Check(std::isinf(black.C2()) && black.C2() < 0, "black neutral coefficient is negative infinity");
  Check(std::isinf(white.C2()) && white.C2() > 0, "white neutral coefficient is positive infinity");
}

void TestColorSpaceAndWrappers(const std::string& assetDirectory) {
  RGBColorSpace colorSpace = LoadSRGBColorSpace(assetDirectory);
  Check(colorSpace.HasSpectralReconstruction(), "loaded sRGB color space has spectral reconstruction");
  Check(colorSpace.illuminant != nullptr, "loaded sRGB color space has illuminant");
  Check(colorSpace.rgbToSpectrumTable != nullptr, "loaded sRGB color space has table");
  Check(Approx(colorSpace.w.x, 0.3127f, static_cast<Float>(8e-4)), "spectral D65 white x");
  Check(Approx(colorSpace.w.y, 0.3290f, static_cast<Float>(8e-4)), "spectral D65 white y");

  NamedSpectrumRegistry registry = NamedSpectrumRegistry::LoadFromDirectory(assetDirectory);
  const Spectrum& x = registry.GetOrThrow("cie-x");
  const Spectrum& y = registry.GetOrThrow("cie-y");
  const Spectrum& z = registry.GetOrThrow("cie-z");

  RGBAlbedoSpectrum albedoWhite(colorSpace, RGB(1, 1, 1));
  Check(albedoWhite.MaxValue() <= 1, "white albedo is bounded");
  for (int lambda = 360; lambda <= 830; lambda += 17) {
    Float value = albedoWhite(static_cast<Float>(lambda));
    Check(value >= 0 && value <= 1, "albedo spectrum stays in [0, 1]");
  }

  RGBIlluminantSpectrum whiteIlluminant(colorSpace, RGB(1, 1, 1));
  Check(
    Approx(whiteIlluminant(560), (*colorSpace.illuminant)(560), static_cast<Float>(1e-6)),
    "white RGB illuminant equals standard illuminant"
  );
  Check(
    !Approx(whiteIlluminant(560), 1, static_cast<Float>(1e-3)),
    "white RGB illuminant is not equal-energy"
  );

  RGBUnboundedSpectrum unbounded(colorSpace, RGB(2, 1, 0.5f));
  Check(Approx(unbounded.Scale(), 4), "unbounded spectrum uses pbrt scale");
  Check(unbounded(500) >= 0, "unbounded spectrum evaluates non-negative selected wavelength");

  const RGB reflectanceSamples[] = {
    RGB(0.05f, 0.1f, 0.2f),
    RGB(0.25f, 0.5f, 0.75f),
    RGB(0.9f, 0.2f, 0.1f)
  };
  for (const RGB& sample : reflectanceSamples) {
    RGBAlbedoSpectrum spectrum(colorSpace, sample);
    XYZ xyz = ReflectanceToXYZ(spectrum, *colorSpace.illuminant, x, y, z);
    RGB reconstructed = colorSpace.ToLinearRGB(xyz);
    CheckRGBApprox(reconstructed, sample, static_cast<Float>(0.025), "albedo RGB reconstruction");
  }

  XYZ whiteIlluminantXYZ = SpectrumToXYZ(whiteIlluminant, x, y, z);
  RGB whiteIlluminantRGB = colorSpace.ToLinearRGB(whiteIlluminantXYZ);
  CheckRGBApprox(whiteIlluminantRGB, RGB(1, 1, 1), static_cast<Float>(0.003), "white illuminant RGB reconstruction");

  CheckThrowsContaining(
    [&]() { (void)RGBAlbedoSpectrum(colorSpace, RGB(1.2f, 0, 0)); },
    "[0, 1]",
    "albedo rejects out-of-range RGB"
  );
}

void TestBinaryDiagnosticsAndCache(const std::string& assetDirectory) {
  std::string sourcePath = detail::JoinPath(assetDirectory, "rgb-to-spectrum-srgb-v1.bin");
  std::vector<unsigned char> valid = ReadBinaryFile(sourcePath);

  CheckThrowsContaining(
    [&]() { (void)RGBToSpectrumTable::LoadFromFile(detail::JoinPath(assetDirectory, "missing-rgb-table.bin")); },
    "Unable to open",
    "missing RGB table fails clearly"
  );

  std::vector<unsigned char> wrongVersion = valid;
  wrongVersion[16] = 2;
  std::string wrongVersionPath = TempPath("wrong-version");
  WriteBinaryFile(wrongVersionPath, wrongVersion);
  CheckThrowsContaining(
    [&]() { (void)RGBToSpectrumTable::LoadFromFile(wrongVersionPath); },
    "version",
    "wrong-version RGB table fails clearly"
  );
  std::remove(wrongVersionPath.c_str());

  std::vector<unsigned char> wrongEndian = valid;
  wrongEndian[20] = 0;
  std::string wrongEndianPath = TempPath("wrong-endian");
  WriteBinaryFile(wrongEndianPath, wrongEndian);
  CheckThrowsContaining(
    [&]() { (void)RGBToSpectrumTable::LoadFromFile(wrongEndianPath); },
    "byte-order",
    "wrong-byte-order RGB table fails clearly"
  );
  std::remove(wrongEndianPath.c_str());

  std::vector<unsigned char> corruptPayload = valid;
  corruptPayload.back() ^= 1u;
  std::string corruptPath = TempPath("checksum");
  WriteBinaryFile(corruptPath, corruptPayload);
  CheckThrowsContaining(
    [&]() { (void)RGBToSpectrumTable::LoadFromFile(corruptPath); },
    "checksum",
    "corrupt RGB table fails clearly"
  );
  std::remove(corruptPath.c_str());

  std::vector<unsigned char> truncated(valid.begin(), valid.begin() + 40);
  std::string truncatedPath = TempPath("truncated");
  WriteBinaryFile(truncatedPath, truncated);
  CheckThrowsContaining(
    [&]() { (void)RGBToSpectrumTable::LoadFromFile(truncatedPath); },
    "truncated",
    "truncated RGB table fails clearly"
  );
  std::remove(truncatedPath.c_str());

  std::string cachePath = TempPath("cache");
  WriteBinaryFile(cachePath, valid);
  std::shared_ptr<const RGBToSpectrumTable> first = LoadRGBToSpectrumTable(cachePath);
  std::vector<unsigned char> mutated = valid;
  mutated.back() ^= 1u;
  WriteBinaryFile(cachePath, mutated);
  std::shared_ptr<const RGBToSpectrumTable> second = LoadRGBToSpectrumTable(cachePath);
  Check(first.get() == second.get(), "RGB table cache loads once per process");
  std::remove(cachePath.c_str());
}

} // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "usage: pr5-rgb-spectrum-tests <spectral-asset-directory>" << std::endl;
    return 1;
  }

  std::string assetDirectory = argv[1];
  TestTransferFunctions();
  TestMatrixConstruction();
  TestTableLayoutAndInterpolation(assetDirectory);
  TestColorSpaceAndWrappers(assetDirectory);
  TestBinaryDiagnosticsAndCache(assetDirectory);

  std::cout << "PR5 RGB spectrum tests passed" << std::endl;
  return 0;
}
