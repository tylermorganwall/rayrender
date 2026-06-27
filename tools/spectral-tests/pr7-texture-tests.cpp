#include "src/materials/spectral_texture.h"

#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <new>
#include <string>
#include <utility>
#include <vector>

using namespace rayrender::base;
using namespace rayrender::materials;

Float random_gen::unif_rand() {
  return std::ldexp(rng(), -32);
}

random_gen::~random_gen() {}

namespace {

bool countAllocations = false;
std::size_t allocationCount = 0;

void CountAllocation() {
  if (countAllocations) {
    ++allocationCount;
  }
}

bool Approx(Float lhs, Float rhs, Float tolerance = static_cast<Float>(1e-5)) {
  return std::fabs(lhs - rhs) <= tolerance;
}

void Check(bool condition, const char* message) {
  if (!condition) {
    std::cerr << "PR7 texture test failed: " << message << std::endl;
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
    std::cerr << "PR7 texture test failed: " << message
              << " (unexpected diagnostic: " << e.what() << ")" << std::endl;
    std::exit(1);
  }
  Check(false, message);
}

TextureEvalContext EvalContext(Float u = static_cast<Float>(0.5), Float v = static_cast<Float>(0.5)) {
  TextureEvalContext ctx;
  ctx.uv = point2f(u, v);
  ctx.p = point3f(static_cast<Float>(0.25), static_cast<Float>(0.5), static_cast<Float>(0.75));
  return ctx;
}

RGBColorSpace ColorSpaceWithConstantIlluminant(Float value) {
  RGBColorSpace colorSpace = RGBColorSpace::SRGB();
  colorSpace.illuminant = std::make_shared<const DenselySampledSpectrum>(
    DenselySampledSpectrum::SampleFunction([&](Float) { return value; })
  );
  return colorSpace;
}

void CheckSpectrumConstant(const SampledSpectrum& spectrum, Float expected, const char* message) {
  for (int i = 0; i < NSpectrumSamples; ++i) {
    if (!Approx(spectrum[i], expected)) {
      std::cerr << "PR7 texture test failed: " << message
                << " component " << i << " expected " << expected
                << " got " << spectrum[i] << std::endl;
      std::exit(1);
    }
  }
}

std::shared_ptr<const ImageTextureData> OnePixelImage(const ImageCacheKey& key, Float value) {
  return ImageTextureData::FromInterleaved(
    key,
    1,
    1,
    4,
    {value, value, value, 1},
    ColorSpaceWithConstantIlluminant(1)
  );
}

void TestAnalyticTextures() {
  TextureEvalContext ctx = EvalContext();
  FloatTexture mixed = FloatTexture::Mix(
    FloatTexture::Constant(0),
    FloatTexture::Constant(2),
    static_cast<Float>(0.25)
  );
  FloatTexture scaled = FloatTexture::Scale(std::move(mixed), 2);
  Check(Approx(scaled.Evaluate(ctx), 1), "float mix/scale value");
  Check(scaled.Kind() == TextureKind::Scale, "float kind reports scale");

  TextureEvalContext oddChecker = EvalContext(static_cast<Float>(0.75), static_cast<Float>(0.25));
  FloatTexture checker = FloatTexture::Checker(
    FloatTexture::Constant(3),
    FloatTexture::Constant(7),
    2,
    2,
    2
  );
  Check(Approx(checker.Evaluate(oddChecker), 7), "float checker odd cell");

  FloatTexture noise = FloatTexture::Noise(2, 4, 0, 0, 0);
  Check(Approx(noise.Evaluate(ctx), 3), "float procedural noise endpoints");

  SampledWavelengths lambda = SampledWavelengths::SampleUniform(static_cast<Float>(0.25));
  SpectrumTexture spectrum = SpectrumTexture::Scale(
    SpectrumTexture::Mix(
      SpectrumTexture::Constant(1),
      SpectrumTexture::Constant(3),
      static_cast<Float>(0.25)
    ),
    2
  );
  CheckSpectrumConstant(spectrum.Evaluate(ctx, lambda), 3, "spectrum mix/scale value");
}

void TestImageDecodeFilterReconstruct() {
  ImageCacheKey key = ImageCacheKey::Color(
    "encoded.png",
    TextureSemanticRole::Albedo,
    TextureImageEncoding::SRGB,
    "sRGB",
    TextureWrapMode::Clamp,
    TextureFilterMode::Bilinear
  );
  std::shared_ptr<const ImageTextureData> image = ImageTextureData::FromInterleaved(
    key,
    2,
    1,
    4,
    {
      static_cast<Float>(0.25), static_cast<Float>(0.25), static_cast<Float>(0.25), 1,
      static_cast<Float>(0.75), static_cast<Float>(0.75), static_cast<Float>(0.75), 1
    }
  );

  TextureEvalContext ctx = EvalContext(static_cast<Float>(0.5), static_cast<Float>(0.5));
  Float expected = static_cast<Float>(0.5) *
                   (RGBColorEncoding::SRGBToLinear(static_cast<Float>(0.25)) +
                    RGBColorEncoding::SRGBToLinear(static_cast<Float>(0.75)));
  RGB filtered = image->SampleRGB(ctx);
  Check(image->WasDecodedToLinear(), "sRGB color image decoded once at load");
  Check(Approx(filtered.r, expected, static_cast<Float>(1e-6)), "linear red filtered after decode");
  Check(Approx(filtered.g, expected, static_cast<Float>(1e-6)), "linear green filtered after decode");
  Check(Approx(filtered.b, expected, static_cast<Float>(1e-6)), "linear blue filtered after decode");

  SampledWavelengths lambda = SampledWavelengths::SampleUniform(static_cast<Float>(0.5));
  SpectrumTexture texture = SpectrumTexture::Image(
    image,
    TextureSemanticRole::Albedo,
    ColorSpaceWithConstantIlluminant(1)
  );
  CheckSpectrumConstant(texture.Evaluate(ctx, lambda), expected, "albedo reconstructs filtered linear RGB");
}

void TestScalarImagesRemainLinear() {
  ImageCacheKey key = ImageCacheKey::Scalar(
    "roughness.png",
    TextureScalarChannel::Red,
    TextureImageEncoding::SRGB,
    TextureWrapMode::Clamp,
    TextureFilterMode::Nearest
  );
  std::shared_ptr<const ImageTextureData> image = ImageTextureData::FromInterleaved(
    key,
    1,
    1,
    4,
    {static_cast<Float>(0.5), static_cast<Float>(0.25), static_cast<Float>(0.75), 1}
  );
  TextureEvalContext ctx = EvalContext();
  Check(!image->WasDecodedToLinear(), "scalar image does not gamma decode");
  Check(Approx(image->SampleScalar(ctx), static_cast<Float>(0.5)), "scalar sample remains raw linear data");
  Check(Approx(FloatTexture::Image(image).Evaluate(ctx), static_cast<Float>(0.5)), "float image texture stays scalar");
}

void TestCacheSemanticSeparation() {
  SpectralImageCache cache;
  int loads = 0;
  auto loader = [&](const ImageCacheKey& key) {
    ++loads;
    return OnePixelImage(key, 1);
  };

  ImageCacheKey albedo = ImageCacheKey::Color(
    "shared.png",
    TextureSemanticRole::Albedo,
    TextureImageEncoding::Linear
  );
  ImageCacheKey illuminant = ImageCacheKey::Color(
    "shared.png",
    TextureSemanticRole::Illuminant,
    TextureImageEncoding::Linear
  );

  std::shared_ptr<const ImageTextureData> a1 = cache.LookupOrCreate(albedo, loader);
  std::shared_ptr<const ImageTextureData> a2 = cache.LookupOrCreate(albedo, loader);
  std::shared_ptr<const ImageTextureData> light = cache.LookupOrCreate(illuminant, loader);

  Check(a1.get() == a2.get(), "same image key reuses cache entry");
  Check(a1.get() != light.get(), "semantic role gets a distinct cache entry");
  Check(cache.Size() == 2, "cache has two semantic entries");
  Check(loads == 2, "cache loader called only for distinct keys");
  Check(albedo.ToString() != illuminant.ToString(), "cache key string includes role");
}

void TestRoleSpecificReconstruction() {
  TextureEvalContext ctx = EvalContext();
  SampledWavelengths lambda = SampledWavelengths::SampleUniform(static_cast<Float>(0.1));
  RGBColorSpace colorSpace = ColorSpaceWithConstantIlluminant(2);

  ImageCacheKey albedoKey = ImageCacheKey::Color(
    "white.png",
    TextureSemanticRole::Albedo,
    TextureImageEncoding::Linear,
    "sRGB",
    TextureWrapMode::Clamp,
    TextureFilterMode::Nearest
  );
  ImageCacheKey illuminantKey = ImageCacheKey::Color(
    "white.png",
    TextureSemanticRole::Illuminant,
    TextureImageEncoding::Linear,
    "sRGB",
    TextureWrapMode::Clamp,
    TextureFilterMode::Nearest
  );

  std::shared_ptr<const ImageTextureData> albedoImage = OnePixelImage(albedoKey, 1);
  std::shared_ptr<const ImageTextureData> illuminantImage = OnePixelImage(illuminantKey, 1);
  SpectrumTexture albedo = SpectrumTexture::Image(albedoImage, TextureSemanticRole::Albedo, colorSpace);
  SpectrumTexture illuminant = SpectrumTexture::Image(illuminantImage, TextureSemanticRole::Illuminant, colorSpace);

  CheckSpectrumConstant(albedo.Evaluate(ctx, lambda), 1, "white albedo reconstructs as bounded reflectance");
  CheckSpectrumConstant(illuminant.Evaluate(ctx, lambda), 2, "white illuminant reconstructs with illuminant");
  CheckThrowsContaining(
    [&]() { (void)SpectrumTexture::Image(albedoImage, TextureSemanticRole::Illuminant, colorSpace); },
    "role",
    "spectrum image role must match cache key"
  );
}

void TestDerivativeFootprintAndHandles() {
  TextureEvalContext ctx = EvalContext();
  ctx.dudx = static_cast<Float>(0.1);
  ctx.dvdx = static_cast<Float>(0.2);
  ctx.dudy = static_cast<Float>(0.3);
  ctx.dvdy = static_cast<Float>(0.4);
  TextureFilterFootprint footprint = TextureFilterFootprintFromContext(ctx, 10, 20);
  Check(Approx(footprint.width, std::sqrt(static_cast<Float>(73))), "derivative footprint width");

  TextureHandleTable table;
  TextureHandle f = table.Add(FloatTexture::Constant(4));
  TextureHandle s = table.Add(SpectrumTexture::Constant(5));
  SampledWavelengths lambda = SampledWavelengths::SampleUniform(static_cast<Float>(0.3));
  Check(table.ValueType(f) == TextureValueType::Float, "float handle type");
  Check(table.ValueType(s) == TextureValueType::Spectrum, "spectrum handle type");
  Check(Approx(table.GetFloat(f).Evaluate(ctx), 4), "float handle evaluate");
  CheckSpectrumConstant(table.GetSpectrum(s).Evaluate(ctx, lambda), 5, "spectrum handle evaluate");
  CheckThrowsContaining(
    [&]() { (void)table.GetFloat(s); },
    "FloatTexture",
    "typed handle rejects wrong dispatch"
  );
}

void TestNoHeapAllocationDuringEvaluation() {
  TextureEvalContext ctx = EvalContext();
  ImageCacheKey scalarKey = ImageCacheKey::Scalar(
    "scalar.png",
    TextureScalarChannel::Red,
    TextureImageEncoding::Linear,
    TextureWrapMode::Clamp,
    TextureFilterMode::Bilinear
  );
  std::shared_ptr<const ImageTextureData> scalarImage = ImageTextureData::FromInterleaved(
    scalarKey,
    2,
    1,
    4,
    {0, 0, 0, 1, 1, 1, 1, 1}
  );
  FloatTexture scalar = FloatTexture::Mix(
    FloatTexture::Image(scalarImage),
    FloatTexture::Checker(FloatTexture::Constant(0), FloatTexture::Constant(1), 2, 2, 2),
    static_cast<Float>(0.25)
  );

  ImageCacheKey albedoKey = ImageCacheKey::Color(
    "albedo.png",
    TextureSemanticRole::Albedo,
    TextureImageEncoding::Linear,
    "sRGB",
    TextureWrapMode::Clamp,
    TextureFilterMode::Bilinear
  );
  std::shared_ptr<const ImageTextureData> colorImage = ImageTextureData::FromInterleaved(
    albedoKey,
    2,
    1,
    4,
    {
      static_cast<Float>(0.25), static_cast<Float>(0.25), static_cast<Float>(0.25), 1,
      static_cast<Float>(0.75), static_cast<Float>(0.75), static_cast<Float>(0.75), 1
    }
  );
  SpectrumTexture spectrum = SpectrumTexture::Mix(
    SpectrumTexture::Image(colorImage, TextureSemanticRole::Albedo, ColorSpaceWithConstantIlluminant(1)),
    SpectrumTexture::Constant(1),
    FloatTexture::Constant(static_cast<Float>(0.5))
  );
  SampledWavelengths lambda = SampledWavelengths::SampleUniform(static_cast<Float>(0.2));

  allocationCount = 0;
  countAllocations = true;
  volatile Float scalarAccumulator = 0;
  volatile Float spectrumAccumulator = 0;
  for (int i = 0; i < 1000; ++i) {
    scalarAccumulator += scalar.Evaluate(ctx);
    SampledSpectrum sample = spectrum.Evaluate(ctx, lambda);
    spectrumAccumulator += sample[0];
  }
  countAllocations = false;
  (void)scalarAccumulator;
  (void)spectrumAccumulator;

  Check(allocationCount == 0, "texture evaluation performs no heap allocations");
}

} // namespace

void* operator new(std::size_t size) {
  CountAllocation();
  if (void* pointer = std::malloc(size)) {
    return pointer;
  }
  throw std::bad_alloc();
}

void* operator new[](std::size_t size) {
  CountAllocation();
  if (void* pointer = std::malloc(size)) {
    return pointer;
  }
  throw std::bad_alloc();
}

void operator delete(void* pointer) noexcept {
  std::free(pointer);
}

void operator delete[](void* pointer) noexcept {
  std::free(pointer);
}

void operator delete(void* pointer, std::size_t) noexcept {
  std::free(pointer);
}

void operator delete[](void* pointer, std::size_t) noexcept {
  std::free(pointer);
}

int main() {
  try {
    TestAnalyticTextures();
    TestImageDecodeFilterReconstruct();
    TestScalarImagesRemainLinear();
    TestCacheSemanticSeparation();
    TestRoleSpecificReconstruction();
    TestDerivativeFootprintAndHandles();
    TestNoHeapAllocationDuringEvaluation();
  } catch (const std::exception& e) {
    std::cerr << "PR7 texture test failed with exception: " << e.what() << std::endl;
    return 1;
  }

  std::cout << "PR7 texture tests passed" << std::endl;
  return 0;
}
