#include "spectral_texture.h"

#ifdef NOT_CRAN
#include "../hitables/hitable.h"

#include <cmath>
#include <testthat.h>

namespace {

using namespace rayrender::base;
using namespace rayrender::materials;

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

void ExpectSpectrumConstant(const SampledSpectrum& spectrum, Float value) {
  for (int i = 0; i < NSpectrumSamples; ++i) {
    expect_true(spectrum[i] == Approx(value).epsilon(1e-5));
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

} // namespace

context("PR7 spectral texture evaluation") {
  test_that("FloatTexture analytic nodes dispatch expected values") {
    TextureEvalContext ctx = EvalContext();

    FloatTexture mixed = FloatTexture::Mix(
      FloatTexture::Constant(0),
      FloatTexture::Constant(2),
      static_cast<Float>(0.25)
    );
    FloatTexture scaled = FloatTexture::Scale(std::move(mixed), 2);
    expect_true(scaled.Evaluate(ctx) == Approx(1));
    expect_true(scaled.Kind() == TextureKind::Scale);

    TextureEvalContext oddChecker = EvalContext(static_cast<Float>(0.75), static_cast<Float>(0.25));
    FloatTexture checker = FloatTexture::Checker(
      FloatTexture::Constant(3),
      FloatTexture::Constant(7),
      2,
      2,
      2
    );
    expect_true(checker.Evaluate(oddChecker) == Approx(7));

    FloatTexture noise = FloatTexture::Noise(2, 4, 0, 0, 0);
    expect_true(noise.Evaluate(ctx) == Approx(3));
  }

  test_that("SpectrumTexture analytic nodes dispatch expected sampled spectra") {
    TextureEvalContext ctx = EvalContext();
    SampledWavelengths lambda = SampledWavelengths::SampleUniform(static_cast<Float>(0.25));

    SpectrumTexture mixed = SpectrumTexture::Mix(
      SpectrumTexture::Constant(1),
      SpectrumTexture::Constant(3),
      static_cast<Float>(0.25)
    );
    SpectrumTexture scaled = SpectrumTexture::Scale(std::move(mixed), 2);
    ExpectSpectrumConstant(scaled.Evaluate(ctx, lambda), 3);

    TextureEvalContext evenChecker = EvalContext(static_cast<Float>(0.25), static_cast<Float>(0.25));
    SpectrumTexture checker = SpectrumTexture::Checker(
      SpectrumTexture::Constant(5),
      SpectrumTexture::Constant(9),
      2,
      2,
      2
    );
    ExpectSpectrumConstant(checker.Evaluate(evenChecker, lambda), 5);
  }

  test_that("sRGB color images decode before linear filtering and reconstruction") {
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
    expect_true(image->WasDecodedToLinear());
    expect_true(filtered.r == Approx(expected).epsilon(1e-6));
    expect_true(filtered.g == Approx(expected).epsilon(1e-6));
    expect_true(filtered.b == Approx(expected).epsilon(1e-6));

    SampledWavelengths lambda = SampledWavelengths::SampleUniform(static_cast<Float>(0.5));
    SpectrumTexture texture = SpectrumTexture::Image(
      image,
      TextureSemanticRole::Albedo,
      ColorSpaceWithConstantIlluminant(1)
    );
    ExpectSpectrumConstant(texture.Evaluate(ctx, lambda), expected);
  }

  test_that("scalar images remain linear even when the key asks for sRGB encoding") {
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
    expect_true(!image->WasDecodedToLinear());
    expect_true(image->SampleScalar(ctx) == Approx(0.5));
    expect_true(FloatTexture::Image(image).Evaluate(ctx) == Approx(0.5));
  }

  test_that("image cache keys separate semantic interpretations") {
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

    expect_true(a1.get() == a2.get());
    expect_true(a1.get() != light.get());
    expect_true(cache.Size() == 2);
    expect_true(loads == 2);
    expect_true(albedo.ToString() != illuminant.ToString());
  }

  test_that("the consuming texture role controls RGB spectral reconstruction") {
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

    ExpectSpectrumConstant(albedo.Evaluate(ctx, lambda), 1);
    ExpectSpectrumConstant(illuminant.Evaluate(ctx, lambda), 2);
    expect_error(
      SpectrumTexture::Image(albedoImage, TextureSemanticRole::Illuminant, colorSpace),
      "role"
    );
  }

  test_that("hit records adapt to texture contexts with derivative footprints") {
    hit_record hit;
    hit.p = point3f(1, 2, 3);
    hit.u = static_cast<Float>(0.2);
    hit.v = static_cast<Float>(0.8);
    hit.dpdu = vec3f(1, 0, 0);
    hit.dpdv = vec3f(0, 1, 0);

    TextureRayDifferentials differentials;
    differentials.dudx = static_cast<Float>(0.1);
    differentials.dvdx = static_cast<Float>(0.2);
    differentials.dudy = static_cast<Float>(0.3);
    differentials.dvdy = static_cast<Float>(0.4);

    TextureEvalContext ctx = TextureEvalContextFromHitRecord(hit, differentials);
    TextureFilterFootprint footprint = TextureFilterFootprintFromContext(ctx, 10, 20);

    expect_true(ctx.p.xyz.x == Approx(1));
    expect_true(ctx.uv[0] == Approx(0.2));
    expect_true(ctx.uv[1] == Approx(0.8));
    expect_true(ctx.dpdu.xyz.x == Approx(1));
    expect_true(footprint.width == Approx(std::sqrt(static_cast<Float>(73))).epsilon(1e-6));
  }

  test_that("texture handles retain typed dispatch") {
    TextureHandleTable table;
    base::TextureHandle f = table.Add(FloatTexture::Constant(4));
    base::TextureHandle s = table.Add(SpectrumTexture::Constant(5));

    TextureEvalContext ctx = EvalContext();
    SampledWavelengths lambda = SampledWavelengths::SampleUniform(static_cast<Float>(0.3));
    expect_true(table.ValueType(f) == TextureValueType::Float);
    expect_true(table.ValueType(s) == TextureValueType::Spectrum);
    expect_true(table.GetFloat(f).Evaluate(ctx) == Approx(4));
    ExpectSpectrumConstant(table.GetSpectrum(s).Evaluate(ctx, lambda), 5);
    expect_error(table.GetFloat(s), "FloatTexture");
  }
}
#endif
