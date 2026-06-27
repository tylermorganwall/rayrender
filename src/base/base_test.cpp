#include "base.h"

#include "../math/vectypes.h"

#include <type_traits>

static_assert(!std::is_convertible<Float, rayrender::base::RGB>::value,
              "RGB scalar construction must remain explicit");
static_assert(!std::is_convertible<Float, rayrender::base::XYZ>::value,
              "XYZ scalar construction must remain explicit");
static_assert(!std::is_constructible<rayrender::base::RGB, point3f>::value,
              "RGB must not be constructible from geometry point types");
static_assert(!std::is_constructible<rayrender::base::XYZ, vec3f>::value,
              "XYZ must not be constructible from geometry vector types");
static_assert(!std::is_same<rayrender::base::SpectrumHandle, rayrender::base::TextureHandle>::value,
              "Tagged handle aliases must remain distinct types");
static_assert(rayrender::base::NSpectrumSamples == 4, "Spectral packets must have four samples");
static_assert(sizeof(rayrender::base::SampledSpectrum) == sizeof(Float) * rayrender::base::NSpectrumSamples,
              "SampledSpectrum must keep fixed inline storage");
static_assert(sizeof(rayrender::base::SampledWavelengths) ==
                sizeof(Float) * rayrender::base::NSpectrumSamples * 2,
              "SampledWavelengths must keep fixed inline storage");

#ifdef NOT_CRAN
#include <testthat.h>

namespace {

using namespace rayrender::base;

template <typename Handle>
void ExpectHandleContract(DispatchTag expectedTag) {
  Handle invalid;
  expect_true(!invalid.IsValid());

  Handle handle = Handle::FromIndex(7, 3);
  expect_true(handle.IsValid());
  expect_true(handle.Index() == 7);
  expect_true(handle.Generation() == 3);
  expect_true(Handle::Tag == expectedTag);
}

struct alignas(32) TrackedScratchObject {
  explicit TrackedScratchObject(int* counter) : counter(counter) {}
  ~TrackedScratchObject() noexcept { ++(*counter); }
  int* counter = nullptr;
};

} // namespace

context("PR2 base spectral foundation") {
  test_that("RGB, XYZ, and color-space math are explicit and finite") {
    RGB rgb(0.25f, 0.5f, 0.75f);
    RGB doubled = rgb * 2.0f;
    expect_true(doubled.r == Approx(0.5f));
    expect_true(doubled.g == Approx(1.0f));
    expect_true(doubled.b == Approx(1.5f));

    RGB safe = SafeDiv(RGB(1, 2, 3), RGB(1, 0, 3), -1);
    expect_true(safe.r == Approx(1));
    expect_true(safe.g == Approx(-1));
    expect_true(safe.b == Approx(1));

    RGBColorSpace srgb = RGBColorSpace::SRGB();
    RGB encoded(0.25f, 0.5f, 0.75f);
    RGB linear = srgb.Decode(encoded);
    RGB encodedRoundTrip = srgb.Encode(linear);
    expect_true(encodedRoundTrip.r == Approx(encoded.r).epsilon(1e-5));
    expect_true(encodedRoundTrip.g == Approx(encoded.g).epsilon(1e-5));
    expect_true(encodedRoundTrip.b == Approx(encoded.b).epsilon(1e-5));

    RGB linearRoundTrip = srgb.ToLinearRGB(srgb.ToXYZ(linear));
    expect_true(linearRoundTrip.r == Approx(linear.r).epsilon(1e-5));
    expect_true(linearRoundTrip.g == Approx(linear.g).epsilon(1e-5));
    expect_true(linearRoundTrip.b == Approx(linear.b).epsilon(1e-5));
  }

  test_that("Tagged handles cover every registered dispatch tag") {
    ExpectHandleContract<SpectrumHandle>(DispatchTag::Spectrum);
    ExpectHandleContract<TextureHandle>(DispatchTag::Texture);
    ExpectHandleContract<MaterialHandle>(DispatchTag::Material);
    ExpectHandleContract<BxDFHandle>(DispatchTag::BxDF);
    ExpectHandleContract<LightHandle>(DispatchTag::Light);
    ExpectHandleContract<MediumHandle>(DispatchTag::Medium);
    ExpectHandleContract<PhaseFunctionHandle>(DispatchTag::PhaseFunction);
  }

  test_that("ScratchBuffer allocation is aligned and resettable") {
    int destructorCount = 0;
    ScratchBuffer scratch(256);
    auto* tracked = scratch.Create<TrackedScratchObject>(&destructorCount);

    expect_true(reinterpret_cast<std::uintptr_t>(tracked) % alignof(TrackedScratchObject) == 0);
    expect_true(scratch.Used() >= sizeof(TrackedScratchObject));
    expect_true(scratch.HighWater() == scratch.Used());
    expect_true(scratch.TrackedObjectCount() == 1);

    scratch.Reset();
    expect_true(destructorCount == 1);
    expect_true(scratch.Used() == 0);
    expect_true(scratch.HighWater() >= sizeof(TrackedScratchObject));
    expect_true(scratch.TrackedObjectCount() == 0);
  }

  test_that("BxDF flags and optional samples use pbrt-style conventions") {
    BxDFFlags flags = BxDFFlags::Reflection | BxDFFlags::Diffuse;
    expect_true(HasFlag(flags, BxDFFlags::Reflection));
    expect_true(HasAny(flags, BxDFFlags::Reflection | BxDFFlags::Transmission));
    expect_true(!HasFlag(flags, BxDFFlags::Specular));

    RGBBSDFSampleResult none;
    expect_true(!none.has_value());

    RGBBSDFSample sample;
    sample.f = RGB(1, 1, 1);
    sample.wi = vec3f(0, 0, 1);
    sample.pdf = 0.5f;
    sample.flags = BxDFFlags::Transmission | BxDFFlags::Specular;
    RGBBSDFSampleResult result = sample;

    expect_true(result.has_value());
    expect_true(result->IsTransmission());
    expect_true(result->IsSpecular());
    expect_true(!result->IsReflection());
  }
}

context("PR3 sampled spectrum and wavelength packets") {
  test_that("SampledSpectrum arithmetic and diagnostics are componentwise") {
    SampledSpectrum a{1, 2, 3, 4};
    SampledSpectrum b{4, 3, 2, 1};
    SampledSpectrum product = a * b;

    expect_true(product[0] == Approx(4));
    expect_true(product[1] == Approx(6));
    expect_true(product[2] == Approx(6));
    expect_true(product[3] == Approx(4));
    expect_true(a.Average() == Approx(2.5));
    expect_true(a.MinComponentValue() == Approx(1));
    expect_true(a.MaxComponentValue() == Approx(4));
    expect_true(a.IsPositive());
    expect_true(!a.HasNaNs());
    expect_true(!a.IsInf());
  }

  test_that("SampledWavelengths follows pbrt shifted sampling and termination") {
    SampledWavelengths uniform = SampledWavelengths::SampleUniform(0.9f, 400.f, 700.f);
    expect_true(uniform[0] == Approx(670));
    expect_true(uniform[1] == Approx(445));
    expect_true(uniform[2] == Approx(520));
    expect_true(uniform[3] == Approx(595));
    expect_true(uniform.PDF(0) == Approx(1.f / 300.f));
    expect_true(uniform.InvariantsHold(400.f, 700.f));

    SampledWavelengths visible = SampledWavelengths::SampleVisible(0.37f);
    Float originalPDF0 = visible.PDF(0);
    expect_true(!visible.SecondaryTerminated());
    visible.TerminateSecondary();
    expect_true(visible.SecondaryTerminated());
    expect_true(visible.PDF(0) == Approx(originalPDF0 / NSpectrumSamples));
    expect_true(visible.PDF(1) == Approx(0));
    visible.TerminateSecondary();
    expect_true(visible.PDF(0) == Approx(originalPDF0 / NSpectrumSamples));
  }
}

context("PR4 spectrum representations and named assets") {
  test_that("spectrum implementations evaluate and validate setup data") {
    ConstantSpectrum constant(2.5f);
    expect_true(constant(410) == Approx(2.5));

    PiecewiseLinearSpectrum piecewise({400, 500, 600}, {1, 3, 7});
    expect_true(piecewise(450) == Approx(2));
    expect_true(piecewise(399) == Approx(0));

    SpectrumDataPolicy sortedPolicy;
    sortedPolicy.order = SpectrumOrderPolicy::Sort;
    PiecewiseLinearSpectrum sorted({600, 400, 500}, {7, 1, 3}, sortedPolicy);
    expect_true(sorted(450) == Approx(2));

    expect_error(PiecewiseLinearSpectrum({500, 400}, {1, 2}));
    expect_error(PiecewiseLinearSpectrum({400, 400}, {1, 2}));
    expect_error(PiecewiseLinearSpectrum({400, 500}, {1, -1}));

    DenselySampledSpectrum dense(400, 402, {10, 20, 30});
    expect_true(dense(401.51f) == Approx(30));

    expect_true(Blackbody(483, 6000) == Approx(3.1849e13).epsilon(0.001));
    BlackbodySpectrum blackbody(5000);
    expect_true(blackbody.MaxValue() == Approx(1));
  }
}
#endif
