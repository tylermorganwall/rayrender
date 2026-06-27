#include "src/base/base.h"
#include "src/math/vectypes.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>

using namespace rayrender::base;

static_assert(!std::is_convertible<Float, RGB>::value,
              "RGB scalar construction must remain explicit");
static_assert(!std::is_convertible<Float, XYZ>::value,
              "XYZ scalar construction must remain explicit");
static_assert(!std::is_constructible<RGB, point3f>::value,
              "RGB must not be constructible from point3f");
static_assert(!std::is_constructible<XYZ, vec3f>::value,
              "XYZ must not be constructible from vec3f");
static_assert(!std::is_same<SpectrumHandle, TextureHandle>::value,
              "Tagged handle aliases must remain distinct types");

namespace {

bool Approx(Float lhs, Float rhs, Float tolerance = static_cast<Float>(1e-5)) {
  return std::fabs(lhs - rhs) <= tolerance;
}

void Check(bool condition, const char* message) {
  if (!condition) {
    std::cerr << "PR2 base test failed: " << message << std::endl;
    std::exit(1);
  }
}

template <typename Handle>
void CheckHandle(DispatchTag expectedTag, const char* expectedName) {
  Handle invalid;
  Check(!invalid.IsValid(), "default handle is invalid");

  Handle handle = Handle::FromIndex(42, 9);
  Check(handle.IsValid(), "constructed handle is valid");
  Check(handle.Index() == 42, "handle index is retained");
  Check(handle.Generation() == 9, "handle generation is retained");
  Check(Handle::Tag == expectedTag, "handle tag is retained");
  Check(std::string(DispatchTagName(expectedTag)) == expectedName, "dispatch tag name is stable");
}

void TestColorTypes() {
  RGB rgb(0.25f, 0.5f, 0.75f);
  RGB doubled = rgb * 2.0f;
  Check(Approx(doubled.r, 0.5f), "RGB scalar multiply red");
  Check(Approx(doubled.g, 1.0f), "RGB scalar multiply green");
  Check(Approx(doubled.b, 1.5f), "RGB scalar multiply blue");

  RGB safe = SafeDiv(RGB(1, 2, 3), RGB(1, 0, 3), -1);
  Check(Approx(safe.r, 1), "SafeDiv finite component");
  Check(Approx(safe.g, -1), "SafeDiv zero denominator fallback");
  Check(Approx(safe.b, 1), "SafeDiv third component");

  RGB nonfinite(std::numeric_limits<Float>::infinity(), 0, 1);
  Check(!nonfinite.IsFinite(), "RGB finite check rejects infinity");

  ColorSpaceMatrix3x3 identity = ColorSpaceMatrix3x3::Identity();
  XYZ xyz = TransformRGBToXYZ(identity, RGB(1, 2, 3));
  Check(Approx(xyz.x, 1), "identity RGB to XYZ x");
  Check(Approx(xyz.y, 2), "identity RGB to XYZ y");
  Check(Approx(xyz.z, 3), "identity RGB to XYZ z");

  RGBColorSpace srgb = RGBColorSpace::SRGB();
  RGB encoded(0.25f, 0.5f, 0.75f);
  RGB linear = srgb.Decode(encoded);
  RGB encodedRoundTrip = srgb.Encode(linear);
  Check(Approx(encodedRoundTrip.r, encoded.r), "sRGB encode/decode red");
  Check(Approx(encodedRoundTrip.g, encoded.g), "sRGB encode/decode green");
  Check(Approx(encodedRoundTrip.b, encoded.b), "sRGB encode/decode blue");

  RGB linearRoundTrip = srgb.ToLinearRGB(srgb.ToXYZ(linear));
  Check(Approx(linearRoundTrip.r, linear.r), "sRGB matrix round-trip red");
  Check(Approx(linearRoundTrip.g, linear.g), "sRGB matrix round-trip green");
  Check(Approx(linearRoundTrip.b, linear.b), "sRGB matrix round-trip blue");
}

void TestHandles() {
  CheckHandle<SpectrumHandle>(DispatchTag::Spectrum, "Spectrum");
  CheckHandle<TextureHandle>(DispatchTag::Texture, "Texture");
  CheckHandle<MaterialHandle>(DispatchTag::Material, "Material");
  CheckHandle<BxDFHandle>(DispatchTag::BxDF, "BxDF");
  CheckHandle<LightHandle>(DispatchTag::Light, "Light");
  CheckHandle<MediumHandle>(DispatchTag::Medium, "Medium");
  CheckHandle<PhaseFunctionHandle>(DispatchTag::PhaseFunction, "PhaseFunction");
}

struct alignas(16) Align16 {
  int value = 0;
};

struct alignas(32) Align32 {
  int value = 0;
};

struct alignas(64) TrackedAlign64 {
  explicit TrackedAlign64(int* counter) : counter(counter) {}
  ~TrackedAlign64() noexcept { ++(*counter); }
  int* counter = nullptr;
};

void TestScratchBuffer() {
  int destructorCount = 0;
  ScratchBuffer scratch(512);

  Align16* align16 = scratch.Create<Align16>();
  Align32* align32 = scratch.Create<Align32>();
  TrackedAlign64* align64 = scratch.Create<TrackedAlign64>(&destructorCount);

  Check(reinterpret_cast<std::uintptr_t>(align16) % alignof(Align16) == 0, "Align16 allocation alignment");
  Check(reinterpret_cast<std::uintptr_t>(align32) % alignof(Align32) == 0, "Align32 allocation alignment");
  Check(reinterpret_cast<std::uintptr_t>(align64) % alignof(TrackedAlign64) == 0, "Align64 allocation alignment");
  Check(scratch.Used() <= scratch.Capacity(), "scratch usage stays within capacity");
  Check(scratch.HighWater() == scratch.Used(), "scratch high-water tracks peak usage");
  Check(scratch.TrackedObjectCount() == 1, "scratch tracks nontrivial destructors");

  scratch.Reset();
  Check(destructorCount == 1, "scratch reset destroys tracked object");
  Check(scratch.Used() == 0, "scratch reset clears used bytes");
  Check(scratch.TrackedObjectCount() == 0, "scratch reset clears destructor records");

  Align32* reused = scratch.Create<Align32>();
  Check(reinterpret_cast<std::uintptr_t>(reused) % alignof(Align32) == 0, "scratch reuse preserves alignment");
  scratch.Reset();
}

void TestFlagsAndSampling() {
  BxDFFlags flags = BxDFFlags::Reflection | BxDFFlags::Diffuse;
  Check(HasFlag(flags, BxDFFlags::Reflection), "BxDF reflection flag");
  Check(HasAny(flags, BxDFFlags::Reflection | BxDFFlags::Transmission), "BxDF any reflection/transmission flag");
  Check(!HasFlag(flags, BxDFFlags::Specular), "BxDF missing specular flag");

  BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::Reflection | BxDFReflTransFlags::Transmission;
  Check(HasFlag(sampleFlags, BxDFReflTransFlags::Reflection), "sample reflection flag");
  Check(HasFlag(sampleFlags, BxDFReflTransFlags::Transmission), "sample transmission flag");

  RGBBSDFSampleResult none;
  Check(!none.has_value(), "optional sample can represent no sample");

  RGBBSDFSample sample;
  sample.f = RGB(1, 1, 1);
  sample.wi = vec3f(0, 0, 1);
  sample.pdf = 0.5f;
  sample.flags = BxDFFlags::Transmission | BxDFFlags::Specular;
  RGBBSDFSampleResult result = sample;
  Check(result.has_value(), "optional sample can hold a value");
  Check(result->IsTransmission(), "sample transmission predicate");
  Check(result->IsSpecular(), "sample specular predicate");
  Check(!result->IsReflection(), "sample reflection predicate");
}

} // namespace

int main() {
  TestColorTypes();
  TestHandles();
  TestScratchBuffer();
  TestFlagsAndSampling();
  std::cout << "PR2 base tests passed" << std::endl;
  return 0;
}
