#include "src/render/spectral_light.h"

#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace rayrender::base;
using namespace rayrender::render;

namespace {

constexpr Float Pi = static_cast<Float>(3.14159265358979323846264338327950288);

bool Approx(Float lhs, Float rhs, Float tolerance = static_cast<Float>(1e-5)) {
  return std::fabs(lhs - rhs) <= tolerance;
}

void Check(bool condition, const char* message) {
  if (!condition) {
    std::cerr << "PR11 Light test failed: " << message << std::endl;
    std::exit(1);
  }
}

void CheckApprox(Float actual, Float expected, Float tolerance, const char* message) {
  if (!Approx(actual, expected, tolerance)) {
    std::cerr << "PR11 Light test failed: " << message
              << " expected " << expected << " got " << actual << std::endl;
    std::exit(1);
  }
}

void CheckSpectrumConstant(const SampledSpectrum& spectrum, Float expected, const char* message) {
  for (int i = 0; i < NSpectrumSamples; ++i) {
    if (!Approx(spectrum[i], expected, static_cast<Float>(1e-4))) {
      std::cerr << "PR11 Light test failed: " << message
                << " component " << i << " expected " << expected
                << " got " << spectrum[i] << std::endl;
      std::exit(1);
    }
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
    std::cerr << "PR11 Light test failed: " << message
              << " (unexpected diagnostic: " << e.what() << ")" << std::endl;
    std::exit(1);
  }
  Check(false, message);
}

LightSampleContext ReferenceContext(point3f p = point3f(0, 0, 0)) {
  Interaction interaction;
  interaction.p = p;
  interaction.n = normal3f(0, 1, 0);
  interaction.pError = vec3f(static_cast<Float>(1e-4));
  interaction.time = static_cast<Float>(0.25);
  return LightSampleContext(interaction);
}

SampledWavelengths TestWavelengths() {
  return SampledWavelengths::SampleUniform(static_cast<Float>(0.31));
}

RGBColorSpace ColorSpaceWithConstantIlluminant(Float value) {
  RGBColorSpace colorSpace = RGBColorSpace::SRGB();
  colorSpace.illuminant = std::make_shared<const DenselySampledSpectrum>(
    DenselySampledSpectrum::SampleFunction([&](Float) { return value; })
  );
  return colorSpace;
}

void TestLightSpectrumRGBIlluminantFallback() {
  RGBColorSpace colorSpace = ColorSpaceWithConstantIlluminant(2);
  LightSpectrum spectrum = LightSpectrum::FromRGBIlluminant(colorSpace, RGB(3, 3, 3));
  CheckSpectrumConstant(spectrum.Sample(TestWavelengths()), 6, "neutral RGB lights multiply by illuminant");
}

void TestPointSpotAndDistantLights() {
  SampledWavelengths lambda = TestWavelengths();
  Light point = Light::Point(point3f(0, 2, 0), LightSpectrum::Constant(8));
  std::optional<LightLiSample> ps = point.SampleLi(
    ReferenceContext(),
    point2f(static_cast<Float>(0.2), static_cast<Float>(0.4)),
    lambda
  );
  Check(ps.has_value(), "point light returns a sample");
  Check(ps->delta, "point light sample is delta");
  CheckApprox(ps->pdf, 1, static_cast<Float>(1e-6), "point light sample pdf");
  CheckApprox(point.PDF_Li(ReferenceContext(), ps->wi), 0, static_cast<Float>(1e-6), "point PDF_Li is zero");
  CheckSpectrumConstant(ps->L, 2, "point light inverse-square radiance");
  SpawnedRay pointShadow = ps->visibility.SpawnRay();
  Check(pointShadow.ray.tMax < 1, "visibility endpoint caps shadow ray at target");

  Light spot = Light::SpotFromCosines(
    point3f(0, 2, 0),
    vec3f(0, -1, 0),
    LightSpectrum::Constant(4),
    static_cast<Float>(0.5),
    static_cast<Float>(0.75)
  );
  std::optional<LightLiSample> ss = spot.SampleLi(
    ReferenceContext(),
    point2f(static_cast<Float>(0.1), static_cast<Float>(0.7)),
    lambda
  );
  Check(ss.has_value(), "spot light samples inside cone");
  CheckSpectrumConstant(ss->L, 1, "spot light center cone inverse-square radiance");

  Light distant = Light::Distant(vec3f(0, 0, 1), LightSpectrum::Constant(5));
  distant.Preprocess(Bounds3f(point3f(-2, -2, -2), point3f(2, 2, 2)));
  std::optional<LightLiSample> ds = distant.SampleLi(
    ReferenceContext(),
    point2f(static_cast<Float>(0.6), static_cast<Float>(0.6)),
    lambda
  );
  Check(ds.has_value(), "distant light returns a sample");
  Check(ds->delta, "distant light is delta direction");
  CheckSpectrumConstant(ds->L, 5, "distant light radiance");
  CheckApprox(distant.PDF_Li(ReferenceContext(), ds->wi), 0, static_cast<Float>(1e-6), "distant PDF_Li is zero");
}

void TestDiffuseAreaLightPDFAndBinding() {
  SampledWavelengths lambda = TestWavelengths();
  Shape sphere = Shape::Sphere(1);
  DiffuseAreaLightOptions twoSided;
  twoSided.twoSided = true;
  Light area = Light::DiffuseArea(sphere, LightSpectrum::Constant(7), twoSided);
  LightSampleContext ctx = ReferenceContext(point3f(0, 0, -4));
  std::optional<LightLiSample> sample = area.SampleLi(
    ctx,
    point2f(static_cast<Float>(0.75), static_cast<Float>(0.5)),
    lambda
  );
  Check(sample.has_value(), "area light can sample shape");
  Check(!sample->delta, "area light is not delta");
  CheckApprox(
    sample->pdf,
    sphere.PDFDirection(ctx.ToInteraction(), sample->wi),
    static_cast<Float>(2e-5),
    "area light sample PDF agrees with shape solid-angle PDF"
  );
  CheckSpectrumConstant(area.L(sample->pLight, -sample->wi, lambda), 7, "sampled area emission matches direct L");
  CheckSpectrumConstant(sample->L, 7, "area SampleLi returns emitted radiance");

  PrimitiveBinding binding = CompileAreaLightPrimitiveBinding(
    MaterialHandle::FromIndex(2, 1),
    LightHandle::FromIndex(3, 1),
    sphere
  );
  Check(binding.material == MaterialHandle::FromIndex(2, 1), "area primitive keeps material");
  Check(binding.areaLight == LightHandle::FromIndex(3, 1), "area primitive stores light handle");
  Check(binding.shapeRequirements.areaLight, "area primitive requests area-light capability");

  CSGShapeOptions options;
  options.bounds = Bounds3f(point3f(-1, -1, -1), point3f(1, 1, 1));
  Shape csg = Shape::CSGSignedDistance(
    [](const point3f& p) {
      return p.length() - static_cast<Float>(0.5);
    },
    options,
    "implicit sphere"
  );
  CheckThrowsContaining(
    [&]() {
      (void)CompileAreaLightPrimitiveBinding(
        MaterialHandle::FromIndex(1, 1),
        LightHandle::FromIndex(1, 1),
        csg
      );
    },
    "csg_mesh_sampler",
    "CSG area lights are rejected without render-mesh sampling"
  );
}

void TestInfiniteLightPDFAndEmission() {
  SampledWavelengths lambda = TestWavelengths();
  Light uniform = Light::UniformInfinite(LightSpectrum::Constant(3));
  uniform.Preprocess(Bounds3f(point3f(-1, -1, -1), point3f(1, 1, 1)));
  std::optional<LightLiSample> us = uniform.SampleLi(
    ReferenceContext(),
    point2f(static_cast<Float>(0.4), static_cast<Float>(0.8)),
    lambda
  );
  Check(us.has_value(), "uniform infinite light samples directions");
  CheckApprox(us->pdf, static_cast<Float>(1) / (static_cast<Float>(4) * Pi), static_cast<Float>(1e-6), "uniform infinite PDF");
  CheckApprox(uniform.PDF_Li(ReferenceContext(), us->wi), us->pdf, static_cast<Float>(1e-6), "uniform infinite PDF_Li agrees");
  CheckSpectrumConstant(
    uniform.Le(Ray(point3f(0, 0, 0), us->wi), lambda),
    3,
    "uniform infinite Le agrees with SampleLi"
  );

  std::vector<LightSpectrum> texels = {
    LightSpectrum::Constant(1),
    LightSpectrum::Constant(3),
    LightSpectrum::Constant(1),
    LightSpectrum::Constant(1)
  };
  Light image = Light::ImageInfinite(2, 2, std::move(texels));
  image.Preprocess(Bounds3f(point3f(-2, -2, -2), point3f(2, 2, 2)));
  std::optional<LightLiSample> is = image.SampleLi(
    ReferenceContext(),
    point2f(static_cast<Float>(0.30), static_cast<Float>(0.50)),
    lambda
  );
  Check(is.has_value(), "image infinite samples distribution");
  Float expectedPdf = static_cast<Float>(1) / (static_cast<Float>(2) * Pi);
  CheckApprox(is->pdf, expectedPdf, static_cast<Float>(1e-5), "image infinite PDF includes equal-area 4*pi Jacobian");
  CheckApprox(image.PDF_Li(ReferenceContext(), is->wi), is->pdf, static_cast<Float>(1e-5), "image infinite PDF_Li agrees with sample");
  CheckSpectrumConstant(is->L, 3, "image infinite selected bright texel");
  CheckSpectrumConstant(
    image.Le(Ray(point3f(0, 0, 0), is->wi), lambda),
    3,
    "image infinite sampled radiance matches Le"
  );
}

void TestRegistriesSamplersAndAdapters() {
  SampledWavelengths lambda = TestWavelengths();
  Scene scene;
  SpectralLightTable table;
  LightHandle point = RegisterLight(scene, table, Light::Point(point3f(0, 1, 0), LightSpectrum::Constant(1)));
  LightHandle area = RegisterLight(scene, table, Light::DiffuseArea(Shape::Sphere(1), LightSpectrum::Constant(2)));
  LightHandle env = RegisterLight(scene, table, Light::UniformInfinite(LightSpectrum::Constant(3)));

  Check(scene.Lights().size() == 2, "finite and area lights are scene finite lights");
  Check(scene.InfiniteLights().size() == 1, "infinite lights are scene infinite lights");
  Check(scene.ShapeCount() == 0, "registering infinite light does not add environment geometry");
  Check(scene.PrimitiveCount() == 0, "registering infinite light does not add environment primitive");
  Check(!ShouldCompileLegacyEnvironmentSphereToSpectralGeometry(), "legacy environment sphere is excluded");

  UniformLightSampler uniform(table, {point, area, env});
  CheckApprox(uniform.PMFSum(), 1, static_cast<Float>(1e-6), "uniform sampler PMFs sum to one");
  std::optional<SampledLight> uniformSample = uniform.Sample(static_cast<Float>(0.7));
  Check(uniformSample.has_value(), "uniform sampler returns a light");
  CheckApprox(uniformSample->pmf, uniform.PMF(uniformSample->handle), static_cast<Float>(1e-6), "uniform sample PMF agrees");

  table.Preprocess(Bounds3f(point3f(-1, -1, -1), point3f(1, 1, 1)));
  PowerLightSampler power(table, {point, area, env});
  CheckApprox(power.PMFSum(), 1, static_cast<Float>(1e-5), "power sampler PMFs sum to one");
  Float pmfBefore = power.PMF(point);
  (void)table.Get(point).Phi(lambda);
  (void)table.Get(point).Phi(SampledWavelengths::SampleUniform(static_cast<Float>(0.81)));
  CheckApprox(power.PMF(point), pmfBefore, static_cast<Float>(1e-7), "power sampler PMF is wavelength independent");

  RGBColorSpace colorSpace = ColorSpaceWithConstantIlluminant(1);
  LegacyEmissiveMaterialDescriptor diffuse;
  diffuse.color = RGB(2, 2, 2);
  diffuse.intensity = 3;
  DiffuseAreaLightOptions legacyAreaOptions;
  legacyAreaOptions.twoSided = true;
  Light legacyArea = ConvertLegacyEmissiveMaterialToAreaLight(
    diffuse,
    Shape::Sphere(1),
    colorSpace,
    legacyAreaOptions
  );
  std::optional<LightLiSample> legacyAreaSample = legacyArea.SampleLi(
    ReferenceContext(point3f(0, 0, -4)),
    point2f(static_cast<Float>(0.75), static_cast<Float>(0.5)),
    lambda
  );
  Check(legacyAreaSample.has_value(), "legacy diffuse_light descriptor converts to area light");
  CheckSpectrumConstant(legacyAreaSample->L, 6, "legacy diffuse_light RGB uses illuminant reconstruction");

  LegacySpotLightMaterialDescriptor legacySpot;
  legacySpot.position = point3f(0, 2, 0);
  legacySpot.direction = vec3f(0, -1, 0);
  legacySpot.color = RGB(2, 2, 2);
  legacySpot.intensity = 2;
  Light convertedSpot = ConvertLegacySpotLightMaterialToLight(legacySpot, colorSpace);
  std::optional<LightLiSample> convertedSpotSample = convertedSpot.SampleLi(
    ReferenceContext(),
    point2f(static_cast<Float>(0.2), static_cast<Float>(0.2)),
    lambda
  );
  Check(convertedSpotSample.has_value(), "legacy spot_light descriptor converts to free SpotLight");
  CheckSpectrumConstant(convertedSpotSample->L, 1, "legacy spot_light RGB uses illuminant reconstruction");

  Check(!SpectralMaterialEmissionPathEnabled(), "spectral material emission path remains disabled");
}

} // namespace

int main() {
  TestLightSpectrumRGBIlluminantFallback();
  TestPointSpotAndDistantLights();
  TestDiffuseAreaLightPDFAndBinding();
  TestInfiniteLightPDFAndEmission();
  TestRegistriesSamplersAndAdapters();
  return 0;
}
