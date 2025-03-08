// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#include <gtest/gtest.h>

#include <pbrt/cameras.h>
#include <pbrt/cpu/aggregates.h>
#include <pbrt/cpu/integrators.h>
#include <pbrt/filters.h>
#include <pbrt/lights.h>
#include <pbrt/materials.h>
#include <pbrt/options.h>
#include <pbrt/pbrt.h>
#include <pbrt/samplers.h>
#include <pbrt/shapes.h>
#include <pbrt/textures.h>
#include <pbrt/util/color.h>
#include <pbrt/util/colorspace.h>
#include <pbrt/util/image.h>
#include <pbrt/util/spectrum.h>
#include <pbrt/util/vecmath.h>

#include <memory>

using namespace pbrt;

static std::string inTestDir(const std::string &path) {
    return path;
}

struct TestScene {
    Primitive aggregate;
    std::vector<Light> lights;
    std::string description;
    float expected;
};

struct TestIntegrator {
    Integrator *integrator;
    const Film film;
    std::string description;
    TestScene scene;
};

void PrintTo(const TestIntegrator &tr, ::std::ostream *os) {
    *os << tr.description;
}

void CheckSceneAverage(const std::string &filename, float expected) {
    pstd::optional<ImageAndMetadata> im = Image::Read(filename);
    ASSERT_TRUE((bool)im);
    ASSERT_EQ(im->image.NChannels(), 3);

    float delta = .025;
    float sum = 0;

    Image &image = im->image;
    for (int t = 0; t < image.Resolution()[1]; ++t)
        for (int s = 0; s < image.Resolution()[0]; ++s)
            for (int c = 0; c < 3; ++c)
                sum += image.GetChannel(Point2i(s, t), c);
    int nPixels = image.Resolution().x * image.Resolution().y * 3;
    EXPECT_NEAR(expected, sum / nPixels, delta);
}

std::vector<TestScene> GetScenes() {
    std::vector<TestScene> scenes;

    Allocator alloc;
    static Transform identity;
    {
        // Unit sphere, Kd = 0.5, point light I = 3.1415 at center
        // -> With GI, should have radiance of 1.
        Shape sphere = new Sphere(&identity, &identity, true /* reverse orientation */, 1,
                                  -1, 1, 360);

        static ConstantSpectrum cs(0.5);
        SpectrumTexture Kd = alloc.new_object<SpectrumConstantTexture>(&cs);
        // FIXME: here and below, Materials leak...
        Material material = new DiffuseMaterial(Kd, nullptr, nullptr);

        MediumInterface mediumInterface;
        std::vector<Primitive> prims;
        prims.push_back(Primitive(
            new GeometricPrimitive(sphere, material, nullptr, mediumInterface)));
        Primitive bvh(new BVHAggregate(std::move(prims)));

        // We have to do this little dance here to make sure the spectrum is
        // properly normalized (this is usually all handled inside *Light::Create())
        ConstantSpectrum I(1);
        Float scale = Pi / SpectrumToPhotometric(&I);
        std::vector<Light> lights;
        lights.push_back(new PointLight(identity, MediumInterface(), &I, scale));

        scenes.push_back({bvh, lights, "Sphere, 1 light, Kd = 0.5", 1.0});
    }

    {
        // Unit sphere, Kd = 0.5, 4 point lights I = 3.1415/4 at center
        // -> With GI, should have radiance of 1.
        Shape sphere = new Sphere(&identity, &identity, true /* reverse orientation */, 1,
                                  -1, 1, 360);

        static ConstantSpectrum cs(0.5);
        SpectrumTexture Kd = alloc.new_object<SpectrumConstantTexture>(&cs);
        const Material material = new DiffuseMaterial(Kd, nullptr, nullptr);

        MediumInterface mediumInterface;
        std::vector<Primitive> prims;
        prims.push_back(Primitive(
            new GeometricPrimitive(sphere, material, nullptr, mediumInterface)));
        Primitive bvh(new BVHAggregate(std::move(prims)));

        // We have to do this little dance here to make sure the spectrum is
        // properly normalized (this is usually all handled inside *Light::Create())
        ConstantSpectrum I(1);
        Float scale = Pi / (4 * SpectrumToPhotometric(&I));
        std::vector<Light> lights;
        lights.push_back(new PointLight(identity, MediumInterface(), &I, scale));
        lights.push_back(new PointLight(identity, MediumInterface(), &I, scale));
        lights.push_back(new PointLight(identity, MediumInterface(), &I, scale));
        lights.push_back(new PointLight(identity, MediumInterface(), &I, scale));

        scenes.push_back({bvh, lights, "Sphere, 1 light, Kd = 0.5", 1.0});
    }

    {
        // Unit sphere, Kd = 0.5, Le = 0.5
        // -> With GI, should have radiance of 1.
        Shape sphere = new Sphere(&identity, &identity, true /* reverse orientation */, 1,
                                  -1, 1, 360);

        static ConstantSpectrum cs(0.5);
        SpectrumTexture Kd = alloc.new_object<SpectrumConstantTexture>(&cs);
        const Material material = new DiffuseMaterial(Kd, nullptr, nullptr);

        // We have to do this little dance here to make sure the spectrum is
        // properly normalized (this is usually all handled inside *Light::Create())
        ConstantSpectrum Le(1);
        Float scale = 0.5 / SpectrumToPhotometric(&Le);
        Light areaLight = new DiffuseAreaLight(identity, MediumInterface(), &Le, scale,
                                               sphere, nullptr, Image(), nullptr, false);

        std::vector<Light> lights;
        lights.push_back(areaLight);

        MediumInterface mediumInterface;
        std::vector<Primitive> prims;
        prims.push_back(Primitive(
            new GeometricPrimitive(sphere, material, lights.back(), mediumInterface)));
        Primitive bvh(new BVHAggregate(std::move(prims)));

        scenes.push_back({bvh, lights, "Sphere, Kd = 0.5, Le = 0.5", 1.0});
    }

#if 0
    {
        // Unit sphere, Kd = 0.25, Kr = .5, point light I = 7.4 at center
        // -> With GI, should have radiance of ~1.
        Shape sphere = new Sphere(
            &identity, &identity, true /* reverse orientation */, 1, -1, 1, 360);

        static ConstantSpectrum cs5(0.5), cs25(0.25);
        SpectrumTexture Kd =
            alloc.new_object<SpectrumConstantTexture>(&cs25);
        SpectrumTexture Kr =
            alloc.new_object<SpectrumConstantTexture>(&cs5);
        SpectrumTexture black =
            alloc.new_object<SpectrumConstantTexture>(Spectra::Zero());
        SpectrumTexture white =
            alloc.new_object<SpectrumConstantTexture>(Spectra::One());
        FloatTexture zero =
            alloc.new_object<FloatConstantTexture>(0.);
        FloatTexture one =
            alloc.new_object<FloatConstantTexture>(1.);
        const Material material = new UberMaterial(
            Kd, black, Kr, black, zero, zero, one, nullptr, false, nullptr);

        MediumInterface mediumInterface;
        std::vector<Primitive> prims;
        prims.push_back(Primitive(new GeometricPrimitive(
            sphere, material, nullptr, mediumInterface)));
        Primitive bvh(new BVHAggregate(std::move(prims)));

        static ConstantSpectrum I(3. * Pi);
        std::vector<Light> lights;
        lights.push_back(std::make_unique<PointLight>(identity,
                                                      nullptr, &I, Allocator()));

        scenes.push_back({bvh, lights, "Sphere, 1 light, Kd = 0.25 Kr = 0.5", 1.0});
    }
#endif

#if 0
  {
    // Unit sphere, Kd = 0.25, Kr = .5, Le .587
    // -> With GI, should have radiance of ~1.
    Shape sphere = new Sphere(
        &identity, &identity, true /* reverse orientation */, 1, -1, 1, 360);

    static ConstantSpectrum cs5(0.5), cs25(0.25);
    SpectrumTexture Kd =
        alloc.new_object<SpectrumConstantTexture>(&cs25);
    SpectrumTexture Kr =
        alloc.new_object<SpectrumConstantTexture>(&cs5);
    SpectrumTexture black =
        alloc.new_object<SpectrumConstantTexture>(Spectra::Zero());
    SpectrumTexture white =
        alloc.new_object<SpectrumConstantTexture>(Spectra::One());
    FloatTexture zero =
        alloc.new_object<FloatConstantTexture>(0.);
    FloatTexture one =
        alloc.new_object<FloatConstantTexture>(1.);
    std::shared_ptr<Material> material = std::make_shared<UberMaterial>(
        Kd, black, Kr, black, zero, zero, zero, white, one, nullptr, false, nullptr);

    static ConstantSpectrum Le(0.587);
    std::shared_ptr<AreaLight> areaLight = std::make_shared<DiffuseAreaLight>(
        identity, nullptr, &Le, 8, sphere, nullptr, true, false,
        std::make_shared<ParameterDictionary>(std::initializer_list<const NamedValues *>{}, nullptr));

    MediumInterface mediumInterface;
    std::vector<std::shared_ptr<Primitive>> prims;
    prims.push_back(Primitive(new GeometricPrimitive(
        sphere, material, areaLight, mediumInterface)));
    Primitive bvh(new BVHAggregate(std::move(prims)));

    std::vector<std::shared_ptr<Light>> lights;
    lights.push_back(std::move(areaLight));

    scenes.push_back({bvh, lights, "Sphere, Kd = 0.25 Kr = 0.5, Le = 0.587", 1.0});
  }
#endif

    return scenes;
}

std::vector<std::pair<Sampler, std::string>> GetSamplers(const Point2i &resolution) {
    static std::vector<std::pair<Sampler, std::string>> samplers;
    if (!samplers.empty())
        return samplers;

    int sqrtSpp = 16;
    if (NSpectrumSamples < 4)
        sqrtSpp *= 2;
    int spp = sqrtSpp * sqrtSpp;

    samplers.push_back(std::make_pair(new HaltonSampler(spp, resolution), "Halton"));
    samplers.push_back(
        std::make_pair(new PaddedSobolSampler(spp, RandomizeStrategy::PermuteDigits),
                       "Padded Sobolspp"));
    samplers.push_back(std::make_pair(
        new ZSobolSampler(spp, Point2i(16, 16), RandomizeStrategy::PermuteDigits),
        "Z Sobol"));
    samplers.push_back(
        std::make_pair(new SobolSampler(spp, resolution, RandomizeStrategy::None),
                       "Sobol Not Randomized"));
    samplers.push_back(std::make_pair(
        new SobolSampler(spp, resolution, RandomizeStrategy::PermuteDigits),
        "Sobol XOR Scramble"));
    samplers.push_back(
        std::make_pair(new SobolSampler(spp, resolution, RandomizeStrategy::Owen),
                       "Sobol Owen Scramble"));
    samplers.push_back(std::make_pair(new IndependentSampler(spp), "Independent"));
    samplers.push_back(
        std::make_pair(new StratifiedSampler(sqrtSpp,sqrtSpp, true), "Stratified"));
    samplers.push_back(std::make_pair(new PMJ02BNSampler(spp), "PMJ02bn"));

    return samplers;
}

std::vector<TestIntegrator> GetIntegrators() {
    std::vector<TestIntegrator> integrators;

    Point2i resolution(10, 10);
    static Transform id;
    AnimatedTransform identity(id, 0, id, 1);
    for (const auto &scene : GetScenes()) {
        // Path tracing integrators
        for (auto &sampler : GetSamplers(resolution)) {
            Filter filter = new BoxFilter(Vector2f(0.5, 0.5));
            FilmBaseParameters fp(resolution, Bounds2i(Point2i(0, 0), resolution), filter,
                                  1., PixelSensor::CreateDefault(),
                                  inTestDir("test.exr"));
            RGBFilm *film = new RGBFilm(fp, RGBColorSpace::sRGB);
            CameraBaseParameters cbp(CameraTransform(identity), film, nullptr, {},
                                     nullptr);
            PerspectiveCamera *camera = new PerspectiveCamera(
                cbp, 45, Bounds2f(Point2f(-1, -1), Point2f(1, 1)), 0., 10.);

            const Film filmp = camera->GetFilm();
            Integrator *integrator = new PathIntegrator(8, camera, sampler.first,
                                                        scene.aggregate, scene.lights);
            integrators.push_back({integrator, filmp,
                                   "Path, depth 8, Perspective, " + sampler.second +
                                       ", " + scene.description,
                                   scene});
        }

        for (auto &sampler : GetSamplers(resolution)) {
            Filter filter = new BoxFilter(Vector2f(0.5, 0.5));
            FilmBaseParameters fp(resolution, Bounds2i(Point2i(0, 0), resolution), filter,
                                  1., PixelSensor::CreateDefault(),
                                  inTestDir("test.exr"));
            RGBFilm *film = new RGBFilm(fp, RGBColorSpace::sRGB);
            CameraBaseParameters cbp(CameraTransform(identity), film, nullptr, {},
                                     nullptr);
            OrthographicCamera *camera = new OrthographicCamera(
                cbp, Bounds2f(Point2f(-.1, -.1), Point2f(.1, .1)), 0., 10.);
            const Film filmp = camera->GetFilm();

            Integrator *integrator = new PathIntegrator(8, camera, sampler.first,
                                                        scene.aggregate, scene.lights);
            integrators.push_back(
                {integrator, filmp,
                 "Path, depth 8, Ortho, " + sampler.second + ", " + scene.description,
                 scene});
        }

        // Volume path tracing integrators
        for (auto &sampler : GetSamplers(resolution)) {
            Filter filter = new BoxFilter(Vector2f(0.5, 0.5));
            FilmBaseParameters fp(resolution, Bounds2i(Point2i(0, 0), resolution), filter,
                                  1., PixelSensor::CreateDefault(),
                                  inTestDir("test.exr"));
            RGBFilm *film = new RGBFilm(fp, RGBColorSpace::sRGB);
            CameraBaseParameters cbp(CameraTransform(identity), film, nullptr, {},
                                     nullptr);
            PerspectiveCamera *camera = new PerspectiveCamera(
                cbp, 45, Bounds2f(Point2f(-1, -1), Point2f(1, 1)), 0., 10.);
            const Film filmp = camera->GetFilm();

            Integrator *integrator = new VolPathIntegrator(8, camera, sampler.first,
                                                           scene.aggregate, scene.lights);
            integrators.push_back({integrator, filmp,
                                   "VolPath, depth 8, Perspective, " + sampler.second +
                                       ", " + scene.description,
                                   scene});
        }
        for (auto &sampler : GetSamplers(resolution)) {
            Filter filter = new BoxFilter(Vector2f(0.5, 0.5));
            FilmBaseParameters fp(resolution, Bounds2i(Point2i(0, 0), resolution), filter,
                                  1., PixelSensor::CreateDefault(),
                                  inTestDir("test.exr"));
            RGBFilm *film = new RGBFilm(fp, RGBColorSpace::sRGB);
            CameraBaseParameters cbp(CameraTransform(identity), film, nullptr, {},
                                     nullptr);
            OrthographicCamera *camera = new OrthographicCamera(
                cbp, Bounds2f(Point2f(-.1, -.1), Point2f(.1, .1)), 0., 10.);
            const Film filmp = camera->GetFilm();

            Integrator *integrator = new VolPathIntegrator(8, camera, sampler.first,
                                                           scene.aggregate, scene.lights);
            integrators.push_back(
                {integrator, filmp,
                 "VolPath, depth 8, Ortho, " + sampler.second + ", " + scene.description,
                 scene});
        }

        // Simple path (perspective only, still sample light and BSDFs). Yolo
        for (auto &sampler : GetSamplers(resolution)) {
            Filter filter = new BoxFilter(Vector2f(0.5, 0.5));
            FilmBaseParameters fp(resolution, Bounds2i(Point2i(0, 0), resolution), filter,
                                  1., PixelSensor::CreateDefault(),
                                  inTestDir("test.exr"));
            RGBFilm *film = new RGBFilm(fp, RGBColorSpace::sRGB);
            CameraBaseParameters cbp(CameraTransform(identity), film, nullptr, {},
                                     nullptr);
            PerspectiveCamera *camera = new PerspectiveCamera(
                cbp, 45, Bounds2f(Point2f(-1, -1), Point2f(1, 1)), 0., 10.);

            const Film filmp = camera->GetFilm();
            Integrator *integrator = new SimplePathIntegrator(
                8, true, true, camera, sampler.first, scene.aggregate, scene.lights);
            integrators.push_back({integrator, filmp,
                                   "SimplePath, depth 8, Perspective, " + sampler.second +
                                       ", " + scene.description,
                                   scene});
        }

        // BDPT
        for (auto &sampler : GetSamplers(resolution)) {
            Filter filter = new BoxFilter(Vector2f(0.5, 0.5));
            FilmBaseParameters fp(resolution, Bounds2i(Point2i(0, 0), resolution), filter,
                                  1., PixelSensor::CreateDefault(),
                                  inTestDir("test.exr"));
            RGBFilm *film = new RGBFilm(fp, RGBColorSpace::sRGB);
            CameraBaseParameters cbp(CameraTransform(identity), film, nullptr, {},
                                     nullptr);
            PerspectiveCamera *camera = new PerspectiveCamera(
                cbp, 45, Bounds2f(Point2f(-1, -1), Point2f(1, 1)), 0., 10.);
            const Film filmp = camera->GetFilm();

            Integrator *integrator =
                new BDPTIntegrator(camera, sampler.first, scene.aggregate, scene.lights,
                                   6, false, false, false);
            integrators.push_back({integrator, filmp,
                                   "BDPT, depth 8, Perspective, " + sampler.second +
                                       ", " + scene.description,
                                   scene});
        }

        // MLT
        {
            Filter filter = new BoxFilter(Vector2f(0.5, 0.5));
            FilmBaseParameters fp(resolution, Bounds2i(Point2i(0, 0), resolution), filter,
                                  1., PixelSensor::CreateDefault(),
                                  inTestDir("test.exr"));
            RGBFilm *film = new RGBFilm(fp, RGBColorSpace::sRGB);
            CameraBaseParameters cbp(CameraTransform(identity), film, nullptr, {},
                                     nullptr);
            PerspectiveCamera *camera = new PerspectiveCamera(
                cbp, 45, Bounds2f(Point2f(-1, -1), Point2f(1, 1)), 0., 10.);
            const Film filmp = camera->GetFilm();

            Integrator *integrator =
                new MLTIntegrator(camera, scene.aggregate, scene.lights, 8 /* depth */,
                                  100000 /* n bootstrap */, 1000 /* nchains */,
                                  1024 /* mutations per pixel */, 0.01 /* sigma */,
                                  0.3 /* large step prob */, false /* regularize */);
            integrators.push_back({integrator, filmp,
                                   "MLT, depth 8, Perspective, " + scene.description,
                                   scene});
        }
    }

    return integrators;
}

struct RenderTest : public testing::TestWithParam<TestIntegrator> {};

TEST_P(RenderTest, RadianceMatches) {
    const TestIntegrator &tr = GetParam();
    tr.integrator->Render();
    CheckSceneAverage(inTestDir("test.exr"), tr.scene.expected);
    // The SpatialLightSampler class keeps a per-thread cache that
    // must be cleared out between test runs. In turn, this means that we
    // must delete the Integrator here in order to make sure that its
    // destructor runs. (This is ugly and should be fixed in a better way.)
    delete tr.integrator;

    EXPECT_EQ(0, remove(inTestDir("test.exr").c_str()));
}

INSTANTIATE_TEST_CASE_P(AnalyticTestScenes, RenderTest,
                        testing::ValuesIn(GetIntegrators()));
