#ifdef NOT_CRAN
#include <testthat.h>
#include "wavefront.h"
#include "hlbvh.h"
#include "camera.h"
#include "adaptivesampler.h"
#include "../hitables/rectangle.h"
#include "../hitables/infinite_area_light.h"
#include "../hitables/sphere.h"
#include "../volumes/boundary.h"
#include <chrono>

context("Metal wavefront session") {
  test_that("a numerically repeated dielectric triangle cannot stall a queue") {
    if (MetalBVHAvailable()) {
      // Captured from a glass/milk render on M1 Max: Metal returns a positive
      // 2.3e-7 hit although this outward ray is already outside the first face.
      // Close that face into a tetrahedron and make its interface index matched.
      WavefrontScene scene;
      point3f a(-.213658512f, .150000006f, .434860498f);
      point3f b(-.220242143f, 2.20000005f, .475794971f);
      point3f c(-.235702395f, .150000006f, .463107109f);
      vec3f g(.788308799f, -.00975271221f, .615202606f);
      point3f d = point3f((a + b + c) / 3) - .2f * g;
      point3f center = (a + b + c + d) / 4;
      auto face = [&](point3f x, point3f y, point3f z) {
        WFTriangle t;
        vec3f n = unit_vector(cross(y - x, z - x));
        if (dot(n, x - center) < 0) n = -n;
        t.geometric = {float(n[0]), float(n[1]), float(n[2]), 0};
        point3f vertices[] = {x, y, z};
        for (int k = 0; k < 3; ++k) {
          t.p[k] = {float(vertices[k][0]), float(vertices[k][1]), float(vertices[k][2]), 0};
          t.n[k] = t.geometric;
        }
        t.boundary = 1;
        scene.triangles.push_back(t);
      };
      face(a, b, c); face(a, c, d); face(c, b, d); face(b, a, d);
      scene.triangles[0].geometric = {float(g[0]), float(g[1]), float(g[2]), 0};
      scene.boundaries.push_back(WFBoundary());
      WFMaterial material;
      material.dielectric = 1;
      scene.materials.push_back(material);
      WFTexture white;
      white.a = {1, 1, 1, 0};
      scene.textures.push_back(white);
      WFLight background;
      background.forward = {1, 0, 0, 0};
      background.right = {0, 1, 0, 0};
      background.up = {0, 0, 1, 0};
      scene.lights.push_back(background);
      auto session = MakeMetalWavefront(scene, 1);
      WFParameters p;
      p.width = p.height = p.active = p.lights = p.boundaries = 1;
      p.environments = 1;
      p.max_depth = 4;
      p.camera.orthographic = 1;
      p.camera.lower_left = {-.222824052f, 1.68612814f, .470959127f, 0};
      p.camera.forward = {.579234183f, .655820966f, .484134972f, 0};
      p.boundary_lo = {-1, 0, 0, 0};
      p.boundary_hi = {0, 3, 1, 0};
      const WFPixel *pixels = nullptr;
      WavefrontReport report;
      const auto start = std::chrono::steady_clock::now();
      bool finished = session->Render(p, {0}, scene.lights, [&] {
        return std::chrono::steady_clock::now() - start > std::chrono::seconds(2);
      }, pixels, report);
      expect_true(finished);
      expect_true(report.discarded_paths == 0);
      if (finished) expect_true(std::abs(pixels[0].radiance.x - 1) < 1e-5);
    }
  }
  test_that("a persistent session follows camera and light changes and cancels atomically") {
    if (MetalBVHAvailable()) {
      for (bool with_subsurface : {false, true}) {
        Transform identity, environment, inverse_environment;
        hitable_list world, lights;
        auto white = std::make_shared<constant_texture>(point3f(.1));
        auto background = std::make_shared<ImageInfiniteLight>(white, 8, 4, 0);
        auto disk = std::make_shared<DiskInfiniteLight>(point3f(4), vec3f(0, 0, 1), 60, 0);
        auto sources = std::make_shared<InfiniteLightMixture>(
            std::vector<std::shared_ptr<InfiniteLight>>{background, disk});
        auto diffuse =
            std::make_shared<diffuse_material>(std::make_shared<constant_texture>(point3f(.5)));
        if (with_subsurface) {
          Rcpp::Function make_body = Rcpp::Environment::namespace_env("rayrender")["subsurface"];
          Rcpp::List descriptor =
              make_body(Rcpp::Named("sigma_a") = .03, Rcpp::Named("sigma_s") = 2);
          Rcpp::List parameters = descriptor[0];
          auto medium = std::make_shared<Medium>(Rcpp::as<Rcpp::List>(parameters["subsurface"]));
          auto surface = std::make_shared<dielectric>(point3f(1), 1.3, point3f(0), 0);
          auto shape =
              std::make_shared<sphere>(1, surface, nullptr, nullptr, &identity, &identity, false);
          world.add(std::make_shared<MediumBoundary>(shape, medium, identity, true, 1));
        } else {
          world.add(std::make_shared<xy_rect>(-2, 2, -2, 2, 0, diffuse, nullptr, nullptr, &identity,
                                              &identity, false));
        }
        auto sky = std::make_shared<InfiniteAreaLight>(sources, 100, point3f(0), &environment,
                                                        &inverse_environment);
        world.add(sky);
        camera cam({0, 0, 4}, {0, 0, 0}, {0, 1, 0}, 30, 1, 0, 4, 0, 1, 1);
        WavefrontContext context;
        auto &scene = context.scene;
        auto &session = context.session;
        WavefrontReport report;
        context.BeginFrame(world, lights, cam, 16 * 16, 0, report, 1);
        expect_true(bool(session));
        expect_true(!report.scene_reused);
        if (session) {
          const auto *initial_session = session.get();
          RayMatrix rgb(16, 16, 3), second(16, 16, 3), normal(16, 16, 3), albedo(16, 16, 3),
              alpha(16, 16, 1), draw(16, 16, 3);
          adaptive_sampler film(1, 16, 16, 32, 0, 0, 1, rgb, second, normal, albedo, alpha, draw,
                                false, false);
          auto render = [&]() {
            context.BeginFrame(world, lights, cam, 16 * 16, 0, report, 1);
            expect_true(session.get() == initial_session);
            expect_true(report.scene_reused);
            expect_true(report.completed_samples == 0);
            expect_true(report.sample_seconds == 0);
            expect_true(report.scattering_events == 0);
            film.reset();
            for (size_t s = 0; s < 16; ++s)
              if (!RenderWavefrontSample(
                      *session, scene, cam, film, s, 0, 791, 3, 3, 1e6, [] { return false; },
                      report))
                return false;
            return true;
          };
          expect_true(render());
          const auto original = rgb.data;
          // These are the same mutable camera/transform objects used by DrawImage.
          cam.update_position_absolute({8, 0, 4});
          cam.update_lookat({8, 0, 0});
          expect_true(render());
          const auto moved = rgb.data;
          expect_true(moved != original);
          cam.update_position_absolute({0, 0, 4});
          cam.update_lookat({0, 0, 0});
          expect_true(render());
          expect_true(rgb.data == original);
          environment = RotateY(120);
          inverse_environment = Inverse(environment);
          expect_true(render());
          expect_true(rgb.data != original);
          environment = Transform();
          inverse_environment = Transform();
          expect_true(render());
          expect_true(rgb.data == original);
          int polls = 0;
          bool completed = RenderWavefrontSample(
              *session, scene, cam, film, 16, 0, 791, 8, 3, 1e6, [&] { return ++polls >= 2; },
              report);
          expect_true(!completed);
          expect_true(rgb.data == original);
          // The same allocated queues can service the low-resolution preview,
          // then return to full resolution without stale per-path state.
          RayMatrix small_rgb(4, 4, 3), small_second(4, 4, 3), small_normal(4, 4, 3),
              small_albedo(4, 4, 3), small_alpha(4, 4, 1), small_draw(4, 4, 3);
          adaptive_sampler small(1, 4, 4, 1, 0, 0, 1, small_rgb, small_second, small_normal,
                                 small_albedo, small_alpha, small_draw, false, false);
          expect_true(RenderWavefrontSample(
              *session, scene, cam, small, 0, 0, 791, 3, 3, 1e6, [] { return false; }, report));
          expect_true(render());
          expect_true(rgb.data == original);
          // Larger queues invalidate the allocation; smaller frames can keep it.
          context.BeginFrame(world, lights, cam, 32 * 32, 0, report, 1);
          expect_true(bool(session));
          expect_true(!report.scene_reused);
          context.BeginFrame(world, lights, cam, 4 * 4, 0, report, 1);
          expect_true(report.scene_reused);
          sky->light = background;
          context.BeginFrame(world, lights, cam, 4 * 4, 0, report, 1);
          expect_true(!report.scene_reused);
          expect_true(bool(session));

          Transform pose = Inverse(LookAt(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0)));
          AnimatedTransform camera_motion(&pose, 0, &pose, 1);
          std::vector<Float> lens = {
            29.475f, 3.76f, 1.67f, 25.2f, 84.83f, .12f, 1.f, 25.2f,
            19.275f, 4.025f, 1.67f, 23.f, 40.77f, 3.275f, 1.699f, 23.f,
            12.75f, 5.705f, 1.f, 18.f, 0.f, 4.5f, 0.f, 17.1f,
            -14.495f, 1.18f, 1.603f, 17.f, 40.77f, 6.065f, 1.658f, 20.f,
            -20.385f, .19f, 1.f, 20.f, 437.065f, 3.22f, 1.717f, 20.f,
            -39.73f, 0.f, 1.f, 20.f};
          RealisticCamera lens_camera(camera_motion, 2, 6, 10, 10, 10, 10,
                                      false, lens, .022, 1, 1, vec3f(0, 1, 0), pose, point3f(0));
          context.BeginFrame(world, lights, lens_camera, 4 * 4, 0, report, 1);
          expect_true(bool(session));
          expect_true(!report.scene_reused);
          expect_true(!scene.camera_data.empty());
          context.BeginFrame(world, lights, lens_camera, 4 * 4, 0, report, 1);
          expect_true(report.scene_reused);
          context.BeginFrame(world, lights, cam, 4 * 4, 0, report, 1);
          expect_true(report.scene_reused);
          expect_true(scene.camera_data.empty());
        }
      }
    }
  }
}
#endif
