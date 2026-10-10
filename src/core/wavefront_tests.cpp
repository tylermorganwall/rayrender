#ifdef NOT_CRAN
#include <testthat.h>
#include "wavefront.h"
#include "hlbvh.h"
#include "camera.h"
#include "adaptivesampler.h"
#include "../hitables/rectangle.h"
#include "../hitables/infinite_area_light.h"

context("Metal wavefront session") {
  test_that("a persistent session follows camera and light changes and cancels atomically") {
    if (MetalBVHAvailable()) {
      Transform identity, environment, inverse_environment;
      hitable_list world, lights;
      auto white = std::make_shared<constant_texture>(point3f(.1));
      auto background = std::make_shared<ImageInfiniteLight>(white, 8, 4, 0);
      auto disk = std::make_shared<DiskInfiniteLight>(point3f(4), vec3f(0, 0, 1), 60, 0);
      auto sources = std::make_shared<InfiniteLightMixture>(
          std::vector<std::shared_ptr<InfiniteLight>>{background, disk});
      auto diffuse =
          std::make_shared<diffuse_material>(std::make_shared<constant_texture>(point3f(.5)));
      world.add(std::make_shared<xy_rect>(-2, 2, -2, 2, 0, diffuse, nullptr, nullptr, &identity,
                                          &identity, false));
      world.add(std::make_shared<InfiniteAreaLight>(sources, 100, point3f(0), &environment,
                                                    &inverse_environment));
      camera cam({0, 0, 4}, {0, 0, 0}, {0, 1, 0}, 30, 1, 0, 4, 0, 1, 1);
      WavefrontScene scene;
      WavefrontReport report;
      auto session = PrepareWavefront(world, lights, cam, 16 * 16, 0, scene, report);
      expect_true(bool(session));
      if (session) {
        RayMatrix rgb(16, 16, 3), second(16, 16, 3), normal(16, 16, 3), albedo(16, 16, 3),
            alpha(16, 16, 1), draw(16, 16, 3);
        adaptive_sampler film(1, 16, 16, 32, 0, 0, 1, rgb, second, normal, albedo, alpha, draw,
                              false, false);
        auto render = [&]() {
          film.reset();
          for (size_t s = 0; s < 16; ++s)
            if (!RenderWavefrontSample(
                    *session, scene, cam, film, s, 0, 791, 3, 3, 1e6, [] { return false; }, report))
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
      }
    }
  }
}
#endif
