#ifdef NOT_CRAN
#include "adaptivesampler.h"
#include <testthat.h>

context("Adaptive denoising guides") {
  test_that("early convergence averages guides once and preserves active sums") {
    RayMatrix rgb(2, 2, 3), second(2, 2, 3), normal(2, 2, 3), albedo(2, 2, 3),
        alpha(2, 2, 1), draw(2, 2, 3), preview_normal(2, 2, 3), preview_albedo(2, 2, 3);
    adaptive_sampler film(2, 2, 2, 8, 0, .001f, 1,
                          rgb, second, normal, albedo, alpha, draw, true);
    // Two sample normals average to a non-unit vector. Do not renormalize its
    // length: OIDN needs the coverage-weighted average, not a unit direction.
    for (size_t x = 0; x < 2; ++x) {
      for (size_t y = 0; y < 2; ++y) {
        film.add_normal(x, y, normal3f(1, 0, 0));
        film.add_normal(x, y, normal3f(0, 1, 0));
        film.add_albedo(x, y, point3f(.2f, .4f, .6f));
        film.add_albedo(x, y, point3f(.6f, .4f, .2f));
        film.add_color_main(x, y, point3f(2, 4, 6));
        film.add_color_main(x, y, point3f(4, 6, 8));
      }
    }
    film.pixel_chunks[0].erase = true;
    film.split_remove_chunks(1); // Zero-based sample index: two samples, not ns=8.
    film.max_s = 2;
    film.copy_denoising_features(preview_normal, preview_albedo);
    expect_true(film.finalized[0]);
    expect_true(film.size() == 3);
    expect_true(std::abs(rgb(0, 0, 0) - 3) < 1e-6);
    expect_true(std::abs(normal(0, 0, 0) - .5) < 1e-6);
    expect_true(std::abs(albedo(0, 0, 0) - .4) < 1e-6);
    expect_true(std::abs(normal(1, 1, 0) - 1) < 1e-6);
    for (size_t x = 0; x < 2; ++x) {
      for (size_t y = 0; y < 2; ++y) {
        expect_true(std::abs(preview_normal(x, y, 0) - .5) < 1e-6);
        expect_true(std::abs(preview_normal(x, y, 1) - .5) < 1e-6);
        expect_true(std::abs(preview_albedo(x, y, 0) - .4) < 1e-6);
      }
    }
    // The remaining pixels continue sampling. Final output must average only
    // those pixels, leaving the already-converged guides untouched.
    for (const auto& block : film.pixel_chunks) {
      for (size_t x = block.startx; x < block.endx; ++x) {
        for (size_t y = block.starty; y < block.endy; ++y) {
          for (int sample = 2; sample < 8; ++sample) {
            film.add_normal(x, y, normal3f(.5f, .5f, 0));
            film.add_albedo(x, y, point3f(.4f));
          }
        }
      }
    }
    film.write_final_pixels();
    for (size_t x = 0; x < 2; ++x) {
      for (size_t y = 0; y < 2; ++y) {
        expect_true(std::abs(normal(x, y, 0) - .5) < 1e-6);
        expect_true(std::abs(albedo(x, y, 0) - .4) < 1e-6);
      }
    }
    film.reset();
    expect_true(film.size() == 4);
    expect_false(film.finalized[0]);
    expect_true(normal(0, 0, 0) == 0);
    expect_true(albedo(0, 0, 0) == 0);
  }
}
#endif
