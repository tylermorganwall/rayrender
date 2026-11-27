#ifdef NOT_CRAN
#include "texture.h"
#include <testthat.h>

context("Alpha texture mapping") {
  test_that("coverage has independent repeat and offset without gamma conversion") {
    unsigned char pixels[16] = {255,255,255,0, 255,255,255,64,
                               255,255,255,128, 255,255,255,255};
    alpha_texture original(pixels, 2, 2, 4);
    expect_true(original.value(.25f, .75f, point3f(0)) == 0);
    expect_true(std::abs(original.value(.25f, .25f, point3f(0)) - 128.f/255) < 1e-6f);
    alpha_texture mapped(pixels, 2, 2, 4, .25f, -.1f, 2, 2);
    expect_true(mapped.value(.2f, .25f, point3f(0)) == 1);
    expect_true(mapped.value(-.3f, .75f, point3f(0)) == 1);
    mapped.opacity = .3f;
    expect_true(std::abs(mapped.value(.2f, .25f, point3f(0)) - .3f) < 1e-6f);
    alpha_texture constant(.37f);
    expect_true(constant.value(.8f, -.1f, point3f(0)) == .37f);
  }
}
#endif
