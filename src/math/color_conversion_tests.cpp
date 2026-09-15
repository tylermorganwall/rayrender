#ifdef NOT_CRAN
#include <Rcpp.h>
#include "mathinline.h"
#include <testthat.h>

context("RGB and HSV conversion") {
  test_that("primary and secondary colors retain hue and brightness") {
    point3f colors[] = {
        point3f(1, 0, 0), point3f(1, 1, 0), point3f(0, 1, 0),
        point3f(0, 1, 1), point3f(0, 0, 1), point3f(1, 0, 1)};
    for (int i = 0; i < 6; ++i) {
      point3f hsv = RGBtoHSV(colors[i]);
      expect_true(std::abs(hsv[0] - 60 * i) < 1e-5);
      expect_true(std::abs(hsv[1] - 1) < 1e-6);
      expect_true(std::abs(hsv[2] - 1) < 1e-6);
      point3f restored = HSVtoRGB(hsv);
      for (int channel = 0; channel < 3; ++channel)
        expect_true(std::abs(restored[channel] - colors[i][channel]) < 1e-6);
    }
  }

  test_that("achromatic and blue-dominant HDR colors round trip") {
    for (Float value : {Float(0), Float(0.25), Float(1), Float(4)}) {
      point3f gray(value);
      point3f hsv = RGBtoHSV(gray);
      expect_true(hsv[0] == 0);
      expect_true(hsv[1] == 0);
      expect_true(hsv[2] == value);
    }
    point3f blue(1, 2, 4);
    point3f hsv = RGBtoHSV(blue);
    expect_true(std::abs(hsv[0] - 220) < 1e-5);
    expect_true(std::abs(hsv[1] - 0.75) < 1e-6);
    expect_true(hsv[2] == 4);
    point3f restored = HSVtoRGB(hsv);
    for (int channel = 0; channel < 3; ++channel)
      expect_true(std::abs(restored[channel] - blue[channel]) < 1e-6);
  }
}
#endif
