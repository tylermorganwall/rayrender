#ifdef NOT_CRAN
#include <Rcpp.h>
#include "transform.h"
#include <cmath>
#include <testthat.h>

// Keep these calls outside transform.cpp: a definition available only through
// inlining in that file cannot satisfy geometry's references to these overloads.
context("Transform error bounds across translation units") {
  test_that("point and vector overloads preserve positions and propagate errors") {
    const Transform transform = Translate(vec3f(1, 2, 3)) * Scale(2, -3, 4);
    const point3f point(1, -2, .5);
    const vec3f vector(1, -2, .5), input_error(.01, .02, .03);
    vec3f point_error, propagated_point_error, vector_error, propagated_vector_error;
    const auto mapped_point = transform(point, &point_error);
    const auto mapped_point_with_error = transform(point, input_error, &propagated_point_error);
    const auto mapped_vector = transform(vector, &vector_error);
    const auto mapped_vector_with_error = transform(vector, input_error, &propagated_vector_error);
    const point3f expected_point(3, 8, 5);
    const vec3f expected_vector(2, 6, 2), minimum_error(.02, .06, .12);

    for (int axis = 0; axis < 3; ++axis) {
      expect_true(std::abs(mapped_point[axis] - expected_point[axis]) < 1e-6);
      expect_true(std::abs(mapped_point_with_error[axis] - expected_point[axis]) < 1e-6);
      expect_true(std::abs(mapped_vector[axis] - expected_vector[axis]) < 1e-6);
      expect_true(std::abs(mapped_vector_with_error[axis] - expected_vector[axis]) < 1e-6);
      expect_true(std::isfinite(point_error[axis]));
      expect_true(point_error[axis] >= 0);
      expect_true(std::isfinite(vector_error[axis]));
      expect_true(vector_error[axis] >= 0);
      expect_true(propagated_point_error[axis] >= minimum_error[axis]);
      expect_true(propagated_vector_error[axis] >= minimum_error[axis]);
    }
  }
}
#endif
