#ifdef NOT_CRAN
#include <Rcpp.h>
#include "transform.h"
#include "transformcache.h"
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

context("Transform cache probing and growth") {
  test_that("cached pointers survive collisions and multiple table growths") {
    TransformCache cache;
    std::vector<Transform> transforms;
    std::vector<Transform *> pointers;
    size_t missed_before_growth = 0, missed_after_growth = 0;
    for (int i = 0; i < 2000; ++i) {
      transforms.push_back(Translate(vec3f(i * .125f, (i % 17) * .25f, -i * .375f)));
      pointers.push_back(cache.Lookup(transforms.back()));
      if (i == 199) {
        for (size_t j = 0; j < pointers.size(); ++j)
          missed_before_growth += cache.Lookup(transforms[j]) != pointers[j];
      }
    }
    for (size_t i = 0; i < pointers.size(); ++i) {
      missed_after_growth += cache.Lookup(transforms[i]) != pointers[i];
      expect_true(*pointers[i] == transforms[i]);
    }
    expect_true(missed_before_growth == 0);
    expect_true(missed_after_growth == 0);
    cache.Clear();
    Transform *fresh = cache.Lookup(transforms.back());
    expect_true(*fresh == transforms.back());
    expect_true(cache.Lookup(transforms.back()) == fresh);
    Float signed_zeros[4][4] = {};
    for (int row = 0; row < 4; ++row)
      for (int column = 0; column < 4; ++column)
        signed_zeros[row][column] = row == column ? Float(1) : -Float(0);
    expect_true(cache.Lookup(Transform()) == cache.Lookup(Transform(signed_zeros)));
  }
}
#endif
