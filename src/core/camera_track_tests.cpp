#ifdef NOT_CRAN
#include "camera_track.h"
#include <testthat.h>

context("Scene camera tracks") {
  test_that("nonuniform camera rows sample without seek history") {
    animation::CameraKey first{2, {0,0,-10, 0,0,0, 0,60,10, 4,4, 0,1,0}};
    auto middle = first, last = first;
    middle.time = 3; middle.values[0] = 4;
    last.time = 9; last.values[0] = -2; last.values[7] = 30;
    animation::CameraTrack track({first,middle,last});
    camera cam(point3f(0,0,-10),point3f(0),vec3f(0,1,0),60,1,0,10,0,0,1);
    for (double time : {6.,9.,2.,6.,3.,-2.,15.,6.}) {
      auto key = track.At(time);
      track.Apply(key,cam);
      double expected = time <= 2 ? 0 : time <= 3 ? 4*(time-2) : time < 9 ? 4-(time-3) : -2;
      expect_true(std::abs(cam.get_origin()[0] - expected) < 1e-6);
      if (time == 6) {
        expect_true(cam.get_fov() == 45);
        auto ray = cam.get_ray(.5,.5,point3f(0),0);
        expect_true((unit_vector(ray.d)-unit_vector(point3f(0)-cam.get_origin())).length() < 1e-6);
      }
    }
    auto invalid = first;
    invalid.values[3] = invalid.values[0]; invalid.values[4] = invalid.values[1];
    invalid.values[5] = invalid.values[2];
    point3f previous = cam.get_origin();
    expect_error(track.Apply(invalid,cam));
    expect_true((cam.get_origin()-previous).length() == 0);
    expect_error(animation::CameraTrack({first,first}));
    auto overflow = first;
    overflow.values[0] = 1e300;
    expect_error(animation::CameraTrack({overflow}));
    middle.values[7] = 0;
    expect_error(animation::CameraTrack({first,middle}));
  }
}
#endif
