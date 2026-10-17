#ifdef NOT_CRAN
#include "camera_track.h"
#include <testthat.h>

context("Scene camera tracks") {
  test_that("physical tracks validate constant optics before rendering") {
    animation::CameraKey a{2,{0,0,-10, 0,0,0, 10,-1,10, 4,4, 0,1,0}};
    auto b=a; b.time=9; b.values[0]=3;
    animation::CameraTrack track({a,b});
    expect_true(track.At(5.5).values[0] == 1.5);
    auto invalid=b; invalid.values[6]=5;
    expect_error(animation::CameraTrack({a,invalid}));
    invalid=b; invalid.values[8]=5;
    expect_error(animation::CameraTrack({a,invalid}));
    invalid=b; invalid.values[7]=60;
    expect_error(animation::CameraTrack({a,invalid}));
    expect_error(animation::CameraTrack({invalid,a}));
    expect_error(track.MakeCamera(a,1,1,0,0,false));
    expect_error(track.SampleRay(2,0,.5,.5,point3f(0),1,1,false));
  }
  test_that("exposures visit interior camera keys and share scene time with transport") {
    animation::CameraKey a{2,{0,0,-10, 0,0,0, 0,60,10, 4,4, 0,1,0}};
    auto b=a, c=a;
    b.time=2.2; b.values[0]=4; b.values[7]=30;
    c.time=2.4;
    auto track=std::make_shared<const animation::CameraTrack>(std::vector<animation::CameraKey>{a,b,c});
    auto clock=std::make_shared<SceneTimeContext>();
    clock->epoch=2; clock->ray_anchor=7;
    camera cam(point3f(0,0,-10),point3f(0),vec3f(0,1,0),60,1,0,10,7,7.4,1);
    cam.set_shutter_speed(1);
    cam.scene_camera_track=track; cam.scene_camera_time=clock;
    for (Float fraction : {0.f,.5f,1.f,.25f,.5f}) {
      auto ray=cam.get_ray(.5,.5,point3f(0),fraction);
      Float time=7+(7.4f-7)*fraction;
      double absolute=2+(double(time)-7);
      auto expected=track->SampleRay(absolute,time,.5,.5,point3f(0),1,1,false);
      expect_true((ray.o-expected.o).length()==0);
      expect_true((ray.d-expected.d).length()==0);
      expect_true(ray.time()==time);
      if(fraction==.5f) expect_true(std::abs(ray.o[0]-4)<1e-4);
    }
    auto invalid=b; invalid.values[11]=0; invalid.values[12]=-1; invalid.values[13]=0;
    expect_error(animation::CameraTrack({a,invalid}));
  }
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
    animation::CameraTrack cut({first,middle});
    expect_true(cut.At(2.9).values == first.values);
    expect_true(cut.At(3).values == middle.values);
    expect_error(cut.Apply(middle,cam));
    auto changed = cut.MakeCamera(middle,1,1,0,1,false);
    expect_true(changed->get_fov() == 0);
  }
}
#endif
