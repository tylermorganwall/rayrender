#ifdef NOT_CRAN
#include "keyframe.h"
#include "animatedtransform.h"
#include <testthat.h>

using animation::Value;
using animation::Track;
using animation::Property;
using animation::Interpolation;

context("Absolute-time animation tracks") {
  test_that("frame anchors preserve exact repeat and STEP boundaries after distant seeks") {
    auto track = std::make_shared<animation::TransformTrack>();
    track->position = std::make_shared<Track>(Property::Position, Interpolation::Step,
      std::vector<double>{0, 1, 2}, std::vector<Value>{{-1,0,0,0},{0,0,0,0},{1,0,0,0}});
    track->clock = std::make_shared<animation::ClipClock>(0, 2, 0, 0, 1, animation::Loop::Repeat);
    auto frame = std::make_shared<SceneTimeContext>();
    Transform identity;
    AnimatedTransform adapted(&identity, 0, &identity, 0, track, -.001, frame);
    for (double requested : {4., 1., 0., 2., -.001, 1.}) {
      frame->epoch = requested;
      frame->ray_anchor = Float(requested + .001);
      Transform actual;
      adapted.Interpolate(frame->ray_anchor, &actual);
      expect_true(actual == track->At(requested));
    }
  }

  test_that("clip clocks define boundaries, negative times, reverse and frozen playback") {
    using animation::ClipClock;
    using animation::Loop;
    ClipClock repeat(2, 4, 10, .5, 2, Loop::Repeat);
    expect_true(repeat.At(10) == 2.5);
    expect_true(repeat.At(9.75) == 2);
    expect_true(repeat.At(10.75) == 2);
    expect_true(repeat.At(9.5) == 3.5);
    ClipClock reverse(2, 4, 0, 2, -1, Loop::Hold);
    expect_true(reverse.At(-10) == 4);
    expect_true(reverse.At(.5) == 3.5);
    expect_true(reverse.At(3) == 2);
    ClipClock pingpong(2, 4, 0, 0, 1, Loop::PingPong);
    for (const auto &sample : std::vector<std::pair<double,double>>{{-4,2},{-2,4},{-1,3},{0,2},{1,3},{2,4},{3,3},{4,2}})
      expect_true(pingpong.At(sample.first) == sample.second);
    ClipClock zero(7, 7, 1e12, .5, -2, Loop::Repeat);
    expect_true(zero.At(-1e12) == 7);
    ClipClock paused(2, 4, 1e12, .5, 0, Loop::Repeat);
    expect_true(paused.At(-1e12) == 2.5);
    expect_error(ClipClock(2, 1, 0, 0, 1, Loop::Hold));
    expect_error(ClipClock(0, 1, 0, INFINITY, 1, Loop::Repeat));
    expect_error(repeat.At(INFINITY));
  }

  test_that("loop interval bounds enclose seams, reflections and many cycles") {
    using animation::ClipClock;
    using animation::Loop;
    animation::TransformTrack track;
    track.position = std::make_shared<Track>(Property::Position, Interpolation::Cubic,
      std::vector<double>{2, 2.3, 3.4, 4},
      std::vector<Value>{{0,0,0,0},{3,2,0,0},{-2,1,0,0},{1,0,0,0}});
    aabb box(point3f(-.1,-.2,-.3), point3f(.1,.2,.3));
    for (auto mode : {Loop::Repeat, Loop::PingPong, Loop::Hold}) {
      for (double speed : {-2., 0., .5, 2.}) {
        track.clock = std::make_shared<ClipClock>(2,4,.2,.7,speed,mode);
        for (auto interval : std::vector<std::pair<double,double>>{{-1,.2},{.2,2.7},{-10000,10000},{2,2}}) {
          auto bounds = track.Bounds(box, interval.first, interval.second);
          for (int i = 0; i <= 100; ++i) {
            double t = interval.first + (interval.second - interval.first) * i / 100.;
            auto pose = track.At(t);
            for (int corner = 0; corner < 8; ++corner) {
              auto p = pose(box.Corner(corner));
              for (int k = 0; k < 3; ++k)
                expect_true((p[k] >= bounds.min()[k] && p[k] <= bounds.max()[k]));
            }
          }
        }
      }
    }
    ClipClock repeat(2,4,0,0,1,Loop::Repeat);
    expect_true(repeat.Intervals(1.9,2.1).size() == 2);
    expect_true(repeat.Intervals(-10000,10000).size() == 1);
    ClipClock pingpong(2,4,0,0,1,Loop::PingPong);
    expect_true(pingpong.Intervals(1.9,2.1).size() == 2);
  }

  test_that("nonuniform smooth rotation has continuous angular velocity") {
    std::vector<double> times{0, .7, 2.8, 4};
    std::vector<Value> values;
    for (double angle : {0., .6, 1.4, 2.5}) values.push_back({0, std::sin(angle / 2), 0, std::cos(angle / 2)});
    Track track(Property::Rotation, Interpolation::Smooth, times, values);
    double h = 1e-5;
    for (int i : {1, 2}) {
      auto before = track.At(times[i] - h), at = track.At(times[i]), after = track.At(times[i] + h);
      for (int k = 0; k < 4; ++k)
        expect_true(std::abs((at[k] - before[k]) / h - (after[k] - at[k]) / h) < 1e-4);
    }
    auto expected = track.At(.3);
    for (double time : {20., -2., 2.8, .3, 1.7}) track.At(time);
    expect_true(track.At(.3) == expected);
  }

  test_that("bounds contain cubic extrema and closed rotational excursions") {
    auto motion = std::make_shared<animation::TransformTrack>();
    motion->position = std::make_shared<Track>(Property::Position, Interpolation::Cubic,
      std::vector<double>{0, 1}, std::vector<Value>{{0, 0, 0, 0}, {0, 0, 0, 0}},
      std::vector<Value>{{0, 0, 0, 0}, {-12, 0, 0, 0}}, std::vector<Value>{{12, 0, 0, 0}, {0, 0, 0, 0}});
    motion->rotation = std::make_shared<Track>(Property::Rotation, Interpolation::Linear,
      std::vector<double>{0, .2, .6, 1},
      std::vector<Value>{{0, 0, 0, 1}, {0, 0, 1, 0}, {0, 0, 0, -1}, {0, 0, 0, -1}});
    motion->scale = std::make_shared<Track>(Property::Scale, Interpolation::Cubic,
      std::vector<double>{0, .3, 1},
      std::vector<Value>{{1, -2, .5, 0}, {2, -3, .7, 0}, {1, -2, .5, 0}});
    aabb box(point3f(1, -1, -.1), point3f(3, 1, .1));
    auto range = motion->position->Range(0, 1);
    expect_true(std::abs(range[1][0] - 3) < 1e-12);
    auto bounds = motion->Bounds(box, 0, 1);
    for (int i = 0; i <= 500; ++i) for (int corner = 0; corner < 8; ++corner) {
      auto p = motion->At(i / 500.)(box.Corner(corner));
      for (int k = 0; k < 3; ++k) expect_true((p[k] >= bounds.min()[k] && p[k] <= bounds.max()[k]));
    }
    auto snapshot = motion->Bounds(box, .5, .5), reference = motion->At(.5)(box);
    for (int k = 0; k < 3; ++k) {
      expect_true(snapshot.min()[k] == reference.min()[k]);
      expect_true(snapshot.max()[k] == reference.max()[k]);
    }
  }

  test_that("frame epoch preserves nearby double keys and legacy shear") {
    auto motion = std::make_shared<animation::TransformTrack>();
    const double epoch = 1e12;
    motion->position = std::make_shared<Track>(Property::Position, Interpolation::Linear,
      std::vector<double>{epoch, epoch + .25, epoch + .5, epoch + 1},
      std::vector<Value>{{0, 0, 0, 0}, {1, 0, 0, 0}, {2, 0, 0, 0}, {0, 0, 0, 0}});
    Transform identity;
    AnimatedTransform adapted(&identity, 0, &identity, 0, motion, epoch + .125);
    Transform result;
    adapted.Interpolate(0, &result);
    expect_true(std::abs(result.GetMatrix().m[0][3] - .5) < 1e-6);
    adapted.Interpolate(.25, &result);
    expect_true(std::abs(result.GetMatrix().m[0][3] - 1.5) < 1e-6);

    Float m[4][4] = {{1, .3, 0, 2}, {0, 2, .2, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
    Transform end(m);
    AnimatedTransform legacy(&identity, 0, &end, 2), shifted(&identity, 0, &end, 2, nullptr, .5);
    aabb box(point3f(-1, -2, -3), point3f(1, 2, 3));
    auto old_bounds = legacy.MotionBounds(box), adapted_bounds = shifted.MotionBounds(box);
    for (int k = 0; k < 3; ++k) {
      expect_true(old_bounds.min()[k] == adapted_bounds.min()[k]);
      expect_true(old_bounds.max()[k] == adapted_bounds.max()[k]);
    }
    for (Float time : {-.5f, 0.f, .2f, .8f, 1.5f, 4.f}) {
      Transform a, b;
      legacy.Interpolate(time + .5f, &a); shifted.Interpolate(time, &b);
      expect_true(a == b);
    }
  }
}
#endif
