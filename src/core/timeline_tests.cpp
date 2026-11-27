#ifdef NOT_CRAN
#include "timeline.h"
#include <cmath>
#include <testthat.h>

context("Scene timeline controller") {
  test_that("playback uses elapsed wall seconds and reanchors on rate, pause and seek") {
    double wall = 100;
    animation::TimelineController timeline(0, 20, 2, 24, [&] { return wall; });
    timeline.Pause(false);
    wall += 3;
    expect_true(timeline.Request().time == 5);
    timeline.SetRate(-2);
    wall += 1;
    expect_true(timeline.Request().time == 3);
    timeline.Seek(9);
    expect_true(timeline.Request().time == 9);
    wall += .5;
    timeline.Pause(true);
    expect_true(timeline.Request().time == 8);
    wall += 100;
    expect_true(timeline.Request().time == 8);
    timeline.Command(animation::TimelineCommand::StepBack);
    expect_true(std::abs(timeline.Request().time - (8 - 1./24)) < 1e-12);
    timeline.Command(animation::TimelineCommand::End);
    expect_true(timeline.Request().time == 20);
    timeline.SetRate(0);
    timeline.Pause(false);
    wall += 10;
    expect_true(timeline.Request().time == 20);
  }
  test_that("seeks coalesce and stale prepared or displayed revisions cannot publish") {
    double wall = 0;
    animation::TimelineController timeline(0, 10, 2, 30, [&] { return wall; });
    auto first = timeline.Request();
    expect_true(timeline.Publish(first));
    expect_true(timeline.Display(first.revision));
    timeline.Seek(9);
    auto stale = timeline.Request();
    timeline.Seek(3);
    timeline.Seek(7);
    auto latest = timeline.Request();
    expect_false(timeline.Publish(stale));
    expect_true(timeline.Displayed().time == 2);
    expect_true(timeline.Publish(latest));
    expect_false(timeline.Display(first.revision));
    expect_true(timeline.Displayed().time == 2);
    expect_true(timeline.Display(latest.revision));
    expect_true(timeline.Displayed().time == 7);
    expect_false(timeline.Pending());
    timeline.Drag(4);
    auto drag = timeline.Request();
    expect_true(drag.dragging);
    expect_true(timeline.Paused());
    timeline.EndDrag();
    expect_false(timeline.Request().dragging);
    expect_false(timeline.Publish(drag));
    timeline.Reject(timeline.Request().revision);
    expect_false(timeline.Pending());
    expect_true(timeline.Published().time == 7);
    timeline.Seek(5);
    expect_true(timeline.Pending());
  }
  test_that("loop policy, zero ranges and dirty dependencies are explicit") {
    double wall = 0;
    animation::TimelineController timeline(-2, 2, 1, 30, [&] { return wall; });
    timeline.SetLoop(true);
    timeline.Pause(false);
    wall += 3;
    expect_true(timeline.Request().time == 0);
    timeline.SetRate(-1);
    wall += 3;
    expect_true(timeline.Request().time == 1);
    timeline.Invalidate(animation::Pose | animation::Shading);
    expect_true((timeline.Request().dirty & animation::Pose) != 0);
    expect_true((timeline.Request().dirty & animation::Shading) != 0);
    animation::TimelineController still(4, 4, 4, 30, [&] { return wall; });
    still.Pause(false);
    still.Drag(100);
    wall += 10;
    expect_true(still.Paused());
    expect_false(still.Request().dragging);
    expect_true(still.Request().time == 4);
    expect_error(still.Seek(INFINITY));
    expect_error(timeline.SetRate(NAN));
    expect_error(animation::TimelineController(2, 1, 1, 30));
  }
}
#endif
