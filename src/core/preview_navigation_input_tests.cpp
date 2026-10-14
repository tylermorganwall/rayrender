#ifdef NOT_CRAN
#include "preview_navigation_input.h"
#include <testthat.h>
#include <limits>

context("Preview navigation input transitions") {
  test_that("separate taps never inherit the preceding render duration") {
    const char* keys = "wasdqz";
    for(int key = 0; key < 6; ++key) {
      for(bool shifted : {false, true}) {
        PreviewNavigationInput input;
        for(double render_seconds : {.016, .2, 2.0, .05, .7}) {
          input.Press(keys[key], shifted);
          auto steps = input.Consume(1u << key, shifted, render_seconds);
          expect_true(steps[shifted][key] == Approx(1.0 / 60));
          expect_true(steps[!shifted][key] == 0);
          // Release and re-press both arrive before the next render poll.
          input.Release(keys[key]);
        }
      }
    }
  }

  test_that("quick taps survive as distinct presses with their own modifiers") {
    PreviewNavigationInput input;
    input.Press('w', false); input.Release('w');
    input.Press('w', false); input.Release('w');
    input.Press('w', true); input.Release('w');
    auto steps = input.Consume(0, false, 5);
    expect_true(steps[0][0] == Approx(2.0 / 60));
    expect_true(steps[1][0] == Approx(1.0 / 60));
    expect_true(input.Consume(0, false, 5) == PreviewNavigationInput::Steps{});
  }

  test_that("only an uninterrupted press accrues held movement") {
    PreviewNavigationInput input;
    input.Press('w', false);
    input.Consume(1, false, 10);
    input.Press('w', false); // Native auto-repeat must not enqueue a tap.
    expect_true(input.Consume(1, false, .1)[0][0] == Approx(.1));
    expect_true(input.Consume(1, false, 5)[0][0] == Approx(.25));
    input.Release('w');
    expect_true(input.Consume(0, false, .2) == PreviewNavigationInput::Steps{});
    input.Press('w', false);
    expect_true(input.Consume(1, false, .2)[0][0] == Approx(1.0 / 60));
  }

  test_that("losing focus or a release cannot leave a held key moving") {
    PreviewNavigationInput input;
    input.Press('w', false);
    input.Consume(1, false, .1);
    input.Reset();
    expect_true(input.Consume(1, false, 10) == PreviewNavigationInput::Steps{});
    input.Press('w', false);
    input.Consume(1, false, .1);
    expect_true(input.Consume(0, false, .1) == PreviewNavigationInput::Steps{});
    expect_true(input.Consume(1, false, .1) == PreviewNavigationInput::Steps{});
    input.Press('w', false);
    expect_true(input.Consume(1, false, .1)[0][0] == Approx(1.0 / 60));
  }

  test_that("a new key does not erase another key's continuous hold") {
    PreviewNavigationInput input;
    input.Press('w', false);
    input.Consume(1, false, .1);
    input.Press('d', false);
    auto steps = input.Consume(1 | (1u << 3), false, .2);
    expect_true(steps[0][0] == Approx(.2));
    expect_true(steps[0][3] == Approx(1.0 / 60));
    expect_true(input.Consume(1, false, -1) == PreviewNavigationInput::Steps{});
    expect_true(input.Consume(1, false, std::numeric_limits<double>::infinity()) ==
                PreviewNavigationInput::Steps{});
    input.Press('\0', false); input.Press('!', true); input.Release('!');
    expect_true(input.Consume(0, false, .1) == PreviewNavigationInput::Steps{});
  }
}
#endif
