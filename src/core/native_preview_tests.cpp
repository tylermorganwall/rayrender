#ifdef NOT_CRAN
#include "PreviewDisplay.h"
#include "preview_sky.h"
#include "preview_sky_controls.h"
#include "../materials/texturecache.h"
#include "integrator.h"
#include "../hitables/infinite_area_light.h"
#include <testthat.h>

namespace {
std::unique_ptr<PreviewDisplay> NativeTestDisplay(RayCamera& cam, Transform& object,
                                                  Transform& world) {
#ifdef HAS_OIDN
  return std::unique_ptr<PreviewDisplay>(new PreviewDisplay(4,
                                                            4,
                                                            false,
                                                            true,
                                                            false,
                                                            10,
                                                            &cam,
                                                            &object,
                                                            &world,
                                                            nullptr,
                                                            nullptr,
                                                            nullptr,
                                                            false,
                                                            false));
#else
  return std::unique_ptr<PreviewDisplay>(
      new PreviewDisplay(4, 4, false, true, false, 10, &cam, &object, &world, false));
#endif
}
}
context("Native preview controls") {
  test_that("camera and quality restarts preserve sample counts and exposure") {
    // A constant HDR environment removes sampling noise and view-dependent light
    // changes. Exercise the real render loop, including its increment after a
    // completed frame commits camera input and clears the accumulation buffers.
    for (bool deferred : {false, true}) {
      for (bool automatic : {false, true}) {
        const size_t size = 16;
        Transform object, world_transform;
        camera cam(
            point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
        RayMatrix rgb(size, size, 3), normal(size, size, 3), albedo(size, size, 3);
        RayMatrix alpha(size, size, 1), filtered(size, size, 3);
#ifdef HAS_OIDN
        PreviewDisplay display(size,
                               size,
                               false,
                               true,
                               deferred,
                               10,
                               &cam,
                               &object,
                               &world_transform,
                               nullptr,
                               nullptr,
                               nullptr,
                               false,
                               automatic);
#else
        PreviewDisplay display(size,
                               size,
                               false,
                               true,
                               deferred,
                               10,
                               &cam,
                               &object,
                               &world_transform,
                               automatic);
#endif
        display.preview_exposure_adjustment = .1;
        RayrenderGui gui;
        gui.width = gui.height = size;
        display.AttachNativeGui(&gui, true, deferred);
        struct Frames {
          RayrenderGui* gui;
          PreviewDisplay* display;
          uint64_t version = 0;
          std::vector<size_t> samples;
          std::vector<Float> exposure, blue;
        } frames{&gui, &display};
        rayimgui_api_v1 api{};
        api.texture_create = [](uint64_t,
                                const rayimgui_image_v1*,
                                uint64_t* handle,
                                rayimgui_error_v1*) -> int32_t {
          *handle = 1;
          return 0;
        };
        api.texture_update = [](uint64_t,
                                uint64_t,
                                const rayimgui_image_v1*,
                                rayimgui_error_v1*) -> int32_t {
          return 0;
        };
        api.step =
            [](uint64_t owner, rayimgui_step_v1* step, rayimgui_error_v1*) -> int32_t {
          auto& state = *reinterpret_cast<Frames*>(uintptr_t(owner));
          auto& gui = *state.gui;
          if (gui.version == state.version) {
            return 0;
          }
          state.version = gui.version;
          state.samples.push_back(gui.samples);
          state.exposure.push_back(state.display->preview_exposure_scale);
          state.blue.push_back(gui.display_rgb[2]);
          if (state.version == 3 || state.version == 9) {
            gui.keys.push_back({RAYIMGUI_KEY_W, 0, 1});
          } else if (state.version == 6 || state.version == 12) {
            gui.fast_preview = state.version == 6;
            gui.fast_pending = true;
          } else if (state.version == 15) {
            if (state.display->deferred_render) {
              gui.request_render = true;
            } else {
              gui.request_reset = true;
            }
          } else if (state.version >= 18) {
            step->close_requested = 1;
          }
          return 0;
        };
        gui.api = &api;
        gui.session = uint64_t(reinterpret_cast<uintptr_t>(&frames));
        auto light = std::make_shared<ImageInfiniteLight>(
            std::make_shared<constant_texture>(point3f(4, 2, 1)), 4, 2, 0);
        hitable_list world, lights;
        world.add(std::make_shared<InfiniteAreaLight>(
            light, 1000, point3f(0), &object, &world_transform));
        random_gen rng(123);
        pathtracer(2,
                   size,
                   size,
                   32,
                   0,
                   0,
                   1,
                   rgb,
                   normal,
                   albedo,
                   alpha,
                   filtered,
                   false,
                   0,
                   1,
                   1,
                   false,
                   &cam,
                   60,
                   world,
                   lights,
                   10,
                   2,
                   2,
                   display,
                   IntegratorType::ShadowRays,
                   &rng);
        expect_true(frames.samples.size() == 18);
        const double expected_exposure =
            automatic ? 1 / (.2126 * 4 + .7152 * 2 + .0722) : 1;
        PreviewColorTransform mapping;
        const double expected_blue = mapping.Apply({0, 0, .1 * expected_exposure})[2];
        bool correct_counts = true, stable_exposure = true, stable_brightness = true;
        for (size_t frame = 0; frame < frames.samples.size(); ++frame) {
          correct_counts = correct_counts && frames.samples[frame] == frame % 3 + 1;
          stable_exposure = stable_exposure &&
                            std::abs(frames.exposure[frame] - expected_exposure) < 1e-6;
          stable_brightness =
              stable_brightness && std::abs(frames.blue[frame] - expected_blue) < 1e-6;
        }
        expect_true(correct_counts);
        expect_true(stable_exposure);
        expect_true(stable_brightness);
      }
    }
  }

  test_that("history coalesces input drags and branches after undo") {
    Transform object, world;
    camera cam(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
    auto display = NativeTestDisplay(cam, object, world);
    RayrenderGui gui;
    display->AttachNativeGui(&gui, true, false);
    display->BeginNativeHistory();
    auto w = RayrenderGui::widget(RAYIMGUI_DOUBLE, 3, "Exposure");
    for (int i = 2; i <= 4; ++i) {
      rayimgui_item_v1 item{sizeof(item),
                            uint32_t(RAYIMGUI_CHANGED | (i == 2 ? RAYIMGUI_BEGIN : 0))};
      gui.track_input(w, item);
      gui.exposure = i;
      gui.exposure_pending = true;
      display->CommitNativeEdits(nullptr);
    }
    expect_true(gui.can_undo);
    expect_false(gui.can_redo);
    gui.history_requests.push_back(-1);
    display->CommitNativeEdits(nullptr);
    expect_true(display->preview_exposure_adjustment == 1);
    expect_false(gui.can_undo);
    expect_true(gui.can_redo);
    gui.history_requests.push_back(1);
    display->CommitNativeEdits(nullptr);
    expect_true(display->preview_exposure_adjustment == 4);
    gui.history_requests.push_back(-1);
    display->CommitNativeEdits(nullptr);
    ++gui.history_epoch;
    gui.exposure = 7;
    gui.exposure_pending = true;
    display->CommitNativeEdits(nullptr);
    expect_false(gui.can_redo);
    gui.history_requests.push_back(-1);
    display->CommitNativeEdits(nullptr);
    expect_true(display->preview_exposure_adjustment == 1);
    expect_false(gui.can_undo);
  }

  test_that("history bounds memory and does not record unchanged inputs") {
    Transform object, world;
    camera cam(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
    auto display = NativeTestDisplay(cam, object, world);
    RayrenderGui gui;
    display->AttachNativeGui(&gui, true, false);
    display->BeginNativeHistory();
    gui.exposure_pending = true;
    display->CommitNativeEdits(nullptr);
    expect_false(gui.can_undo);
    for (int i = 2; i <= 102; ++i) {
      ++gui.history_epoch;
      gui.exposure = i;
      gui.exposure_pending = true;
      display->CommitNativeEdits(nullptr);
    }
    gui.history_requests.assign(100, -1);
    display->CommitNativeEdits(nullptr);
    expect_true(display->preview_exposure_adjustment == 2);
    expect_false(gui.can_undo);
    expect_true(gui.can_redo);
  }

  test_that("history restores camera, preview controls, and keyframes with snapshots") {
    Transform object, world;
    camera cam(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
    auto display = NativeTestDisplay(cam, object, world);
    RayrenderGui gui;
    display->AttachNativeGui(&gui, true, false);
    display->BeginNativeHistory();
    auto commit = [&] {
      return display->CommitNativeEdits(nullptr);
    };
    gui.keys.push_back({RAYIMGUI_KEY_W, 0, 1});
    ++gui.history_epoch;
    commit();
    expect_true(cam.get_origin()[2] > -10);
    gui.history_requests.push_back(-1);
    commit();
    expect_true(cam.get_origin()[2] == -10);
    gui.history_requests.push_back(1);
    commit();
    const auto moved = cam.get_origin();
    ++gui.history_epoch;
    gui.animation_closed = gui.animation_blur = 1;
    gui.animation_shutter = .25;
    gui.animation_settings_pending = true;
    commit();
    gui.history_requests.push_back(-1);
    commit();
    expect_false(display->CameraMotionBlurEnabled());
    expect_false(display->KeyframeMotionClosed());
    expect_true(display->GetShutterSpeed() == 2);
    gui.history_requests.push_back(1);
    commit();
    expect_true(display->CameraMotionBlurEnabled());
    expect_true(display->KeyframeMotionClosed());
    expect_true(display->GetShutterSpeed() == 4);
    expect_true((cam.get_origin() - moved).length() == 0);
    gui.width = gui.height = 4;
    gui.pixels.assign(4 * 4 * 4, 123);
    auto state = display->CreateCurrentKeyframe(0);
    for (R_xlen_t i = 0; i < state.size(); ++i) {
      gui.snapshot_camera.push_back(Rcpp::as<double>(state[i]));
    }
    ++gui.history_epoch;
    gui.animation_action = RayrenderGui::AnimationAction::Save;
    commit();
    const auto id = gui.keyframe_snapshots[0].id;
    expect_true(gui.keyframe_snapshots[0].pixels[0] == 123);
    ++gui.history_epoch;
    gui.animation_action = RayrenderGui::AnimationAction::Delete;
    commit();
    expect_true(display->Keyframes.empty());
    gui.history_requests.push_back(-1);
    commit();
    expect_true(display->Keyframes.size() == 1);
    expect_true(gui.keyframe_snapshots[0].id == id);
    expect_true(gui.keyframe_snapshots[0].pixels[0] == 123);
    gui.history_requests.push_back(1);
    commit();
    expect_true(gui.keyframe_snapshots.empty());
  }

  test_that("playback poses do not enter history but input edits during playback do") {
    Transform object, world;
    camera cam(point3f(0,0,-10),point3f(0),vec3f(0,1,0),60,1,0,5,0,1,1);
    auto display = NativeTestDisplay(cam, object, world);
    RayrenderGui gui;
    display->AttachNativeGui(&gui,true,false);
    display->SetKeyframeMotionArgs(Rcpp::List::create(Rcpp::_["frames"]=20));
    display->BeginNativeHistory();
    auto action = [&](RayrenderGui::AnimationAction action) {
      ++gui.history_epoch;
      gui.animation_action = action;
      display->CommitNativeEdits(nullptr);
    };
    action(RayrenderGui::AnimationAction::Save);
    gui.keys.push_back({RAYIMGUI_KEY_W, 0, 1});
    ++gui.history_epoch;
    display->CommitNativeEdits(nullptr);
    action(RayrenderGui::AnimationAction::Save);
    const auto original = cam.get_origin();
    action(RayrenderGui::AnimationAction::Play);
    for (int frame = 0; frame < 3; ++frame) {
      display->CommitNativeEdits(nullptr);
    }
    ++gui.history_epoch;
    gui.orbit = 0;
    gui.history_dirty = true;
    display->CommitNativeEdits(nullptr);
    expect_true(display->IsPreviewMotionActive());
    gui.history_requests.push_back(-1);
    display->CommitNativeEdits(nullptr);
    expect_false(display->IsPreviewMotionActive());
    expect_true(gui.orbit == 1);
    expect_true((cam.get_origin() - original).length() == 0);
    gui.history_requests.push_back(-1);
    display->CommitNativeEdits(nullptr);
    expect_true(display->Keyframes.size() == 1);
  }

  test_that("failed sky undo preserves the live camera and history position") {
    Transform object, world;
    camera cam(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
    auto display = NativeTestDisplay(cam, object, world);
    RayrenderGui gui;
    display->AttachNativeGui(&gui, true, false);
    Rcpp::List sky = Rcpp::List::create(Rcpp::_["model"] = 0,
                                        Rcpp::_["elevation"] = 10.0,
                                        Rcpp::_["azimuth"] = 90.0,
                                        Rcpp::_["haze"] = false,
                                        Rcpp::_["altitude"] = false);
    display->export_sky = [&] {
      return Rcpp::clone(sky);
    };
    display->SetSunControls(10, 90, [&](double e, double a) {
      sky["elevation"] = e;
      sky["azimuth"] = a;
    });
    bool fail = true;
    display->prepare_sky_restore =
        [&](const Rcpp::List& saved) -> std::function<void()> {
      if (fail) {
        throw std::runtime_error("Sky file unavailable");
      }
      auto next = Rcpp::clone(saved);
      return [&, next] {
        sky = next;
      };
    };
    display->BeginNativeHistory();
    ++gui.history_epoch;
    gui.sun_elevation = 40;
    gui.sun_pending = true;
    gui.keys.push_back({RAYIMGUI_KEY_W, 0, 1});
    display->CommitNativeEdits(nullptr);
    const auto moved = cam.get_origin();
    gui.history_requests.push_back(-1);
    display->CommitNativeEdits(nullptr);
    expect_true((cam.get_origin() - moved).length() == 0);
    expect_true(Rcpp::as<double>(sky["elevation"]) == 40);
    expect_true(gui.can_undo);
    expect_false(gui.can_redo);
    expect_true(gui.history_message == "Undo failed: Sky file unavailable");
    fail = false;
    gui.history_requests.push_back(-1);
    display->CommitNativeEdits(nullptr);
    expect_true(cam.get_origin()[2] == -10);
    expect_true(Rcpp::as<double>(sky["elevation"]) == 10);
    expect_false(gui.can_undo);
    expect_true(gui.can_redo);
  }

  test_that("animation buttons and keyboard share keyframes and settings") {
    Transform object, world;
    camera cam(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
    auto display = NativeTestDisplay(cam, object, world);
    RayrenderGui gui;
    display->AttachNativeGui(&gui, true, false);
    using Action = RayrenderGui::AnimationAction;
    auto action = [&](Action value) {
      gui.animation_action = value;
      return display->ApplyNativeControls(nullptr);
    };
    auto key = [&](unsigned value, unsigned modifiers = 0) {
      gui.keys.push_back({value, modifiers, 1});
      return display->ApplyNativeControls(nullptr);
    };
    expect_false(action(Action::Save));
    expect_true((display->Keyframes.size() == 1 && gui.keyframe_snapshots.size() == 1));
    expect_true(gui.animation_current == 0);
    key(RAYIMGUI_KEY_W);
    const auto second = cam.get_origin();
    key(RAYIMGUI_KEY_K);
    expect_true((display->Keyframes.size() == 2 && gui.keyframe_snapshots.size() == 2));
    expect_true(gui.animation_current == 1);
    expect_true(action(Action::Previous));
    expect_true(cam.get_origin()[2] == -10);
    expect_true(gui.animation_current == 0);
    expect_true(action(Action::Next));
    expect_true((cam.get_origin() - second).length() == 0);
    gui.animation_target = gui.keyframe_snapshots[0].id;
    expect_true(action(Action::Select));
    expect_true((gui.animation_current == 0 && cam.get_origin()[2] == -10));

    gui.animation_blur = 1;
    gui.animation_closed = 1;
    gui.animation_shutter = .25;
    gui.animation_settings_pending = true;
    expect_true(display->ApplyNativeControls(nullptr));
    expect_true(display->CameraMotionBlurEnabled());
    expect_true(display->KeyframeMotionClosed());
    expect_true(display->GetShutterSpeed() == 4);
    key(RAYIMGUI_KEY_B);
    key(RAYIMGUI_KEY_L, RAYIMGUI_SHIFT);
    expect_false(gui.animation_blur);
    expect_false(gui.animation_closed);
    gui.animation_shutter = 0;
    gui.animation_settings_pending = true;
    display->ApplyNativeControls(nullptr);
    expect_true(std::isinf(display->GetShutterSpeed()));
    expect_true(gui.animation_shutter == 0);

    // Navigation wraps consistently with the existing keyboard controls.
    action(Action::Previous);
    expect_true(gui.animation_current == 1);
    expect_true(action(Action::Delete));
    expect_true((display->Keyframes.size() == 1 && gui.keyframe_snapshots.size() == 1));
    expect_true(gui.animation_current == 0);
    key(RAYIMGUI_KEY_SLASH);
    expect_true((display->Keyframes.empty() && gui.keyframe_snapshots.empty()));
    expect_true(gui.animation_current == -1);
  }

  test_that("duration boxes edit the transitions between snapshots") {
    RayrenderGui gui;
    gui.can_edit = true;
    struct Input {
      uint64_t target = 0;
      int32_t value = 0;
      int disabled = 0;
      std::vector<bool> scopes;
      std::vector<uint64_t> order;
      std::vector<float> widths;
    } input;
    rayimgui_api_v1 api{};
    api.widget = [](uint64_t owner,
                    rayimgui_widget_v1* widget,
                    rayimgui_item_v1* item,
                    rayimgui_error_v1*) -> int32_t {
      auto& input = *reinterpret_cast<Input*>(uintptr_t(owner));
      item->flags = RAYIMGUI_VISIBLE;
      if (widget->kind == RAYIMGUI_DISABLED_BEGIN) {
        input.scopes.push_back(*widget->integers != 0);
        input.disabled += input.scopes.back();
      } else if (widget->kind == RAYIMGUI_DISABLED_END) {
        input.disabled -= input.scopes.back();
        input.scopes.pop_back();
      }
      if (widget->kind == RAYIMGUI_INT ||
          ((widget->kind == RAYIMGUI_BUTTON || widget->kind == RAYIMGUI_IMAGE_BUTTON) &&
           widget->id >= 100000)) {
        input.order.push_back(widget->id);
      }
      if (widget->kind == RAYIMGUI_INT) {
        input.widths.push_back(widget->width);
        if (widget->id == input.target && !input.disabled) {
          *widget->integers = input.value;
          item->flags |= RAYIMGUI_BEGIN | RAYIMGUI_CHANGED;
        }
      }
      return 0;
    };
    auto draw = [&] {
      input.order.clear();
      input.widths.clear();
      expect_true(gui.draw_animation(&api,
                                     uint64_t(reinterpret_cast<uintptr_t>(&input)),
                                     nullptr) == 0);
    };
    draw();
    expect_true(input.widths.empty());
    for (int i = 0; i < 3; ++i) {
      gui.save_keyframe_snapshot({double(i)});
    }
    input.target = (uint64_t(1) << 62) | gui.keyframe_snapshots[0].id;
    input.value = 12;
    draw();
    expect_true(input.widths == std::vector<float>({30, 30}));
    expect_true((input.order ==
                 std::vector<uint64_t>{
                     100001, input.target, 100002, (uint64_t(1) << 62) | 2, 100003}));
    expect_true((gui.animation_timing_custom && gui.history_dirty));
    expect_true(gui.keyframe_snapshots[0].frames_to_next == 12);
    input.value = 0;
    draw();
    expect_true(gui.keyframe_snapshots[0].frames_to_next == 1);
    input.value = 500000;
    draw();
    expect_true(gui.keyframe_snapshots[0].frames_to_next ==
                RayrenderGui::MaxSegmentFrames);
    gui.animation_closed = 1;
    input.target = (uint64_t(1) << 62) | gui.keyframe_snapshots[2].id;
    input.value = 21;
    draw();
    expect_true(input.widths.size() == 3);
    expect_true(gui.keyframe_snapshots[2].frames_to_next == 21);
    gui.animation_playing = true;
    input.value = 22;
    draw();
    expect_true(gui.keyframe_snapshots[2].frames_to_next == 21);
    gui.animation_playing = false;
    gui.animation_closed = 0;
    draw();
    expect_true(input.widths.size() == 2);
    expect_true(gui.keyframe_snapshots[2].frames_to_next == 21);
  }

  test_that("edited keyframe durations control playback and its selected thumbnail") {
    Transform object, world;
    camera cam(point3f(0, 1, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
    auto display = NativeTestDisplay(cam, object, world);
    RayrenderGui gui;
    display->AttachNativeGui(&gui, true, false);
    for (int i = 0; i < 3; ++i) {
      cam.update_position_absolute(point3f(i * 3, 1, -10));
      display->SaveCurrentKeyframe(0);
    }
    // Untouched controls and playback use 30 intervals per transition,
    // independently of the number of saved views and open/closed mode.
    Float default_angle = 0;
    for (bool closed : {false, true}) {
      display->keyframe_motion_closed = closed;
      expect_true(display->StartPreviewMotion(default_angle));
      expect_true(display->preview_motion.nrows() == (closed ? 91 : 61));
      for (const auto& snapshot : gui.keyframe_snapshots) {
        expect_true(snapshot.frames_to_next == 30);
      }
      display->CancelPreviewMotion(&default_angle);
    }
    display->SetKeyframeMotionArgs(Rcpp::List::create(Rcpp::_["frames"] = 20));
    display->SyncNativeAnimationState();
    expect_true(gui.keyframe_snapshots[0].frames_to_next == 10);
    expect_true(gui.keyframe_snapshots[1].frames_to_next == 9);
    gui.animation_timing_custom = true;
    gui.keyframe_snapshots[0].frames_to_next = 2;
    gui.keyframe_snapshots[1].frames_to_next = 5;
    gui.keyframe_snapshots[2].frames_to_next = 3;
    Float angle = 0;
    for (bool closed : {false, true}) {
      display->keyframe_motion_closed = closed;
      display->SyncNativeAnimationState();
      expect_true(display->StartPreviewMotion(angle));
      expect_true(display->preview_motion.nrows() == (closed ? 11 : 8));
      const int count = display->preview_motion.nrows();
      for (int frame = 0; frame < count; ++frame) {
        expect_true(display->AdvancePreviewMotion(&angle));
        const int selected = frame < 2 ? 0 : frame < 7 ? 1 : frame < 10 ? 2 : 0;
        expect_true(display->current_keyframe == selected);
        if (frame == 0 || frame == 2 || frame == 7 || frame == 10) {
          expect_true(std::abs(cam.get_origin()[0] - selected * 3) < 1e-5);
        }
      }
      display->AdvancePreviewMotion(&angle);
      expect_false(display->IsPreviewMotionActive());
      expect_true(cam.get_origin()[0] == 6);
    }
  }

  test_that("keyframe timing survives undo, deletion and loop changes") {
    Transform object, world;
    camera cam(point3f(0, 1, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
    auto display = NativeTestDisplay(cam, object, world);
    RayrenderGui gui;
    display->AttachNativeGui(&gui, true, false);
    for (int i = 0; i < 3; ++i) {
      cam.update_position_absolute(point3f(i * 3, 1, -10));
      display->SaveCurrentKeyframe(0);
    }
    display->SyncNativeAnimationState();
    display->BeginNativeHistory();
    const int original = gui.keyframe_snapshots[0].frames_to_next;
    const auto image = gui.keyframe_snapshots[0].image;
    ++gui.history_epoch;
    gui.keyframe_snapshots[0].frames_to_next = 12;
    gui.animation_timing_custom = gui.history_dirty = true;
    display->CommitNativeEdits(nullptr);
    gui.history_requests.push_back(-1);
    display->CommitNativeEdits(nullptr);
    expect_false(gui.animation_timing_custom);
    expect_true(gui.keyframe_snapshots[0].frames_to_next == original);
    gui.history_requests.push_back(1);
    display->CommitNativeEdits(nullptr);
    expect_true(gui.animation_timing_custom);
    expect_true(gui.keyframe_snapshots[0].frames_to_next == 12);
    expect_true(gui.keyframe_snapshots[0].image == image);
    ++gui.history_epoch;
    gui.animation_closed = 1;
    gui.animation_settings_pending = true;
    display->CommitNativeEdits(nullptr);
    expect_true(gui.keyframe_snapshots[0].frames_to_next == 12);
    const auto last = gui.keyframe_snapshots[2].id;
    display->current_keyframe = 1;
    display->SyncNativeAnimationState();
    ++gui.history_epoch;
    gui.animation_action = RayrenderGui::AnimationAction::Delete;
    display->CommitNativeEdits(nullptr);
    expect_true(gui.keyframe_snapshots.size() == 2);
    expect_true(gui.keyframe_snapshots[1].id == last);
    expect_true(gui.keyframe_snapshots[0].frames_to_next == 12);
    gui.history_requests.push_back(-1);
    display->CommitNativeEdits(nullptr);
    expect_true(gui.keyframe_snapshots.size() == 3);
    expect_true(gui.keyframe_snapshots[0].frames_to_next == 12);
    expect_true(gui.keyframe_snapshots[2].id == last);
  }

  test_that("populated animation panels draw through the real provider") {
    Rcpp::Function require = Rcpp::Environment::base_env()["requireNamespace"];
    if (!Rcpp::as<bool>(require("rayimgui", Rcpp::_["quietly"] = true))) {
      return;
    }
    Rcpp::Function acquire =
        Rcpp::Environment::namespace_env("rayimgui")["acquire_api"];
    Rcpp::RObject handle = acquire();
    RayrenderGui gui;
    expect_true(rayimgui_api_from_R_v1(
                    handle, sizeof(rayimgui_api_v1), RAYIMGUI_CAP_HEADLESS, &gui.api) ==
                0);
    if (gui.api->header.abi_minor < 9) {
      return;
    }
    struct CloseGui {
      RayrenderGui& gui;
      ~CloseGui() {
        gui.close();
      }
    } cleanup{gui};
    rayimgui_session_desc_v1 desc{
        sizeof(desc), RAYIMGUI_HEADLESS, "Animation test", 14, 1440, 900};
    expect_true(gui.api->open(&desc, &gui.session, &gui.error) == 0);
    rayimgui_callback_v1 callback{sizeof(callback), 1, RayrenderGui::draw, &gui};
    expect_true(gui.api->register_callback(
                    gui.session, &callback, &gui.callback, &gui.error) == 0);
    Transform object, world;
    camera cam(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
    auto display = NativeTestDisplay(cam, object, world);
    display->AttachNativeGui(&gui, true, false);
    gui.width = gui.height = 16;
    gui.pixels.assign(16 * 16 * 4, 160);
    for (int i = 0; i < 12; ++i) {
      cam.update_position_absolute(point3f(i, 0, -10));
      auto state = display->CreateCurrentKeyframe(0);
      gui.snapshot_camera.clear();
      for (R_xlen_t j = 0; j < state.size(); ++j) {
        gui.snapshot_camera.push_back(Rcpp::as<double>(state[j]));
      }
      display->SaveCurrentKeyframe(0);
    }
    display->SyncNativeAnimationState();
    gui.publish();
    for (int frame = 0; frame < 4; ++frame) {
      expect_false(gui.poll(true));
    }
    expect_true(gui.keyframe_snapshots.size() == 12);
    for (const auto& snapshot : gui.keyframe_snapshots) {
      expect_true(snapshot.texture != 0);
    }
    gui.close();
    expect_true(gui.keyframe_snapshots.empty());
    expect_true(gui.session == 0);
  }

  test_that(
      "replacing a keyframe preserves timing and restores camera and thumbnail versions") {
    Transform object, world;
    camera cam(point3f(0, 1, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
    auto display = NativeTestDisplay(cam, object, world);
    RayrenderGui gui;
    display->AttachNativeGui(&gui, true, false);
    struct Provider {
      uint64_t next = 1, context = 0;
      bool replace = false;
      std::map<uint64_t, uint8_t> textures;
      std::string label;
    } provider;
    rayimgui_api_v1 api{};
    api.texture_create = [](uint64_t owner,
                            const rayimgui_image_v1* image,
                            uint64_t* texture,
                            rayimgui_error_v1*) -> int32_t {
      auto& provider = *reinterpret_cast<Provider*>(uintptr_t(owner));
      *texture = provider.next++;
      provider.textures[*texture] = static_cast<const uint8_t*>(image->pixels)[0];
      return 0;
    };
    api.texture_destroy =
        [](uint64_t owner, uint64_t* texture, rayimgui_error_v1*) -> int32_t {
      auto& provider = *reinterpret_cast<Provider*>(uintptr_t(owner));
      provider.textures.erase(*texture);
      *texture = 0;
      return 0;
    };
    api.widget = [](uint64_t owner,
                    rayimgui_widget_v1* widget,
                    rayimgui_item_v1* item,
                    rayimgui_error_v1*) -> int32_t {
      auto& provider = *reinterpret_cast<Provider*>(uintptr_t(owner));
      item->flags = RAYIMGUI_VISIBLE;
      if (widget->kind == RAYIMGUI_CONTEXT_BEGIN) {
        item->flags = widget->id == ((uint64_t(1) << 61) | provider.context)
                          ? RAYIMGUI_VISIBLE
                          : 0;
      } else if (widget->kind == RAYIMGUI_MENU_ITEM) {
        provider.label.assign(widget->label, widget->label_length);
        if (provider.replace) {
          item->flags |= RAYIMGUI_CHANGED | RAYIMGUI_COMMIT;
        }
      }
      return 0;
    };
    gui.api = &api;
    gui.session = uint64_t(reinterpret_cast<uintptr_t>(&provider));
    gui.width = gui.height = 4;
    auto rendered = [&](uint8_t value) {
      gui.pixels.assign(4 * 4 * 4, value);
      auto keyframe = display->CreateCurrentKeyframe(0);
      gui.snapshot_camera.clear();
      for (R_xlen_t i = 0; i < keyframe.size(); ++i) {
        gui.snapshot_camera.push_back(Rcpp::as<double>(keyframe[i]));
      }
      gui.capture_keyframe_snapshots();
    };
    rendered(20);
    display->SaveCurrentKeyframe(0);
    cam.update_position_absolute(point3f(3, 1, -10));
    rendered(80);
    display->SaveCurrentKeyframe(0);
    gui.animation_timing_custom = true;
    gui.keyframe_snapshots[0].frames_to_next = 42;
    const auto old_id = gui.keyframe_snapshots[0].id;
    const auto next_id = gui.keyframe_snapshots[1].id;
    const auto old_image = gui.keyframe_snapshots[0].image;
    const auto old_texture = gui.keyframe_snapshots[0].texture;
    // Navigate away without publishing a completed render of the new camera yet.
    cam.update_position_absolute(point3f(9, 2, -7));
    cam.update_fov_absolute(35);
    const auto current = cam.get_origin();
    display->SyncNativeAnimationState();
    display->SyncNativeCameraControls(true);
    display->BeginNativeHistory();
    provider.context = old_id;
    expect_true(gui.draw_animation(&api, gui.session, nullptr) == 0);
    expect_true(provider.label == "Replace with current view");
    expect_true(gui.animation_action == RayrenderGui::AnimationAction::None);
    expect_true((cam.get_origin() - current).length() == 0);
    provider.replace = true;
    expect_true(gui.draw_animation(&api, gui.session, nullptr) == 0);
    expect_true(gui.animation_action == RayrenderGui::AnimationAction::Replace);
    display->CommitNativeEdits(nullptr);
    expect_true(gui.can_undo);
    expect_true((cam.get_origin() - current).length() == 0);
    expect_true(Rcpp::as<double>(display->Keyframes[0]["x"]) == 9);
    expect_true(Rcpp::as<double>(display->Keyframes[0]["fov"]) == 35);
    expect_true((display->Keyframes.size() == 2 && gui.keyframe_snapshots.size() == 2));
    expect_true(gui.keyframe_snapshots[0].frames_to_next == 42);
    expect_true(gui.keyframe_snapshots[1].id == next_id);
    expect_true(gui.keyframe_snapshots[0].pixels.empty());
    expect_true(provider.textures.count(old_texture) == 0);
    const auto replacement_id = gui.keyframe_snapshots[0].id;
    expect_true(replacement_id != old_id);
    // Delayed capture updates only the replacement image, never the undo image.
    rendered(190);
    expect_true(gui.keyframe_snapshots[0].pixels[0] == 190);
    expect_true(old_image->pixels[0] == 20);
    expect_true(provider.textures.size() == 2);
    gui.history_requests.push_back(-1);
    display->CommitNativeEdits(nullptr);
    expect_true(gui.keyframe_snapshots[0].id == old_id);
    expect_true(Rcpp::as<double>(display->Keyframes[0]["x"]) == 0);
    expect_true(gui.keyframe_snapshots[0].pixels[0] == 20);
    expect_true(provider.textures.at(gui.keyframe_snapshots[0].texture) == 20);
    expect_true(gui.keyframe_snapshots[0].frames_to_next == 42);
    gui.history_requests.push_back(1);
    display->CommitNativeEdits(nullptr);
    expect_true(gui.keyframe_snapshots[0].id == replacement_id);
    expect_true(Rcpp::as<double>(display->Keyframes[0]["x"]) == 9);
    expect_true(gui.keyframe_snapshots[0].pixels[0] == 190);
    expect_true(provider.textures.at(gui.keyframe_snapshots[0].texture) == 190);
    expect_true(provider.textures.size() == 2);
    // A menu request for a version that has since been replaced is harmless.
    gui.animation_target = old_id;
    gui.animation_action = RayrenderGui::AnimationAction::Replace;
    display->CommitNativeEdits(nullptr);
    expect_true(gui.keyframe_snapshots[0].id == replacement_id);
    gui.delete_keyframe_snapshot(1);
    gui.delete_keyframe_snapshot(0);
    expect_true(provider.textures.empty());
    gui.api = nullptr;
    gui.session = 0;
  }

  test_that("native keyframe limits preserve room for viewport textures") {
    Transform object, world;
    camera cam(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
    auto display = NativeTestDisplay(cam, object, world);
    RayrenderGui gui;
    display->AttachNativeGui(&gui, true, false);
    for (size_t i = 0; i <= RayrenderGui::MaxKeyframeSnapshots; ++i) {
      display->SaveCurrentKeyframe(0);
    }
    expect_true(display->Keyframes.size() == RayrenderGui::MaxKeyframeSnapshots);
    expect_true(gui.keyframe_snapshots.size() == display->Keyframes.size());
    expect_false(gui.animation_message.empty());
    display->DeleteCurrentKeyframe(nullptr);
    display->SaveCurrentKeyframe(0);
    expect_true(gui.animation_message.empty());
    expect_true(gui.keyframe_snapshots.size() == RayrenderGui::MaxKeyframeSnapshots);
  }

  test_that("keyframe snapshots match rendered cameras and are immutable") {
    RayrenderGui gui;
    gui.width = 288;
    gui.height = 180;
    gui.pixels.assign(size_t(gui.width) * gui.height * 4, 60);
    gui.snapshot_camera = {1, 2, 3};
    int uploads = 0, destroys = 0;
    struct Textures {
      int* uploads;
      int* destroys;
    } textures{&uploads, &destroys};
    rayimgui_api_v1 api{};
    api.texture_create = [](uint64_t owner,
                            const rayimgui_image_v1* image,
                            uint64_t* handle,
                            rayimgui_error_v1*) -> int32_t {
      auto& state = *reinterpret_cast<Textures*>(uintptr_t(owner));
      ++*state.uploads;
      *handle = *state.uploads;
      return image->width == 144 && image->height == 90 ? 0 : RAYIMGUI_INVALID;
    };
    api.texture_destroy =
        [](uint64_t owner, uint64_t* handle, rayimgui_error_v1*) -> int32_t {
      auto& state = *reinterpret_cast<Textures*>(uintptr_t(owner));
      ++*state.destroys;
      *handle = 0;
      return 0;
    };
    gui.api = &api;
    gui.session = uint64_t(reinterpret_cast<uintptr_t>(&textures));
    gui.save_keyframe_snapshot({1, 2, 3});
    expect_true(uploads == 1);
    expect_true(gui.keyframe_snapshots[0].pixels[0] == 60);
    expect_true(gui.keyframe_snapshots[0].pixels[3] == 255);
    gui.save_keyframe_snapshot({4, 5, 6});
    expect_true((uploads == 1 && gui.keyframe_snapshots[1].pixels.empty()));
    const auto id = gui.keyframe_snapshots[1].id;
    gui.pixels.assign(gui.pixels.size(), 180);
    gui.snapshot_camera = {4, 5, 6};
    gui.capture_keyframe_snapshots();
    expect_true(uploads == 2);
    expect_true(gui.keyframe_snapshots[0].pixels[0] == 60);
    expect_true(gui.keyframe_snapshots[1].pixels[0] == 180);
    gui.capture_keyframe_snapshots();
    expect_true(uploads == 2);
    gui.delete_keyframe_snapshot(0);
    expect_true((destroys == 1 && gui.keyframe_snapshots[0].id == id));
    gui.delete_keyframe_snapshot(0);
    expect_true((destroys == 2 && gui.keyframe_snapshots.empty()));
    gui.api = nullptr;
    gui.session = 0;
  }

  test_that("animation playback locks editing and restores the original camera") {
    Transform object, world;
    camera cam(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
    auto display = NativeTestDisplay(cam, object, world);
    RayrenderGui gui;
    display->AttachNativeGui(&gui, true, false);
    using Action = RayrenderGui::AnimationAction;
    auto action = [&](Action value) {
      gui.animation_action = value;
      return display->ApplyNativeControls(nullptr);
    };
    display->SetKeyframeMotionArgs(Rcpp::List::create(Rcpp::Named("frames") = 6));
    action(Action::Save);
    cam.update_position_absolute(point3f(5, 1, -10));
    action(Action::Save);
    const auto original = cam.get_origin();
    expect_true(action(Action::Play));
    expect_true(gui.animation_playing);
    action(Action::Delete);
    action(Action::Save);
    expect_true((gui.keyframe_snapshots.size() == 2 && display->Keyframes.size() == 2));
    expect_true(action(Action::Stop));
    expect_false(gui.animation_playing);
    expect_true((cam.get_origin() - original).length() == 0);
    expect_true(gui.animation_current == 1);
    action(Action::Play);
    for (int i = 0; i < 10; ++i) {
      display->ApplyNativeControls(nullptr);
    }
    expect_false(gui.animation_playing);
    expect_true((cam.get_origin() - original).length() == 0);
    display->AttachNativeGui(&gui, false, false);
    action(Action::Save);
    expect_true(display->Keyframes.size() == 2);
  }
  test_that(
      "viewport clicks stop playing or paused animation without picking the old image") {
    Transform object, world;
    camera cam(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
    auto display = NativeTestDisplay(cam, object, world);
    RayrenderGui gui;
    display->AttachNativeGui(&gui, true, false);
    display->SaveCurrentKeyframe(0);
    cam.update_position_absolute(point3f(5, 1, -10));
    display->SaveCurrentKeyframe(0);
    cam.update_position_absolute(point3f(9, 2, -10));
    const auto original = cam.get_origin();
    gui.animation_loop = 1;
    gui.object.enabled = true;
    using Action = RayrenderGui::AnimationAction;
    auto action = [&](Action value) {
      gui.animation_action = value;
      display->ApplyNativeControls(nullptr);
    };
    rayimgui_viewport_v1 viewport{sizeof(viewport)};
    viewport.width = viewport.height = 100;
    for (unsigned button :
         {RAYIMGUI_MOUSE_LEFT, RAYIMGUI_MOUSE_RIGHT, RAYIMGUI_MOUSE_MIDDLE}) {
      action(Action::Play);
      if (button == RAYIMGUI_MOUSE_RIGHT) {
        action(Action::Play); // Pause.
      }
      rayimgui_input_v1 input{sizeof(input)};
      input.keyboard_available = input.mouse_available = 1;
      gui.collect_input(input, viewport); // Hovering is not an interaction.
      display->ApplyNativeControls(nullptr);
      expect_true(gui.animation_playing);
      input.mouse_available = 0;
      input.mouse_clicked = button;
      gui.collect_input(input, viewport); // A widget owns this click.
      expect_true(gui.animation_action == Action::None);
      input.mouse_available = 1;
      input.modifiers = RAYIMGUI_SHIFT;
      gui.collect_input(input, viewport);
      expect_true(gui.animation_action == Action::Stop);
      expect_false((gui.pick_pending || gui.object.pick_pending));
      // A later unfocused poll cannot discard Stop; a queued M cannot undo it.
      input.keyboard_available = input.mouse_available = 0;
      gui.collect_input(input, viewport);
      gui.keys.push_back({RAYIMGUI_KEY_M, 0, 1});
      display->CommitNativeEdits(nullptr);
      expect_false((gui.animation_playing || gui.animation_paused));
      expect_true((cam.get_origin() - original).length() == 0);
      expect_true(gui.animation_loop == 1);
    }
    // Once restored, the next Shift-click uses ordinary object selection again.
    rayimgui_input_v1 input{sizeof(input)};
    input.mouse_available = 1;
    input.mouse_clicked = RAYIMGUI_MOUSE_LEFT;
    input.modifiers = RAYIMGUI_SHIFT;
    gui.collect_input(input, viewport);
    expect_true(gui.object.pick_pending);
    expect_true(gui.animation_action == Action::None);
  }

  test_that("viewport navigation stops playback then moves from the restored camera") {
    for (bool paused : {false, true}) {
      Transform object, world;
      camera cam(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
      auto display = NativeTestDisplay(cam, object, world);
      RayrenderGui gui;
      display->AttachNativeGui(&gui, true, false);
      display->SaveCurrentKeyframe(0);
      cam.update_position_absolute(point3f(5, 1, -10));
      display->SaveCurrentKeyframe(0);
      cam.update_position_absolute(point3f(9, 2, -10));
      const auto original = cam.get_origin();
      const auto expected = original + cam.get_w() * Float(.5);
      gui.animation_loop = 1;
      gui.orbit = 0;
      display->SyncNativeCameraControls(true);
      display->BeginNativeHistory();
      gui.animation_action = RayrenderGui::AnimationAction::Play;
      display->CommitNativeEdits(nullptr);
      if (paused) {
        gui.animation_action = RayrenderGui::AnimationAction::Play;
        display->CommitNativeEdits(nullptr);
      }
      rayimgui_viewport_v1 viewport{sizeof(viewport)};
      rayimgui_input_v1 input{sizeof(input)};
      input.keyboard_available = 1;
      input.keys_pressed = RAYIMGUI_KEY_BIT(RAYIMGUI_KEY_W);
      gui.collect_input(input, viewport);
      display->CommitNativeEdits(nullptr);
      expect_false((gui.animation_playing || gui.animation_paused));
      expect_true((cam.get_origin() - expected).length() < 1e-5);
      gui.history_requests.push_back(-1);
      display->CommitNativeEdits(nullptr);
      expect_false(gui.animation_playing);
      expect_true((cam.get_origin() - original).length() < 1e-5);
    }
  }

  test_that("loop playback repeats paths and preserves pause, stop and undo") {
    for (bool closed : {false, true}) {
      Transform object, world;
      camera cam(point3f(0, 1, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
      auto display = NativeTestDisplay(cam, object, world);
      RayrenderGui gui;
      display->AttachNativeGui(&gui, true, false);
      display->SaveCurrentKeyframe(0);
      cam.update_position_absolute(point3f(5, 1, -10));
      display->SaveCurrentKeyframe(0);
      display->keyframe_motion_closed = closed;
      gui.animation_timing_custom = true;
      for (auto& snapshot : gui.keyframe_snapshots) {
        snapshot.frames_to_next = 2;
      }
      cam.update_position_absolute(point3f(9, 2, -10));
      const auto original = cam.get_origin();
      display->SetCameraMotionBlur(true);
      display->SetShutterSpeed(1);
      display->SyncNativeCameraControls(true);
      display->BeginNativeHistory();
      expect_true(gui.animation_loop == 0);
      gui.animation_loop = 1;
      gui.history_dirty = true;
      ++gui.history_epoch;
      display->CommitNativeEdits(nullptr);
      gui.history_requests.push_back(-1);
      display->CommitNativeEdits(nullptr);
      expect_true(gui.animation_loop == 0);
      gui.history_requests.push_back(1);
      display->CommitNativeEdits(nullptr);
      expect_true(gui.animation_loop == 1);
      using Action = RayrenderGui::AnimationAction;
      auto action = [&](Action command) {
        gui.animation_action = command;
        display->ApplyNativeControls(nullptr);
      };
      action(Action::Play);
      const int count = display->preview_motion.nrows();
      const Rcpp::NumericVector xs = display->preview_motion["x"];
      bool repeated = true, sharp_restart = true;
      // Exercise two full passes through the real playback checkpoint.
      for (int step = 1; step <= count * 2; ++step) {
        display->ApplyNativeControls(nullptr);
        repeated &= gui.animation_playing &&
                    display->preview_motion_frame == step % count + 1 &&
                    std::abs(cam.get_origin().xyz.x - xs[step % count]) < 1e-5;
        if (step % count == 0) {
          const auto ray = cam.get_ray(.5, .5, point3f(0), 0);
          sharp_restart &= (ray.origin() - cam.get_origin()).length() < 1e-5;
        }
      }
      expect_true(repeated);
      expect_true(sharp_restart);
      action(Action::Play); // Pause the repeating path.
      const int paused_frame = display->preview_motion_frame;
      gui.animation_loop = 0; // The paused view must stay put until resumed.
      display->ApplyNativeControls(nullptr);
      expect_true(
          (gui.animation_paused && display->preview_motion_frame == paused_frame));
      action(Action::Play);
      for (int i = 0; i < count; ++i) {
        display->ApplyNativeControls(nullptr);
      }
      expect_false(gui.animation_playing);
      expect_true((cam.get_origin() - original).length() == 0);
      // Stop must restore the original camera even after another complete pass.
      gui.animation_loop = 1;
      action(Action::Play);
      for (int i = 0; i < count; ++i) {
        display->ApplyNativeControls(nullptr);
      }
      action(Action::Stop);
      expect_false(gui.animation_playing);
      expect_true((cam.get_origin() - original).length() == 0);
    }
  }

  test_that("animation panel M and buttons share pause resume and stop") {
    Transform object, world;
    camera cam(point3f(0, 1, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
    auto display = NativeTestDisplay(cam, object, world);
    RayrenderGui gui;
    display->AttachNativeGui(&gui, true, false);
    display->SaveCurrentKeyframe(0);
    cam.update_position_absolute(point3f(5, 1, -10));
    display->SaveCurrentKeyframe(0);
    display->SyncNativeAnimationState();
    struct Panel {
      rayimgui_input_v1 input{sizeof(input)};
      uint64_t click = 0;
      std::vector<std::string> buttons;
    } panel;
    rayimgui_api_v1 api{};
    api.widget = [](uint64_t owner,
                    rayimgui_widget_v1* widget,
                    rayimgui_item_v1* item,
                    rayimgui_error_v1*) -> int32_t {
      auto& panel = *reinterpret_cast<Panel*>(uintptr_t(owner));
      item->flags = RAYIMGUI_VISIBLE;
      if (widget->kind == RAYIMGUI_BUTTON) {
        panel.buttons.emplace_back(widget->label, widget->label_length);
        if (widget->id == panel.click) {
          item->flags |= RAYIMGUI_CHANGED;
        }
      } else if (widget->kind == RAYIMGUI_CHECKBOX && widget->id == panel.click) {
        *widget->integers = !*widget->integers;
        item->flags |= RAYIMGUI_CHANGED | RAYIMGUI_BEGIN | RAYIMGUI_COMMIT;
      }
      return 0;
    };
    api.window_input =
        [](uint64_t owner, rayimgui_input_v1* input, rayimgui_error_v1*) -> int32_t {
      *input = reinterpret_cast<Panel*>(uintptr_t(owner))->input;
      return 0;
    };
    auto draw = [&] {
      panel.buttons.clear();
      expect_true(gui.draw_animation(&api,
                                     uint64_t(reinterpret_cast<uintptr_t>(&panel)),
                                     nullptr) == 0);
    };
    // Both clicks select the view; the following M is read in the same panel.
    panel.click = 100000 + gui.keyframe_snapshots[0].id;
    for (int click = 0; click < 2; ++click) {
      draw();
      display->ApplyNativeControls(nullptr);
    }
    const auto original = cam.get_origin();
    expect_true(gui.animation_current == 0);
    panel.click = 0;
    panel.input.keyboard_available = 1;
    panel.input.keys_pressed = RAYIMGUI_KEY_BIT(RAYIMGUI_KEY_M);
    draw();
    expect_true(panel.buttons.front() == "Play (M)");
    // The unfocused viewport discards navigation keys without discarding the
    // separate playback command already queued by the Animation panel.
    rayimgui_input_v1 unfocused{sizeof(unfocused)};
    rayimgui_viewport_v1 view{sizeof(view)};
    gui.collect_input(unfocused, view);
    display->ApplyNativeControls(nullptr);
    expect_true(gui.animation_playing);
    expect_false(gui.animation_paused);
    // Loop playback remains editable while the path is running.
    panel.input.keys_pressed = 0;
    panel.click = 51;
    draw();
    expect_true((gui.animation_loop == 1 && gui.history_dirty));
    const int before_pause = display->preview_motion_frame;
    panel.input.keys_pressed = 0;
    panel.click = 45;
    draw();
    expect_true(panel.buttons.front() == "Pause (M)");
    display->ApplyNativeControls(nullptr);
    expect_true((gui.animation_playing && gui.animation_paused));
    const auto paused = cam.get_origin();
    panel.click = 0;
    for (int frame = 0; frame < 3; ++frame) {
      display->ApplyNativeControls(nullptr);
      expect_true(display->preview_motion_frame == before_pause);
      expect_true((cam.get_origin() - paused).length() == 0);
    }
    draw();
    expect_true(panel.buttons.front() == "Resume (M)");
    // The existing viewport shortcut reaches the same resume implementation.
    gui.keys.push_back({RAYIMGUI_KEY_M, 0, 1});
    display->ApplyNativeControls(nullptr);
    expect_false(gui.animation_paused);
    expect_true(display->preview_motion_frame == before_pause + 1);
    panel.click = 50;
    draw();
    display->ApplyNativeControls(nullptr);
    expect_false(gui.animation_playing);
    expect_false(gui.animation_paused);
    expect_true((cam.get_origin() - original).length() == 0);
    panel.click = 0;
    panel.input.keys_pressed = RAYIMGUI_KEY_BIT(RAYIMGUI_KEY_M);
    panel.input.modifiers = RAYIMGUI_CTRL;
    draw();
    expect_true(gui.animation_action == RayrenderGui::AnimationAction::None);
    panel.input.modifiers = 0;
    panel.input.keyboard_available = 0;
    draw();
    expect_true(gui.animation_action == RayrenderGui::AnimationAction::None);
  }

  test_that("export captures committed renderer state and reports write failures") {
    Transform object, world;
    camera cam(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
    auto display = NativeTestDisplay(cam, object, world);
    RayrenderGui gui;
    display->AttachNativeGui(&gui, true, false);
    display->preview_exposure_scale = .25;
    display->preview_exposure_adjustment = 2;
    gui.exposure = 99; // Uncommitted drafts must not enter exported settings.
    int calls = 0;
    display->export_scene = [&](const Rcpp::List& state, const std::string& filename) {
      ++calls;
      expect_true(Rcpp::as<double>(state["exposure"]) == .5);
      Rcpp::List camera = state["camera"];
      expect_true(Rcpp::as<double>(camera["z"]) == -10);
      expect_true(filename == "rayrender_scene.R");
      return std::string("saved-scene.R");
    };
    display->ApplyNativeExport();
    expect_true(calls == 0);
    gui.export_pending = true;
    display->ApplyNativeExport();
    expect_true(calls == 1);
    expect_false(gui.export_pending);
    expect_true(gui.export_message == "Exported: saved-scene.R");
    display->export_scene = [](const Rcpp::List&, const std::string&) -> std::string {
      throw std::runtime_error("No such export directory");
    };
    gui.export_pending = true;
    display->ApplyNativeExport();
    expect_true(gui.export_message == "Export failed: No such export directory");
    expect_false(display->terminate);
  }
  test_that("F toggles fast quality once from focused panes during playback too") {
    Transform object, world;
    camera cam(point3f(0, 1, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
    auto display = NativeTestDisplay(cam, object, world);
    RayrenderGui gui;
    display->AttachNativeGui(&gui, true, false);
    display->SaveCurrentKeyframe(0);
    cam.update_position_absolute(point3f(5, 1, -10));
    display->SaveCurrentKeyframe(0);
    display->SyncNativeCameraControls(true);
    display->BeginNativeHistory();
    struct Panel {
      uint64_t current = 0, focused = 40;
      uint32_t modifiers = 0;
      bool pressed = true, available = true;
    } panel;
    rayimgui_api_v1 api{};
    api.widget = [](uint64_t owner,
                    rayimgui_widget_v1* widget,
                    rayimgui_item_v1* item,
                    rayimgui_error_v1*) -> int32_t {
      auto& panel = *reinterpret_cast<Panel*>(uintptr_t(owner));
      item->flags = widget->kind == RAYIMGUI_CONTEXT_BEGIN ? 0 : RAYIMGUI_VISIBLE;
      if (widget->kind == RAYIMGUI_WINDOW_BEGIN) {
        panel.current = widget->id;
      }
      return 0;
    };
    api.window_input =
        [](uint64_t owner, rayimgui_input_v1* input, rayimgui_error_v1*) -> int32_t {
      const auto& panel = *reinterpret_cast<Panel*>(uintptr_t(owner));
      *input = {sizeof(*input)};
      input->keyboard_available = panel.available && panel.current == panel.focused;
      input->modifiers = panel.modifiers;
      if (input->keyboard_available) {
        input->keys_pressed = panel.pressed ? RAYIMGUI_KEY_BIT(RAYIMGUI_KEY_F) : 0;
        input->keys_repeated = RAYIMGUI_KEY_BIT(RAYIMGUI_KEY_F);
      }
      return 0;
    };
    api.input = api.window_input;
    api.viewport =
        [](uint64_t, rayimgui_viewport_v1* viewport, rayimgui_error_v1*) -> int32_t {
      *viewport = {sizeof(*viewport)};
      viewport->width = viewport->height = 100;
      return 0;
    };
    // Include the image input path to catch double toggles in the viewport.
    gui.texture = 1;
    auto draw = [&] {
      expect_true(gui.draw_impl(&api,
                                uint64_t(reinterpret_cast<uintptr_t>(&panel)),
                                nullptr) == 0);
    };
    for (uint64_t pane : {40, 1, 30, 31, 9}) {
      panel.focused = pane;
      panel.pressed = true;
      const bool before = display->write_fast_output;
      draw();
      expect_true((gui.fast_pending && bool(gui.fast_preview) != before));
      // Poll again before rendering completes. Holding F and losing viewport
      // keyboard ownership must neither repeat nor discard the pending toggle.
      panel.pressed = false;
      draw();
      display->CommitNativeEdits(nullptr);
      expect_true((display->write_fast_output != before &&
                   bool(gui.fast_preview) == display->write_fast_output));
    }
    const bool toggled = display->write_fast_output;
    gui.history_requests.push_back(-1);
    display->CommitNativeEdits(nullptr);
    expect_true(display->write_fast_output != toggled);
    gui.history_requests.push_back(1);
    display->CommitNativeEdits(nullptr);
    expect_true(display->write_fast_output == toggled);
    panel.pressed = true;
    panel.available = false; // Provider blocks text/numeric edits and focus loss.
    draw();
    expect_false(gui.fast_pending);
    panel.available = true;
    for (uint32_t modifier : {RAYIMGUI_CTRL, RAYIMGUI_ALT, RAYIMGUI_SUPER}) {
      panel.modifiers = modifier;
      draw();
      expect_false(gui.fast_pending);
    }
    panel.modifiers = 0;
    gui.animation_loop = 1;
    gui.animation_action = RayrenderGui::AnimationAction::Play;
    display->CommitNativeEdits(nullptr);
    panel.focused = 40;
    draw();
    display->CommitNativeEdits(nullptr);
    expect_true((gui.animation_playing && display->write_fast_output != toggled));
    gui.animation_action = RayrenderGui::AnimationAction::Play;
    display->CommitNativeEdits(nullptr);
    panel.focused = 1;
    draw();
    display->CommitNativeEdits(nullptr);
    expect_true((gui.animation_paused && display->write_fast_output == toggled));
    gui.animation_action = RayrenderGui::AnimationAction::Stop;
    display->CommitNativeEdits(nullptr);
  }

  test_that("viewport input queues bounded actions and respects focus and modifiers") {
    RayrenderGui gui;
    gui.can_edit = true;
    rayimgui_viewport_v1 view{sizeof(view)};
    view.width = 100;
    view.height = 200;
    view.mouse_x = 25;
    view.mouse_y = 50;
    rayimgui_input_v1 input{sizeof(input)};
    input.keyboard_available = 1;
    input.keys_repeated = RAYIMGUI_KEY_BIT(RAYIMGUI_KEY_W);
    for (int i = 0; i < 1000; ++i) {
      gui.collect_input(input, view);
    }
    expect_true((gui.keys.size() == 1 && gui.keys.front().count == 8));
    input.keyboard_available = 0;
    gui.collect_input(input, view);
    expect_true((gui.keys.empty()));
    input.keyboard_available = 1;
    input.modifiers = RAYIMGUI_CTRL;
    gui.collect_input(input, view);
    expect_true((gui.keys.empty()));
    input.mouse_available = 1;
    input.mouse_clicked = RAYIMGUI_MOUSE_RIGHT;
    gui.collect_input(input, view);
    expect_true((gui.pick_pending && !gui.pick_focus));
    expect_true((gui.pick_u == .75f && gui.pick_v == .75f));
    gui.pick_pending = false;
    input.mouse_available = 0;
    gui.collect_input(input, view);
    expect_false(gui.pick_pending);
  }
  test_that(
      "standard camera movement, orbit, pitch and reset use the renderer camera") {
    Transform object, world;
    camera cam(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
    auto display = NativeTestDisplay(cam, object, world);
    RayrenderGui gui;
    display->AttachNativeGui(&gui, true, true);
    auto key = [&](unsigned code, unsigned modifiers = 0) {
      gui.keys.push_back({code, modifiers, 1});
      return display->ApplyNativeControls(nullptr);
    };
    expect_true((key(RAYIMGUI_KEY_W)));
    expect_true((std::abs(cam.get_origin()[2] + 9.5f) < 1e-6));
    expect_true((key(RAYIMGUI_KEY_S)));
    expect_true((std::abs(cam.get_origin()[2] + 10) < 1e-6));
    expect_false(key(RAYIMGUI_KEY_E));
    expect_true((gui.movement_speed == 2));
    expect_true((key(RAYIMGUI_KEY_Q)));
    expect_true((std::abs(cam.get_origin()[1] - 10 / std::sqrt(101.0)) < 1e-6 &&
                 std::abs(cam.get_origin().length() - 10) < 1e-6));
    expect_true((key(RAYIMGUI_KEY_R)));
    expect_true((gui.movement_speed == 1 && cam.get_origin()[1] == 0));
    const auto original_direction = cam.get_w();
    expect_true((key(RAYIMGUI_KEY_W, RAYIMGUI_SHIFT)));
    expect_true(((cam.get_w() - original_direction).length() > 0));
    expect_true((cam.get_origin()[2] == -10));
    key(RAYIMGUI_KEY_R);
    key(RAYIMGUI_KEY_TAB);
    expect_false(gui.orbit);
    const auto direction = cam.get_w();
    key(RAYIMGUI_KEY_A);
    expect_true(((cam.get_w() - direction).length() < 1e-6));
    expect_true((std::abs(cam.get_origin()[0]) > .1));
    key(RAYIMGUI_KEY_R);
    key(RAYIMGUI_KEY_TAB);
    gui.movement_speed = 128;
    key(RAYIMGUI_KEY_W);
    expect_true((cam.get_origin()[2] == -10));
    key(RAYIMGUI_KEY_F);
    expect_true((display->write_fast_output && gui.fast_preview));
    gui.fast_preview = 0;
    gui.fast_pending = true;
    expect_true((display->ApplyNativeControls(nullptr)));
    expect_false(display->write_fast_output);
    expect_false(gui.fast_pending);
    gui.fast_preview = 1;
    gui.fast_pending = true;
    expect_true((display->ApplyNativeControls(nullptr)));
    expect_true((display->write_fast_output));
    key(RAYIMGUI_KEY_F);
    expect_false(display->write_fast_output);
    expect_false(gui.fast_preview);
    expect_false(display->render_requested);
    key(RAYIMGUI_KEY_ENTER);
    expect_true((display->render_requested));
    display->interactive = false;
    const auto before = cam.get_origin();
    key(RAYIMGUI_KEY_Q);
    expect_true(((cam.get_origin() - before).length() == 0));
  }
  test_that("date components retain typed values and clamp to real calendar dates") {
    RayrenderGui gui;
    gui.set_datetime("2024-01-31 21:42:13");
    expect_true((gui.date[0] == 2024 && gui.date[1] == 1 && gui.date[2] == 31));
    expect_true((gui.time[0] == 21 && gui.time[1] == 42 && gui.time[2] == 13));
    gui.date[1] = 2;
    gui.format_datetime();
    expect_true((std::string(gui.datetime) == "2024-02-29 21:42:13"));
    gui.date[0] = 2025;
    gui.format_datetime();
    expect_true((std::string(gui.datetime) == "2025-02-28 21:42:13"));
    gui.date[0] = 1900;
    gui.date[2] = 29;
    gui.format_datetime();
    expect_true((gui.date[2] == 28));
    gui.date[0] = 2000;
    gui.date[2] = 29;
    gui.format_datetime();
    expect_true((gui.date[2] == 29));
    gui.time[0] = 99;
    gui.time[1] = -1;
    gui.time[2] = 60;
    gui.format_datetime();
    expect_true((std::string(gui.datetime) == "2000-02-29 23:00:59"));
    gui.date[1] = 13;
    gui.date[2] = 0;
    gui.format_datetime();
    expect_true((std::string(gui.datetime) == "2000-12-01 23:00:59"));
  }
  test_that(
      "sky changes apply at a safe checkpoint while dragging; invalid dates keep rendering") {
    Transform object, world;
    camera cam(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
    auto display = NativeTestDisplay(cam, object, world);
    RayrenderGui gui;
    display->AttachNativeGui(&gui, true, false);
    int updates = 0;
    double elevation = 0;
    display->SetSunControls(30, 180, [&](double e, double) {
      ++updates;
      elevation = e;
    });
    gui.sun_pending = true;
    gui.sun_editing = true;
    gui.sun_elevation = 45;
    gui.sun_azimuth = 90;
    display->ApplyNativeSkyControls();
    expect_true((updates == 1 && elevation == 45));
    expect_true(display->ConsumeAtmosphereChange());
    gui.sun_pending = true;
    gui.sun_elevation = 50;
    gui.sun_editing = false;
    display->ApplyNativeSkyControls();
    expect_true((updates == 2 && elevation == 50));
    expect_true((display->ConsumeAtmosphereChange()));
    display->ApplyNativeSkyControls();
    expect_true((updates == 2));
    display->SetSkyControls(0, 0, "bad", [](double, double, const std::string&) {
      return std::string("Invalid date");
    });
    gui.location_pending = true;
    display->ApplyNativeSkyControls();
    expect_true((gui.sky_error == "Invalid date"));
    expect_false(display->ConsumeAtmosphereChange());
    display->SetSkyModelControls(0, [](int model) {
      return model == 1 ? std::string() : std::string("Unavailable sky");
    });
    gui.sky_model = 1;
    gui.sky_model_pending = true;
    display->ApplyNativeSkyControls();
    expect_true((gui.sky_model == 1 && gui.sky_error.empty()));
    expect_true(display->ConsumeAtmosphereChange());
    gui.sky_model = 0;
    gui.sky_model_pending = true;
    display->ApplyNativeSkyControls();
    expect_true((gui.sky_model == 1 && gui.sky_error == "Unavailable sky"));
    expect_false(display->ConsumeAtmosphereChange());
    display->SetSunControls(50, 90, [](double, double) {
      throw std::runtime_error("Sky failure");
    });
    gui.sun_pending = true;
    display->ApplyNativeSkyControls();
    expect_true((gui.sky_error == "Sky failure"));
    expect_false(gui.sun_pending);
    expect_false(display->ConsumeAtmosphereChange());
  }
  test_that("sun drags update live at Fast quality and undo as one gesture") {
    for (bool preferred_fast : {false, true}) {
      Transform object, world;
      camera cam(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
      auto display = NativeTestDisplay(cam, object, world);
      RayrenderGui gui;
      display->write_fast_output = preferred_fast;
      display->AttachNativeGui(&gui, true, false);
      Rcpp::List sky = Rcpp::List::create(Rcpp::_["model"] = 0,
                                          Rcpp::_["elevation"] = 10.,
                                          Rcpp::_["azimuth"] = 90.,
                                          Rcpp::_["haze"] = false,
                                          Rcpp::_["altitude"] = false);
      display->export_sky = [&] {
        return Rcpp::clone(sky);
      };
      display->prepare_sky_restore = [&](const Rcpp::List& saved) {
        auto next = Rcpp::clone(saved);
        return [&, next] {
          sky = next;
        };
      };
      int updates = 0;
      bool fail = false;
      display->SetSunControls(10, 90, [&](double elevation, double azimuth) {
        if (gui.sun_editing) {
          expect_true(display->write_fast_output);
        }
        if (fail) {
          throw std::runtime_error("Cannot rebuild sky");
        }
        ++updates;
        sky["elevation"] = elevation;
        sky["azimuth"] = azimuth;
      });
      display->BeginNativeHistory();
      ++gui.history_epoch;
      gui.sun_editing = true;
      for (double elevation : {20., 30., 40.}) {
        gui.sun_elevation = elevation;
        gui.sun_azimuth = 135;
        gui.sun_pending = true;
        display->CommitNativeEdits(nullptr);
        expect_true(display->ConsumeAtmosphereChange());
        expect_true(display->write_fast_output);
        expect_true(Rcpp::as<double>(sky["elevation"]) == elevation);
        expect_false(gui.sun_pending);
      }
      expect_true(updates == 3);
      // A failed live update retains the prior sky and the gesture still ends.
      fail = true;
      gui.sun_elevation = 50;
      gui.sun_pending = true;
      display->CommitNativeEdits(nullptr);
      expect_true(gui.sky_error == "Cannot rebuild sky");
      expect_false(display->ConsumeAtmosphereChange());
      expect_true(gui.sun_elevation == 40);
      fail = false;
      gui.sun_editing = false;
      display->CommitNativeEdits(nullptr);
      expect_true(display->write_fast_output == preferred_fast);
      expect_true(updates == 4);
      expect_false(display->native_sun_preview);
      expect_true(display->ConsumeAtmosphereChange());
      gui.history_requests.push_back(-1);
      display->CommitNativeEdits(nullptr);
      expect_true(gui.sun_elevation == 10);
      expect_true(Rcpp::as<double>(sky["elevation"]) == 10);
      expect_false(gui.can_undo);
      gui.history_requests.push_back(1);
      display->CommitNativeEdits(nullptr);
      expect_true(gui.sun_elevation == 40);
      expect_true(display->write_fast_output == preferred_fast);
    }
  }

  test_that(
      "camera inputs follow navigation, validate complete poses and export edits") {
    using Inputs = PreviewCameraInputs;
    Transform object, world;
    camera cam(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
    auto display = NativeTestDisplay(cam, object, world);
    RayrenderGui gui;
    display->AttachNativeGui(&gui, true, false);
    auto& input = gui.camera;
    expect_true((input.available && input.editable && input.Valid()));
    expect_true(input.values[Inputs::Position][2] == -10);
    expect_true(input.values[Inputs::Fov][0] == 60);
    display->BeginNativeHistory();
    ++gui.history_epoch;
    input.editing = true;
    for (int step = 1; step <= 3; ++step) {
      input.values[Inputs::Position][0] = step;
      input.values[Inputs::Target][0] = .5 * step;
      input.values[Inputs::Fov][0] = 60 + 5 * step;
      input.values[Inputs::Aperture][0] = .1 * step;
      input.values[Inputs::Focus][0] = 5 + step;
      input.pending = true;
      expect_true(display->CommitNativeEdits(nullptr));
      expect_true(display->write_fast_output);
      expect_true(cam.get_origin()[0] == step);
      expect_true(cam.get_lookat()[0] == .5 * step);
      expect_true(cam.get_fov() == 60 + 5 * step);
      expect_true(std::abs(cam.get_aperture() - .1 * step) < 1e-6);
      expect_true(cam.get_focal_distance() == 5 + step);
    }
    input.editing = false;
    expect_true(display->CommitNativeEdits(nullptr));
    expect_false(display->write_fast_output);
    Rcpp::List exported = display->NativeEditorState()["camera"];
    expect_true(Rcpp::as<double>(exported["x"]) == 3);
    expect_true(Rcpp::as<double>(exported["dx"]) == 1.5);
    expect_true(Rcpp::as<double>(exported["fov"]) == 75);
    expect_true(Rcpp::as<double>(exported["focal"]) == 8);
    gui.history_requests.push_back(-1);
    display->CommitNativeEdits(nullptr);
    expect_true(cam.get_origin()[0] == 0);
    expect_true(input.values[Inputs::Fov][0] == 60);
    expect_false(gui.can_undo);
    gui.history_requests.push_back(1);
    display->CommitNativeEdits(nullptr);
    expect_true(cam.get_origin()[0] == 3);
    expect_true(cam.get_lookat()[0] == 1.5);
    expect_true(input.values[Inputs::Focus][0] == 8);

    const auto valid_position = cam.get_origin();
    input.values[Inputs::Position] = input.values[Inputs::Target];
    input.pending = true;
    ++gui.history_epoch;
    expect_false(display->CommitNativeEdits(nullptr));
    expect_true((cam.get_origin() - valid_position).length() == 0);
    expect_false(input.errors[Inputs::Position].empty());
    expect_false(input.errors[Inputs::Target].empty());
    gui.history_requests.push_back(-1);
    display->CommitNativeEdits(nullptr);
    expect_true(input.Valid());
    expect_true(input.values[Inputs::Position][0] == 3);
    input.values[Inputs::Up] = {0, 0, 0};
    input.pending = true;
    expect_false(display->CommitNativeEdits(nullptr));
    expect_false(input.errors[Inputs::Up].empty());
    display->SyncNativeCameraControls(true);
    for (int axis = 0; axis < 3; ++axis) {
      input.values[Inputs::Up][axis] =
          input.values[Inputs::Target][axis] - input.values[Inputs::Position][axis];
    }
    input.pending = true;
    expect_false(display->CommitNativeEdits(nullptr));
    expect_false(input.errors[Inputs::Up].empty());
    display->SyncNativeCameraControls(true);
    input.values[Inputs::Fov][0] = 180;
    input.pending = true;
    expect_false(display->CommitNativeEdits(nullptr));
    expect_true(cam.get_fov() == 75);
    expect_false(input.errors[Inputs::Fov].empty());
    display->SyncNativeCameraControls(true);
    input.values[Inputs::Focus][0] = 0;
    input.pending = true;
    expect_false(display->CommitNativeEdits(nullptr));
    expect_true(cam.get_focal_distance() == 8);
    display->SyncNativeCameraControls(true);

    // Navigation and reset must update the pane, rather than letting its old
    // numbers overwrite the camera on a later unrelated GUI action.
    gui.keys.push_back({RAYIMGUI_KEY_W, 0, 1});
    expect_true(display->CommitNativeEdits(nullptr));
    expect_true(input.values[Inputs::Position][2] == double(cam.get_origin()[2]));
    expect_true(input.values[Inputs::Position][2] != -10);
    gui.request_reset = true;
    expect_true(display->CommitNativeEdits(nullptr));
    expect_true(input.values[Inputs::Position][2] == -10);
    expect_true(input.values[Inputs::Fov][0] == 60);
  }

  test_that("orthographic camera sizes and keyframes update the camera pane") {
    using Inputs = PreviewCameraInputs;
    Transform object, world;
    ortho_camera cam(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0), 8, 6, 0, 1, 1);
    auto display = NativeTestDisplay(cam, object, world);
    RayrenderGui gui;
    display->AttachNativeGui(&gui, true, false);
    expect_true(gui.camera.projection == Inputs::Orthographic);
    expect_true(gui.camera.Valid());
    expect_false(gui.camera.Visible(Inputs::Focus));
    auto action = [&](RayrenderGui::AnimationAction action) {
      gui.animation_action = action;
      display->CommitNativeEdits(nullptr);
    };
    action(RayrenderGui::AnimationAction::Save);
    gui.camera.values[Inputs::Ortho] = {16, 12, 0};
    gui.camera.values[Inputs::Target] = {1, 0, 0};
    gui.camera.pending = true;
    expect_true(display->CommitNativeEdits(nullptr));
    expect_true((cam.get_ortho()[0] == 16 && cam.get_ortho()[1] == 12));
    expect_true(cam.get_fov() == 0);
    action(RayrenderGui::AnimationAction::Save);
    action(RayrenderGui::AnimationAction::Previous);
    expect_true(gui.camera.values[Inputs::Ortho][0] == 8);
    expect_true(gui.camera.values[Inputs::Target][0] == 0);
    action(RayrenderGui::AnimationAction::Next);
    expect_true(gui.camera.values[Inputs::Ortho][0] == 16);
    expect_true(gui.camera.values[Inputs::Target][0] == 1);
  }

  test_that("camera input validation follows the active projection") {
    PreviewCameraInputs input;
    input.values[PreviewCameraInputs::Position] = {0, 0, -10};
    input.values[PreviewCameraInputs::Up] = {0, 1, 0};
    input.values[PreviewCameraInputs::Focus][0] = 10;
    input.values[PreviewCameraInputs::Ortho] = {8, 6, 0};
    input.projection = PreviewCameraInputs::Orthographic;
    expect_true(input.Validate());
    expect_false(input.Visible(PreviewCameraInputs::Fov));
    expect_true(input.Visible(PreviewCameraInputs::Ortho));
    input.values[PreviewCameraInputs::Ortho][1] = 0;
    expect_false(input.Validate());
    input.projection = PreviewCameraInputs::Environment;
    expect_true(input.Validate());
    expect_false(input.Visible(PreviewCameraInputs::Focus));
    input.values[PreviewCameraInputs::Position][0] =
        std::numeric_limits<double>::infinity();
    expect_false(input.Validate());
  }

  test_that("image sky callbacks replace lighting atomically after setup returns") {
    Rcpp::Function parse = Rcpp::Environment::base_env()["parse"];
    Rcpp::Function eval = Rcpp::Environment::base_env()["eval"];
    Rcpp::List fixture = eval(parse(Rcpp::_["text"] = R"(
      local({
        file = tempfile(fileext=".png")
        png::writePNG(array(.5,c(4,8,3)), file)
        lights = list(infinite_light(file))
        list(lights=lights, controls=list(index=0L,model=0L,
          latitude=40,longitude=-74,datetime="2026-06-21 16:00:00",
          elevation=45,azimuth=180,
          update=function(latitude,longitude,datetime,model,elevation,azimuth,...) {
            if(latitude > 90) return(list(error="Invalid latitude"))
            list(error="",lights=lights,
              elevation=if(is.null(elevation)) 45 else elevation,
              azimuth=if(is.null(azimuth)) 180 else azimuth)
          }))
      })
    )"),
                              Rcpp::Environment::namespace_env("rayrender"));
    Transform object, world;
    camera cam(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
    auto display = NativeTestDisplay(cam, object, world);
    RayrenderGui gui;
    display->AttachNativeGui(&gui, true, false);
    TextureCache textures;
    Rcpp::List lights = fixture["lights"];
    auto environment = std::make_shared<InfiniteAreaLight>(
        BuildInfiniteLights(lights, textures), 100, point3f(0), &object, &world);
    auto volume = std::make_shared<VolumeScene>();
    ConfigureNativeSky(
        *display, environment, volume, textures, fixture["controls"], lights);
    expect_true((gui.has_sky_model && gui.has_sun && gui.has_location));
    expect_false(gui.has_atmosphere);
    auto original = environment->light;
    gui.sky_model = 1;
    gui.sky_model_pending = true;
    display->ApplyNativeSkyControls();
    expect_true((environment->light != original && gui.sky_error.empty()));
    expect_true(display->ConsumeAtmosphereChange());
    original = environment->light;
    gui.latitude = 100;
    gui.location_pending = true;
    display->ApplyNativeSkyControls();
    expect_true(
        (environment->light == original && gui.sky_error == "Invalid latitude"));
    expect_false(display->ConsumeAtmosphereChange());
    gui.sun_elevation = 25;
    gui.sun_azimuth = 130;
    gui.sun_pending = true;
    display->ApplyNativeSkyControls();
    expect_true((gui.sky_error.empty() && gui.manual_sun));
    expect_true((gui.sun_elevation == 25 && gui.sun_azimuth == 130));
    expect_true(display->ConsumeAtmosphereChange());
  }

  test_that("Prague selection from an image sky enables transport and retains haze edits") {
    // Opt in with installed full-altitude data; ordinary checks need no download.
    if (!std::getenv("RAYRENDER_EDITOR_SKY_TEST")) {
      return;
    }
    Rcpp::Function parse = Rcpp::Environment::base_env()["parse"];
    Rcpp::Function eval = Rcpp::Environment::base_env()["eval"];
    Rcpp::List fixture = eval(parse(Rcpp::_["text"] = R"(
      local({
        sky = sky_light_image(40, -74,
          as.POSIXct("2026-06-21 16:00:00", tz="UTC"),
          resolution=16, sun=FALSE, moon=FALSE)
        controls = native_sky_controls(list(sky))
        update = controls$update
        controls$update = function(...) {
          result = update(...)
          if(identical(result$error, "") && result$lights[[1]]$type == "prague") {
            result$lights[[1]]$resolution = 16L
            result$lights[[1]]$transmission_table = FALSE
          }
          result
        }
        list(lights=prepare_scene_infinite_lights(list(sky)), controls=controls)
      })
    )"),
                              Rcpp::Environment::namespace_env("rayrender"));
    Transform object, world;
    camera cam(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
    auto display = NativeTestDisplay(cam, object, world);
    RayrenderGui gui;
    display->AttachNativeGui(&gui, true, false);
    TextureCache textures;
    Rcpp::List lights = fixture["lights"];
    auto environment = std::make_shared<InfiniteAreaLight>(
        BuildInfiniteLights(lights, textures), 100, point3f(0), &object, &world);
    auto volume = std::make_shared<VolumeScene>();
    ConfigureNativeSky(
        *display, environment, volume, textures, fixture["controls"], lights);
    expect_false(gui.has_atmosphere);
    expect_true(volume->atmosphere == nullptr);

    auto choose = [&](int model) {
      gui.sky_model = model;
      gui.sky_model_pending = true;
      display->ApplyNativeSkyControls();
      expect_true(gui.sky_error.empty());
      expect_true(display->ConsumeAtmosphereChange());
    };
    choose(1);
    expect_true((gui.has_atmosphere && gui.haze && gui.altitude));
    expect_true(environment->light->GetAtmosphere() != nullptr);
    expect_true(volume->atmosphere != nullptr);

    // Numeric edits wait for release, rebuild the light, and form one undo step.
    // Export stores committed values, never a still-active numeric draft.
    display->BeginNativeHistory();
    auto before_parameters = environment->light;
    ++gui.history_epoch;
    gui.meters_per_unit = 100;
    gui.base_altitude = 500;
    gui.atmosphere_parameters_editing = true;
    gui.atmosphere_parameters_pending = gui.history_dirty = true;
    display->CommitNativeEdits(nullptr);
    expect_true(environment->light == before_parameters);
    expect_true(Rcpp::as<double>(display->export_sky()["base_altitude"]) == 0);
    gui.meters_per_unit = 200;
    gui.history_dirty = true;
    display->CommitNativeEdits(nullptr);
    gui.atmosphere_parameters_editing = false;
    display->CommitNativeEdits(nullptr);
    expect_true(display->ConsumeAtmosphereChange());
    expect_true(environment->light != before_parameters);
    auto saved_parameters = display->export_sky();
    expect_true(Rcpp::as<double>(saved_parameters["base_altitude"]) == 500);
    expect_true(Rcpp::as<double>(saved_parameters["meters_per_unit"]) == 200);
    auto scaled_transmission = environment->light->GetAtmosphere()->Transmission(
        point3f(0), vec3f(0, 0, 1), 10);
    auto before_transmission = before_parameters->GetAtmosphere()->Transmission(
        point3f(0), vec3f(0, 0, 1), 10);
    expect_true((scaled_transmission - before_transmission).length() > 1e-5);
    gui.history_requests.push_back(-1);
    display->CommitNativeEdits(nullptr);
    expect_true((gui.base_altitude == 0 && gui.meters_per_unit == 1));
    expect_false(gui.can_undo);
    gui.history_requests.push_back(1);
    display->CommitNativeEdits(nullptr);
    expect_true((gui.base_altitude == 500 && gui.meters_per_unit == 200));
    auto redone_transmission = environment->light->GetAtmosphere()->Transmission(
        point3f(0), vec3f(0, 0, 1), 10);
    expect_true((scaled_transmission - redone_transmission).length() < 1e-6);
    display->ConsumeAtmosphereChange();

    // Reject invalid native requests without replacing any part of the light.
    const auto valid_light = environment->light;
    gui.meters_per_unit = 0;
    gui.atmosphere_parameters_pending = true;
    display->ApplyNativeSkyControls();
    expect_false(gui.sky_error.empty());
    expect_false(display->ConsumeAtmosphereChange());
    expect_true(environment->light == valid_light);
    expect_true(Rcpp::as<double>(display->export_sky()["meters_per_unit"]) == 200);
    gui.meters_per_unit = 200;
    gui.base_altitude = 15001;
    gui.atmosphere_parameters_pending = true;
    display->ApplyNativeSkyControls();
    expect_false(gui.sky_error.empty());
    expect_true(environment->light == valid_light);
    gui.base_altitude = 500;

    gui.haze = gui.altitude = 0;
    gui.atmosphere_pending = true;
    display->ApplyNativeSkyControls();
    expect_true(display->ConsumeAtmosphereChange());
    expect_true(volume->atmosphere == nullptr);
    expect_true(environment->light->GetAtmosphere() != nullptr);
    choose(0);
    expect_false(gui.has_atmosphere);
    expect_true(environment->light->GetAtmosphere() == nullptr);
    choose(1);
    expect_true((gui.has_atmosphere && !gui.haze && !gui.altitude));
    expect_true((gui.base_altitude == 500 && gui.meters_per_unit == 200));
    expect_true(volume->atmosphere == nullptr);

    gui.haze = gui.altitude = 1;
    gui.atmosphere_pending = true;
    display->ApplyNativeSkyControls();
    expect_true(display->ConsumeAtmosphereChange());
    expect_true(volume->atmosphere != nullptr);
    gui.sun_elevation = 25;
    gui.sun_azimuth = 135;
    gui.sun_pending = true;
    display->ApplyNativeSkyControls();
    expect_true((gui.sky_error.empty() && gui.manual_sun));
    choose(0);
    choose(1);
    expect_true((gui.has_atmosphere && gui.haze && gui.altitude));
    expect_true((gui.sun_elevation == 25 && gui.sun_azimuth == 135));
    expect_true(volume->atmosphere != nullptr);

    // Both model switches and repeated endpoint edits retain the actual angle.
    // Exercise the real image/native rebuild and then inspect renderer sampling.
    auto check_lighting = [&] {
      const auto light = environment->light;
      expect_true(std::isfinite(light->SamplingWeight()));
      expect_true(light->SamplingWeight() > 0);
      const auto sky = light->Radiance(point3f(0), vec3f(0, 1, 0), 0);
      expect_true(sky[0] + sky[1] + sky[2] > 0);
      for (int i = 0; i < 16; ++i) {
        const auto direction = light->Sample(point3f(0), vec2f((i + .5) / 16, .37), 0);
        const auto radiance = light->Radiance(point3f(0), direction, 0);
        const auto pdf = light->Pdf(point3f(0), direction, 0);
        expect_true((std::isfinite(pdf) && pdf > 0));
        for (int c = 0; c < 3; ++c) {
          expect_true((std::isfinite(direction[c]) && std::isfinite(radiance[c])));
        }
      }
    };
    gui.sun_elevation = 90;
    gui.sun_pending = true;
    display->ApplyNativeSkyControls();
    expect_true(gui.sky_error.empty());
    expect_true(gui.sun_elevation == MaxSunElevationDegrees);
    check_lighting();
    choose(0);
    expect_true(gui.sun_elevation == MaxSunElevationDegrees);
    check_lighting();
    choose(1);
    expect_true(gui.sun_elevation == MaxSunElevationDegrees);
    check_lighting();

    // Direct C++ atmosphere descriptions bypass the GUI and R callback. Their
    // Sun cone, radiance and importance map must still match the guarded angle.
    Rcpp::Function update = Rcpp::as<Rcpp::List>(fixture["controls"])["update"];
    Rcpp::List result = update(40, -74, "2026-06-21 16:00:00", 1, 90, 135);
    Rcpp::List native = Rcpp::as<Rcpp::List>(result["lights"])[0];
    native["haze"] = false;
    native["query_altitude"] = false;
    native["include_sun"] = true;
    native["elevation"] = 90.;
    PragueInfiniteLight zenith(native);
    native["elevation"] = MaxSunElevationDegrees;
    PragueInfiniteLight capped(native);
    expect_true(zenith.SamplingWeight() == capped.SamplingWeight());
    auto zenith_sun = zenith.Sample(point3f(0), vec2f(0, .5), 0);
    auto capped_sun = capped.Sample(point3f(0), vec2f(0, .5), 0);
    expect_true((zenith_sun - capped_sun).length() == 0);
    // A nonzero horizontal component also guards against Float rounding to up.
    expect_true(std::hypot(zenith_sun[0], zenith_sun[2]) > .001);
    expect_true(zenith.Pdf(point3f(0), zenith_sun, 0) > 0);

    // Opening an existing Prague image must also expose the atmosphere panel.
    Rcpp::List controls = Rcpp::clone(Rcpp::as<Rcpp::List>(fixture["controls"]));
    controls["model"] = 1L;
    ConfigureNativeSky(*display, environment, volume, textures, controls, lights);
    expect_true((gui.has_atmosphere && gui.sky_error.empty()));
    expect_true(volume->atmosphere != nullptr);
  }

  test_that(
      "manual sun direction updates the sky and solar disk without changing other lights") {
    Rcpp::List original = Rcpp::List::create(
        Rcpp::List::create(Rcpp::Named("type") = "prague",
                           Rcpp::Named("elevation") = 30.,
                           Rcpp::Named("azimuth") = 180.,
                           Rcpp::Named("rotation") = 25.),
        Rcpp::List::create(Rcpp::Named("type") = "disk",
                           Rcpp::Named("radiance_spectrum") = "sun",
                           Rcpp::Named("direction") =
                               Rcpp::NumericVector::create(0, 1, 0),
                           Rcpp::Named("rotation") = 0.),
        Rcpp::List::create(Rcpp::Named("type") = "disk",
                           Rcpp::Named("radiance_spectrum") = "moon",
                           Rcpp::Named("direction") =
                               Rcpp::NumericVector::create(0, 1, 0)));
    Rcpp::List updated = PreviewSunDescriptions(original, 0, 0, 90);
    Rcpp::List sky = updated[0], sun = updated[1], moon = updated[2],
               before = original[0];
    auto direction = Rcpp::as<Rcpp::NumericVector>(sun["direction"]);
    expect_true((std::abs(direction[0] + 1) < 1e-12 && std::abs(direction[1]) < 1e-12 &&
                 std::abs(direction[2]) < 1e-12));
    expect_true((Rcpp::as<double>(sky["azimuth"]) == 90 &&
                 Rcpp::as<double>(sun["rotation"]) == 25));
    expect_true((Rcpp::as<Rcpp::NumericVector>(moon["direction"])[1] == 1));
    expect_true((Rcpp::as<double>(before["azimuth"]) == 180));
    updated = PreviewSunDescriptions(original, 0, 90, 135);
    sky = updated[0];
    sun = updated[1];
    direction = Rcpp::as<Rcpp::NumericVector>(sun["direction"]);
    expect_true(Rcpp::as<double>(sky["elevation"]) == MaxSunElevationDegrees);
    const double radians = 3.14159265358979323846 / 180;
    expect_true(std::abs(direction[1] - std::sin(MaxSunElevationDegrees * radians)) <
                1e-12);
    expect_true(std::hypot(direction[0], direction[2]) > .001);
    expect_true(Rcpp::as<double>(before["elevation"]) == 30);
  }

  test_that("the left denoise control is disabled without support and queues edits") {
    struct Controls {
      bool disabled = false, saw_checkbox = false, on_left = false;
      int disabled_depth = 0;
    } controls;
    auto owner = uint64_t(reinterpret_cast<uintptr_t>(&controls));
    rayimgui_api_v1 api{};
    api.widget = [](uint64_t owner,
                    rayimgui_widget_v1* w,
                    rayimgui_item_v1* item,
                    rayimgui_error_v1*) -> int32_t {
      auto& state = *reinterpret_cast<Controls*>(uintptr_t(owner));
      item->flags = 0;
      if (w->kind == RAYIMGUI_WINDOW_BEGIN) {
        state.on_left = w->options == RAYIMGUI_DOCK_LEFT;
        if (state.on_left) {
          item->flags = RAYIMGUI_VISIBLE;
        }
      } else if (w->kind == RAYIMGUI_DISABLED_BEGIN) {
        state.disabled = *w->integers != 0;
        ++state.disabled_depth;
      } else if (w->kind == RAYIMGUI_DISABLED_END) {
        --state.disabled_depth;
      } else if (w->kind == RAYIMGUI_CHECKBOX && std::string(w->label) == "Denoise") {
        state.saw_checkbox = state.on_left;
        if (!state.disabled) {
          *w->integers = !*w->integers;
          item->flags = RAYIMGUI_CHANGED;
        }
      }
      return 0;
    };
    RayrenderGui gui;
    expect_true((gui.draw_impl(&api, owner, &gui.error) == 0));
    expect_true(controls.saw_checkbox);
    expect_true(controls.disabled);
    expect_true((controls.disabled_depth == 0));
    expect_false(gui.denoise_pending);
    expect_false(gui.denoise_enabled);
    gui.denoise_available = true;
    expect_true((gui.draw_impl(&api, owner, &gui.error) == 0));
    expect_false(controls.disabled);
    expect_true(gui.denoise_enabled);
    expect_true(gui.denoise_pending);
  }

  test_that("denoise requests cannot enable unsupported displays") {
    Transform object, world;
    camera cam(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
    auto display = NativeTestDisplay(cam, object, world);
    RayrenderGui gui;
    display->AttachNativeGui(&gui, false, false);
    expect_false(gui.denoise_available);
    gui.denoise_enabled = 1;
    gui.denoise_pending = true;
    expect_false(display->ApplyNativeControls(nullptr));
    expect_false(gui.denoise_enabled);
    expect_false(gui.denoise_pending);
  }

#ifdef HAS_OIDN
  test_that(
      "denoising follows the initial setting and toggles without camera editing") {
    Transform object, world;
    camera cam(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
    RayMatrix normal(4, 4, 3), albedo(4, 4, 3);
    RayOidnDenoiser denoiser;
    auto display = NativeTestDisplay(cam, object, world);
    display->SetDenoiser(&denoiser, &albedo, &normal, true);
    RayrenderGui gui;
    display->AttachNativeGui(&gui, false, false);
    expect_true(gui.denoise_available);
    expect_true(gui.denoise_enabled);
    display->MarkDenoisedPreviewReady(5);
    gui.selection_visible = true;
    gui.selection_revision = 7;
    gui.denoise_enabled = 0;
    gui.denoise_pending = true;
    expect_false(display->ApplyNativeControls(nullptr));
    expect_false(display->denoise);
    expect_false(display->HasDenoisedPreview());
    expect_false(gui.denoise_pending);
    gui.denoise_enabled = 1;
    gui.denoise_pending = true;
    expect_false(display->ApplyNativeControls(nullptr));
    expect_true(display->denoise);
    expect_true(gui.selection_visible);
    expect_true((gui.selection_revision == 7));
    expect_true(((cam.get_origin() - point3f(0, 0, -10)).length() == 0));
  }

  test_that(
      "live denoising initializes both filters and follows toggles in fast preview") {
    // Drive completed native frames through a mock provider while exercising the
    // real two-worker renderer and OIDN filters, including startup with denoise off.
    const size_t size = 16;
    Transform object, world_transform;
    camera cam(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 5, 0, 1, 1);
    RayMatrix rgb(size, size, 3), normal(size, size, 3), albedo(size, size, 3);
    RayMatrix alpha(size, size, 1), filtered(size, size, 3);
    RayMatrix oidn_normal(size, size, 3), oidn_albedo(size, size, 3);
    RayOidnDenoiser denoiser;
    PreviewDisplay display(size,
                           size,
                           false,
                           true,
                           true,
                           10,
                           &cam,
                           &object,
                           &world_transform,
                           &denoiser,
                           &oidn_albedo,
                           &oidn_normal,
                           false,
                           false);
    RayrenderGui gui;
    gui.width = gui.height = size;
    display.AttachNativeGui(&gui, true, true);
    expect_true(gui.denoise_available);
    expect_false(gui.denoise_enabled);
    struct Frames {
      RayrenderGui* gui;
      PreviewDisplay* display;
      RayOidnDenoiser* denoiser;
      uint64_t version = 0;
      std::vector<bool> enabled, cached, full_ready, fast;
      std::vector<size_t> samples;
    } frames{&gui, &display, &denoiser};
    rayimgui_api_v1 api{};
    api.texture_create = [](uint64_t,
                            const rayimgui_image_v1*,
                            uint64_t* handle,
                            rayimgui_error_v1*) -> int32_t {
      *handle = 1;
      return 0;
    };
    api.texture_update = [](uint64_t,
                            uint64_t,
                            const rayimgui_image_v1*,
                            rayimgui_error_v1*) -> int32_t {
      return 0;
    };
    api.step =
        [](uint64_t owner, rayimgui_step_v1* step, rayimgui_error_v1*) -> int32_t {
      auto& state = *reinterpret_cast<Frames*>(uintptr_t(owner));
      auto& gui = *state.gui;
      if (gui.version == state.version) {
        return 0;
      }
      state.version = gui.version;
      state.enabled.push_back(state.display->denoise);
      state.cached.push_back(state.display->HasDenoisedPreview());
      state.full_ready.push_back(state.denoiser->Ready());
      state.fast.push_back(state.display->write_fast_output);
      state.samples.push_back(gui.samples);
      if (state.version == 1 || state.version == 4) {
        gui.fast_preview = state.version == 1;
        gui.fast_pending = true;
      }
      if (state.version == 1 || state.version == 3 || state.version == 6) {
        gui.denoise_enabled = 1;
        gui.denoise_pending = true;
      } else if (state.version == 2 || state.version == 5) {
        gui.denoise_enabled = 0;
        gui.denoise_pending = true;
      } else if (state.version >= 7) {
        step->close_requested = 1;
      }
      return 0;
    };
    gui.api = &api;
    gui.session = uint64_t(reinterpret_cast<uintptr_t>(&frames));
    hitable_list world, lights;
    random_gen rng(123);
    pathtracer(2,
               size,
               size,
               8,
               0,
               0,
               1,
               rgb,
               normal,
               albedo,
               alpha,
               filtered,
               false,
               0,
               1,
               1,
               false,
               &cam,
               60,
               world,
               lights,
               10,
               2,
               2,
               display,
               IntegratorType::ShadowRays,
               &rng);
    expect_true((frames.enabled ==
                 std::vector<bool>{false, true, false, true, true, false, true}));
    expect_true((frames.cached == frames.enabled));
    expect_true((frames.full_ready ==
                 std::vector<bool>{false, false, false, false, true, true, true}));
    expect_true((frames.fast ==
                 std::vector<bool>{false, true, true, true, false, false, false}));
    expect_true((frames.samples.size() == 7));
    expect_true((frames.samples[2] == frames.samples[1] + 1));
    expect_true((frames.samples[3] == frames.samples[2] + 1));
    expect_true((frames.samples[5] == frames.samples[4] + 1));
    expect_true((frames.samples[6] == frames.samples[5] + 1));
    expect_true(
        (std::all_of(filtered.data.begin(), filtered.data.end(), [](float value) {
          return std::isfinite(value);
        })));
  }
#endif
}

namespace {
Rcpp::List CameraRigFixture() {
  Rcpp::Function parse = Rcpp::Environment::base_env()["parse"];
  Rcpp::Function eval = Rcpp::Environment::base_env()["eval"];
  return eval(parse(Rcpp::_["text"] = R"(
    local({
      path = system.file('extdata', 'dgauss.50mm.txt', package='rayrender')
      list(nx=16,ny=16,fov=60,lookfrom=c(0,0,-10),lookat=c(0,0,0),
           camera_up=c(0,1,0),aperture=0,focal_distance=10,
           shutteropen=0,shutterclose=1,iso=1,film_size=.022,camera_scale=1,
           ortho_dimensions=c(8,6),real_camera_info=matrix(numeric(),0,4),
           preview_lenses=list(list(name='50 mm',source=path,
             data=as.matrix(utils::read.delim(path,header=FALSE,comment.char='#')))))
    })
  )"));
}
}

context("Editor camera projections") {
  test_that("camera poles clamp by default and free rotation transports the whole frame") {
    PreviewCameraRig rig(CameraRigFixture());
    Transform object, world;
    auto display = NativeTestDisplay(*rig.Active(), object, world);
    display->camera_rig = &rig;
    RayrenderGui gui;
    display->AttachNativeGui(&gui, true, false);
    bool clamped = true, reversible = true, continuous = true, local_axes = true;
    const point3f origin(0, 0, -10), target(0);
    const vec3f up(0, 1, 0);
    for (int model = 0; model < 4; ++model) {
      gui.camera.Select(model);
      display->CommitNativeEdits(nullptr);
      auto &cam = *display->cam;
      cam.set_free_rotation(false);
      for (Float sign : {Float(-1), Float(1)}) {
        cam.update_pose_absolute(origin, target, up);
        // Hold Q/Z past both poles, then keep holding to check that the stop
        // does not oscillate or accumulate drift. Reversing must move away.
        for (int step = 0; step < 80; ++step) {
          const auto right = unit_vector(cam.get_u());
          cam.update_position(sign * .5f * cam.get_v(), true);
          clamped &= dot(right, unit_vector(cam.get_u())) > .999f;
          clamped &= std::abs(cam.get_u().length() - 1) < 1e-5;
          clamped &= std::abs(cam.get_v().length() - 1) < 1e-5;
          clamped &= std::abs((cam.get_origin() - target).length() - 10) < 1e-4;
        }
        const auto stopped = cam.get_origin();
        clamped &= std::abs(dot(cam.get_w(), up)) > .99999f;
        for (int step = 0; step < 30; ++step) {
          cam.update_position(sign * .5f * cam.get_v(), true);
        }
        clamped &= (cam.get_origin() - stopped).length() < 1e-5;
        cam.update_position(-sign * .5f * cam.get_v(), true);
        reversible &= std::abs(dot(cam.get_w(), up)) < .999f;

        cam.update_pose_absolute(origin, target, up);
        cam.rotate_forward(sign * 120);
        const auto pitch_stop = cam.get_w();
        cam.rotate_forward(sign * 5);
        clamped &= (cam.get_w() - pitch_stop).length() < 1e-5;
        cam.rotate_forward(-sign * 5);
        reversible &= std::abs(dot(cam.get_w(), up)) < .999f;
      }

      cam.set_free_rotation(true);
      cam.update_pose_absolute(origin, target, up);
      const Float step_size = 10 * std::tan(Float(5 * M_PI / 180));
      for (int step = 0; step < 72; ++step) {
        const auto previous_up = unit_vector(cam.get_v());
        cam.update_position(step_size * cam.get_v(), true);
        continuous &= dot(previous_up, unit_vector(cam.get_v())) > .99f;
        if (step == 35) {
          continuous &= cam.get_origin()[2] > 9.99f && cam.get_v()[1] < -.99f;
        }
      }
      continuous &= (cam.get_origin() - origin).length() < 5e-4;
      continuous &= (unit_vector(cam.get_v()) - up).length() < 5e-5;

      // Rotating the entire starting scene must rotate the result identically.
      // Mixing orbit, pitch, and roll catches a hidden world-up reset after Shift.
      const Transform global = Rotate(37, unit_vector(vec3f(1, 2, 3)));
      auto gesture = [&] {
        for (int step = 0; step < 20; ++step) {
          cam.update_position(.5f * unit_vector(cam.get_v()), true);
          cam.update_position(.3f * unit_vector(cam.get_u()), true);
          cam.rotate_forward(7);
          cam.rotate_up(13);
        }
      };
      cam.update_pose_absolute(origin, target, up);
      gesture();
      const auto expected_origin = global(cam.get_origin());
      const auto expected_forward = global(cam.get_w());
      const auto expected_up = global(unit_vector(cam.get_v()));
      cam.update_pose_absolute(global(origin), global(target), global(up));
      gesture();
      local_axes &= (cam.get_origin() - expected_origin).length() < 1e-3;
      local_axes &= (cam.get_w() - expected_forward).length() < 1e-4;
      local_axes &= (unit_vector(cam.get_v()) - expected_up).length() < 1e-4;
    }
    expect_true(clamped);
    expect_true(reversible);
    expect_true(continuous);
    expect_true(local_axes);
  }

  test_that("free rotation survives camera changes, history, and export") {
    Rcpp::List info = CameraRigFixture();
    info["free_rotation"] = true;
    PreviewCameraRig rig(info);
    Transform object, world;
    auto display = NativeTestDisplay(*rig.Active(), object, world);
    display->camera_rig = &rig;
    RayrenderGui gui;
    display->AttachNativeGui(&gui, true, false);
    display->BeginNativeHistory();
    expect_true((gui.free_rotation && display->cam->get_free_rotation()));
    gui.free_rotation = 0;
    gui.history_dirty = true;
    ++gui.history_epoch;
    display->CommitNativeEdits(nullptr);
    expect_false(display->cam->get_free_rotation());
    gui.history_requests.push_back(-1);
    display->CommitNativeEdits(nullptr);
    expect_true((gui.free_rotation && display->cam->get_free_rotation()));
    gui.history_requests.push_back(1);
    display->CommitNativeEdits(nullptr);
    expect_false((gui.free_rotation || display->cam->get_free_rotation()));
    gui.free_rotation = 1;
    gui.history_dirty = true;
    ++gui.history_epoch;
    display->CommitNativeEdits(nullptr);
    for (int model : {1, 2, 3, 0}) {
      gui.camera.Select(model);
      display->CommitNativeEdits(nullptr);
      expect_true(display->cam->get_free_rotation());
    }
    expect_true(Rcpp::as<std::string>(display->NativeEditorState()["camera_rotation"]) == "free");
  }

  test_that(
      "projection changes, FOV zero keyframes and depth survive undo and export") {
    using Inputs = PreviewCameraInputs;
    PreviewCameraRig rig(CameraRigFixture());
    Transform object, world;
    auto display = NativeTestDisplay(*rig.Active(), object, world);
    display->camera_rig = &rig;
    RayrenderGui gui;
    display->AttachNativeGui(&gui, true, false);
    display->BeginNativeHistory();
    auto* perspective = display->cam;
    display->SaveCurrentKeyframe(0);
    gui.camera.values[Inputs::Fov][0] = 0;
    gui.camera.pending = true;
    ++gui.history_epoch;
    expect_true(display->CommitNativeEdits(nullptr));
    auto* orthographic = display->cam;
    expect_true((orthographic != perspective && orthographic->get_fov() == 0));
    expect_true((gui.camera.model == 1 && gui.camera.Visible(Inputs::Ortho)));
    display->SaveCurrentKeyframe(0);
    gui.camera.Select(2);
    ++gui.history_epoch;
    expect_true(display->CommitNativeEdits(nullptr));
    expect_true(display->cam->get_fov() == 360);
    gui.history_requests.push_back(-1);
    display->CommitNativeEdits(nullptr);
    expect_true(display->cam == orthographic);
    gui.history_requests.push_back(1);
    display->CommitNativeEdits(nullptr);
    expect_true(display->cam->get_fov() == 360);
    Float angle = 0;
    display->ApplyKeyframe(0, &angle);
    expect_true(display->cam == perspective);
    display->ApplyKeyframe(1, &angle);
    expect_true(display->cam == orthographic);
    display->SyncNativeCameraControls(true);
    gui.camera.Select(0);
    expect_true(
        gui.camera
            .Validate()); // A zero-focus projection supplies a new focus distance.
    display->CommitNativeEdits(nullptr);
    ++gui.history_epoch;
    gui.max_depth = 7;
    gui.depth_pending = true;
    expect_true(display->CommitNativeEdits(nullptr));
    expect_true(display->max_depth == 7);
    expect_true(Rcpp::as<int>(display->NativeEditorState()["max_depth"]) == 7);
    gui.history_requests.push_back(-1);
    display->CommitNativeEdits(nullptr);
    expect_true((display->max_depth == 50 && gui.max_depth == 50));
    gui.history_requests.push_back(1);
    display->CommitNativeEdits(nullptr);
    expect_true((display->max_depth == 7 && gui.max_depth == 7));
  }

  test_that("realistic image orientation and navigation agree with perspective") {
    PreviewCameraRig rig(CameraRigFixture());
    Transform object, world;
    auto display = NativeTestDisplay(*rig.Active(), object, world);
    display->camera_rig = &rig;
    RayrenderGui gui;
    display->AttachNativeGui(&gui, true, false);
    gui.camera.Select(3);
    expect_true(display->CommitNativeEdits(nullptr));
    auto *lens = display->cam;
    bool axes_match = true, image_matches = true, movement_matches = true;
    for (point3f origin : {point3f(0, 0, -10), point3f(4, 3, -8)}) {
      const point3f target(0);
      const vec3f up(0, 1, 0);
      camera reference(origin, target, up, 60, 1, 0, 10, 0, 1, 1);
      lens->update_pose_absolute(origin, target, up);
      axes_match &= dot(lens->get_w(), reference.get_w()) > 0;
      axes_match &= dot(lens->get_u(), reference.get_u()) > 0;
      axes_match &= dot(lens->get_v(), reference.get_v()) > 0;

      // Exercise the film coordinates used by rendering and viewport picking.
      // An optical lens changes the field of view, but each quadrant must stay
      // on the same side of the forward axis, including during motion blur.
      lens->set_camera_motion_blur_range(origin, target, up, 10, origin, target, up, 10);
      for (bool moving : {false, true}) {
        lens->set_camera_motion_blur(moving);
        for (Float u : {.35f, .65f}) {
          for (Float v : {.35f, .65f}) {
            Ray optical;
            const Float weight = lens->GenerateRay(
                CameraSample(point2f(1 - u, 1 - v), point2f(.5, .5), .5), &optical);
            const Ray pinhole = reference.get_ray(u, v, point3f(0), .5);
            image_matches &= weight > 0;
            if (weight > 0) {
              image_matches &= dot(optical.d, reference.get_w()) > 0;
              image_matches &=
                  dot(optical.d, reference.get_u()) * dot(pinhole.d, reference.get_u()) > 0;
              image_matches &=
                  dot(optical.d, reference.get_v()) * dot(pinhole.d, reference.get_v()) > 0;
            }
          }
        }
      }
      lens->set_camera_motion_blur(false);

      // Drive the actual GUI key handler rather than invoking camera setters.
      // Reset the pose for each key so a backwards W cannot be hidden by S.
      gui.orbit = false;
      for (const auto &key : {std::make_pair(RAYIMGUI_KEY_W, reference.get_w()),
                              std::make_pair(RAYIMGUI_KEY_S, -reference.get_w()),
                              std::make_pair(RAYIMGUI_KEY_A, reference.get_u()),
                              std::make_pair(RAYIMGUI_KEY_D, -reference.get_u())}) {
        lens->update_pose_absolute(origin, target, up);
        gui.keys.push_back({key.first, 0, 1});
        display->ApplyNativeControls(nullptr);
        movement_matches &= dot(lens->get_origin() - origin, key.second) > 0;
      }
      lens->update_pose_absolute(origin, target, up);
      gui.orbit = true;
      gui.keys.push_back({RAYIMGUI_KEY_W, 0, 1});
      display->ApplyNativeControls(nullptr);
      movement_matches &= (lens->get_origin() - target).length() < (origin - target).length();
      movement_matches &= dot(lens->get_w(), reference.get_w()) > 0;
    }
    expect_true(axes_match);
    expect_true(image_matches);
    expect_true(movement_matches);
  }

  test_that("realistic lenses have editable optics and atomic pose restoration") {
    using Inputs = PreviewCameraInputs;
    PreviewCameraRig rig(CameraRigFixture());
    Transform object, world;
    auto display = NativeTestDisplay(*rig.Active(), object, world);
    display->camera_rig = &rig;
    RayrenderGui gui;
    display->AttachNativeGui(&gui, true, false);
    display->BeginNativeHistory();
    gui.camera.Select(3);
    ++gui.history_epoch;
    expect_true(display->CommitNativeEdits(nullptr));
    expect_true((display->cam->get_fov() == -1 && gui.camera.editable));
    expect_true((gui.camera.Visible(Inputs::Film) && !gui.camera.Visible(Inputs::Fov)));
    auto saved = display->CreateCurrentKeyframe(0);
    const double aperture = display->cam->get_aperture();
    gui.camera.values[Inputs::Aperture][0] = aperture / 2;
    gui.camera.values[Inputs::Focus][0] = 12;
    gui.camera.values[Inputs::Position] = {
        0, 0, 0}; // Previous target: no singular intermediate LookAt.
    gui.camera.values[Inputs::Target] = {0, 0, 10};
    gui.camera.pending = true;
    ++gui.history_epoch;
    expect_true(display->CommitNativeEdits(nullptr));
    expect_true(display->cam->get_aperture() == Approx(aperture / 2));
    expect_true(display->cam->get_focal_distance() == 12);
    expect_true(display->cam->get_origin()[2] == 0);
    CameraSample sample(point2f(.5, .5), point2f(.5, .5), .5);
    Ray ray;
    expect_true(display->cam->GenerateRay(sample, &ray) > 0);
    expect_true((std::isfinite(ray.direction()[2]) && ray.direction()[2] > 0));
    gui.history_requests.push_back(-1);
    display->CommitNativeEdits(nullptr);
    expect_true(display->cam->get_origin()[2] == -10);
    expect_true(display->cam->get_aperture() == Approx(aperture));
    // A lens cannot focus arbitrarily close. Rejected optics retain the live view.
    auto* previous = display->cam;
    gui.camera.values[Inputs::Focus][0] = .001;
    gui.camera.pending = true;
    display->CommitNativeEdits(nullptr);
    expect_true(display->cam == previous);
    expect_false(gui.camera.optical_error.empty());
    display->SyncNativeCameraControls(true);
    expect_true(gui.camera.Valid());
    gui.camera.Select(2);
    display->CommitNativeEdits(nullptr);
    Float angle = 0;
    display->ApplyCameraState(saved, &angle);
    expect_true(display->cam->get_fov() == -1);
    expect_true(
        Rcpp::as<std::string>(display->NativeEditorState()["camera_description_file"])
            .find("dgauss.50mm") != std::string::npos);
  }
}

namespace {
// Deterministic optical throughput separates GenerateRay from get_ray without
// making the render-loop regression depend on statistical lens acceptance.
class CameraSwitchProbe : public camera {
public:
  explicit CameraSwitchProbe(bool physical)
      : camera(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0), 60, 1, 0, 10, 0, 1, 1),
        physical(physical) {
  }
  Float get_fov() override {
    return physical ? -1 : 60;
  }
  Float GenerateRay(const CameraSample&, Ray* ray) const override {
    *ray = Ray(point3f(0, 0, -10), vec3f(0, 0, 1));
    return 2;
  }

private:
  bool physical;
};

// Emit one unit per bounce and continue specularly. The pixel value is exactly
// max_depth, so this checks the value reaching the integrator, not just the UI.
class DepthProbeMaterial : public material {
public:
  bool scatter(const Ray& ray, const hit_record&, scatter_record& out,
               Sampler*) override {
    out.is_specular = true;
    out.specular_ray = ray;
    out.attenuation = point3f(1);
    return true;
  }
  point3f emitted(const Ray&, const hit_record&, Float, Float, const point3f&,
                  bool&) override {
    return point3f(1);
  }
  const std::string GetName() override {
    return "depth probe";
  }
  size_t GetSize() override {
    return sizeof(*this);
  }
};
class DepthProbeSurface : public hitable {
public:
  DepthProbeMaterial probe;
  const bool Hit(const Ray& ray, hit_record& rec) const {
    rec = hit_record{};
    rec.t = 1;
    rec.p = ray(1);
    rec.normal = normal3f(0, 0, 1);
    rec.mat_ptr = const_cast<DepthProbeMaterial*>(&probe);
    return true;
  }
  const bool hit(const Ray& ray, Float, Float, hit_record& rec,
                 random_gen&) const override {
    return Hit(ray, rec);
  }
  const bool hit(const Ray& ray, Float, Float, hit_record& rec,
                 Sampler*) const override {
    return Hit(ray, rec);
  }
  bool bounding_box(Float, Float, aabb&) const override {
    return false;
  }
  std::string GetName() const override {
    return "depth probe";
  }
  size_t GetSize() override {
    return sizeof(*this);
  }
  void hitable_info_bounds(Float, Float) const override {
  }
};
}

context("Live camera and depth rendering") {
  test_that("full and Fast workers use the current optical camera and path depth") {
    const size_t size = 16;
    CameraSwitchProbe perspective(false), physical(true);
    Transform object, transform;
    auto display = NativeTestDisplay(perspective, object, transform);
    RayrenderGui gui;
    gui.width = gui.height = size;
    display->AttachNativeGui(&gui, true, true);
    struct Frames {
      PreviewDisplay* display;
      RayrenderGui* gui;
      RayCamera* physical;
      uint64_t version = 0;
      std::vector<Float> pixels;
    } frames{display.get(), &gui, &physical};
    rayimgui_api_v1 api{};
    api.texture_create = [](uint64_t,
                            const rayimgui_image_v1*,
                            uint64_t* texture,
                            rayimgui_error_v1*) -> int32_t {
      *texture = 1;
      return 0;
    };
    api.texture_update = [](uint64_t,
                            uint64_t,
                            const rayimgui_image_v1*,
                            rayimgui_error_v1*) -> int32_t {
      return 0;
    };
    api.step =
        [](uint64_t owner, rayimgui_step_v1* step, rayimgui_error_v1*) -> int32_t {
      auto& state = *reinterpret_cast<Frames*>(uintptr_t(owner));
      if (state.version == state.gui->version) {
        return 0;
      }
      state.version = state.gui->version;
      state.pixels.push_back(state.gui->display_rgb[0]);
      if (state.version == 2) {
        state.display->SetCamera(state.physical);
        state.gui->max_depth = 3;
        state.gui->depth_pending = true;
      } else if (state.version == 4) {
        state.gui->fast_preview = 1;
        state.gui->fast_pending = true;
      } else if (state.version == 6) {
        step->close_requested = 1;
      }
      return 0;
    };
    gui.api = &api;
    gui.session = uint64_t(reinterpret_cast<uintptr_t>(&frames));
    hitable_list world, lights;
    world.add(std::make_shared<DepthProbeSurface>());
    RayMatrix rgb(size, size, 3), normal(size, size, 3), albedo(size, size, 3);
    RayMatrix alpha(size, size, 1), filtered(size, size, 3);
    random_gen rng(123);
    pathtracer(2,
               size,
               size,
               8,
               0,
               0,
               1,
               rgb,
               normal,
               albedo,
               alpha,
               filtered,
               false,
               0,
               1,
               1,
               false,
               &perspective,
               60,
               world,
               lights,
               100,
               1,
               10,
               *display,
               IntegratorType::Basic,
               &rng);
    expect_true(frames.pixels.size() == 6);
    PreviewColorTransform mapping;
    bool correct = frames.pixels.size() == 6;
    for (size_t i = 0; i < frames.pixels.size(); ++i) {
      const double value = i < 2 ? 1 : 6;
      const double expected = mapping.Apply({value, value, value})[0];
      correct &= std::abs(frames.pixels[i] - expected) < 1e-6;
    }
    expect_true(correct);
  }
}
#endif
