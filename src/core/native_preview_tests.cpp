#ifdef NOT_CRAN
#include "PreviewDisplay.h"
#include "preview_sky.h"
#include "preview_sky_controls.h"
#include "../materials/texturecache.h"
#include "integrator.h"
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
    if (gui.api->header.abi_minor < 8) {
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
    expect_true(action(Action::Play));
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
      "sky changes wait for release and a safe checkpoint; invalid dates keep rendering") {
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
    expect_true((updates == 0));
    expect_false(display->ConsumeAtmosphereChange());
    gui.sun_elevation = 50;
    gui.sun_editing = false;
    display->ApplyNativeSkyControls();
    expect_true((updates == 1 && elevation == 50));
    expect_true((display->ConsumeAtmosphereChange()));
    display->ApplyNativeSkyControls();
    expect_true((updates == 1));
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
          update=function(latitude,longitude,datetime,model,elevation,azimuth) {
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
#endif
