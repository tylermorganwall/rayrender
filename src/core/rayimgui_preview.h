/* Copyright (c) 2026 Tyler Morgan-Wall. MIT; LICENSE.protocol.
 * Included once, at the end of PreviewDisplay.cpp. These checkpoints run on
 * R's main thread after wait_for_render_jobs has drained all sample tasks.
 */
#include "preview_scene.h"
#include "preview_history.h"

void PreviewDisplay::AttachNativeGui(RayrenderGui* gui, bool edit, bool deferred) {
  gui->max_depth = int32_t(max_depth);
  gui->free_rotation = cam->get_free_rotation();
  native_gui = gui;
  native_history.reset();
  native_drag_fast = false;
  native_fast_saved = write_fast_output;
  native_sun_preview = false;
  gui->animation_paused = false;
  gui->can_undo = gui->can_redo = gui->history_dirty = false;
  gui->history_requests.clear();
  preview = true;
  interactive = edit;
  deferred_render = edit && deferred;
  render_requested = !deferred_render;
  gui->can_edit = edit;
  gui->deferred = deferred_render;
  gui->exposure = preview_exposure_adjustment;
  SyncNativeCameraControls(true);
  SyncNativeAnimationState();
#ifdef HAS_OIDN
  gui->denoise_available = denoiser != nullptr && oidn_albedo_output != nullptr &&
                           oidn_normal_output != nullptr;
  gui->denoise_enabled = denoise;
#endif
  native_base_step =
      std::max(Float(.001), (cam->get_origin() - cam->get_lookat()).length() / 20);
}

void PreviewDisplay::SetSunPosition(double elevation, double azimuth) {
  elevation = ClampSunElevation(elevation);
  sun_elevation = elevation;
  sun_azimuth = azimuth;
  if (native_gui) {
    native_gui->sun_elevation = elevation;
    native_gui->sun_azimuth = azimuth;
  }
}

void PreviewDisplay::SetSkyModelControls(int model,
                                         std::function<std::string(int)> update) {
  sky_model = model;
  update_sky_model = std::move(update);
  if (native_gui) {
    native_gui->has_sky_model = bool(update_sky_model);
    native_gui->sky_model = model;
  }
}

void PreviewDisplay::SetSunControls(double elevation, double azimuth,
                                    std::function<void(double, double)> update) {
  SetSunPosition(elevation, azimuth);
  update_sun = std::move(update);
}

void PreviewDisplay::SetSkyControls(
    double latitude, double longitude, const std::string& datetime,
    std::function<std::string(double, double, const std::string&)> update) {
  update_sky = std::move(update);
  if (native_gui) {
    native_gui->has_location = true;
    native_gui->latitude = latitude;
    native_gui->longitude = longitude;
    native_gui->set_datetime(datetime);
  }
}

// Publish committed animation settings without overwriting a pending widget edit.
void PreviewDisplay::SyncNativeAnimationState() {
  if (!native_gui) {
    return;
  }
  auto& gui = *native_gui;
  gui.animation_current = current_keyframe;
  gui.animation_playing = IsPreviewMotionActive();
  if (!gui.animation_playing) {
    gui.animation_paused = false;
  }
  const std::string type = keyframe_motion_args.containsElementNamed("type")
                               ? Rcpp::as<std::string>(keyframe_motion_args["type"])
                               : "spline";
  gui.animation_timing_available = type == "spline" || type == "linear" ||
                                   type == "quad" || type == "cubic" || type == "exp";
  // Each transition starts at 30 frame intervals. An explicitly supplied total
  // can override that allocation until the first individual duration edit.
  if (!gui.animation_timing_custom && Keyframes.size() > 1 &&
      gui.keyframe_snapshots.size() == Keyframes.size()) {
    const size_t segments = Keyframes.size() - (KeyframeMotionClosed() ? 0 : 1);
    double total = segments * 30 + 1;
    if (keyframe_motion_args.containsElementNamed("frames")) {
      const double supplied = Rcpp::as<double>(keyframe_motion_args["frames"]);
      if (std::isfinite(supplied)) {
        total = supplied;
      }
    }
    const size_t intervals =
        size_t(std::clamp(total - 1,
                          double(segments),
                          double(segments * RayrenderGui::MaxSegmentFrames)));
    for (size_t i = 0; i < segments; ++i) {
      gui.keyframe_snapshots[i].frames_to_next =
          int32_t(intervals / segments + (i < intervals % segments));
    }
    if (gui.animation_timing_available &&
        keyframe_motion_args.containsElementNamed("segment_frames")) {
      Rcpp::NumericVector supplied = keyframe_motion_args["segment_frames"];
      const bool valid =
          supplied.size() == R_xlen_t(segments) &&
          std::all_of(supplied.begin(), supplied.end(), [](double value) {
            return std::isfinite(value) && value == std::floor(value) && value >= 1 &&
                   value <= RayrenderGui::MaxSegmentFrames;
          });
      if (valid) {
        for (size_t i = 0; i < segments; ++i) {
          gui.keyframe_snapshots[i].frames_to_next = int32_t(supplied[i]);
        }
        gui.animation_timing_custom = true;
      }
    }
  }
  if (!gui.animation_settings_pending) {
    gui.animation_blur = CameraMotionBlurEnabled();
    gui.animation_closed = KeyframeMotionClosed();
    gui.animation_shutter =
        std::isfinite(GetShutterSpeed()) ? 1.0 / GetShutterSpeed() : 0;
  }
}

bool PreviewDisplay::ApplyNativeAnimationControls() {
  auto& gui = *native_gui;
  using Action = RayrenderGui::AnimationAction;
  const Action action = gui.animation_action;
  gui.animation_action = Action::None;
  bool reset = false;
  if (!interactive) {
    gui.animation_settings_pending = false;
    SyncNativeAnimationState();
    return false;
  }
  if (gui.animation_settings_pending) {
    if (!IsPreviewMotionActive()) {
      const bool blur = gui.animation_blur != 0;
      const double interval = std::clamp(gui.animation_shutter, 0.0, 1.0);
      const Float speed = interval > 0 ? Float(1.0 / interval) : Infinity;
      reset = blur != CameraMotionBlurEnabled() || speed != GetShutterSpeed();
      SetCameraMotionBlur(blur);
      SetShutterSpeed(speed);
      keyframe_motion_closed = gui.animation_closed != 0;
    }
    gui.animation_settings_pending = false;
  }
  if (action == Action::Play) {
    gui.animation_message.clear();
    if (IsPreviewMotionActive()) {
      // Pausing keeps both the next frame index and the original restore pose.
      // The normal reset below freezes the shutter at the displayed camera.
      gui.animation_paused = !gui.animation_paused;
      reset = true;
    } else {
      gui.animation_paused = false;
      const bool started = StartPreviewMotion(native_env_angle);
      if (!started) {
        gui.animation_message =
            "Could not play camera path; check keyframe motion settings in the R console.";
      }
      reset = started || reset;
    }
  } else if (action == Action::Stop) {
    reset = CancelPreviewMotion(&native_env_angle) || reset;
    gui.animation_paused = false;
  } else if (!IsPreviewMotionActive()) {
    switch (action) {
    case Action::Save:
      SaveCurrentKeyframe(native_env_angle);
      break;
    case Action::Previous:
      reset = JumpKeyframe(-1, &native_env_angle) || reset;
      break;
    case Action::Next:
      reset = JumpKeyframe(1, &native_env_angle) || reset;
      break;
    case Action::Delete:
      reset = DeleteCurrentKeyframe(&native_env_angle) || reset;
      break;
    case Action::Replace:
      // Opening a context menu must not apply its saved camera first. Capture
      // the live camera here, after renderer workers have finished their sample.
      for (size_t i = 0; i < gui.keyframe_snapshots.size(); ++i) {
        if (gui.keyframe_snapshots[i].id == gui.animation_target &&
            i < Keyframes.size()) {
          auto keyframe = CreateCurrentKeyframe(native_env_angle);
          gui.replace_keyframe_snapshot(i, NativeKeyframeCamera(keyframe));
          Keyframes[i] = keyframe;
          current_keyframe = int(i);
          gui.animation_message.clear();
          break;
        }
      }
      break;
    case Action::Select:
      // Stable IDs prevent a queued thumbnail click selecting a different frame
      // if a deletion changes vector indices before the request is consumed.
      for (size_t i = 0; i < gui.keyframe_snapshots.size(); ++i) {
        if (gui.keyframe_snapshots[i].id == gui.animation_target) {
          reset = ApplyKeyframe(int(i), &native_env_angle) || reset;
          break;
        }
      }
      break;
    default:
      break;
    }
  }
  return reset;
}

// Synchronize the pane after navigation and keyframe changes without overwriting
// an active text edit or a rejected draft the user still needs to correct.
void PreviewDisplay::SyncNativeCameraControls(bool force) {
  if (!native_gui || !cam) {
    return;
  }
  auto& input = native_gui->camera;
  const double fov = cam->get_fov();
  input.available = true;
  input.editable = interactive && !IsPreviewMotionActive();
  if (!force && (input.pending || input.editing || !input.Valid())) {
    return;
  }
  input.model = camera_rig ? camera_rig->model : PreviewCameraRig::ModelForFov(fov);
  input.models.resize(3);
  input.lens_apertures.clear();
  if (camera_rig) {
    for (const auto& lens : camera_rig->lenses) {
      input.models.push_back(lens.name);
      input.lens_apertures.push_back(lens.aperture);
    }
  }
  input.projection = fov < 0      ? PreviewCameraInputs::Realistic
                     : fov == 0   ? PreviewCameraInputs::Orthographic
                     : fov >= 180 ? PreviewCameraInputs::Environment
                                  : PreviewCameraInputs::Perspective;

  input.Read(CreateCurrentKeyframe(native_env_angle));
  input.input_errors.fill({});
  if (force) {
    input.optical_error.clear();
  }
  input.Validate();
  if (force) {
    input.pending = input.editing = false;
  }
}

bool PreviewDisplay::ApplyNativeCameraControls() {
  auto& input = native_gui->camera;
  if (!input.pending) {
    return false;
  }
  // Film size and scene-to-lens scale rebuild the sampled optical bounds.
  // Retain the last draft through a drag and prepare it once on release.
  if (input.editing && camera_rig &&
      input.projection == PreviewCameraInputs::Realistic &&
      (input.values[PreviewCameraInputs::Film][0] != camera_rig->film_size ||
       input.values[PreviewCameraInputs::Scale][0] != camera_rig->camera_scale)) {
    return false;
  }
  input.pending = false;
  if (!interactive || IsPreviewMotionActive() || !input.editable || !input.Validate()) {
    return false;
  }
  // Preserve exposure/environment and parameters unsupported by this camera.
  // Validate the complete pose before any setter can change the live camera.
  auto state = CreateCurrentKeyframe(native_env_angle);
  input.Write(state);
  try {
    ApplyCameraState(state, &native_env_angle);
  } catch (const std::exception& error) {
    input.optical_error = error.what();
    input.Validate();
    return false;
  }
  ApplyStaticPreviewCameraMotionRange(cam);
  // FOV zero switches the rig immediately; the next pane shows orthographic size.
  if (input.model != (camera_rig ? camera_rig->model : input.model)) {
    SyncNativeCameraControls(true);
  }
  return true;
}

// Replay queued input through the existing camera controls. Returning true asks
// the caller to discard accumulated samples after the camera or preview mode changes.
bool PreviewDisplay::ApplyNativeControls(hitable* world) {
  RayrenderGui& gui = *native_gui;
  bool stopped = gui.animation_action == RayrenderGui::AnimationAction::Stop;
  bool reset = ApplyNativeAnimationControls();
  // Denoising changes display/output processing, so it is available even when
  // camera editing is off and does not discard accumulated samples or masks.
  if (gui.denoise_pending) {
#ifdef HAS_OIDN
    const bool enabled = gui.denoise_enabled != 0 && gui.denoise_available;
    if (denoise != enabled) {
      denoise = enabled;
      has_denoised_preview = false;
      denoised_preview_sample_count = 0;
    }
    gui.denoise_enabled = denoise;
#else
    gui.denoise_enabled = 0;
#endif
    gui.denoise_pending = false;
  }
  if (gui.depth_pending) {
    if (gui.max_depth >= 1 && gui.max_depth <= 10000) {
      reset = max_depth != size_t(gui.max_depth) || reset;
      max_depth = size_t(gui.max_depth);
    }
    gui.depth_pending = false;
  }
  if (!interactive) {
    gui.keys.clear();
    gui.pick_pending = false;
    return reset;
  }

  if (gui.fast_pending) {
    if (native_drag_fast) {
      native_fast_saved = gui.fast_preview != 0;
    } else {
      reset = write_fast_output != (gui.fast_preview != 0) || reset;
      write_fast_output = gui.fast_preview != 0;
    }
    gui.fast_pending = false;
  }

  // This navigation preference changes no pixels until the camera moves.
  cam->set_free_rotation(gui.free_rotation != 0);

  for (const auto& key : gui.keys) {
    for (unsigned count = 0; count < key.count; ++count) {
      const bool shift = (key.modifiers & RAYIMGUI_SHIFT) != 0;
      const Float speed =
          static_cast<Float>(std::clamp(gui.movement_speed, .001, 128.0));
      // Manual viewport input takes control from playback. Restore the same
      // camera as Stop before applying navigation; M and F retain their own
      // playback/quality behavior. This also covers a paused or looping path.
      if (IsPreviewMotionActive() && key.code != RAYIMGUI_KEY_M &&
          key.code != RAYIMGUI_KEY_F) {
        gui.animation_action = RayrenderGui::AnimationAction::Stop;
        reset = ApplyNativeAnimationControls() || reset;
        stopped = true;
      }
      switch (key.code) {
      case RAYIMGUI_KEY_W:
        if (shift) {
          cam->rotate_forward(speed);
        } else if (!gui.orbit || (cam->get_origin() - cam->get_lookat()).length() >
                                     speed * native_base_step) {
          cam->update_position(
              speed * cam->get_w() * native_base_step, gui.orbit, false);
        }
        reset = true;
        break;
      case RAYIMGUI_KEY_S:
        if (shift) {
          cam->rotate_forward(-speed);
        } else {
          cam->update_position(
              -speed * cam->get_w() * native_base_step, gui.orbit, false);
        }
        reset = true;
        break;
      case RAYIMGUI_KEY_A:
        if (shift) {
          cam->rotate_up(-speed);
        } else {
          cam->update_position(speed * cam->get_u() * native_base_step, gui.orbit);
        }
        reset = true;
        break;
      case RAYIMGUI_KEY_D:
        if (shift) {
          cam->rotate_up(speed);
        } else {
          cam->update_position(-speed * cam->get_u() * native_base_step, gui.orbit);
        }
        reset = true;
        break;
      case RAYIMGUI_KEY_Q:
        cam->update_position(speed * cam->get_v() * native_base_step, gui.orbit);
        reset = true;
        break;
      case RAYIMGUI_KEY_Z:
        cam->update_position(-speed * cam->get_v() * native_base_step, gui.orbit);
        reset = true;
        break;
      case RAYIMGUI_KEY_E:
        gui.movement_speed = std::min(gui.movement_speed * 2, 128.0);
        break;
      case RAYIMGUI_KEY_C:
        gui.movement_speed = std::max(gui.movement_speed / 2, .001);
        break;
      case RAYIMGUI_KEY_TAB:
        gui.orbit = !gui.orbit;
        break;
      case RAYIMGUI_KEY_F:
        if (native_drag_fast) {
          native_fast_saved = !native_fast_saved;
        } else {
          write_fast_output = !write_fast_output;
        }
        gui.fast_preview = native_drag_fast ? native_fast_saved : write_fast_output;
        reset = true;
        break;
      case RAYIMGUI_KEY_UP:
        cam->update_fov(-speed);
        reset = true;
        break;
      case RAYIMGUI_KEY_DOWN:
        cam->update_fov(speed);
        reset = true;
        break;
      case RAYIMGUI_KEY_LEFT:
        cam->update_aperture(-speed * .1f);
        reset = true;
        break;
      case RAYIMGUI_KEY_RIGHT:
        cam->update_aperture(speed * .1f);
        reset = true;
        break;
      case RAYIMGUI_KEY_1:
        cam->update_focal_distance(-speed);
        reset = true;
        break;
      case RAYIMGUI_KEY_2:
        cam->update_focal_distance(speed);
        reset = true;
        break;
      case RAYIMGUI_KEY_3:
      case RAYIMGUI_KEY_4: {
        Float angle = (key.code == RAYIMGUI_KEY_3 ? 1 : -1) * speed * 8;
        *EnvObjectToWorld = RotateY(angle) * (*EnvObjectToWorld);
        *EnvWorldToObject = RotateY(-angle) * (*EnvWorldToObject);
        native_env_angle -= angle;
        reset = true;
        break;
      }
      case RAYIMGUI_KEY_LEFT_BRACKET:
      case RAYIMGUI_KEY_RIGHT_BRACKET:
        if (shift) {
          AdjustShutterSpeedStops(key.code == RAYIMGUI_KEY_LEFT_BRACKET ? -1.f / 3
                                                                        : 1.f / 3);
          reset = true;
        } else {
          if (key.code == RAYIMGUI_KEY_LEFT_BRACKET) {
            DecreasePreviewExposure();
          } else {
            IncreasePreviewExposure();
          }
          gui.exposure = preview_exposure_adjustment;
        }
        break;
      case RAYIMGUI_KEY_R:
        gui.request_reset = true;
        break;
      case RAYIMGUI_KEY_ENTER:
      case RAYIMGUI_KEY_KEYPAD_ENTER:
        if (shift) {
          SavePreviewSnapshot();
        } else if (deferred_render) {
          render_requested = !render_requested;
        }
        break;
      case RAYIMGUI_KEY_P:
        PrintCameraInfo(native_env_angle);
        break;
      case RAYIMGUI_KEY_K:
        SaveCurrentKeyframe(native_env_angle);
        PrintCameraInfo(native_env_angle);
        break;
      case RAYIMGUI_KEY_L:
        if (shift) {
          ToggleKeyframeMotionClosed();
        } else if (!Keyframes.empty()) {
          reset = ApplyKeyframe(static_cast<int>(Keyframes.size()) - 1,
                                &native_env_angle) ||
                  reset;
        }
        break;
      case RAYIMGUI_KEY_COMMA:
        if (shift) {
          reset = JumpKeyframe(-1, &native_env_angle) || reset;
        }
        break;
      case RAYIMGUI_KEY_PERIOD:
        if (shift) {
          reset = JumpKeyframe(1, &native_env_angle) || reset;
        }
        break;
      case RAYIMGUI_KEY_SLASH:
        reset = DeleteCurrentKeyframe(&native_env_angle) || reset;
        break;
      case RAYIMGUI_KEY_M:
        // A viewport click or navigation in this batch takes priority over M.
        if (stopped) {
          break;
        }
        gui.animation_action = RayrenderGui::AnimationAction::Play;
        reset = ApplyNativeAnimationControls() || reset;
        break;
      case RAYIMGUI_KEY_B:
        ToggleCameraMotionBlur();
        reset = true;
        break;
      case RAYIMGUI_KEY_H:
        if (update_atmosphere) {
          gui.haze = !gui.haze;
          if (gui.haze) {
            gui.altitude = 1;
          }
          gui.atmosphere_pending = true;
        }
        break;
      case RAYIMGUI_KEY_Y:
        if (update_atmosphere) {
          gui.altitude = !gui.altitude;
          if (!gui.altitude) {
            gui.haze = 0;
          }
          gui.atmosphere_pending = true;
        }
        break;
      default:
        break;
      }
    }
  }

  gui.keys.clear();
  if (gui.pick_pending) {
    if (world && !IsPreviewMotionActive()) {
      reset = PickCameraTarget(gui.pick_u, gui.pick_v, gui.pick_focus, world) || reset;
    }
    gui.pick_pending = false;
  }

  reset = ApplyNativeCameraControls() || reset;
  if (gui.request_reset) {
    cam->reset();
    SyncNativeCameraControls(true);
    gui.movement_speed = 1;
    gui.request_reset = false;
    *EnvObjectToWorld = Start_EnvObjectToWorld;
    *EnvWorldToObject = Start_EnvWorldToObject;
    native_env_angle = 0;
    reset = true;
  }

  if (reset) {
    ApplyStaticPreviewCameraMotionRange(cam);
  }

  if (!gui.animation_paused && AdvancePreviewMotion(&native_env_angle)) {
    reset = true;
  }
  SyncNativeAnimationState();
  gui.exposure = preview_exposure_adjustment;
  SyncNativeCameraControls();

  return reset;
}

// One checkpoint records input intent, applies ordinary edits, then consumes
// history requests. Renderer rebuilds and all R callbacks remain outside drawing.
bool PreviewDisplay::CommitNativeEdits(hitable* world) {
  auto& gui = *native_gui;
  BeginNativeHistory();
  const bool was_playing = IsPreviewMotionActive();
  bool reset = false;
  if (was_playing && !gui.history_requests.empty()) {
    reset = CancelPreviewMotion(&native_env_angle);
  }
  const bool edited =
      gui.history_dirty || gui.exposure_pending || gui.denoise_pending ||
      gui.fast_pending || gui.depth_pending || gui.request_reset || gui.pick_pending ||
      !gui.keys.empty() || gui.object.apply_transform || gui.object.apply_material ||
      gui.object.revert || gui.object.cancel_transform || gui.sky_model_pending ||
      gui.atmosphere_pending || gui.sun_pending || gui.location_pending ||
      gui.camera.pending ||
      (gui.atmosphere_parameters_pending && !gui.atmosphere_parameters_editing) ||
      gui.animation_settings_pending ||
      (gui.animation_action != RayrenderGui::AnimationAction::None &&
       gui.animation_action != RayrenderGui::AnimationAction::Play);
  if (gui.exposure_pending) {
    preview_exposure_adjustment = static_cast<Float>(gui.exposure);
    gui.exposure_pending = false;
  }
  reset = ApplyNativeControls(world) || reset;
  reset = ApplyNativeObjectControls() || reset;
  const bool dragging = (gui.object.selected && gui.object.transform_active) ||
                        (gui.sun_editing && bool(update_sun)) ||
                        (gui.camera.editing && gui.camera.editable);
  if (dragging != native_drag_fast) {
    if (dragging) {
      native_fast_saved = write_fast_output;
      reset = !write_fast_output || reset;
      write_fast_output = true;
    } else {
      reset = write_fast_output != native_fast_saved || reset;
      write_fast_output = native_fast_saved;
    }
    native_drag_fast = dragging;
  }
  ApplyNativeSkyControls();
  FinishNativeHistory(edited);
  reset = ApplyNativeHistory() || reset;
  if (gui.request_render) {
    render_requested = true;
    gui.request_render = false;
  }
  ApplyNativeExport();
  // A restart changes the scene being sampled, not the image already displayed.
  // Keep that image's outline visible through worker-wait polls; the next
  // completed sample replaces both layers before presenting another frame.
  return reset;
}

// Convert a completed sample into display pixels, publish the UI snapshot, and
// consume queued edits while no sampling worker can read mutable renderer state.
bool PreviewDisplay::DrawNativeGui(adaptive_sampler& sampler, size_t ns,
                                   Float percent_done, hitable* world,
                                   random_gen& rng) {
  RayrenderGui& gui = *native_gui;
#ifdef HAS_OIDN
  // Fast preview writes the small filter's result into the display buffer; its
  // full-size filter may still be uninitialized after enabling denoising there.
  bool use_denoised = denoise && denoiser != nullptr &&
                      (write_fast_output ? HasDenoisedPreview() : denoiser->Ready());
  if (use_denoised && !write_fast_output) {
    denoiser->Execute();
    if (!denoiser->ReportError()) {
      MarkDenoisedPreviewReady(ns + 1);
    }
  }

  RayMatrix& rgb = use_denoised ? sampler.draw_rgb_output : sampler.rgb;
#else
  RayMatrix& rgb = sampler.rgb;
#endif
  PreparePreviewColor(sampler, rgb, ns);
  const size_t width = gui.width, height = gui.height, count = width * height;
  gui.display_rgb.resize(count * 3);
  gui.pixels.resize(count * 4);
  snapshot_width = gui.width;
  snapshot_height = gui.height;
  snapshot_pixels.resize(count * 3);
  // Film coordinates match the legacy preview's horizontal and vertical flip.
  for (size_t y = 0; y < height; ++y) {
    for (size_t x = 0; x < width; ++x) {
      const size_t sx = width - 1 - x, sy = height - 1 - y, index = sx + width * sy;
      const Float samples =
          sampler.finalized[index] ? (interactive ? 1.f : 4.f) : Float(ns + 1);
      const point3f mapped = ApplyPreviewColor(rgb, sx, sy, samples);
      for (size_t channel = 0; channel < 3; ++channel) {
        gui.display_rgb[3 * (x + width * y) + channel] =
            clamp(mapped[channel], 0.f, 1.f);
      }
      sampler.just_finalized[index] = false;
    }
  }
  // Renderer-owned overlays use the same display transform as the old preview.
  std::vector<Float> overlay(gui.display_rgb.begin(), gui.display_rgb.end());
  CompositeTextOverlaysToFloatBuffer(overlay, world, rng);
  CompositeLineOverlaysToFloatBuffer(overlay, world, rng);
  for (size_t i = 0; i < count; ++i) {
    for (size_t channel = 0; channel < 3; ++channel) {
      const uint8_t value =
          static_cast<uint8_t>(255 * clamp(overlay[3 * i + channel], 0.f, 1.f));
      gui.pixels[4 * i + channel] = value;
      snapshot_pixels[3 * i + channel] = value;
    }
    gui.pixels[4 * i + 3] = 255;
  }

  CaptureVolumeSnapshot(sampler, rgb, ns + 1, world, rng);
  gui.samples = ns + 1;
  gui.progress = percent_done;
  // Refresh from renderer state only when doing so would not overwrite a user
  // change that has been queued by an intervening GUI poll.
  if (!gui.fast_pending) {
    gui.fast_preview = native_drag_fast ? native_fast_saved : write_fast_output;
  }
#ifdef HAS_OIDN
  if (!gui.denoise_pending) {
    gui.denoise_enabled = denoise;
  }
#endif

  gui.has_atmosphere = static_cast<bool>(update_atmosphere);
  gui.has_sun = static_cast<bool>(update_sun);
  if (!gui.sun_pending && !gui.sun_editing) {
    gui.sun_elevation = sun_elevation;
    gui.sun_azimuth = sun_azimuth;
  }

  if (!gui.atmosphere_pending) {
    gui.haze = atmosphere_haze;
    gui.altitude = atmosphere_query_altitude;
  }

  SyncNativeCameraControls();
  UpdateNativeObjectCamera();
  UpdateNativeSelectionMask();
  gui.snapshot_camera = NativeKeyframeCamera(CreateCurrentKeyframe(native_env_angle));
  gui.capture_keyframe_snapshots();
  SyncNativeAnimationState();
  BeginNativeHistory();
  gui.publish(); // Synchronous copy; no worker can be writing this snapshot.
  terminate = gui.poll(true);
  if (terminate) {
    return false;
  }

  bool reset = CommitNativeEdits(world);
  // After a scene/camera edit, retain the published image and its cached outline
  // as a pair. Updating or hiding only the outline here makes GUI-only frames
  // disagree with the render frames while the replacement image is sampled.
  if (!reset && UpdateNativeSelectionMask()) {
    // Selection alone does not restart rendering; show its feedback immediately.
    terminate = gui.poll(true);
  }
  return reset;
}

// Apply the latest queued sky edit at a worker-drained renderer checkpoint.
// Sun drags update each preview sample; location/time still requires validation.
void PreviewDisplay::ApplyNativeSkyControls() {
  RayrenderGui& gui = *native_gui;
  if (gui.sky_model_pending && update_sky_model) {
    gui.sky_error = update_sky_model(gui.sky_model);
    gui.sky_model_pending = false;
    if (gui.sky_error.empty()) {
      sky_model = gui.sky_model;
      native_sun_preview = false;
      atmosphere_changed = true;
    } else {
      gui.sky_model = sky_model;
    }
    return;
  }
  if (gui.atmosphere_pending && update_atmosphere) {
    // Unlike the legacy OS callback helper, this safe checkpoint propagates
    // scene/sky errors and R interrupts. They never trigger GUI fallback.
    update_atmosphere(gui.haze != 0, gui.altitude != 0);
    atmosphere_haze = gui.haze != 0;
    atmosphere_query_altitude = gui.altitude != 0;
    atmosphere_changed = true;
    gui.atmosphere_pending = false;
  }

  if (gui.atmosphere_parameters_pending && !gui.atmosphere_parameters_editing &&
      update_atmosphere_parameters) {
    // Rebuild all sky components together, retaining the active light if any
    // validation, metadata query or texture preparation fails.
    gui.sky_error =
        update_atmosphere_parameters(gui.base_altitude, gui.meters_per_unit);
    gui.atmosphere_parameters_pending = false;
    if (gui.sky_error.empty()) {
      atmosphere_changed = true;
    }
  }

  // A drag uses a cheaper sky map/proposal as well as Fast scene sampling.
  // Refine that sky once on release, even if the last angle did not change.
  if ((gui.sun_pending || (native_sun_preview && !gui.sun_editing)) && update_sun) {
    const bool previous_preview = native_sun_preview;
    native_sun_preview = false;
    gui.sun_elevation = std::clamp(gui.sun_elevation, -90.0, MaxSunElevationDegrees);
    gui.sun_azimuth = std::clamp(gui.sun_azimuth, 0.0, 360.0);
    try {
      update_sun(gui.sun_elevation, gui.sun_azimuth);
      gui.sky_error.clear();
    } catch (const std::exception& error) {
      native_sun_preview = gui.sun_editing && previous_preview;
      gui.sky_error = error.what();
      gui.sun_pending = false;
      SetSunPosition(sun_elevation, sun_azimuth);
      return;
    }
    native_sun_preview = gui.sun_editing;
    sun_elevation = gui.sun_elevation;
    sun_azimuth = gui.sun_azimuth;
    atmosphere_changed = true;
    gui.sun_pending = false;
    gui.manual_sun = true;
  }

  if (gui.location_pending && update_sky) {
    gui.sky_error = update_sky(gui.latitude, gui.longitude, gui.datetime);
    gui.location_pending = false;
    if (gui.sky_error.empty()) {
      atmosphere_changed = true;
      gui.manual_sun = false;
      gui.sun_pending = false;
    }
  }
}

// Resolve selection/cancel requests before applying edits, so each operation
// uses the current selected root and the renderer's committed transform baseline.
bool PreviewDisplay::ApplyNativeObjectControls() {
  if (!scene_editor || !native_gui || !interactive) {
    return false;
  }

  auto& ui = native_gui->object;
  bool reset = false;
  if (ui.begin_transform && ui.selected) {
    scene_editor->BeginTransform();
    ui.begin_transform = false;
  }
  if (ui.clear_pending || ui.select_pending || ui.pick_pending) {
    scene_editor->EndTransform();
  }
  if (ui.clear_pending) {
    scene_editor->Describe(0, ui);
    ui.clear_pending = false;
  }

  if (ui.select_pending) {
    const auto id = ui.select_id;
    ui.select_pending = false;
    ui.pick_pending = false;
    scene_editor->Describe(id, ui);
  }

  if (ui.cancel_transform) {
    try {
      reset = scene_editor->CancelTransform(ui);
    } catch (const Rcpp::internal::InterruptedException&) {
      throw;
    } catch (const std::exception& error) {
      ui.error = std::string("Cannot cancel transform: ") + error.what();
      ui.cancel_transform = ui.transform_active = false;
      scene_editor->EndTransform();
    }
  }

  if (ui.pick_pending) {
    ui.pick_pending = false;
    if (!IsPreviewMotionActive()) {
      // Use a centered lens sample and midpoint shutter sample for repeatable
      // picking. RealisticCamera uses the opposite film-coordinate convention
      // from get_ray(), so its normalized coordinates are flipped here.
      Ray ray;
      if (cam->get_fov() < 0) {
        CameraSample sample(
            point2f(1 - ui.pick_u, 1 - ui.pick_v), point2f(.5f, .5f), .5f);
        if (!(cam->GenerateRay(sample, &ray) > 0)) {
          return false;
        }
      } else {
        ray = cam->get_ray(ui.pick_u, ui.pick_v, point3f(0), .5f);
      }
      scene_editor->Pick(ray, ui, [this]() {
        return PollCloseEvent();
      });
    }
  }

  // Keep numeric fields in sync with gizmo motion without overwriting typed edits.
  if (ui.transform_pending && !ui.numeric_transform) {
    PreviewDecompose(ui);
  }

  reset = scene_editor->Apply(ui) || reset;
  if (ui.end_transform && !ui.transform_active) {
    scene_editor->EndTransform();
    ui.end_transform = false;
  }
  return reset;
}

// Visibility is refreshed only when selection, geometry, camera or film size
// changes. Cap the query resolution so large output images stay interactive.
bool PreviewDisplay::UpdateNativeSelectionMask() {
  auto& gui = *native_gui;
  if (!scene_editor || !interactive || !gui.object.selected) {
    const bool changed = gui.selection_visible;
    gui.selection_id = 0;
    gui.selection_visible = false;
    return changed;
  }
  if (IsPreviewMotionActive()) {
    const bool changed = gui.selection_visible;
    gui.selection_visible = false;
    return changed;
  }
  // Keep silhouette feedback inexpensive during live transforms; refine it on release.
  const double mask_limit = gui.object.transform_active ? 192.0 : 768.0;
  const double reduction = std::min(1.0, mask_limit / std::max(gui.width, gui.height));
  const uint32_t width = std::max(1u, uint32_t(gui.width * reduction));
  const uint32_t height = std::max(1u, uint32_t(gui.height * reduction));
  std::vector<double> camera_state;
  auto append = [&](const auto& value) {
    for (int i = 0; i < 3; ++i) {
      camera_state.push_back(value[i]);
    }
  };
  append(cam->get_origin());
  append(cam->get_lookat());
  append(cam->get_u());
  append(cam->get_v());
  append(cam->get_w());
  camera_state.push_back(cam->get_fov());
  camera_state.push_back(cam->get_aperture());
  camera_state.push_back(cam->get_focal_distance());
  camera_state.push_back(cam->get_shutter_speed());
  const auto ortho = cam->get_ortho();
  camera_state.push_back(ortho[0]);
  camera_state.push_back(ortho[1]);
  const uint64_t id = gui.object.id, revision = scene_editor->Revision();
  if (gui.selection_id == id && gui.selection_revision == revision &&
      gui.selection_width == width && gui.selection_height == height &&
      gui.selection_camera == camera_state) {
    const bool changed = !gui.selection_visible;
    gui.selection_visible = true;
    if (changed) {
      gui.publish_selection();
    }
    return changed;
  }
  // Use the same centered lens and midpoint shutter as picking. Coordinates
  // here are top-left display UVs, accounting for the renderer's film flip.
  auto mask = scene_editor->SelectionMask(
      id,
      width,
      height,
      [this](float u, float v, Ray& ray) {
        if (cam->get_fov() < 0) {
          return cam->GenerateRay(CameraSample(point2f(u, v), point2f(.5f, .5f), .5f),
                                  &ray) > 0;
        }
        ray = cam->get_ray(1 - u, 1 - v, point3f(0), .5f);
        return true;
      },
      [this]() {
        return PollCloseEvent();
      });
  if (mask.empty()) {
    gui.selection_visible = false;
    return false;
  }
  // Coverage and outline share one lifetime. Updating the render's colors alone
  // never re-traces visibility or recomputes which pixels form the boundary.
  auto outline = PreviewSelectionOverlay::Outline(mask, width, height);
  gui.selection_mask = std::move(mask);
  gui.selection_outline = std::move(outline);
  gui.selection_width = width;
  gui.selection_height = height;
  gui.selection_camera = std::move(camera_state);
  gui.selection_revision = revision;
  gui.selection_id = id;
  gui.selection_visible = true;
  gui.publish_selection();
  return true;
}

// Build gizmo matrices that project onto the displayed renderer image. Cameras
// without a matching linear projection retain numeric transform editing only.
void PreviewDisplay::UpdateNativeObjectCamera() {
  if (!native_gui) {
    return;
  }

  auto& ui = native_gui->object;
  ui.enabled = scene_editor && interactive;
  const double fov = cam->get_fov();
  ui.projection_valid =
      ui.enabled && !IsPreviewMotionActive() && (fov == 0 || (fov > 0 && fov < 180));
  if (!ui.projection_valid) {
    return;
  }

  // Display pixels reverse the horizontal film axis; OpenGL looks down -Z.
  // Negate camera right and forward while retaining camera up for screen Y.
  auto u = -cam->get_u(), v = cam->get_v(), w = -cam->get_w();
  auto origin = cam->get_origin();
  for (int i = 0; i < 3; ++i) {
    ui.view[4 * i] = u[i];
    ui.view[4 * i + 1] = v[i];
    ui.view[4 * i + 2] = w[i];
    ui.view[4 * i + 3] = 0;
  }

  const auto o = convert_to_vec3(origin);
  ui.view[12] = -dot(u, o);
  ui.view[13] = -dot(v, o);
  ui.view[14] = -dot(w, o);
  ui.view[15] = 1;
  ui.projection.fill(0);
  const double near = .001, far = 1e8;
  if (fov == 0) {
    auto dimensions = cam->get_ortho();
    ui.projection[0] = 2 / dimensions[0];
    ui.projection[5] = 2 / dimensions[1];
    ui.projection[10] = -2 / (far - near);
    ui.projection[14] = -(far + near) / (far - near);
    ui.projection[15] = 1;
    ui.projection_kind = RAYIMGUI_ORTHOGRAPHIC;
  } else {
    const double height = std::tan(fov * M_PI / 360),
                 aspect = double(native_gui->width) / native_gui->height;
    ui.projection[0] = 1 / (height * aspect);
    ui.projection[5] = 1 / height;
    ui.projection[10] = -(far + near) / (far - near);
    ui.projection[11] = -1;
    ui.projection[14] = -2 * far * near / (far - near);
    ui.projection_kind = RAYIMGUI_PERSPECTIVE;
  }
}

// Capture committed renderer values, not widget drafts. This is also returned
// with the image so callers can inspect the state after closing the editor.
Rcpp::List PreviewDisplay::NativeEditorState() const {
  Rcpp::List camera = CreateCurrentKeyframe(native_env_angle);
  camera["fov"] = cam->get_fov();
  Rcpp::RObject sky = R_NilValue;
  if (export_sky) {
    sky = export_sky();
  }
  bool use_denoising = false;
#ifdef HAS_OIDN
  use_denoising = denoise;
#endif
  return Rcpp::List::create(
      Rcpp::_["version"] = 1,
      Rcpp::_["objects"] = scene_editor ? scene_editor->ExportEdits() : Rcpp::List(),
      Rcpp::_["camera"] = camera,
      Rcpp::_["sky"] = sky,
      Rcpp::_["denoise"] = use_denoising,
      Rcpp::_["exposure"] =
          double(preview_exposure_scale) * preview_exposure_adjustment,
      Rcpp::_["camera_motion_blur"] = CameraMotionBlurEnabled(),
      Rcpp::_["shutter_speed"] = GetShutterSpeed(),
      Rcpp::_["integrator_type"] = native_integrator,
      Rcpp::_["max_depth"] = double(max_depth),
      Rcpp::_["camera_rotation"] = cam->get_free_rotation() ? "free" : "clamped",
      Rcpp::_["camera_description_file"] = camera_rig ? camera_rig->Source() : "");
}

void PreviewDisplay::ApplyNativeExport() {
  if (!native_gui || !native_gui->export_pending) {
    return;
  }
  auto& gui = *native_gui;
  gui.export_pending = false;
  if (!export_scene) {
    gui.export_message = "Scene export is unavailable.";
    return;
  }
  try {
    const auto filename = export_scene(NativeEditorState(), gui.export_filename);
    gui.export_message = "Exported: " + filename;
    Rprintf("Exported rayrender scene: %s\n", filename.c_str());
  } catch (const Rcpp::internal::InterruptedException&) {
    throw;
  } catch (const std::exception& error) {
    // File/path errors belong in the panel and must not stop the live renderer.
    gui.export_message = std::string("Export failed: ") + error.what();
  }
}
