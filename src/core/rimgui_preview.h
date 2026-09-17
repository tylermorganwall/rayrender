/* Copyright (c) 2026 Tyler Morgan-Wall. MIT; LICENSE.protocol.
 * Included once, at the end of PreviewDisplay.cpp. These checkpoints run on
 * R's main thread after wait_for_render_jobs has drained all sample tasks.
 */
#include "preview_scene.h"

void PreviewDisplay::AttachNativeGui(RayrenderGui* gui, bool edit, bool deferred) {
  native_gui = gui;
  preview = true;
  interactive = edit;
  deferred_render = edit && deferred;
  render_requested = !deferred_render;
  gui->can_edit = edit;
  gui->deferred = deferred_render;
  gui->exposure = preview_exposure_adjustment;
  native_base_step =
      std::max(Float(.001), (cam->get_origin() - cam->get_lookat()).length() / 20);
}

void PreviewDisplay::SetSunControls(double elevation, double azimuth,
                                    std::function<void(double, double)> update) {
  sun_elevation = elevation;
  sun_azimuth = azimuth;
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

// Replay queued input through the existing camera controls. Returning true asks
// the caller to discard accumulated samples after the camera or preview mode changes.
bool PreviewDisplay::ApplyNativeControls(hitable* world) {
  RayrenderGui& gui = *native_gui;
  bool reset = false;
  if (!interactive) {
    gui.keys.clear();
    gui.pick_pending = false;
    return false;
  }

  if (gui.fast_pending) {
    reset = write_fast_output != (gui.fast_preview != 0);
    write_fast_output = gui.fast_preview != 0;
    gui.fast_pending = false;
  }

  for (const auto& key : gui.keys) {
    for (unsigned count = 0; count < key.count; ++count) {
      const bool shift = (key.modifiers & RIMGUI_SHIFT) != 0;
      const Float speed =
          static_cast<Float>(std::clamp(gui.movement_speed, .001, 128.0));
      if (IsPreviewMotionActive() && key.code != RIMGUI_KEY_M) {
        continue;
      }
      switch (key.code) {
      case RIMGUI_KEY_W:
        if (shift) {
          cam->rotate_forward(speed);
        } else if (!gui.orbit || (cam->get_origin() - cam->get_lookat()).length() >
                                     speed * native_base_step) {
          cam->update_position(
              speed * cam->get_w() * native_base_step, gui.orbit, false);
        }
        reset = true;
        break;
      case RIMGUI_KEY_S:
        if (shift) {
          cam->rotate_forward(-speed);
        } else {
          cam->update_position(
              -speed * cam->get_w() * native_base_step, gui.orbit, false);
        }
        reset = true;
        break;
      case RIMGUI_KEY_A:
        if (shift) {
          cam->rotate_up(-speed);
        } else {
          cam->update_position(speed * cam->get_u() * native_base_step, gui.orbit);
        }
        reset = true;
        break;
      case RIMGUI_KEY_D:
        if (shift) {
          cam->rotate_up(speed);
        } else {
          cam->update_position(-speed * cam->get_u() * native_base_step, gui.orbit);
        }
        reset = true;
        break;
      case RIMGUI_KEY_Q:
        cam->update_position(speed * cam->get_v() * native_base_step, gui.orbit);
        reset = true;
        break;
      case RIMGUI_KEY_Z:
        cam->update_position(-speed * cam->get_v() * native_base_step, gui.orbit);
        reset = true;
        break;
      case RIMGUI_KEY_E:
        gui.movement_speed = std::min(gui.movement_speed * 2, 128.0);
        break;
      case RIMGUI_KEY_C:
        gui.movement_speed = std::max(gui.movement_speed / 2, .001);
        break;
      case RIMGUI_KEY_TAB:
        gui.orbit = !gui.orbit;
        break;
      case RIMGUI_KEY_F:
        write_fast_output = !write_fast_output;
        gui.fast_preview = write_fast_output;
        reset = true;
        break;
      case RIMGUI_KEY_UP:
        cam->update_fov(-speed);
        reset = true;
        break;
      case RIMGUI_KEY_DOWN:
        cam->update_fov(speed);
        reset = true;
        break;
      case RIMGUI_KEY_LEFT:
        cam->update_aperture(-speed * .1f);
        reset = true;
        break;
      case RIMGUI_KEY_RIGHT:
        cam->update_aperture(speed * .1f);
        reset = true;
        break;
      case RIMGUI_KEY_1:
        cam->update_focal_distance(-speed);
        reset = true;
        break;
      case RIMGUI_KEY_2:
        cam->update_focal_distance(speed);
        reset = true;
        break;
      case RIMGUI_KEY_3:
      case RIMGUI_KEY_4: {
        Float angle = (key.code == RIMGUI_KEY_3 ? 1 : -1) * speed * 8;
        *EnvObjectToWorld = RotateY(angle) * (*EnvObjectToWorld);
        *EnvWorldToObject = RotateY(-angle) * (*EnvWorldToObject);
        native_env_angle -= angle;
        reset = true;
        break;
      }
      case RIMGUI_KEY_LEFT_BRACKET:
      case RIMGUI_KEY_RIGHT_BRACKET:
        if (shift) {
          AdjustShutterSpeedStops(key.code == RIMGUI_KEY_LEFT_BRACKET ? -1.f / 3
                                                                      : 1.f / 3);
          reset = true;
        } else {
          if (key.code == RIMGUI_KEY_LEFT_BRACKET) {
            DecreasePreviewExposure();
          } else {
            IncreasePreviewExposure();
          }
          gui.exposure = preview_exposure_adjustment;
        }
        break;
      case RIMGUI_KEY_R:
        gui.request_reset = true;
        break;
      case RIMGUI_KEY_ENTER:
      case RIMGUI_KEY_KEYPAD_ENTER:
        if (shift) {
          SavePreviewSnapshot();
        } else if (deferred_render) {
          render_requested = !render_requested;
        }
        break;
      case RIMGUI_KEY_P:
        PrintCameraInfo(native_env_angle);
        break;
      case RIMGUI_KEY_K:
        SaveCurrentKeyframe(native_env_angle);
        PrintCameraInfo(native_env_angle);
        break;
      case RIMGUI_KEY_L:
        if (shift) {
          ToggleKeyframeMotionClosed();
        } else if (!Keyframes.empty()) {
          reset = ApplyKeyframe(static_cast<int>(Keyframes.size()) - 1,
                                &native_env_angle) ||
                  reset;
        }
        break;
      case RIMGUI_KEY_COMMA:
        if (shift) {
          reset = JumpKeyframe(-1, &native_env_angle) || reset;
        }
        break;
      case RIMGUI_KEY_PERIOD:
        if (shift) {
          reset = JumpKeyframe(1, &native_env_angle) || reset;
        }
        break;
      case RIMGUI_KEY_SLASH:
        reset = DeleteCurrentKeyframe(&native_env_angle) || reset;
        break;
      case RIMGUI_KEY_M:
        reset = (IsPreviewMotionActive() ? CancelPreviewMotion(&native_env_angle)
                                         : StartPreviewMotion(native_env_angle)) ||
                reset;
        break;
      case RIMGUI_KEY_B:
        ToggleCameraMotionBlur();
        reset = true;
        break;
      case RIMGUI_KEY_H:
        if (update_atmosphere) {
          gui.haze = !gui.haze;
          if (gui.haze) {
            gui.altitude = 1;
          }
          gui.atmosphere_pending = true;
        }
        break;
      case RIMGUI_KEY_Y:
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

  if (gui.request_reset) {
    cam->reset();
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

  if (AdvancePreviewMotion(&native_env_angle)) {
    reset = true;
  }

  return reset;
}

// Convert a completed sample into display pixels, publish the UI snapshot, and
// consume queued edits while no sampling worker can read mutable renderer state.
bool PreviewDisplay::DrawNativeGui(adaptive_sampler& sampler, size_t ns,
                                   Float percent_done, hitable* world,
                                   random_gen& rng) {
  RayrenderGui& gui = *native_gui;
#ifdef HAS_OIDN
  bool use_denoised = denoise && denoiser != nullptr && denoiser->Ready();
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
    gui.fast_preview = write_fast_output;
  }

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

  UpdateNativeObjectCamera();
  gui.publish(); // Synchronous copy; no worker can be writing this snapshot.
  terminate = gui.poll(true);
  if (terminate) {
    return false;
  }

  if (gui.exposure_pending) {
    preview_exposure_adjustment = static_cast<Float>(gui.exposure);
    gui.exposure_pending = false;
  }

  bool reset = ApplyNativeControls(world);
  reset = ApplyNativeObjectControls() || reset;
  if (gui.request_render) {
    render_requested = true;
    gui.request_render = false;
  }

  ApplyNativeSkyControls();
  return reset;
}

// Apply expensive sky rebuilds at the renderer checkpoint. Sun drags wait for
// release; location/time updates commit only after the R callback validates them.
void PreviewDisplay::ApplyNativeSkyControls() {
  RayrenderGui& gui = *native_gui;
  if (gui.atmosphere_pending) {
    // Unlike the legacy OS callback helper, this safe checkpoint propagates
    // scene/sky errors and R interrupts. They never trigger GUI fallback.
    update_atmosphere(gui.haze != 0, gui.altitude != 0);
    atmosphere_haze = gui.haze != 0;
    atmosphere_query_altitude = gui.altitude != 0;
    atmosphere_changed = true;
    gui.atmosphere_pending = false;
  }

  if (gui.sun_pending && !gui.sun_editing && update_sun) {
    gui.sun_elevation = std::clamp(gui.sun_elevation, -90.0, 90.0);
    gui.sun_azimuth = std::clamp(gui.sun_azimuth, 0.0, 360.0);
    update_sun(gui.sun_elevation, gui.sun_azimuth);
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
  if (ui.clear_pending) {
    scene_editor->Describe(0, ui);
    ui.clear_pending = false;
  }

  if (ui.cancel_transform) {
    scene_editor->CancelTransform(ui);
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

  return scene_editor->Apply(ui);
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
    ui.projection_kind = RIMGUI_ORTHOGRAPHIC;
  } else {
    const double height = std::tan(fov * M_PI / 360),
                 aspect = double(native_gui->width) / native_gui->height;
    ui.projection[0] = 1 / (height * aspect);
    ui.projection[5] = 1 / height;
    ui.projection[10] = -(far + near) / (far - near);
    ui.projection[11] = -1;
    ui.projection[14] = -2 * far * near / (far - near);
    ui.projection_kind = RIMGUI_PERSPECTIVE;
  }
}
