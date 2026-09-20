#ifndef RAYRENDER_PREVIEW_HISTORY_H
#define RAYRENDER_PREVIEW_HISTORY_H

// Included once by rayimgui_preview.h. History belongs to R's main thread and
// stores values, never live geometry pointers or borrowed provider resources.
struct PreviewHistoryAnimation {
  std::vector<Rcpp::List> keyframes;
  std::vector<RayrenderGui::KeyframeSnapshot> thumbnails;
};
struct PreviewHistoryState {
  Rcpp::List camera, sky, inputs, object_values;
  PreviewObjectState object;
  std::shared_ptr<const PreviewSceneSettings> scene;
  std::shared_ptr<const PreviewHistoryAnimation> animation;
  bool denoise = false, fast = false, blur = false, closed = false;
  Float shutter = 2;
  int keyframe = -1;
};
struct PreviewHistory {
  struct Entry {
    PreviewHistoryState before, after;
    uint64_t epoch = 0;
  };
  std::vector<Entry> entries;
  size_t cursor = 0;
  uint64_t scene_revision = 0;
  PreviewHistoryState current;
};

static bool HistoryEqual(SEXP a, SEXP b) {
  return R_compute_identical(a, b, 0);
}

// Only editable values enter comparisons. Projection matrices, hierarchy layout,
// sampling progress and transient hover flags must not create undo steps.
static Rcpp::List HistoryObjectValues(const PreviewObjectState& object) {
  Rcpp::List materials;
  for (const auto& material : object.materials) {
    Rcpp::List fields;
    for (const auto& field : material.fields) {
      fields.push_back(Rcpp::List::create(Rcpp::_["values"] = Rcpp::wrap(field.values),
                                          Rcpp::_["text"] = field.text));
    }
    materials.push_back(fields);
  }
  return Rcpp::List::create(Rcpp::_["id"] = double(object.id),
                            Rcpp::_["position"] = Rcpp::wrap(object.translation),
                            Rcpp::_["rotation"] = Rcpp::wrap(object.rotation),
                            Rcpp::_["scale"] = Rcpp::wrap(object.scale),
                            Rcpp::_["model"] = Rcpp::wrap(object.model),
                            Rcpp::_["tool"] = object.operation,
                            Rcpp::_["axes"] = object.mode,
                            Rcpp::_["slot"] = object.material_slot,
                            Rcpp::_["materials"] = materials);
}

static Rcpp::List HistoryInputs(const RayrenderGui& gui) {
  return Rcpp::List::create(Rcpp::_["exposure"] = gui.exposure,
                            Rcpp::_["orbit"] = gui.orbit,
                            Rcpp::_["speed"] = gui.movement_speed,
                            Rcpp::_["sky_model"] = gui.sky_model,
                            Rcpp::_["latitude"] = gui.latitude,
                            Rcpp::_["longitude"] = gui.longitude,
                            Rcpp::_["datetime"] = std::string(gui.datetime),
                            Rcpp::_["elevation"] = gui.sun_elevation,
                            Rcpp::_["azimuth"] = gui.sun_azimuth,
                            Rcpp::_["manual"] = gui.manual_sun,
                            Rcpp::_["haze"] = gui.haze,
                            Rcpp::_["altitude"] = gui.altitude,
                            Rcpp::_["export_file"] = std::string(gui.export_filename));
}

static bool HistoryEqual(const PreviewHistoryState& a, const PreviewHistoryState& b) {
  return a.scene == b.scene && a.animation == b.animation && a.denoise == b.denoise &&
         a.fast == b.fast && a.blur == b.blur && a.closed == b.closed &&
         a.shutter == b.shutter && a.keyframe == b.keyframe &&
         HistoryEqual(a.camera, b.camera) && HistoryEqual(a.sky, b.sky) &&
         HistoryEqual(a.inputs, b.inputs) &&
         HistoryEqual(a.object_values, b.object_values);
}

// Reuse immutable scene/keyframe snapshots across unrelated input changes. Images
// are shared once per saved keyframe, so a long drag cannot duplicate thumbnails.
static PreviewHistoryState CaptureHistory(PreviewDisplay& display,
                                          PreviewHistory& history,
                                          const Rcpp::List& camera) {
  const auto& gui = *display.native_gui;
  PreviewHistoryState state;
  // Normalize temporary playback poses to the camera restored by Stop. Inputs
  // can still be edited during playback without recording every movie frame.
  state.camera = display.IsPreviewMotionActive()
                     ? Rcpp::clone(display.preview_motion_restore_state)
                     : camera;
  state.sky = display.export_sky ? display.export_sky() : Rcpp::List();
  state.inputs = HistoryInputs(gui);
  if (display.IsPreviewMotionActive()) {
    state.inputs["exposure"] = state.camera["exposure"];
  }
  state.object = gui.object;
  state.object_values = HistoryObjectValues(gui.object);
  state.scene = history.current.scene;
  if (display.scene_editor &&
      (!state.scene || history.scene_revision != display.scene_editor->Revision())) {
    state.scene =
        std::make_shared<PreviewSceneSettings>(display.scene_editor->settings);
    history.scene_revision = display.scene_editor->Revision();
  }
  state.animation = history.current.animation;
  bool same_keyframes = state.animation && state.animation->thumbnails.size() ==
                                               gui.keyframe_snapshots.size();
  if (same_keyframes) {
    for (size_t i = 0; i < gui.keyframe_snapshots.size(); ++i) {
      same_keyframes &=
          state.animation->thumbnails[i].id == gui.keyframe_snapshots[i].id;
    }
  }
  if (!same_keyframes) {
    auto animation = std::make_shared<PreviewHistoryAnimation>();
    animation->keyframes = display.Keyframes;
    for (const auto& snapshot : gui.keyframe_snapshots) {
      RayrenderGui::KeyframeSnapshot saved;
      saved.id = snapshot.id;
      saved.camera = snapshot.camera;
      saved.image = snapshot.image;
      animation->thumbnails.push_back(std::move(saved));
    }
    state.animation = std::move(animation);
  }
#ifdef HAS_OIDN
  state.denoise = display.denoise;
#endif
  state.fast = display.write_fast_output;
  state.blur = display.CameraMotionBlurEnabled();
  state.closed = display.KeyframeMotionClosed();
  state.shutter = display.GetShutterSpeed();
  state.keyframe = display.IsPreviewMotionActive()
                       ? display.preview_motion_restore_keyframe
                       : display.current_keyframe;
  return state;
}

void PreviewDisplay::BeginNativeHistory() {
  if (!native_history) {
    native_history = std::make_shared<PreviewHistory>();
    native_history->current =
        CaptureHistory(*this, *native_history, CreateCurrentKeyframe(native_env_angle));
  }
}

void PreviewDisplay::FinishNativeHistory(bool edited) {
  auto& gui = *native_gui;
  BeginNativeHistory();
  auto& history = *native_history;
  auto next = CaptureHistory(*this, history, CreateCurrentKeyframe(native_env_angle));
  if (edited && !HistoryEqual(history.current, next)) {
    history.entries.erase(history.entries.begin() + history.cursor,
                          history.entries.end());
    if (!history.entries.empty() && gui.history_epoch != 0 &&
        history.entries.back().epoch == gui.history_epoch) {
      history.entries.back().after = next;
      // Escaping a gizmo or dragging back to the starting value is a no-op.
      if (HistoryEqual(history.entries.back().before, next)) {
        history.entries.pop_back();
      }
    } else {
      history.entries.push_back({history.current, next, gui.history_epoch});
      if (history.entries.size() > 100) {
        history.entries.erase(history.entries.begin());
      }
    }
    history.cursor = history.entries.size();
    gui.history_message.clear();
  }
  history.current = std::move(next);
  gui.history_dirty = false;
  gui.can_undo = history.cursor > 0;
  gui.can_redo = history.cursor < history.entries.size();
}

// Clear requests, but preserve unapplied material/transform drafts saved in a
// snapshot. Restoring a snapshot must never replay an Apply, pick, or export click.
static void ClearHistoryRequests(RayrenderGui& gui) {
  gui.keys.clear();
  gui.pick_pending = gui.request_reset = gui.request_render = false;
  gui.fast_pending = gui.denoise_pending = gui.exposure_pending = false;
  gui.sky_model_pending = gui.sun_pending = gui.sun_editing = false;
  gui.location_pending = gui.atmosphere_pending = gui.export_pending = false;
  gui.animation_settings_pending = false;
  gui.animation_action = RayrenderGui::AnimationAction::None;
  gui.object.apply_transform = gui.object.apply_material = gui.object.revert = false;
  gui.object.pick_pending = gui.object.select_pending = gui.object.clear_pending =
      false;
  gui.object.cancel_transform = gui.object.transform_active = false;
}

static void RestoreHistoryInputs(RayrenderGui& gui, const PreviewHistoryState& target) {
  const auto& input = target.inputs;
  gui.exposure = Rcpp::as<double>(input["exposure"]);
  gui.orbit = Rcpp::as<int>(input["orbit"]);
  gui.movement_speed = Rcpp::as<double>(input["speed"]);
  gui.sky_model = Rcpp::as<int>(input["sky_model"]);
  gui.latitude = Rcpp::as<double>(input["latitude"]);
  gui.longitude = Rcpp::as<double>(input["longitude"]);
  gui.set_datetime(Rcpp::as<std::string>(input["datetime"]));
  gui.sun_elevation = Rcpp::as<double>(input["elevation"]);
  gui.sun_azimuth = Rcpp::as<double>(input["azimuth"]);
  gui.manual_sun = Rcpp::as<bool>(input["manual"]);
  gui.haze = Rcpp::as<int>(input["haze"]);
  gui.altitude = Rcpp::as<int>(input["altitude"]);
  const auto filename = Rcpp::as<std::string>(input["export_file"]);
  std::snprintf(
      gui.export_filename, sizeof(gui.export_filename), "%s", filename.c_str());
  gui.fast_preview = target.fast;
  gui.denoise_enabled = target.denoise;
  gui.object = target.object;
  gui.object.error.clear();
  gui.sky_error.clear();
  ClearHistoryRequests(gui);
}

bool PreviewDisplay::ApplyNativeHistory() {
  auto& gui = *native_gui;
  if (gui.history_requests.empty()) {
    return false;
  }
  BeginNativeHistory();
  auto& history = *native_history;
  bool reset = CancelPreviewMotion(&native_env_angle);
  const auto requests = std::move(gui.history_requests);
  gui.history_requests.clear();
  for (const int direction : requests) {
    if ((direction < 0 && history.cursor == 0) ||
        (direction > 0 && history.cursor == history.entries.size())) {
      continue;
    }
    const size_t index = direction < 0 ? history.cursor - 1 : history.cursor;
    const auto& target =
        direction < 0 ? history.entries[index].before : history.entries[index].after;
    std::vector<rayimgui_texture_handle> created;
    try {
      const bool scene_changed = target.scene != history.current.scene;
      const bool sky_changed = !HistoryEqual(target.sky, history.current.sky);
      const bool animation_changed = target.animation != history.current.animation;
      const bool camera_changed = !HistoryEqual(target.camera, history.current.camera);
      std::function<void()> commit_scene, commit_sky;
      if (scene_changed && scene_editor) {
        commit_scene = scene_editor->PrepareRestore(*target.scene, target.object.id);
      }
      if (sky_changed) {
        if (!prepare_sky_restore) {
          throw std::runtime_error("Sky history is unavailable.");
        }
        commit_sky = prepare_sky_restore(target.sky);
      }
      // Prepare new thumbnails before replacing existing ones. Stable keyframe
      // IDs reuse their textures; adding/deleting a frame needs at most one upload.
      std::vector<RayrenderGui::KeyframeSnapshot> thumbnails;
      std::vector<Rcpp::List> keyframes;
      if (animation_changed) {
        thumbnails = target.animation->thumbnails;
        keyframes = target.animation->keyframes;
        for (auto& snapshot : thumbnails) {
          const auto found = std::find_if(gui.keyframe_snapshots.begin(),
                                          gui.keyframe_snapshots.end(),
                                          [&](const auto& value) {
                                            return value.id == snapshot.id;
                                          });
          snapshot.width = snapshot.image->width;
          snapshot.height = snapshot.image->height;
          snapshot.pixels = snapshot.image->pixels;
          if (found != gui.keyframe_snapshots.end()) {
            snapshot.texture = found->texture;
          } else if (gui.api && gui.session && !snapshot.pixels.empty()) {
            rayimgui_image_v1 image{sizeof(image),
                                    snapshot.width,
                                    snapshot.height,
                                    RAYIMGUI_TOP_LEFT,
                                    RAYIMGUI_RGBA,
                                    RAYIMGUI_OPAQUE,
                                    RAYIMGUI_DISPLAY_ENCODED,
                                    uint64_t(snapshot.width) * 4,
                                    uint64_t(snapshot.pixels.size()),
                                    1,
                                    snapshot.pixels.data()};
            gui.check(gui.api->texture_create(
                gui.session, &image, &snapshot.texture, &gui.error));
            created.push_back(snapshot.texture);
          }
        }
      }
      // Both expensive scene/sky preparations succeeded. Publish the restored
      // values at the same checkpoint before sampling can resume.
      if (commit_scene) {
        commit_scene();
      }
      if (commit_sky) {
        commit_sky();
      }
      // A draft-only undo can switch the inspected object without rebuilding
      // geometry. Rebind its committed pivot before the next transform edit.
      if (scene_editor && !scene_changed && gui.object.id != target.object.id) {
        PreviewObjectState inspector;
        scene_editor->Describe(target.object.id, inspector);
      }
      ApplyCameraState(target.camera, &native_env_angle);
      if (animation_changed) {
        for (auto& snapshot : gui.keyframe_snapshots) {
          const bool retained =
              std::any_of(thumbnails.begin(), thumbnails.end(), [&](const auto& value) {
                return value.id == snapshot.id;
              });
          if (!retained && snapshot.texture && gui.api && gui.session) {
            gui.api->texture_destroy(gui.session, &snapshot.texture, nullptr);
          }
        }
        gui.keyframe_snapshots.swap(thumbnails);
        Keyframes.swap(keyframes);
      }
      current_keyframe = target.keyframe;
      keyframe_motion_closed = target.closed;
      SetCameraMotionBlur(target.blur);
      SetShutterSpeed(target.shutter);
      ApplyStaticPreviewCameraMotionRange(cam);
#ifdef HAS_OIDN
      if (denoise != target.denoise) {
        denoise = target.denoise;
        has_denoised_preview = false;
        denoised_preview_sample_count = 0;
      }
#endif
      reset = scene_changed || sky_changed || camera_changed ||
              write_fast_output != target.fast || history.current.blur != target.blur ||
              history.current.shutter != target.shutter || reset;
      write_fast_output = target.fast;
      RestoreHistoryInputs(gui, target);
      if (sky_changed) {
        sky_model = Rcpp::as<int>(target.sky["model"]);
        sun_elevation = Rcpp::as<double>(target.sky["elevation"]);
        sun_azimuth = Rcpp::as<double>(target.sky["azimuth"]);
        atmosphere_haze = gui.has_atmosphere && Rcpp::as<bool>(target.sky["haze"]);
        atmosphere_query_altitude =
            gui.has_atmosphere && Rcpp::as<bool>(target.sky["altitude"]);
        atmosphere_changed = true;
      }
      history.current = target;
      history.scene_revision = scene_editor ? scene_editor->Revision() : 0;
      history.cursor = direction < 0 ? index : index + 1;
      ++gui.history_epoch; // A new edit after undo starts a new branch.
      gui.history_message.clear();
      SyncNativeAnimationState();
    } catch (const Rcpp::internal::InterruptedException&) {
      for (auto texture : created) {
        gui.api->texture_destroy(gui.session, &texture, nullptr);
      }
      throw;
    } catch (const std::exception& error) {
      for (auto texture : created) {
        gui.api->texture_destroy(gui.session, &texture, nullptr);
      }
      gui.history_message =
          std::string(direction < 0 ? "Undo failed: " : "Redo failed: ") + error.what();
      break; // Keep the cursor and live state intact when preparation fails.
    }
  }
  gui.can_undo = history.cursor > 0;
  gui.can_redo = history.cursor < history.entries.size();
  gui.history_dirty = false;
  return reset;
}
#endif
