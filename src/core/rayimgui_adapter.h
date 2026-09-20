/* Copyright (c) 2026 Tyler Morgan-Wall. MIT; LICENSE.protocol.
 * Consumer-owned native state. No provider or GUI link dependency.
 */
#ifndef RAYRENDER_RAYIMGUI_ADAPTER_H
#define RAYRENDER_RAYIMGUI_ADAPTER_H
#include "rayimgui/rayimgui_r.h"
#include "preview_object_state.h"
#include "preview_selection_overlay.h"
#include "../lights/sun_direction.h"
#include <algorithm>
#include <chrono>
#include <cstring>
#include <cstdio>
#include <string>
#include <stdexcept>
#include <vector>
#include <map>
#include <functional>
#include <cmath>
#include <memory>

struct RayrenderGui {
  const rayimgui_api_v1* api = nullptr;
  rayimgui_session_handle session = 0;
  rayimgui_callback_handle callback = 0;
  rayimgui_texture_handle texture = 0, selection_texture = 0;
  uint64_t selection_id = 0, selection_revision = 0;
  uint32_t selection_width = 0, selection_height = 0;
  std::vector<double> selection_camera;
  std::vector<uint8_t> selection_mask, selection_outline, selection_pixels;
  // Visibility is separate from the cache key: restarting accumulation or
  // temporarily hiding an overlay must not retrace unchanged geometry.
  bool selection_visible = false;
  char message[512]{};
  rayimgui_error_v1 error{sizeof(error), 0, message, sizeof(message)};
  std::vector<uint8_t> pixels;
  std::vector<float> display_rgb;
  uint32_t width = 0, height = 0;
  uint64_t version = 0;
  size_t samples = 0, objects = 0;
  double render_fps = 0;
  size_t fps_frames = 0;
  std::chrono::steady_clock::time_point fps_window_start{};
  double progress = 0, exposure = 1;
  int32_t haze = 0, altitude = 0;
  bool can_edit = false, has_atmosphere = false, deferred = false;
  bool close_requested = false, request_render = false, request_reset = false;
  bool atmosphere_pending = false, exposure_pending = false;
  bool has_sun = false, sun_pending = false, sun_editing = false;
  bool has_sky_model = false, sky_model_pending = false;
  int32_t sky_model = 0;
  double sun_elevation = 0, sun_azimuth = 0;
  bool has_location = false, location_pending = false, manual_sun = false;
  double latitude = 0, longitude = 0;
  char datetime[32]{};
  int32_t date[3]{2026, 1, 1}, time[3]{0, 0, 0};
  int32_t fast_preview = 0;
  bool fast_pending = false;
  int32_t denoise_enabled = 0;
  bool denoise_available = false, denoise_pending = false;
  std::string sky_error;
  bool export_available = false, export_pending = false;
  char export_filename[4096]{"rayrender_scene.R"};
  std::string export_message;
  // These are plain GUI snapshots and queued requests. Renderer/R operations
  // happen only after sampling workers finish, never inside draw callbacks.
  enum class AnimationAction { None, Save, Previous, Next, Delete, Play, Select };
  AnimationAction animation_action = AnimationAction::None;
  uint64_t animation_target = 0, next_keyframe_id = 1;
  int animation_current = -1;
  bool animation_playing = false, animation_settings_pending = false;
  int32_t animation_blur = 0, animation_closed = 0;
  double animation_shutter = 1;
  std::string animation_message;
  // History shares immutable thumbnail pixels; it never owns provider handles.
  struct KeyframeImage {
    std::vector<uint8_t> pixels;
    uint32_t width = 0, height = 0;
  };
  struct KeyframeSnapshot {
    uint64_t id = 0;
    std::vector<double> camera;
    std::vector<uint8_t> pixels;
    uint32_t width = 0, height = 0;
    rayimgui_texture_handle texture = 0;
    std::shared_ptr<KeyframeImage> image = std::make_shared<KeyframeImage>();
  };
  // The provider allows 128 live textures; reserve space for viewport/overlays.
  static constexpr size_t MaxKeyframeSnapshots = 120;
  std::vector<KeyframeSnapshot> keyframe_snapshots;
  std::vector<double> snapshot_camera;

  bool history_dirty = false, can_undo = false, can_redo = false;
  uint64_t history_epoch = 0;
  std::vector<int> history_requests;
  std::string history_message;

  // Widget activation identifies a gesture, so repeated drag frames coalesce.
  // All notifications remain value-only; history snapshots live in the renderer.
  void track_input(const rayimgui_widget_v1& widget, const rayimgui_item_v1& item) {
    const bool editable =
        widget.kind == RAYIMGUI_DOUBLE || widget.kind == RAYIMGUI_INT ||
        widget.kind == RAYIMGUI_COLOR || widget.kind == RAYIMGUI_INPUT_TEXT ||
        widget.kind == RAYIMGUI_CHECKBOX || widget.kind == RAYIMGUI_COMBO ||
        widget.kind == RAYIMGUI_ANGLE;
    if ((editable || widget.kind == RAYIMGUI_BUTTON) && (item.flags & RAYIMGUI_BEGIN)) {
      ++history_epoch;
    }
    if (editable && (item.flags & RAYIMGUI_CHANGED)) {
      history_dirty = true;
    }
  }

  int32_t orbit = 1;
  double movement_speed = 1;
  struct Key {
    unsigned code, modifiers, count;
  };
  std::vector<Key> keys;
  bool pick_pending = false, pick_focus = false;
  PreviewObjectState object;
  float pick_u = 0, pick_v = 0;
  std::chrono::steady_clock::time_point last_step{};

  static bool same_camera(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size() || a.empty()) {
      return false;
    }
    for (size_t i = 0; i < a.size(); ++i) {
      if (std::abs(a[i] - b[i]) > 1e-6 * std::max(1.0, std::abs(a[i]))) {
        return false;
      }
    }
    return true;
  }

  // A keyframe keeps its first matching completed preview, without the selection
  // outline or gizmo. Camera movement in the same input batch defers capture
  // until that camera has actually rendered, avoiding a thumbnail of the old view.
  void capture_keyframe_snapshots() {
    if (!width || !height || pixels.size() != size_t(width) * height * 4) {
      return;
    }
    for (auto& snapshot : keyframe_snapshots) {
      if (!snapshot.pixels.empty() || !same_camera(snapshot.camera, snapshot_camera)) {
        continue;
      }
      const double scale = std::min({1.0, 144.0 / width, 90.0 / height});
      snapshot.width = std::max(1u, uint32_t(width * scale));
      snapshot.height = std::max(1u, uint32_t(height * scale));
      snapshot.pixels.resize(size_t(snapshot.width) * snapshot.height * 4);
      // Box averaging suppresses single-sample noise and thin-edge aliasing.
      // The small immutable buffers also bound CPU/GPU memory per keyframe.
      for (uint32_t y = 0; y < snapshot.height; ++y) {
        const uint32_t y0 = uint64_t(y) * height / snapshot.height;
        const uint32_t y1 = uint64_t(y + 1) * height / snapshot.height;
        for (uint32_t x = 0; x < snapshot.width; ++x) {
          const uint32_t x0 = uint64_t(x) * width / snapshot.width;
          const uint32_t x1 = uint64_t(x + 1) * width / snapshot.width;
          uint64_t sum[3]{};
          for (uint32_t sy = y0; sy < y1; ++sy) {
            for (uint32_t sx = x0; sx < x1; ++sx) {
              for (size_t c = 0; c < 3; ++c) {
                sum[c] += pixels[4 * (size_t(sy) * width + sx) + c];
              }
            }
          }
          const size_t offset = 4 * (size_t(y) * snapshot.width + x);
          for (size_t c = 0; c < 3; ++c) {
            snapshot.pixels[offset + c] = sum[c] / (uint64_t(x1 - x0) * (y1 - y0));
          }
          snapshot.pixels[offset + 3] = 255;
        }
      }
      snapshot.image->pixels = snapshot.pixels;
      snapshot.image->width = snapshot.width;
      snapshot.image->height = snapshot.height;
      if (api && session) {
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
        check(api->texture_create(session, &image, &snapshot.texture, &error));
      }
    }
  }

  void save_keyframe_snapshot(std::vector<double> camera) {
    KeyframeSnapshot snapshot;
    snapshot.id = next_keyframe_id++;
    snapshot.camera = std::move(camera);
    keyframe_snapshots.push_back(std::move(snapshot));
    capture_keyframe_snapshots();
  }

  void delete_keyframe_snapshot(size_t index) {
    if (index >= keyframe_snapshots.size()) {
      return;
    }
    auto& snapshot = keyframe_snapshots[index];
    if (snapshot.texture && api && session) {
      check(api->texture_destroy(session, &snapshot.texture, &error));
    }
    keyframe_snapshots.erase(keyframe_snapshots.begin() + index);
  }

  void set_datetime(const std::string& value) {
    std::snprintf(datetime, sizeof(datetime), "%s", value.c_str());
    int year, month, day, hour, minute, second;
    if (std::sscanf(datetime,
                    "%d-%d-%d %d:%d:%d",
                    &year,
                    &month,
                    &day,
                    &hour,
                    &minute,
                    &second) == 6) {
      date[0] = year;
      date[1] = month;
      date[2] = day;
      time[0] = hour;
      time[1] = minute;
      time[2] = second;
    }
  }

  // Numeric date fields can temporarily describe impossible dates while dragging.
  // Clamp them to a real calendar date before formatting the UTC update request.
  void format_datetime() {
    date[0] = std::clamp(date[0], 1, 9999);
    date[1] = std::clamp(date[1], 1, 12);
    const bool leap = date[0] % 4 == 0 && (date[0] % 100 != 0 || date[0] % 400 == 0);
    const int days[] = {31, leap ? 29 : 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    date[2] = std::clamp(date[2], 1, days[date[1] - 1]);
    time[0] = std::clamp(time[0], 0, 23);
    time[1] = std::clamp(time[1], 0, 59);
    time[2] = std::clamp(time[2], 0, 59);
    std::snprintf(datetime,
                  sizeof(datetime),
                  "%04d-%02d-%02d %02d:%02d:%02d",
                  date[0],
                  date[1],
                  date[2],
                  time[0],
                  time[1],
                  time[2]);
  }

  int32_t acquire(SEXP handle) {
    int32_t status = rayimgui_api_from_R_v1(handle,
                                            sizeof(rayimgui_api_v1),
                                            RAYIMGUI_CAP_RGBA8 | RAYIMGUI_CAP_WIDGETS |
                                                RAYIMGUI_CAP_INPUT | RAYIMGUI_CAP_GIZMO,
                                            &api);
    if (status || !api->input || api->header.abi_minor < 8) {
      std::snprintf(message,
                    sizeof(message),
                    "Editor panels require rayimgui 0.0.12 or later (ABI 1.8).");
      return status ? status : RAYIMGUI_ABI;
    }
    if (!(api->header.capabilities & RAYIMGUI_CAP_NATIVE)) {
      std::snprintf(message,
                    sizeof(message),
                    "rayimgui was built without a native window backend.");
      return RAYIMGUI_BACKEND_NOT_BUILT;
    }
    return RAYIMGUI_OK;
  }

  static bool expected_unavailability(int32_t status) {
    return status == RAYIMGUI_ABI || status == RAYIMGUI_BACKEND_NOT_BUILT ||
           status == RAYIMGUI_BACKEND_INIT;
  }

  int32_t open(uint32_t w, uint32_t h) {
    width = w;
    height = h;
    // Window-level Escape also works while a side panel, text field or gizmo
    // owns input. The provider reports a close; renderer workers finish first.
    rayimgui_session_desc_v1 desc{sizeof(desc),
                                  RAYIMGUI_HISTORY_SHORTCUTS | RAYIMGUI_ESCAPE_CLOSE,
                                  "rayrender",
                                  9,
                                  1440,
                                  900};
    int32_t status = api->open(&desc, &session, &error);
    if (status) {
      return status;
    }
    rayimgui_callback_v1 client{sizeof(client), 1, draw, this};
    status = api->register_callback(session, &client, &callback, &error);
    if (status) {
      close();
    }
    return status;
  }

  void check(int32_t status) {
    if (status) {
      throw std::runtime_error(std::string("rayrender native GUI: ") +
                               (message[0] ? message : "provider operation failed"));
    }
  }

  // Upload a complete display snapshot; the provider copies the pixel buffer
  // synchronously, so later frames may safely resize or reuse this vector.
  void publish() {
    // Count completed preview images, excluding GUI polls and selection-overlay
    // uploads. Average over half a second so short frame-time variations settle.
    const auto now = std::chrono::steady_clock::now();
    if (fps_window_start == std::chrono::steady_clock::time_point{}) {
      fps_window_start = now;
    } else {
      ++fps_frames;
      const double elapsed =
          std::chrono::duration<double>(now - fps_window_start).count();
      if (elapsed >= 0.5) {
        render_fps = fps_frames / elapsed;
        fps_frames = 0;
        fps_window_start = now;
      }
    }

    rayimgui_image_v1 image{sizeof(image),
                            width,
                            height,
                            RAYIMGUI_TOP_LEFT,
                            RAYIMGUI_RGBA,
                            RAYIMGUI_OPAQUE,
                            RAYIMGUI_DISPLAY_ENCODED,
                            static_cast<uint64_t>(width) * 4,
                            static_cast<uint64_t>(pixels.size()),
                            ++version,
                            pixels.data()};
    check(texture ? api->texture_update(session, texture, &image, &error)
                  : api->texture_create(session, &image, &texture, &error));
    // The selection outline is independent of render colors. Its texture is
    // refreshed only when UpdateNativeSelectionMask updates or restores it.
  }

  // Upload the cached outline as an independent layer with a transparent interior.
  // The provider remains generic: it only composites this RGBA texture.
  void publish_selection() {
    PreviewSelectionOverlay::Compose(width,
                                     height,
                                     selection_outline,
                                     selection_width,
                                     selection_height,
                                     selection_pixels);
    rayimgui_image_v1 image{sizeof(image),
                            width,
                            height,
                            RAYIMGUI_TOP_LEFT,
                            RAYIMGUI_RGBA,
                            RAYIMGUI_STRAIGHT,
                            RAYIMGUI_DISPLAY_ENCODED,
                            uint64_t(width) * 4,
                            uint64_t(selection_pixels.size()),
                            ++version,
                            selection_pixels.data()};
    check(selection_texture
              ? api->texture_update(session, selection_texture, &image, &error)
              : api->texture_create(session, &image, &selection_texture, &error));
  }

  // Worker-wait loops may poll frequently. Limit GUI frame frequency there,
  // while force allows a completed renderer snapshot to be shown immediately.
  bool poll(bool force = false) {
    if (!session || close_requested) {
      return close_requested;
    }
    auto now = std::chrono::steady_clock::now();
    if (!force && now - last_step < std::chrono::milliseconds(16)) {
      return false;
    }
    last_step = now;
    rayimgui_step_v1 result{};
    result.struct_size = sizeof(result);
    check(api->step(session, &result, &error));
    // Shortcut events apply after worker jobs finish. Do not mutate the scene
    // from the provider callback, including when polling while workers run.
    if (api->event_poll) {
      rayimgui_event_v1 event{sizeof(event)};
      int32_t status;
      while ((status = api->event_poll(session, &event, &error)) == RAYIMGUI_OK) {
        if (event.application_id == 0 && event.phase == RAYIMGUI_COMMIT &&
            history_requests.size() < 32) {
          if (event.command == RAYIMGUI_COMMAND_UNDO) {
            history_requests.push_back(-1);
          }
          if (event.command == RAYIMGUI_COMMAND_REDO) {
            history_requests.push_back(1);
          }
        }
      }
      if (status != RAYIMGUI_EMPTY) {
        check(status);
      }
    }
    close_requested = result.close_requested != 0;
    return close_requested;
  }

  void close() noexcept {
    // Also called by explicit R_UnwindProtect cleanup, never a finalizer.
    if (session && api) {
      api->close(&session, nullptr);
    }
    // Closing the provider session destroys every thumbnail texture as well.
    keyframe_snapshots.clear();
    snapshot_camera.clear();
    callback = 0;
    texture = selection_texture = 0;
    api = nullptr;
  }

  // Respect provider focus/ownership decisions and queue intent only. These
  // callbacks can run while sampling workers still read the active scene.
  void collect_input(const rayimgui_input_v1& input,
                     const rayimgui_viewport_v1& viewport) {
    if (input.keyboard_available &&
        !(input.modifiers & (RAYIMGUI_CTRL | RAYIMGUI_ALT | RAYIMGUI_SUPER))) {
      // Limit queued motion while a slow sample renders. No renderer or R
      // state is touched from a GUI callback or worker-wait poll.
      const uint64_t repeat =
          RAYIMGUI_KEY_BIT(RAYIMGUI_KEY_W) | RAYIMGUI_KEY_BIT(RAYIMGUI_KEY_A) |
          RAYIMGUI_KEY_BIT(RAYIMGUI_KEY_S) | RAYIMGUI_KEY_BIT(RAYIMGUI_KEY_D) |
          RAYIMGUI_KEY_BIT(RAYIMGUI_KEY_Q) | RAYIMGUI_KEY_BIT(RAYIMGUI_KEY_Z) |
          RAYIMGUI_KEY_BIT(RAYIMGUI_KEY_UP) | RAYIMGUI_KEY_BIT(RAYIMGUI_KEY_DOWN) |
          RAYIMGUI_KEY_BIT(RAYIMGUI_KEY_LEFT) | RAYIMGUI_KEY_BIT(RAYIMGUI_KEY_RIGHT) |
          RAYIMGUI_KEY_BIT(RAYIMGUI_KEY_1) | RAYIMGUI_KEY_BIT(RAYIMGUI_KEY_2) |
          RAYIMGUI_KEY_BIT(RAYIMGUI_KEY_3) | RAYIMGUI_KEY_BIT(RAYIMGUI_KEY_4);
      uint64_t pressed = input.keys_pressed | (input.keys_repeated & repeat);
      if (can_edit && pressed) {
        if (input.keys_pressed) {
          ++history_epoch;
        }
        history_dirty = true;
      }
      if (can_edit) {
        for (unsigned code = 0; code < RAYIMGUI_KEY_COUNT; ++code) {
          if (pressed & RAYIMGUI_KEY_BIT(code)) {
            if (!keys.empty() && keys.back().code == code &&
                keys.back().modifiers == input.modifiers) {
              keys.back().count = std::min(keys.back().count + 1, 8u);
            } else if (keys.size() < 64) {
              keys.push_back({code, input.modifiers, 1});
            }
          }
        }
      }
    } else {
      keys.clear();
    }
    if (can_edit && input.mouse_available &&
        (input.mouse_clicked & (RAYIMGUI_MOUSE_LEFT | RAYIMGUI_MOUSE_RIGHT))) {
      // Undo the displayed image's two-axis flip to obtain renderer film coordinates.
      const float u = 1 - viewport.mouse_x / viewport.width,
                  v = 1 - viewport.mouse_y / viewport.height;
      if (object.enabled && (input.modifiers & RAYIMGUI_SHIFT) &&
          !(input.modifiers & RAYIMGUI_ALT) &&
          (input.mouse_clicked & RAYIMGUI_MOUSE_LEFT)) {
        object.pick_pending = true;
        object.pick_u = u;
        object.pick_v = v;
        pick_pending = false;
      } else {
        ++history_epoch;
        history_dirty = true;
        pick_pending = true;
        pick_focus = (input.mouse_clicked & RAYIMGUI_MOUSE_LEFT) != 0;
        pick_u = u;
        pick_v = v;
      }
    }
  }

  static rayimgui_widget_v1 widget(uint32_t kind, uint64_t id, const char* label) {
    rayimgui_widget_v1 value{};
    value.struct_size = sizeof(value);
    value.kind = kind;
    value.id = id;
    value.label = label;
    value.label_length = static_cast<uint32_t>(std::strlen(label));
    value.speed = .01;
    return value;
  }

  static int32_t draw(const rayimgui_api_v1* api, rayimgui_session_handle session,
                      void* userdata, rayimgui_error_v1* error) noexcept {
    try {
      return static_cast<RayrenderGui*>(userdata)->draw_impl(api, session, error);
    } catch (...) {
      return RAYIMGUI_CALLBACK_ERROR;
    }
  }

  // Widgets edit value snapshots and request flags. Scene rebuilding and material
  // setters run later at the renderer checkpoint, outside this provider callback.
  int32_t draw_object_panel(const rayimgui_api_v1* table, rayimgui_session_handle owner,
                            rayimgui_error_v1* err) {
    auto& o = object;
    rayimgui_item_v1 item{sizeof(item)};
    rayimgui_widget_v1 w{};
    int32_t status;
#define OBJECT_DRAW()                                                                  \
  do {                                                                                 \
    status = table->widget(owner, &w, &item, err);                                     \
    track_input(w, item);                                                              \
    if (status)                                                                        \
      return status;                                                                   \
  } while (0)
#define OBJECT_SIMPLE(kind, id, label)                                                 \
  do {                                                                                 \
    w = widget(kind, id, label);                                                       \
    OBJECT_DRAW();                                                                     \
  } while (0)
    OBJECT_SIMPLE(RAYIMGUI_TEXT, 0, "Shift-click viewport or select in tree.");
    if (o.selected) {
      OBJECT_SIMPLE(RAYIMGUI_TEXT, 0, o.label.c_str());
      OBJECT_SIMPLE(RAYIMGUI_BUTTON, 7000, "Deselect");
      o.clear_pending = o.clear_pending || (item.flags & RAYIMGUI_CHANGED);
      static const char* operations[] = {"Translate", "Rotate", "Scale"};
      static const uint32_t lengths[] = {9, 6, 5};
      w = widget(RAYIMGUI_COMBO, 7001, "Tool");
      w.integers = &o.operation;
      w.count = 1;
      w.choices = operations;
      w.choice_lengths = lengths;
      w.choice_count = 3;
      OBJECT_DRAW();
      static const char* modes[] = {"Local", "World"};
      static const uint32_t mode_lengths[] = {5, 5};
      w = widget(RAYIMGUI_COMBO, 7002, "Axes");
      w.integers = &o.mode;
      w.count = 1;
      w.choices = modes;
      w.choice_lengths = mode_lengths;
      w.choice_count = 2;
      OBJECT_DRAW();
      const char* labels[] = {"Position", "Rotation (degrees)", "Scale"};
      double* values[] = {o.translation.data(), o.rotation.data(), o.scale.data()};
      for (int i = 0; i < 3; ++i) {
        w = widget(RAYIMGUI_DOUBLE, 7003 + i, labels[i]);
        w.numbers = values[i];
        w.count = 3;
        w.speed = i == 1 ? .25 : .01;
        if (i == 2) {
          w.minimum = -10000;
          w.maximum = 10000;
        }
        OBJECT_DRAW();
        if (item.flags & RAYIMGUI_CHANGED) {
          o.transform_pending = true;
          o.numeric_transform = true;
        }
      }
      OBJECT_SIMPLE(RAYIMGUI_BUTTON, 7006, "Apply transform");
      o.apply_transform = o.apply_transform || (item.flags & RAYIMGUI_CHANGED);
      if (o.projection_valid) {
        OBJECT_SIMPLE(RAYIMGUI_TEXT, 0, "Drag a handle; release to apply.");
      } else {
        OBJECT_SIMPLE(RAYIMGUI_TEXT, 0, "Use numeric transforms for this camera.");
      }
      OBJECT_SIMPLE(RAYIMGUI_SEPARATOR, 0, "");
      OBJECT_SIMPLE(RAYIMGUI_TEXT, 0, "Material");
      if (!o.materials.empty()) {
        if (o.materials.size() > 1 && o.materials.size() <= 1024) {
          std::vector<const char*> labels;
          std::vector<uint32_t> sizes;
          for (const auto& panel : o.materials) {
            labels.push_back(panel.label.c_str());
            sizes.push_back(panel.label.size());
          }
          w = widget(RAYIMGUI_COMBO, 7010, "Slot");
          w.integers = &o.material_slot;
          w.count = 1;
          w.choices = labels.data();
          w.choice_lengths = sizes.data();
          w.choice_count = labels.size();
          OBJECT_DRAW();
        } else if (o.materials.size() > 1024) {
          w = widget(RAYIMGUI_INT, 7010, "Slot (zero based)");
          w.integers = &o.material_slot;
          w.count = 1;
          w.minimum = 0;
          w.maximum = o.materials.size() - 1;
          w.speed = 1;
          OBJECT_DRAW();
        }
        o.material_slot = std::clamp(
            o.material_slot, 0, static_cast<int32_t>(o.materials.size() - 1));
        auto& panel = o.materials[o.material_slot];
        OBJECT_SIMPLE(RAYIMGUI_TEXT, 0, panel.type.c_str());
        if (panel.targets.size() > 1) {
          const std::string scope =
              "Applies to " + std::to_string(panel.targets.size()) + " instances";
          OBJECT_SIMPLE(RAYIMGUI_TEXT, 0, scope.c_str());
        }
        std::string section;
        for (size_t i = 0; i < panel.fields.size(); ++i) {
          auto& field = panel.fields[i];
          if (!field.condition.empty()) {
            const auto control = std::find_if(
                panel.fields.begin(), panel.fields.end(), [&](const auto& value) {
                  return value.name == field.condition;
                });
            if (control == panel.fields.end() ||
                std::find(field.visible_choices.begin(),
                          field.visible_choices.end(),
                          static_cast<int>(control->values[0])) ==
                    field.visible_choices.end()) {
              continue;
            }
          }
          if (field.section != section) {
            section = field.section;
            OBJECT_SIMPLE(RAYIMGUI_SEPARATOR, 0, "");
            OBJECT_SIMPLE(RAYIMGUI_TEXT, 0, section.c_str());
          }
          const std::string label = field.name + (field.mixed ? " (mixed)" : "");
          uint32_t kind = RAYIMGUI_DOUBLE;
          if (field.text_input) {
            kind = RAYIMGUI_INPUT_TEXT;
          } else if (!field.choices.empty()) {
            kind = RAYIMGUI_COMBO;
          } else if (field.boolean) {
            kind = RAYIMGUI_CHECKBOX;
          } else if (field.integer) {
            kind = RAYIMGUI_INT;
          } else if (field.color) {
            kind = RAYIMGUI_COLOR;
          }
          // Long material names stay readable when the inspector is narrow.
          if (!field.boolean) {
            OBJECT_SIMPLE(RAYIMGUI_TEXT, 0, label.c_str());
          }
          w = widget(kind,
                     100000 + 1000 * uint64_t(o.material_slot) + i,
                     field.boolean ? label.c_str() : "");
          int32_t integer_value = 0;
          if (field.boolean || field.integer || !field.choices.empty()) {
            integer_value = static_cast<int32_t>(
                std::clamp(field.values[0], double(INT32_MIN), double(INT32_MAX)));
          }
          std::array<char, 4096> text{};
          std::snprintf(text.data(), text.size(), "%s", field.text.c_str());
          std::vector<const char*> choices;
          std::vector<uint32_t> lengths;
          for (const auto& choice : field.choices) {
            choices.push_back(choice.c_str());
            lengths.push_back(choice.size());
          }
          w.numbers = field.values.data();
          w.integers = &integer_value;
          w.text = text.data();
          w.text_capacity = text.size();
          w.choices = choices.data();
          w.choice_lengths = lengths.data();
          w.choice_count = choices.size();
          w.count = field.count;
          w.minimum = field.minimum;
          w.maximum = field.maximum;
          w.speed = field.speed;
          if (field.color) {
            w.options = RAYIMGUI_COLOR_RGB;
            w.fraction = 1;
          }
          OBJECT_DRAW();
          if (item.flags & RAYIMGUI_CHANGED) {
            if (field.text_input) {
              field.text = text.data();
            } else if (field.boolean || field.integer || !field.choices.empty()) {
              field.values[0] = integer_value;
            }
            field.changed = true;
            field.mixed = false;
            panel.changed = true;
            o.material_pending = true;
          }
          if (field.text_input) {
            OBJECT_SIMPLE(
                RAYIMGUI_TOOLTIP,
                0,
                "Enter a file path. Blank uses the scene's original or embedded map.");
          }
        }
        if (panel.fields.empty()) {
          OBJECT_SIMPLE(RAYIMGUI_TEXT, 0, "No editable parameters for this material.");
        }
        OBJECT_SIMPLE(RAYIMGUI_BUTTON, 7011, "Apply material");
        o.apply_material =
            o.apply_material || ((item.flags & RAYIMGUI_CHANGED) && o.material_pending);
      }
      OBJECT_SIMPLE(RAYIMGUI_BUTTON, 7012, "Reset object edits");
      o.revert = o.revert || (item.flags & RAYIMGUI_CHANGED);
      if (!o.error.empty()) {
        OBJECT_SIMPLE(RAYIMGUI_TEXT, 0, o.error.c_str());
      }
    }

#undef OBJECT_SIMPLE
#undef OBJECT_DRAW
    return RAYIMGUI_OK;
  }

  // Draw a value-only scene snapshot. Selecting a node queues the same stable
  // object ID used by viewport picking; renderer mutations wait for a checkpoint.
  int32_t draw_hierarchy(const rayimgui_api_v1* table, rayimgui_session_handle owner,
                         rayimgui_error_v1* err) {
    if (object.hierarchy.empty()) {
      auto w = widget(RAYIMGUI_TEXT, 0, "No scene objects");
      rayimgui_item_v1 item{sizeof(item)};
      return table->widget(owner, &w, &item, err);
    }
    std::map<uint64_t, std::vector<const PreviewHierarchyNode*>> children;
    for (const auto& node : object.hierarchy) {
      children[node.parent].push_back(&node);
    }
    std::function<int32_t(uint64_t)> draw_children = [&](uint64_t parent) -> int32_t {
      for (const auto* node : children[parent]) {
        auto w = widget(RAYIMGUI_TREE_BEGIN, node->id, node->label.c_str());
        const bool leaf = children.find(node->id) == children.end();
        w.options = RAYIMGUI_TREE_DEFAULT_OPEN;
        if (leaf) {
          w.options |= RAYIMGUI_TREE_LEAF;
        }
        if (object.selected && object.id == node->id) {
          w.options |= RAYIMGUI_TREE_SELECTED;
        }
        rayimgui_item_v1 item{sizeof(item)};
        int32_t status = table->widget(owner, &w, &item, err);
        track_input(w, item);
        if (status) {
          return status;
        }
        if ((item.flags & RAYIMGUI_CHANGED) && can_edit && object.enabled) {
          object.select_id = node->id;
          object.select_pending = true;
        }
        if (!leaf && (item.flags & RAYIMGUI_VISIBLE)) {
          status = draw_children(node->id);
          if (status) {
            return status;
          }
        }
        w = widget(RAYIMGUI_TREE_END, 0, "");
        status = table->widget(owner, &w, &item, err);
        track_input(w, item);
        if (status) {
          return status;
        }
      }
      return RAYIMGUI_OK;
    };
    return draw_children(0);
  }

  int32_t draw_animation(const rayimgui_api_v1* table, rayimgui_session_handle owner,
                         rayimgui_error_v1* err) {
    rayimgui_item_v1 item{sizeof(item), 0};
    auto submit = [&](rayimgui_widget_v1 w) {
      const auto status = table->widget(owner, &w, &item, err);
      track_input(w, item);
      return status;
    };
    int32_t status = 0;
#define ANIMATION_DRAW(w)                                                              \
  do {                                                                                 \
    status = submit(w);                                                                \
    if (status)                                                                        \
      return status;                                                                   \
  } while (0)
    auto simple = [&](uint32_t kind, uint64_t id, const char* label) {
      return widget(kind, id, label);
    };
    auto button = [&](uint64_t id,
                      const char* label,
                      AnimationAction action,
                      bool disabled) -> int32_t {
      auto scope = widget(RAYIMGUI_DISABLED_BEGIN, 0, "");
      int32_t value = disabled;
      scope.integers = &value;
      scope.count = 1;
      ANIMATION_DRAW(scope);
      ANIMATION_DRAW(simple(RAYIMGUI_BUTTON, id, label));
      if (item.flags & RAYIMGUI_CHANGED) {
        animation_action = action;
      }
      ANIMATION_DRAW(simple(RAYIMGUI_DISABLED_END, 0, ""));
      return RAYIMGUI_OK;
    };
    auto w = widget(RAYIMGUI_WINDOW_BEGIN, 40, "Animation");
    w.options = RAYIMGUI_DOCK_BOTTOM;
    ANIMATION_DRAW(w);
    if (item.flags & RAYIMGUI_VISIBLE) {
      const bool locked = !can_edit || animation_playing;
      const bool empty = keyframe_snapshots.empty();
      status = button(41, "Previous", AnimationAction::Previous, locked || empty);
      if (status) {
        return status;
      }
      ANIMATION_DRAW(simple(RAYIMGUI_SAME_LINE, 0, ""));
      status = button(42, "Next", AnimationAction::Next, locked || empty);
      if (status) {
        return status;
      }
      ANIMATION_DRAW(simple(RAYIMGUI_SAME_LINE, 0, ""));
      status = button(43,
                      "Save keyframe",
                      AnimationAction::Save,
                      locked || keyframe_snapshots.size() >= MaxKeyframeSnapshots);
      if (status) {
        return status;
      }
      ANIMATION_DRAW(simple(RAYIMGUI_SAME_LINE, 0, ""));
      status = button(44, "Delete", AnimationAction::Delete, locked || empty);
      if (status) {
        return status;
      }
      ANIMATION_DRAW(simple(RAYIMGUI_SAME_LINE, 0, ""));
      status =
          button(45,
                 animation_playing ? "Stop" : "Play",
                 AnimationAction::Play,
                 !can_edit || (!animation_playing && keyframe_snapshots.size() < 2));
      if (status) {
        return status;
      }
      ANIMATION_DRAW(simple(RAYIMGUI_SAME_LINE, 0, ""));
      char detail[96];
      std::snprintf(detail,
                    sizeof(detail),
                    "Keyframe %d / %zu%s",
                    animation_current + 1,
                    keyframe_snapshots.size(),
                    animation_playing ? " | Playing" : "");
      ANIMATION_DRAW(simple(RAYIMGUI_TEXT, 0, detail));

      w = widget(RAYIMGUI_DISABLED_BEGIN, 0, "");
      int32_t disabled = locked;
      w.integers = &disabled;
      w.count = 1;
      ANIMATION_DRAW(w);
      w = widget(RAYIMGUI_CHECKBOX, 46, "Camera motion blur");
      w.integers = &animation_blur;
      w.count = 1;
      ANIMATION_DRAW(w);
      animation_settings_pending |= (item.flags & RAYIMGUI_CHANGED) != 0;
      ANIMATION_DRAW(simple(RAYIMGUI_SAME_LINE, 0, ""));
      const char* loops[] = {"Open path", "Closed loop"};
      const uint32_t loop_lengths[] = {9, 11};
      w = widget(RAYIMGUI_COMBO, 47, "Path");
      w.integers = &animation_closed;
      w.count = 1;
      w.choices = loops;
      w.choice_lengths = loop_lengths;
      w.choice_count = 2;
      w.width = 130;
      ANIMATION_DRAW(w);
      animation_settings_pending |= (item.flags & RAYIMGUI_CHANGED) != 0;
      ANIMATION_DRAW(simple(RAYIMGUI_SAME_LINE, 0, ""));
      w = widget(RAYIMGUI_DOUBLE, 48, "Shutter (frame fraction)");
      w.numbers = &animation_shutter;
      w.count = 1;
      w.minimum = 0;
      w.maximum = 1;
      w.speed = .001;
      w.width = 120;
      ANIMATION_DRAW(w);
      animation_settings_pending |= (item.flags & RAYIMGUI_CHANGED) != 0;
      ANIMATION_DRAW(simple(RAYIMGUI_DISABLED_END, 0, ""));
      ANIMATION_DRAW(simple(
          RAYIMGUI_TOOLTIP,
          0,
          "Shutter: 0 freezes motion; 1 exposes the full frame interval. Motion blur applies during playback. Closed loop joins the last view to the first."));
      if (!animation_message.empty()) {
        ANIMATION_DRAW(simple(RAYIMGUI_TEXT, 0, animation_message.c_str()));
      }

      w = widget(RAYIMGUI_CHILD_BEGIN, 49, "Keyframe snapshots");
      w.options = RAYIMGUI_CHILD_HORIZONTAL_SCROLL;
      ANIMATION_DRAW(w);
      if (item.flags & RAYIMGUI_VISIBLE) {
        if (empty) {
          ANIMATION_DRAW(simple(
              RAYIMGUI_TEXT,
              0,
              "Move the camera and save a keyframe (K). Save two views to play a path."));
        }
        w = widget(RAYIMGUI_DISABLED_BEGIN, 0, "");
        w.integers = &disabled;
        w.count = 1;
        ANIMATION_DRAW(w);
        for (size_t i = 0; i < keyframe_snapshots.size(); ++i) {
          if (i) {
            ANIMATION_DRAW(simple(RAYIMGUI_SAME_LINE, 0, ""));
          }
          const auto& snapshot = keyframe_snapshots[i];
          std::snprintf(detail,
                        sizeof(detail),
                        "Keyframe %zu%s",
                        i + 1,
                        snapshot.texture ? "" : " (pending)");
          w = widget(snapshot.texture ? RAYIMGUI_IMAGE_BUTTON : RAYIMGUI_BUTTON,
                     100000 + snapshot.id,
                     detail);
          w.texture = snapshot.texture;
          w.width = 144;
          w.height = 90;
          if (int(i) == animation_current) {
            w.options = RAYIMGUI_IMAGE_SELECTED;
          }
          ANIMATION_DRAW(w);
          if (item.flags & RAYIMGUI_CHANGED) {
            animation_target = snapshot.id;
            animation_action = AnimationAction::Select;
          }
        }
        ANIMATION_DRAW(simple(RAYIMGUI_DISABLED_END, 0, ""));
      }
      ANIMATION_DRAW(simple(RAYIMGUI_CHILD_END, 0, ""));
    }
    ANIMATION_DRAW(simple(RAYIMGUI_WINDOW_END, 0, ""));
#undef ANIMATION_DRAW
    return RAYIMGUI_OK;
  }

  int32_t draw_impl(const rayimgui_api_v1* table, rayimgui_session_handle owner,
                    rayimgui_error_v1* err) {
    rayimgui_widget_v1 w{};
    rayimgui_item_v1 item{sizeof(item), 0};
    int32_t status = draw_animation(table, owner, err);
    if (status) {
      return status;
    }
    bool viewport_visible = false;
    sun_editing = false;
#define DRAW()                                                                         \
  do {                                                                                 \
    status = table->widget(owner, &w, &item, err);                                     \
    track_input(w, item);                                                              \
    if (status)                                                                        \
      return status;                                                                   \
  } while (0)
#define SIMPLE(kind, id, label)                                                        \
  do {                                                                                 \
    w = widget(kind, id, label);                                                       \
    DRAW();                                                                            \
  } while (0)
    w = widget(RAYIMGUI_WINDOW_BEGIN, 1, "Render controls");
    w.width = 300;
    w.height = 700;
    w.x = 12;
    w.y = 12;
    w.options = RAYIMGUI_DOCK_LEFT;
    DRAW();
    if (item.flags & RAYIMGUI_VISIBLE) {
      SIMPLE(RAYIMGUI_TEXT, 0, "CPU progressive render");
      int32_t history_disabled = !can_undo;
      w = widget(RAYIMGUI_DISABLED_BEGIN, 0, "");
      w.integers = &history_disabled;
      w.count = 1;
      DRAW();
      SIMPLE(RAYIMGUI_BUTTON, 32, "Undo");
      if (item.flags & RAYIMGUI_CHANGED) {
        history_requests.push_back(-1);
      }
      SIMPLE(RAYIMGUI_DISABLED_END, 0, "");
      SIMPLE(RAYIMGUI_SAME_LINE, 0, "");
      history_disabled = !can_redo;
      w = widget(RAYIMGUI_DISABLED_BEGIN, 0, "");
      w.integers = &history_disabled;
      w.count = 1;
      DRAW();
      SIMPLE(RAYIMGUI_BUTTON, 33, "Redo");
      if (item.flags & RAYIMGUI_CHANGED) {
        history_requests.push_back(1);
      }
      SIMPLE(RAYIMGUI_DISABLED_END, 0, "");
      SIMPLE(RAYIMGUI_TOOLTIP, 0, "Control/Command-Z: undo. Add Shift to redo.");
      if (!history_message.empty()) {
        SIMPLE(RAYIMGUI_TEXT, 0, history_message.c_str());
      }
      char detail[128];
      std::snprintf(detail,
                    sizeof(detail),
                    "%zu objects | %zu samples | %.1f FPS",
                    objects,
                    samples,
                    render_fps);
      SIMPLE(RAYIMGUI_TEXT, 0, detail);
      w = widget(RAYIMGUI_PROGRESS, 2, "Progress");
      w.fraction = progress;
      w.width = -1;
      w.height = 18;
      DRAW();
      w = widget(RAYIMGUI_DOUBLE, 3, "Exposure");
      w.numbers = &exposure;
      w.count = 1;
      w.minimum = .001;
      w.maximum = 1024;
      w.speed = .01;
      DRAW();
      exposure_pending = exposure_pending || (item.flags & RAYIMGUI_CHANGED);
      // Keep the setting visible on builds without denoising support.
      int32_t denoise_disabled = !denoise_available;
      w = widget(RAYIMGUI_DISABLED_BEGIN, 0, "");
      w.integers = &denoise_disabled;
      w.count = 1;
      DRAW();
      w = widget(RAYIMGUI_CHECKBOX, 27, "Denoise");
      w.integers = &denoise_enabled;
      w.count = 1;
      DRAW();
      denoise_pending =
          denoise_pending || (denoise_available && (item.flags & RAYIMGUI_CHANGED));
      SIMPLE(RAYIMGUI_DISABLED_END, 0, "");
      if (!denoise_available) {
        SIMPLE(RAYIMGUI_TEXT, 0, "Denoising unavailable in this build.");
      }
      if (export_available) {
        w = widget(RAYIMGUI_INPUT_TEXT, 30, "Export file");
        w.text = export_filename;
        w.text_capacity = sizeof(export_filename);
        DRAW();
        SIMPLE(RAYIMGUI_BUTTON, 31, "Export R code");
        export_pending = export_pending || (item.flags & RAYIMGUI_CHANGED);
        SIMPLE(RAYIMGUI_TOOLTIP, 0, "Save the committed scene and current view as a runnable R script.");
        if (!export_message.empty()) {
          SIMPLE(RAYIMGUI_TEXT, 0, export_message.c_str());
        }
      }
      if (can_edit) {
        w = widget(RAYIMGUI_CHECKBOX, 26, "Fast preview (F)");
        w.integers = &fast_preview;
        w.count = 1;
        DRAW();
        fast_pending = fast_pending || (item.flags & RAYIMGUI_CHANGED);
        w = widget(RAYIMGUI_CHECKBOX, 11, "Orbit camera");
        w.integers = &orbit;
        w.count = 1;
        DRAW();
        w = widget(RAYIMGUI_DOUBLE, 12, "Movement speed");
        w.numbers = &movement_speed;
        w.count = 1;
        w.minimum = .001;
        w.maximum = 128;
        w.speed = .05;
        DRAW();
        SIMPLE(RAYIMGUI_BUTTON, 4, "Reset camera");
        request_reset = request_reset || (item.flags & RAYIMGUI_CHANGED);
        if (deferred) {
          SIMPLE(RAYIMGUI_BUTTON, 5, "Start final render");
          request_render = request_render || (item.flags & RAYIMGUI_CHANGED);
        }
        if (has_atmosphere || has_sun || has_location || has_sky_model) {
          SIMPLE(RAYIMGUI_SEPARATOR, 0, "");
          SIMPLE(RAYIMGUI_TEXT, 0, "Sky");
          if (has_sky_model) {
            static const char* models[] = {"Hosek", "Prague"};
            static const uint32_t lengths[] = {5, 6};
            w = widget(RAYIMGUI_COMBO, 28, "Sky model");
            w.integers = &sky_model;
            w.count = 1;
            w.choices = models;
            w.choice_lengths = lengths;
            w.choice_count = 2;
            DRAW();
            sky_model_pending = sky_model_pending || (item.flags & RAYIMGUI_CHANGED);
          }
          if (has_sun) {
            w = widget(RAYIMGUI_ANGLE, 13, "Elevation (deg)");
            w.options = RAYIMGUI_ANGLE_QUARTER;
            w.numbers = &sun_elevation;
            w.count = 1;
            w.minimum = -90;
            w.maximum = 90;
            w.speed = .25;
            DRAW();
            // The provider defines a full quarter arc; rayrender keeps its Sun
            // just below that endpoint for both drags and typed values.
            sun_elevation = ClampSunElevation(sun_elevation);
            sun_pending = sun_pending || (item.flags & RAYIMGUI_CHANGED);
            sun_editing = sun_editing || (item.flags & RAYIMGUI_ACTIVE);
            SIMPLE(
                RAYIMGUI_TOOLTIP,
                0,
                "Drag the handle from horizon (0) toward zenith (maximum 89.9). Type negative degrees for below the horizon.");
            // Each angle widget groups its label, dial, and numeric input.
            SIMPLE(RAYIMGUI_SAME_LINE, 0, "");
            w = widget(RAYIMGUI_ANGLE, 14, "Azimuth (deg)");
            w.numbers = &sun_azimuth;
            w.count = 1;
            w.minimum = 0;
            w.maximum = 360;
            w.speed = .5;
            DRAW();
            sun_pending = sun_pending || (item.flags & RAYIMGUI_CHANGED);
            sun_editing = sun_editing || (item.flags & RAYIMGUI_ACTIVE);
            SIMPLE(RAYIMGUI_TOOLTIP,
                   0,
                   "North 0, east 90. Release to update the sky and Sun.");
            if (manual_sun) {
              SIMPLE(RAYIMGUI_TEXT, 0, "Sun direction manually adjusted");
            }
          }
          if (has_location) {
            // Keep the sun dials visible while optional ephemeris inputs are folded.
            SIMPLE(RAYIMGUI_TREE_BEGIN, 29, "Location and date/time");
            if (item.flags & RAYIMGUI_VISIBLE) {
              w = widget(RAYIMGUI_DOUBLE, 16, "Latitude");
              w.numbers = &latitude;
              w.count = 1;
              w.minimum = -90;
              w.maximum = 90;
              w.speed = .1;
              DRAW();
              w = widget(RAYIMGUI_DOUBLE, 17, "Longitude");
              w.numbers = &longitude;
              w.count = 1;
              w.minimum = -180;
              w.maximum = 180;
              w.speed = .1;
              DRAW();
              SIMPLE(RAYIMGUI_TEXT, 0, "Date/time (UTC)");
              const char* labels[] = {
                  "Year", "Month", "Day", "Hour", "Minute", "Second"};
              const int maximum[] = {9999, 12, 31, 23, 59, 59};
              for (unsigned component = 0; component < 6; ++component) {
                w = widget(RAYIMGUI_INT, 20 + component, labels[component]);
                w.integers = component < 3 ? &date[component] : &time[component - 3];
                w.count = 1;
                w.minimum = component < 3 ? 1 : 0;
                w.maximum = maximum[component];
                w.speed = .1;
                DRAW();
                if (item.flags & RAYIMGUI_CHANGED) {
                  format_datetime();
                }
                SIMPLE(
                    RAYIMGUI_TOOLTIP,
                    0,
                    "Drag to change; double-click to type. Then Apply location/time.");
              }
              SIMPLE(RAYIMGUI_TEXT, 0, datetime);
              SIMPLE(RAYIMGUI_BUTTON, 19, "Apply location/time");
              location_pending = location_pending || (item.flags & RAYIMGUI_CHANGED);
            }
            SIMPLE(RAYIMGUI_TREE_END, 0, "");
          }
          if (!sky_error.empty()) {
            SIMPLE(RAYIMGUI_TEXT, 0, sky_error.c_str());
          }
          if (has_atmosphere) {
            w = widget(RAYIMGUI_CHECKBOX, 6, "Haze");
            w.integers = &haze;
            w.count = 1;
            DRAW();
            if (item.flags & RAYIMGUI_CHANGED) {
              if (haze) {
                altitude = 1;
              }
              atmosphere_pending = true;
            }
            w = widget(RAYIMGUI_CHECKBOX, 7, "Altitude queries");
            w.integers = &altitude;
            w.count = 1;
            DRAW();
            if (item.flags & RAYIMGUI_CHANGED) {
              if (!altitude) {
                haze = 0;
              }
              atmosphere_pending = true;
            }
          }
        }
        SIMPLE(RAYIMGUI_TREE_BEGIN, 15, "Movement controls");
        if (item.flags & RAYIMGUI_VISIBLE) {
          SIMPLE(RAYIMGUI_TEXT, 0, "Click the image to focus controls.");
          SIMPLE(RAYIMGUI_TEXT, 0, "W/A/S/D move | Q/Z up/down");
          SIMPLE(RAYIMGUI_TEXT, 0, "Shift-W/S pitch | Shift-A/D roll");
          SIMPLE(RAYIMGUI_TEXT, 0, "Tab orbit | E/C speed | F fast preview");
          SIMPLE(RAYIMGUI_TEXT, 0, "Arrows lens | 1/2 focus | 3/4 environment");
          SIMPLE(RAYIMGUI_TEXT, 0, "Shift-click select / drill down");
          SIMPLE(RAYIMGUI_TEXT, 0, "Click target + focus | Right-click target");
          SIMPLE(RAYIMGUI_TEXT, 0, "R reset | [/] exposure | Enter render");
        }
        SIMPLE(RAYIMGUI_TREE_END, 0, "");
      }
      SIMPLE(RAYIMGUI_BUTTON, 8, "Cancel render");
      if (item.flags & RAYIMGUI_CHANGED) {
        status = table->request_close(owner, err);
        if (status) {
          return status;
        }
      }
    }
    SIMPLE(RAYIMGUI_WINDOW_END, 0, "");
    w = widget(RAYIMGUI_WINDOW_BEGIN, 30, "Scene hierarchy");
    w.options = RAYIMGUI_DOCK_RIGHT_TOP;
    DRAW();
    if (item.flags & RAYIMGUI_VISIBLE) {
      status = draw_hierarchy(table, owner, err);
      if (status) {
        return status;
      }
    }
    SIMPLE(RAYIMGUI_WINDOW_END, 0, "");
    w = widget(RAYIMGUI_WINDOW_BEGIN, 31, "Object inspector");
    w.options = RAYIMGUI_DOCK_RIGHT_BOTTOM;
    DRAW();
    if (item.flags & RAYIMGUI_VISIBLE) {
      if (can_edit && object.enabled) {
        status = draw_object_panel(table, owner, err);
        if (status) {
          return status;
        }
      } else {
        SIMPLE(RAYIMGUI_TEXT, 0, "Enable interactive rendering to edit objects.");
      }
    }
    SIMPLE(RAYIMGUI_WINDOW_END, 0, "");
    w = widget(RAYIMGUI_WINDOW_BEGIN, 9, "Viewport");
    w.width = 760;
    w.height = 700;
    w.x = 325;
    w.y = 12;
    w.options = RAYIMGUI_DOCK_CENTER | RAYIMGUI_WINDOW_VIEWPORT;
    DRAW();
    if (item.flags & RAYIMGUI_VISIBLE) {
      if (texture) {
        w = widget(RAYIMGUI_IMAGE, 10, "");
        w.texture = texture;
        // Zero dimensions fit the image into the current docked panel.
        w.height = 0;
        DRAW();
        if (item.flags & RAYIMGUI_VISIBLE) {
          viewport_visible = true;
          rayimgui_viewport_v1 viewport{sizeof(viewport)};
          rayimgui_input_v1 input{sizeof(input)};
          status = table->viewport(owner, &viewport, err);
          if (status) {
            return status;
          }
          if (object.selected && selection_visible && selection_texture &&
              selection_id == object.id && !object.clear_pending &&
              !object.select_pending && !object.pick_pending) {
            w = widget(RAYIMGUI_IMAGE_OVERLAY, 11, "");
            w.texture = selection_texture;
            w.fraction = 1;
            DRAW();
          }
          // Reserve Shift-click for picking before the gizmo can capture it.
          // An already-active drag keeps ownership until it commits or cancels.
          status = table->input(owner, &input, err);
          if (status) {
            return status;
          }
          const bool selecting = (input.modifiers & RAYIMGUI_SHIFT) &&
                                 !(input.modifiers & RAYIMGUI_ALT) &&
                                 !object.transform_active;
          if (object.selected && object.projection_valid && !object.clear_pending &&
              !object.select_pending && !selecting) {
            rayimgui_gizmo_v1 gizmo{sizeof(gizmo)};
            gizmo.id = object.id;
            gizmo.operation = object.operation + 1;
            gizmo.mode = object.mode + 1;
            gizmo.projection = object.projection_kind;
            gizmo.view = object.view.data();
            gizmo.projection_matrix = object.projection.data();
            gizmo.model = object.model.data();
            gizmo.viewport = &viewport;
            status = table->gizmo(owner, &gizmo, &item, err);
            if (status) {
              return status;
            }
            if (item.flags & RAYIMGUI_BEGIN) {
              ++history_epoch;
            }
            if (item.flags & (RAYIMGUI_CHANGED | RAYIMGUI_COMMIT)) {
              history_dirty = true;
            }
            object.transform_active = (item.flags & RAYIMGUI_ACTIVE) != 0;
            if (item.flags & RAYIMGUI_CHANGED) {
              object.transform_pending = true;
              object.numeric_transform = false;
            }
            if (item.flags & RAYIMGUI_COMMIT) {
              object.apply_transform = object.transform_pending;
            }
            if (item.flags & RAYIMGUI_CANCEL) {
              object.cancel_transform = true;
            }
          }
          status = table->input(owner, &input, err);
          if (status) {
            return status;
          }
          collect_input(input, viewport);
        }
      } else {
        SIMPLE(RAYIMGUI_TEXT, 0, "Waiting for the first completed sample");
      }
    }
    SIMPLE(RAYIMGUI_WINDOW_END, 0, "");
    // Hidden images cannot own movement or picking. Cancel an in-progress gizmo
    // drag so reopening the panel does not apply a stale transform.
    if (!viewport_visible) {
      keys.clear();
      pick_pending = false;
      object.pick_pending = false;
      if (object.transform_active) {
        object.cancel_transform = true;
      }
      object.transform_active = false;
    }
#undef SIMPLE
#undef DRAW
    return RAYIMGUI_OK;
  }
};
#endif
