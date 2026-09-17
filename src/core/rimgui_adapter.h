/* Copyright (c) 2026 Tyler Morgan-Wall. MIT; LICENSE.protocol.
 * Consumer-owned native state. No provider or GUI link dependency.
 */
#ifndef RAYRENDER_RIMGUI_ADAPTER_H
#define RAYRENDER_RIMGUI_ADAPTER_H
#include "rimgui/rimgui_r.h"
#include "preview_object_state.h"
#include <algorithm>
#include <chrono>
#include <cstring>
#include <cstdio>
#include <string>
#include <stdexcept>
#include <vector>

struct RayrenderGui {
  const rimgui_api_v1* api = nullptr;
  rimgui_session_handle session = 0;
  rimgui_callback_handle callback = 0;
  rimgui_texture_handle texture = 0;
  char message[512]{};
  rimgui_error_v1 error{sizeof(error), 0, message, sizeof(message)};
  std::vector<uint8_t> pixels;
  std::vector<float> display_rgb;
  uint32_t width = 0, height = 0;
  uint64_t version = 0;
  size_t samples = 0, objects = 0;
  double progress = 0, exposure = 1;
  int32_t haze = 0, altitude = 0;
  bool can_edit = false, has_atmosphere = false, deferred = false;
  bool close_requested = false, request_render = false, request_reset = false;
  bool atmosphere_pending = false, exposure_pending = false;
  bool has_sun = false, sun_pending = false, sun_editing = false;
  double sun_elevation = 0, sun_azimuth = 0;
  bool has_location = false, location_pending = false, manual_sun = false;
  double latitude = 0, longitude = 0;
  char datetime[32]{};
  int32_t date[3]{2026, 1, 1}, time[3]{0, 0, 0};
  int32_t fast_preview = 0;
  bool fast_pending = false;
  std::string sky_error;
  int32_t orbit = 1;
  double movement_speed = 1;
  struct Key {
    unsigned code, modifiers, count;
  };
  std::vector<Key> keys;
  bool pick_pending = false, pick_focus = false, escape_requested = false;
  PreviewObjectState object;
  float pick_u = 0, pick_v = 0;
  std::chrono::steady_clock::time_point last_step{};

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
    int32_t status = rimgui_api_from_R_v1(handle,
                                          sizeof(rimgui_api_v1),
                                          RIMGUI_CAP_RGBA8 | RIMGUI_CAP_WIDGETS |
                                              RIMGUI_CAP_INPUT | RIMGUI_CAP_GIZMO,
                                          &api);
    if (status || !api->input || api->header.abi_minor < 2) {
      std::snprintf(message,
                    sizeof(message),
                    "Object controls require rimgui 0.0.6 or later (ABI 1.2).");
      return status ? status : RIMGUI_ABI;
    }
    if (!(api->header.capabilities & RIMGUI_CAP_NATIVE)) {
      std::snprintf(message,
                    sizeof(message),
                    "rimgui was built without a native window backend.");
      return RIMGUI_BACKEND_NOT_BUILT;
    }
    return RIMGUI_OK;
  }

  static bool expected_unavailability(int32_t status) {
    return status == RIMGUI_ABI || status == RIMGUI_BACKEND_NOT_BUILT ||
           status == RIMGUI_BACKEND_INIT;
  }

  int32_t open(uint32_t w, uint32_t h) {
    width = w;
    height = h;
    rimgui_session_desc_v1 desc{sizeof(desc), 0, "rayrender", 9, 1100, 760};
    int32_t status = api->open(&desc, &session, &error);
    if (status) {
      return status;
    }
    rimgui_callback_v1 client{sizeof(client), 1, draw, this};
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
    rimgui_image_v1 image{sizeof(image),
                          width,
                          height,
                          RIMGUI_TOP_LEFT,
                          RIMGUI_RGBA,
                          RIMGUI_OPAQUE,
                          RIMGUI_DISPLAY_ENCODED,
                          static_cast<uint64_t>(width) * 4,
                          static_cast<uint64_t>(pixels.size()),
                          ++version,
                          pixels.data()};
    check(texture ? api->texture_update(session, texture, &image, &error)
                  : api->texture_create(session, &image, &texture, &error));
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
    rimgui_step_v1 result{};
    result.struct_size = sizeof(result);
    check(api->step(session, &result, &error));
    close_requested = result.close_requested != 0 || escape_requested;
    return close_requested;
  }

  void close() noexcept {
    // Also called by explicit R_UnwindProtect cleanup, never a finalizer.
    if (session && api) {
      api->close(&session, nullptr);
    }
    callback = 0;
    texture = 0;
    api = nullptr;
  }

  // Respect provider focus/ownership decisions and queue intent only. These
  // callbacks can run while sampling workers still read the active scene.
  void collect_input(const rimgui_input_v1& input, const rimgui_viewport_v1& viewport) {
    if (input.keyboard_available &&
        !(input.modifiers & (RIMGUI_CTRL | RIMGUI_ALT | RIMGUI_SUPER))) {
      escape_requested =
          escape_requested || (input.keys_pressed & RIMGUI_KEY_BIT(RIMGUI_KEY_ESCAPE));
      // Limit queued motion while a slow sample renders. No renderer or R
      // state is touched from a GUI callback or worker-wait poll.
      const uint64_t repeat =
          RIMGUI_KEY_BIT(RIMGUI_KEY_W) | RIMGUI_KEY_BIT(RIMGUI_KEY_A) |
          RIMGUI_KEY_BIT(RIMGUI_KEY_S) | RIMGUI_KEY_BIT(RIMGUI_KEY_D) |
          RIMGUI_KEY_BIT(RIMGUI_KEY_Q) | RIMGUI_KEY_BIT(RIMGUI_KEY_Z) |
          RIMGUI_KEY_BIT(RIMGUI_KEY_UP) | RIMGUI_KEY_BIT(RIMGUI_KEY_DOWN) |
          RIMGUI_KEY_BIT(RIMGUI_KEY_LEFT) | RIMGUI_KEY_BIT(RIMGUI_KEY_RIGHT) |
          RIMGUI_KEY_BIT(RIMGUI_KEY_1) | RIMGUI_KEY_BIT(RIMGUI_KEY_2) |
          RIMGUI_KEY_BIT(RIMGUI_KEY_3) | RIMGUI_KEY_BIT(RIMGUI_KEY_4);
      uint64_t pressed = input.keys_pressed | (input.keys_repeated & repeat);
      if (can_edit) {
        for (unsigned code = 0; code < RIMGUI_KEY_COUNT; ++code) {
          if (pressed & RIMGUI_KEY_BIT(code)) {
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
        (input.mouse_clicked & (RIMGUI_MOUSE_LEFT | RIMGUI_MOUSE_RIGHT))) {
      // Undo the displayed image's two-axis flip to obtain renderer film coordinates.
      const float u = 1 - viewport.mouse_x / viewport.width,
                  v = 1 - viewport.mouse_y / viewport.height;
      if (object.enabled && (input.modifiers & RIMGUI_SHIFT) &&
          (input.mouse_clicked & RIMGUI_MOUSE_LEFT)) {
        object.pick_pending = true;
        object.pick_u = u;
        object.pick_v = v;
        pick_pending = false;
      } else {
        pick_pending = true;
        pick_focus = (input.mouse_clicked & RIMGUI_MOUSE_LEFT) != 0;
        pick_u = u;
        pick_v = v;
      }
    }
  }

  static rimgui_widget_v1 widget(uint32_t kind, uint64_t id, const char* label) {
    rimgui_widget_v1 value{};
    value.struct_size = sizeof(value);
    value.kind = kind;
    value.id = id;
    value.label = label;
    value.label_length = static_cast<uint32_t>(std::strlen(label));
    value.speed = .01;
    return value;
  }

  static int32_t draw(const rimgui_api_v1* api, rimgui_session_handle session,
                      void* userdata, rimgui_error_v1* error) noexcept {
    try {
      return static_cast<RayrenderGui*>(userdata)->draw_impl(api, session, error);
    } catch (...) {
      return RIMGUI_CALLBACK_ERROR;
    }
  }

  // Widgets edit value snapshots and request flags. Scene rebuilding and material
  // setters run later at the renderer checkpoint, outside this provider callback.
  int32_t draw_object_panel(const rimgui_api_v1* table, rimgui_session_handle owner,
                            rimgui_error_v1* err) {
    auto& o = object;
    rimgui_item_v1 item{sizeof(item)};
    rimgui_widget_v1 w{};
    int32_t status;
#define OBJECT_DRAW()                                                                  \
  do {                                                                                 \
    status = table->widget(owner, &w, &item, err);                                     \
    if (status)                                                                        \
      return status;                                                                   \
  } while (0)
#define OBJECT_SIMPLE(kind, id, label)                                                 \
  do {                                                                                 \
    w = widget(kind, id, label);                                                       \
    OBJECT_DRAW();                                                                     \
  } while (0)
    OBJECT_SIMPLE(RIMGUI_SEPARATOR, 0, "");
    OBJECT_SIMPLE(RIMGUI_TEXT, 0, "Shift-click an object to select it.");
    if (o.selected) {
      OBJECT_SIMPLE(RIMGUI_TEXT, 0, o.label.c_str());
      OBJECT_SIMPLE(RIMGUI_BUTTON, 7000, "Deselect");
      o.clear_pending = o.clear_pending || (item.flags & RIMGUI_CHANGED);
      static const char* operations[] = {"Translate", "Rotate", "Scale"};
      static const uint32_t lengths[] = {9, 6, 5};
      w = widget(RIMGUI_COMBO, 7001, "Tool");
      w.integers = &o.operation;
      w.count = 1;
      w.choices = operations;
      w.choice_lengths = lengths;
      w.choice_count = 3;
      OBJECT_DRAW();
      static const char* modes[] = {"Local", "World"};
      static const uint32_t mode_lengths[] = {5, 5};
      w = widget(RIMGUI_COMBO, 7002, "Axes");
      w.integers = &o.mode;
      w.count = 1;
      w.choices = modes;
      w.choice_lengths = mode_lengths;
      w.choice_count = 2;
      OBJECT_DRAW();
      const char* labels[] = {"Position", "Rotation (degrees)", "Scale"};
      double* values[] = {o.translation.data(), o.rotation.data(), o.scale.data()};
      for (int i = 0; i < 3; ++i) {
        w = widget(RIMGUI_DOUBLE, 7003 + i, labels[i]);
        w.numbers = values[i];
        w.count = 3;
        w.speed = i == 1 ? .25 : .01;
        if (i == 2) {
          w.minimum = -10000;
          w.maximum = 10000;
        }
        OBJECT_DRAW();
        if (item.flags & RIMGUI_CHANGED) {
          o.transform_pending = true;
          o.numeric_transform = true;
        }
      }
      OBJECT_SIMPLE(RIMGUI_BUTTON, 7006, "Apply transform");
      o.apply_transform = o.apply_transform || (item.flags & RIMGUI_CHANGED);
      if (o.projection_valid) {
        OBJECT_SIMPLE(RIMGUI_TEXT, 0, "Drag a handle; release to apply.");
      } else {
        OBJECT_SIMPLE(RIMGUI_TEXT, 0, "Use numeric transforms for this camera.");
      }
      OBJECT_SIMPLE(RIMGUI_SEPARATOR, 0, "");
      OBJECT_SIMPLE(RIMGUI_TEXT, 0, "Material");
      if (!o.materials.empty()) {
        if (o.materials.size() > 1 && o.materials.size() <= 1024) {
          std::vector<const char*> labels;
          std::vector<uint32_t> sizes;
          for (const auto& panel : o.materials) {
            labels.push_back(panel.label.c_str());
            sizes.push_back(panel.label.size());
          }
          w = widget(RIMGUI_COMBO, 7010, "Slot");
          w.integers = &o.material_slot;
          w.count = 1;
          w.choices = labels.data();
          w.choice_lengths = sizes.data();
          w.choice_count = labels.size();
          OBJECT_DRAW();
        } else if (o.materials.size() > 1024) {
          w = widget(RIMGUI_INT, 7010, "Slot (zero based)");
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
        OBJECT_SIMPLE(RIMGUI_TEXT, 0, panel.type.c_str());
        for (size_t i = 0; i < panel.fields.size(); ++i) {
          auto& field = panel.fields[i];
          w = widget(field.color ? RIMGUI_COLOR : RIMGUI_DOUBLE,
                     100000 + 1000 * uint64_t(o.material_slot) + i,
                     field.name.c_str());
          w.numbers = field.values.data();
          w.count = field.count;
          w.minimum = field.minimum;
          w.maximum = field.maximum;
          w.speed = field.speed;
          if (field.color) {
            w.options = RIMGUI_COLOR_RGB;
            w.fraction = 1;
          }
          OBJECT_DRAW();
          if (item.flags & RIMGUI_CHANGED) {
            panel.changed = true;
            o.material_pending = true;
          }
        }
        if (panel.fields.empty()) {
          OBJECT_SIMPLE(RIMGUI_TEXT, 0, "No editable parameters for this material.");
        }
        OBJECT_SIMPLE(RIMGUI_BUTTON, 7011, "Apply material");
        o.apply_material =
            o.apply_material || ((item.flags & RIMGUI_CHANGED) && o.material_pending);
      }
      OBJECT_SIMPLE(RIMGUI_BUTTON, 7012, "Reset object edits");
      o.revert = o.revert || (item.flags & RIMGUI_CHANGED);
      if (!o.error.empty()) {
        OBJECT_SIMPLE(RIMGUI_TEXT, 0, o.error.c_str());
      }
    }
    OBJECT_SIMPLE(RIMGUI_SEPARATOR, 0, "");
    OBJECT_SIMPLE(RIMGUI_TEXT, 0, "Render and camera");
#undef OBJECT_SIMPLE
#undef OBJECT_DRAW
    return RIMGUI_OK;
  }

  int32_t draw_impl(const rimgui_api_v1* table, rimgui_session_handle owner,
                    rimgui_error_v1* err) {
    rimgui_widget_v1 w{};
    rimgui_item_v1 item{sizeof(item), 0};
    int32_t status;
    bool viewport_visible = false;
    sun_editing = false;
#define DRAW()                                                                         \
  do {                                                                                 \
    status = table->widget(owner, &w, &item, err);                                     \
    if (status)                                                                        \
      return status;                                                                   \
  } while (0)
#define SIMPLE(kind, id, label)                                                        \
  do {                                                                                 \
    w = widget(kind, id, label);                                                       \
    DRAW();                                                                            \
  } while (0)
    w = widget(RIMGUI_WINDOW_BEGIN, 1, "Render controls");
    w.width = 300;
    w.height = 700;
    w.x = 12;
    w.y = 12;
    w.options = RIMGUI_WINDOW_POSITION;
    DRAW();
    if (item.flags & RIMGUI_VISIBLE) {
      SIMPLE(RIMGUI_TEXT, 0, "CPU progressive render");
      char detail[128];
      std::snprintf(
          detail, sizeof(detail), "%zu scene objects | %zu samples", objects, samples);
      SIMPLE(RIMGUI_TEXT, 0, detail);
      w = widget(RIMGUI_PROGRESS, 2, "Progress");
      w.fraction = progress;
      w.width = 260;
      w.height = 18;
      DRAW();
      if (can_edit && object.enabled) {
        status = draw_object_panel(table, owner, err);
        if (status) {
          return status;
        }
      }
      w = widget(RIMGUI_DOUBLE, 3, "Exposure");
      w.numbers = &exposure;
      w.count = 1;
      w.minimum = .001;
      w.maximum = 1024;
      w.speed = .01;
      DRAW();
      exposure_pending = exposure_pending || (item.flags & RIMGUI_CHANGED);
      if (can_edit) {
        w = widget(RIMGUI_CHECKBOX, 26, "Fast preview (F)");
        w.integers = &fast_preview;
        w.count = 1;
        DRAW();
        fast_pending = fast_pending || (item.flags & RIMGUI_CHANGED);
        w = widget(RIMGUI_CHECKBOX, 11, "Orbit camera");
        w.integers = &orbit;
        w.count = 1;
        DRAW();
        w = widget(RIMGUI_DOUBLE, 12, "Movement speed");
        w.numbers = &movement_speed;
        w.count = 1;
        w.minimum = .001;
        w.maximum = 128;
        w.speed = .05;
        DRAW();
        SIMPLE(RIMGUI_BUTTON, 4, "Reset camera");
        request_reset = request_reset || (item.flags & RIMGUI_CHANGED);
        if (deferred) {
          SIMPLE(RIMGUI_BUTTON, 5, "Start final render");
          request_render = request_render || (item.flags & RIMGUI_CHANGED);
        }
        if (has_atmosphere) {
          SIMPLE(RIMGUI_SEPARATOR, 0, "");
          SIMPLE(RIMGUI_TEXT, 0, "Atmosphere");
          if (has_sun) {
            w = widget(RIMGUI_DOUBLE, 13, "Sun elevation (degrees)");
            w.numbers = &sun_elevation;
            w.count = 1;
            w.minimum = -90;
            w.maximum = 90;
            w.speed = .25;
            DRAW();
            sun_pending = sun_pending || (item.flags & RIMGUI_CHANGED);
            sun_editing = sun_editing || (item.flags & RIMGUI_ACTIVE);
            w = widget(RIMGUI_DOUBLE, 14, "Sun azimuth (degrees)");
            w.numbers = &sun_azimuth;
            w.count = 1;
            w.minimum = 0;
            w.maximum = 360;
            w.speed = .5;
            DRAW();
            sun_pending = sun_pending || (item.flags & RIMGUI_CHANGED);
            sun_editing = sun_editing || (item.flags & RIMGUI_ACTIVE);
            SIMPLE(RIMGUI_TOOLTIP,
                   0,
                   "North 0, east 90. Release to update the sky and Sun.");
            if (manual_sun) {
              SIMPLE(RIMGUI_TEXT, 0, "Sun direction manually adjusted");
            }
          }
          if (has_location) {
            w = widget(RIMGUI_DOUBLE, 16, "Latitude");
            w.numbers = &latitude;
            w.count = 1;
            w.minimum = -90;
            w.maximum = 90;
            w.speed = .1;
            DRAW();
            w = widget(RIMGUI_DOUBLE, 17, "Longitude");
            w.numbers = &longitude;
            w.count = 1;
            w.minimum = -180;
            w.maximum = 180;
            w.speed = .1;
            DRAW();
            SIMPLE(RIMGUI_TEXT, 0, "Date/time (UTC)");
            const char* labels[] = {"Year", "Month", "Day", "Hour", "Minute", "Second"};
            const int maximum[] = {9999, 12, 31, 23, 59, 59};
            for (unsigned component = 0; component < 6; ++component) {
              w = widget(RIMGUI_INT, 20 + component, labels[component]);
              w.integers = component < 3 ? &date[component] : &time[component - 3];
              w.count = 1;
              w.minimum = component < 3 ? 1 : 0;
              w.maximum = maximum[component];
              w.speed = .1;
              DRAW();
              if (item.flags & RIMGUI_CHANGED) {
                format_datetime();
              }
              SIMPLE(RIMGUI_TOOLTIP,
                     0,
                     "Drag to change; double-click to type. Then Apply location/time.");
            }
            SIMPLE(RIMGUI_TEXT, 0, datetime);
            SIMPLE(RIMGUI_BUTTON, 19, "Apply location/time");
            location_pending = location_pending || (item.flags & RIMGUI_CHANGED);
            if (!sky_error.empty()) {
              SIMPLE(RIMGUI_TEXT, 0, sky_error.c_str());
            }
          }
          w = widget(RIMGUI_CHECKBOX, 6, "Haze");
          w.integers = &haze;
          w.count = 1;
          DRAW();
          if (item.flags & RIMGUI_CHANGED) {
            if (haze) {
              altitude = 1;
            }
            atmosphere_pending = true;
          }
          w = widget(RIMGUI_CHECKBOX, 7, "Altitude queries");
          w.integers = &altitude;
          w.count = 1;
          DRAW();
          if (item.flags & RIMGUI_CHANGED) {
            if (!altitude) {
              haze = 0;
            }
            atmosphere_pending = true;
          }
        }
        SIMPLE(RIMGUI_TREE_BEGIN, 15, "Movement controls");
        if (item.flags & RIMGUI_VISIBLE) {
          SIMPLE(RIMGUI_TEXT, 0, "Click the image to focus controls.");
          SIMPLE(RIMGUI_TEXT, 0, "W/A/S/D move | Q/Z up/down");
          SIMPLE(RIMGUI_TEXT, 0, "Shift-W/S pitch | Shift-A/D roll");
          SIMPLE(RIMGUI_TEXT, 0, "Tab orbit | E/C speed | F fast preview");
          SIMPLE(RIMGUI_TEXT, 0, "Arrows lens | 1/2 focus | 3/4 environment");
          SIMPLE(RIMGUI_TEXT, 0, "Left click target + focus | Right click target");
          SIMPLE(RIMGUI_TEXT, 0, "R reset | [/] exposure | Enter render");
        }
        SIMPLE(RIMGUI_TREE_END, 0, "");
      }
      SIMPLE(RIMGUI_BUTTON, 8, "Cancel render");
      if (item.flags & RIMGUI_CHANGED) {
        status = table->request_close(owner, err);
        if (status) {
          return status;
        }
      }
    }
    SIMPLE(RIMGUI_WINDOW_END, 0, "");
    w = widget(RIMGUI_WINDOW_BEGIN, 9, "Rendered image");
    w.width = 760;
    w.height = 700;
    w.x = 325;
    w.y = 12;
    w.options = RIMGUI_WINDOW_POSITION | RIMGUI_WINDOW_VIEWPORT;
    DRAW();
    if (item.flags & RIMGUI_VISIBLE) {
      if (texture) {
        w = widget(RIMGUI_IMAGE, 10, "");
        w.texture = texture;
        w.height = 620;
        DRAW();
        if (item.flags & RIMGUI_VISIBLE) {
          viewport_visible = true;
          rimgui_viewport_v1 viewport{sizeof(viewport)};
          rimgui_input_v1 input{sizeof(input)};
          status = table->viewport(owner, &viewport, err);
          if (status) {
            return status;
          }
          if (object.selected && object.projection_valid && !object.clear_pending) {
            rimgui_gizmo_v1 gizmo{sizeof(gizmo)};
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
            object.transform_active = (item.flags & RIMGUI_ACTIVE) != 0;
            if (item.flags & RIMGUI_CHANGED) {
              object.transform_pending = true;
              object.numeric_transform = false;
            }
            if (item.flags & RIMGUI_COMMIT) {
              object.apply_transform = object.transform_pending;
            }
            if (item.flags & RIMGUI_CANCEL) {
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
        SIMPLE(RIMGUI_TEXT, 0, "Waiting for the first completed sample");
      }
    }
    SIMPLE(RIMGUI_WINDOW_END, 0, "");
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
    return RIMGUI_OK;
  }
};
#endif
