/* Copyright (c) 2026 Tyler Morgan-Wall. MIT license; see LICENSE.protocol.
 * Pure C protocol. This development ABI is not yet frozen for production.
 * Complete phase/ownership rules: installed developer/ABI.md.
 */
#ifndef RAYIMGUI_API_H
#define RAYIMGUI_API_H
#include <stddef.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

#define RAYIMGUI_ABI_MAJOR 1u
#define RAYIMGUI_ABI_MINOR 8u
#define RAYIMGUI_CAP_HEADLESS UINT64_C(1)
#define RAYIMGUI_CAP_NATIVE UINT64_C(2)
#define RAYIMGUI_CAP_DOCKING UINT64_C(4)
#define RAYIMGUI_CAP_GIZMO UINT64_C(8)
#define RAYIMGUI_CAP_RGBA8 UINT64_C(16)
#define RAYIMGUI_CAP_WIDGETS UINT64_C(32)
#define RAYIMGUI_CAP_EVENTS UINT64_C(64)
#define RAYIMGUI_CAP_INPUT UINT64_C(128)

typedef uint64_t rayimgui_session_handle;
typedef uint64_t rayimgui_texture_handle;
typedef uint64_t rayimgui_callback_handle;

enum rayimgui_status {
  RAYIMGUI_OK = 0,
  RAYIMGUI_INVALID = 1,
  RAYIMGUI_ABI = 2,
  RAYIMGUI_BACKEND_NOT_BUILT = 3,
  RAYIMGUI_BACKEND_INIT = 4,
  RAYIMGUI_WRONG_THREAD = 5,
  RAYIMGUI_PHASE = 6,
  RAYIMGUI_STALE = 7,
  RAYIMGUI_BUSY = 8,
  RAYIMGUI_CALLBACK_ERROR = 9,
  RAYIMGUI_MEMORY = 10,
  RAYIMGUI_INTERNAL = 11,
  RAYIMGUI_EMPTY = 12,
  RAYIMGUI_LIMIT = 13
};
enum rayimgui_session_flags {
  RAYIMGUI_HEADLESS = 1, RAYIMGUI_HIDDEN = 2,
  /* ABI 1.7: emit application undo/redo events outside active text/widgets. */
  RAYIMGUI_HISTORY_SHORTCUTS = 4,
  /* ABI 1.8: Escape requests a close from any panel, including active inputs. */
  RAYIMGUI_ESCAPE_CLOSE = 8
};
/* Reserved provider commands; macros keep these uint32 values valid in C11. */
#define RAYIMGUI_COMMAND_UNDO UINT32_C(0xffff0001)
#define RAYIMGUI_COMMAND_REDO UINT32_C(0xffff0002)
enum rayimgui_item_flags {
  RAYIMGUI_CHANGED = 1,
  RAYIMGUI_ACTIVE = 2,
  RAYIMGUI_HOVERED = 4,
  RAYIMGUI_BEGIN = 8,
  RAYIMGUI_UPDATE = 16,
  RAYIMGUI_COMMIT = 32,
  RAYIMGUI_CANCEL = 64,
  RAYIMGUI_VISIBLE = 128,
  RAYIMGUI_DEACTIVATED = 256
};
enum rayimgui_widget_kind {
  RAYIMGUI_WINDOW_BEGIN = 1,
  RAYIMGUI_WINDOW_END = 2,
  RAYIMGUI_CHILD_BEGIN = 3,
  RAYIMGUI_CHILD_END = 4,
  RAYIMGUI_ID_PUSH = 5,
  RAYIMGUI_ID_POP = 6,
  RAYIMGUI_TEXT = 7,
  RAYIMGUI_BUTTON = 8,
  RAYIMGUI_CHECKBOX = 9,
  RAYIMGUI_INT = 10,
  RAYIMGUI_DOUBLE = 11,
  RAYIMGUI_COLOR = 12,
  RAYIMGUI_INPUT_TEXT = 13,
  RAYIMGUI_COMBO = 14,
  RAYIMGUI_TREE_BEGIN = 15,
  RAYIMGUI_TREE_END = 16,
  RAYIMGUI_MENU_BAR_BEGIN = 17,
  RAYIMGUI_MENU_BAR_END = 18,
  RAYIMGUI_MENU_BEGIN = 19,
  RAYIMGUI_MENU_END = 20,
  RAYIMGUI_MENU_ITEM = 21,
  RAYIMGUI_SEPARATOR = 22,
  RAYIMGUI_SAME_LINE = 23,
  RAYIMGUI_DISABLED_BEGIN = 24,
  RAYIMGUI_DISABLED_END = 25,
  RAYIMGUI_PROGRESS = 26,
  RAYIMGUI_IMAGE = 27,
  RAYIMGUI_TOOLTIP = 28,
  /* Draw an RGBA texture over the current image viewport without adding an
   * input item or changing layout. fraction is opacity in [0,1]. ABI >= 1.4. */
  RAYIMGUI_IMAGE_OVERLAY = 29,
  /* ABI 1.5: one angle in degrees, draggable arc and embedded numeric input. */
  RAYIMGUI_ANGLE = 30,
  /* ABI 1.6: clickable aspect-fitted texture with a caption; does not register
   * a viewport. Positive width/height describe the image box, excluding caption. */
  RAYIMGUI_IMAGE_BUTTON = 31
};
enum rayimgui_widget_options {
  RAYIMGUI_WINDOW_MENU_BAR = 1,
  RAYIMGUI_WINDOW_POSITION = 2,
  RAYIMGUI_WINDOW_VIEWPORT = 4, /* Route navigation keys to the image consumer. */
  RAYIMGUI_COLOR_RGB =
      8, /* ABI 1.2: COLOR edits bounded RGB without an intensity field. */
  /* ABI 1.3: initial docking slots in a resizable three-column workspace.
   * Use one slot per window. Users may subsequently move or resize the panels. */
  RAYIMGUI_DOCK_LEFT = 16,
  RAYIMGUI_DOCK_CENTER = 32,
  RAYIMGUI_DOCK_RIGHT_TOP = 64,
  RAYIMGUI_DOCK_RIGHT_BOTTOM = 128,
  /* TREE_BEGIN returns CHANGED for selection clicks, independently of VISIBLE. */
  RAYIMGUI_TREE_SELECTED = 256,
  RAYIMGUI_TREE_LEAF = 512,
  RAYIMGUI_TREE_DEFAULT_OPEN = 1024,
  /* ANGLE uses a 0..90 quarter arc instead of a clockwise 0..360 circle. */
  RAYIMGUI_ANGLE_QUARTER = 2048,
  /* ABI 1.6: optional full-width bottom workspace panel. */
  RAYIMGUI_DOCK_BOTTOM = 4096,
  RAYIMGUI_CHILD_HORIZONTAL_SCROLL = 8192,
  RAYIMGUI_IMAGE_SELECTED = 16384
};
enum rayimgui_key {
  RAYIMGUI_KEY_A,
  RAYIMGUI_KEY_B,
  RAYIMGUI_KEY_C,
  RAYIMGUI_KEY_D,
  RAYIMGUI_KEY_E,
  RAYIMGUI_KEY_F,
  RAYIMGUI_KEY_G,
  RAYIMGUI_KEY_H,
  RAYIMGUI_KEY_I,
  RAYIMGUI_KEY_J,
  RAYIMGUI_KEY_K,
  RAYIMGUI_KEY_L,
  RAYIMGUI_KEY_M,
  RAYIMGUI_KEY_N,
  RAYIMGUI_KEY_O,
  RAYIMGUI_KEY_P,
  RAYIMGUI_KEY_Q,
  RAYIMGUI_KEY_R,
  RAYIMGUI_KEY_S,
  RAYIMGUI_KEY_T,
  RAYIMGUI_KEY_U,
  RAYIMGUI_KEY_V,
  RAYIMGUI_KEY_W,
  RAYIMGUI_KEY_X,
  RAYIMGUI_KEY_Y,
  RAYIMGUI_KEY_Z,
  RAYIMGUI_KEY_0,
  RAYIMGUI_KEY_1,
  RAYIMGUI_KEY_2,
  RAYIMGUI_KEY_3,
  RAYIMGUI_KEY_4,
  RAYIMGUI_KEY_5,
  RAYIMGUI_KEY_6,
  RAYIMGUI_KEY_7,
  RAYIMGUI_KEY_8,
  RAYIMGUI_KEY_9,
  RAYIMGUI_KEY_TAB,
  RAYIMGUI_KEY_LEFT,
  RAYIMGUI_KEY_RIGHT,
  RAYIMGUI_KEY_UP,
  RAYIMGUI_KEY_DOWN,
  RAYIMGUI_KEY_ENTER,
  RAYIMGUI_KEY_ESCAPE,
  RAYIMGUI_KEY_LEFT_BRACKET,
  RAYIMGUI_KEY_RIGHT_BRACKET,
  RAYIMGUI_KEY_COMMA,
  RAYIMGUI_KEY_PERIOD,
  RAYIMGUI_KEY_SLASH,
  RAYIMGUI_KEY_SPACE,
  RAYIMGUI_KEY_KEYPAD_ENTER,
  RAYIMGUI_KEY_COUNT
};
#define RAYIMGUI_KEY_BIT(key) (UINT64_C(1) << (key))
enum rayimgui_modifiers {
  RAYIMGUI_SHIFT = 1,
  RAYIMGUI_CTRL = 2,
  RAYIMGUI_ALT = 4,
  RAYIMGUI_SUPER = 8
};
enum rayimgui_mouse_button {
  RAYIMGUI_MOUSE_LEFT = 1,
  RAYIMGUI_MOUSE_RIGHT = 2,
  RAYIMGUI_MOUSE_MIDDLE = 4
};
enum rayimgui_operation {
  RAYIMGUI_TRANSLATE = 1,
  RAYIMGUI_ROTATE = 2,
  RAYIMGUI_SCALE = 3
};
enum rayimgui_mode { RAYIMGUI_LOCAL = 1, RAYIMGUI_WORLD = 2 };
enum rayimgui_projection { RAYIMGUI_PERSPECTIVE = 1, RAYIMGUI_ORTHOGRAPHIC = 2 };
enum rayimgui_origin { RAYIMGUI_TOP_LEFT = 1, RAYIMGUI_BOTTOM_LEFT = 2 };
enum rayimgui_channels { RAYIMGUI_RGBA = 1, RAYIMGUI_BGRA = 2 };
enum rayimgui_alpha { RAYIMGUI_OPAQUE = 1, RAYIMGUI_STRAIGHT = 2 };
enum rayimgui_color_space { RAYIMGUI_DISPLAY_ENCODED = 1 };

typedef struct rayimgui_error_v1 {
  uint32_t struct_size;
  int32_t code;
  char* message;
  uint32_t capacity; /* Includes trailing NUL; truncation always terminates. */
} rayimgui_error_v1;

typedef struct rayimgui_session_desc_v1 {
  uint32_t struct_size;
  uint32_t flags;
  const char* title;
  uint32_t title_length;
  uint32_t width;
  uint32_t height;
} rayimgui_session_desc_v1;

typedef struct rayimgui_step_v1 {
  uint32_t struct_size;
  uint32_t close_requested;
  uint64_t frame;
  double delta_seconds;
  float width, height;
  float framebuffer_scale_x, framebuffer_scale_y;
  uint32_t focused;
  uint32_t want_keyboard;
} rayimgui_step_v1;

typedef struct rayimgui_image_v1 {
  uint32_t struct_size;
  uint32_t width, height;
  uint32_t origin, channels, alpha, color_space;
  uint64_t row_stride;
  uint64_t buffer_length;
  uint64_t version; /* Unchanged version means unchanged pixels. */
  const uint8_t* pixels;
} rayimgui_image_v1;

typedef struct rayimgui_viewport_v1 {
  uint32_t struct_size;
  float x, y, width, height; /* Image bounds in window logical coordinates. */
  float letterbox_x, letterbox_y;
  float mouse_x, mouse_y; /* Relative to the actual image, not its panel. */
  float framebuffer_scale_x, framebuffer_scale_y;
  uint32_t inside, focused, mouse_down, mouse_clicked, mouse_released;
  uint32_t widget_owns_mouse, gizmo_hovered, gizmo_active;
  uint32_t background_click; /* Valid after calling gizmo for this viewport. */
} rayimgui_viewport_v1;

/* ABI 1.1: read during a draw callback, after the current window's image/gizmo.
 * Unavailable input has zero key/button masks. Positions come from viewport().
 * pressed is an edge; repeated includes the first press and keyboard repeat.
 */
typedef struct rayimgui_input_v1 {
  uint32_t struct_size;
  uint32_t keyboard_available, mouse_available, modifiers;
  uint64_t keys_down, keys_pressed, keys_repeated;
  uint32_t mouse_down, mouse_clicked, mouse_released;
} rayimgui_input_v1;

typedef struct rayimgui_widget_v1 {
  uint32_t struct_size;
  uint32_t kind;
  uint64_t id; /* Stable application ID; zero allowed for noninteractive layout. */
  const char* label;
  uint32_t label_length;
  uint32_t options;
  int32_t* integers; /* CHECKBOX/INT/COMBO; count 1..4, COMBO uses one. */
  double* numbers;   /* DOUBLE: 1..4; COLOR: RGB in [0,1], separate intensity. */
  uint32_t count;
  double minimum, maximum, speed;
  char* text; /* INPUT_TEXT: caller-owned writable UTF-8 with terminating NUL. */
  uint32_t text_capacity;
  const char* const* choices;
  const uint32_t* choice_lengths;
  uint32_t choice_count;
  float width, height, x, y; /* WINDOW_POSITION uses x/y on first appearance. */
  double fraction;           /* PROGRESS; or separate COLOR intensity (input/output). */
  rayimgui_texture_handle texture;
} rayimgui_widget_v1;

typedef struct rayimgui_item_v1 {
  uint32_t struct_size;
  uint32_t flags;
} rayimgui_item_v1;

typedef struct rayimgui_gizmo_v1 {
  uint32_t struct_size;
  uint32_t operation, mode, projection;
  uint64_t id;
  const float* view;              /* Borrowed 16 floats, column major, right handed. */
  const float* projection_matrix; /* OpenGL clip Z [-1,1]. */
  float* model; /* 16 floats, column vectors: clip = P * V * M * point. */
  const rayimgui_viewport_v1* viewport;
  uint32_t snap_enabled;
  float snap[3]; /* Translation world units, rotation degrees, scale units. */
} rayimgui_gizmo_v1;

typedef struct rayimgui_event_v1 {
  uint32_t struct_size;
  uint32_t phase; /* Exactly one of BEGIN/UPDATE/COMMIT/CANCEL. */
  uint64_t application_id;
  uint32_t command;
  uint32_t count;
  double values[4]; /* Copied by event_push; interpreted only by consumer. */
} rayimgui_event_v1;

struct rayimgui_api_v1;
typedef int32_t (*rayimgui_draw_ui_v1)(const struct rayimgui_api_v1*,
                                       rayimgui_session_handle, void*,
                                       rayimgui_error_v1*);
typedef struct rayimgui_callback_v1 {
  uint32_t struct_size;
  uint32_t abi_major;
  rayimgui_draw_ui_v1 draw_ui;
  void* userdata;
} rayimgui_callback_v1;

/* Fixed negotiation prefix. Read ONLY this prefix before validating table size.
 * All v1 functions are mandatory. Minor versions may append, never reorder.
 */
typedef struct rayimgui_api_header_v1 {
  uint32_t struct_size;
  uint32_t abi_major;
  uint32_t abi_minor;
  uint32_t reserved;
  uint64_t capabilities;
} rayimgui_api_header_v1;

typedef struct rayimgui_api_v1 {
  rayimgui_api_header_v1 header;
  int32_t (*open)(const rayimgui_session_desc_v1*, rayimgui_session_handle*,
                  rayimgui_error_v1*);
  int32_t (*close)(rayimgui_session_handle*, rayimgui_error_v1*);
  int32_t (*register_callback)(rayimgui_session_handle, const rayimgui_callback_v1*,
                               rayimgui_callback_handle*, rayimgui_error_v1*);
  int32_t (*unregister_callback)(rayimgui_session_handle, rayimgui_callback_handle*,
                                 rayimgui_error_v1*);
  int32_t (*step)(rayimgui_session_handle, rayimgui_step_v1*, rayimgui_error_v1*);
  int32_t (*request_close)(rayimgui_session_handle, rayimgui_error_v1*);
  int32_t (*texture_create)(rayimgui_session_handle, const rayimgui_image_v1*,
                            rayimgui_texture_handle*, rayimgui_error_v1*);
  int32_t (*texture_update)(rayimgui_session_handle, rayimgui_texture_handle,
                            const rayimgui_image_v1*, rayimgui_error_v1*);
  int32_t (*texture_destroy)(rayimgui_session_handle, rayimgui_texture_handle*,
                             rayimgui_error_v1*);
  int32_t (*widget)(rayimgui_session_handle, rayimgui_widget_v1*, rayimgui_item_v1*,
                    rayimgui_error_v1*);
  int32_t (*gizmo)(rayimgui_session_handle, const rayimgui_gizmo_v1*, rayimgui_item_v1*,
                   rayimgui_error_v1*);
  int32_t (*viewport)(rayimgui_session_handle, rayimgui_viewport_v1*,
                      rayimgui_error_v1*);
  int32_t (*event_push)(rayimgui_session_handle, const rayimgui_event_v1*,
                        rayimgui_error_v1*);
  int32_t (*event_poll)(rayimgui_session_handle, rayimgui_event_v1*,
                        rayimgui_error_v1*);
  int32_t (*input)(rayimgui_session_handle, rayimgui_input_v1*, rayimgui_error_v1*);
} rayimgui_api_v1;

typedef const rayimgui_api_v1* (*rayimgui_get_api_v1_fn)(uint32_t major);
#define RAYIMGUI_MEMBER_SIZE(member)                                                   \
  (offsetof(rayimgui_api_v1, member) + sizeof(((rayimgui_api_v1*)0)->member))

/* The caller guarantees that candidate points to readable prefix_bytes bytes.
 * A native forged pointer is outside the contract. No functions are called.
 */
static inline int32_t rayimgui_validate_api_v1(const void* candidate,
                                               uint32_t prefix_bytes,
                                               uint32_t required_size,
                                               uint64_t capabilities) {
  const rayimgui_api_header_v1* header;
  if (!candidate || prefix_bytes < sizeof(rayimgui_api_header_v1)) {
    return RAYIMGUI_ABI;
  }
  header = (const rayimgui_api_header_v1*)candidate;
  if (header->abi_major != 1 || header->struct_size < sizeof(*header) ||
      header->struct_size < required_size ||
      (header->capabilities & capabilities) != capabilities) {
    return RAYIMGUI_ABI;
  }
  return RAYIMGUI_OK;
}

#ifdef __cplusplus
}
#endif
#endif
