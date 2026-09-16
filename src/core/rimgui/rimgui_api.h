/* Copyright (c) 2026 Tyler Morgan-Wall. MIT license; see LICENSE.protocol.
 * Pure C protocol. This development ABI is not yet frozen for production.
 * Complete phase/ownership rules: installed developer/ABI.md.
 */
#ifndef RIMGUI_API_H
#define RIMGUI_API_H
#include <stddef.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

#define RIMGUI_ABI_MAJOR 1u
#define RIMGUI_ABI_MINOR 2u
#define RIMGUI_CAP_HEADLESS UINT64_C(1)
#define RIMGUI_CAP_NATIVE UINT64_C(2)
#define RIMGUI_CAP_DOCKING UINT64_C(4)
#define RIMGUI_CAP_GIZMO UINT64_C(8)
#define RIMGUI_CAP_RGBA8 UINT64_C(16)
#define RIMGUI_CAP_WIDGETS UINT64_C(32)
#define RIMGUI_CAP_EVENTS UINT64_C(64)
#define RIMGUI_CAP_INPUT UINT64_C(128)

typedef uint64_t rimgui_session_handle;
typedef uint64_t rimgui_texture_handle;
typedef uint64_t rimgui_callback_handle;

enum rimgui_status {
    RIMGUI_OK = 0, RIMGUI_INVALID = 1, RIMGUI_ABI = 2,
    RIMGUI_BACKEND_NOT_BUILT = 3, RIMGUI_BACKEND_INIT = 4,
    RIMGUI_WRONG_THREAD = 5, RIMGUI_PHASE = 6, RIMGUI_STALE = 7,
    RIMGUI_BUSY = 8, RIMGUI_CALLBACK_ERROR = 9, RIMGUI_MEMORY = 10,
    RIMGUI_INTERNAL = 11, RIMGUI_EMPTY = 12, RIMGUI_LIMIT = 13
};
enum rimgui_session_flags { RIMGUI_HEADLESS = 1, RIMGUI_HIDDEN = 2 };
enum rimgui_item_flags {
    RIMGUI_CHANGED = 1, RIMGUI_ACTIVE = 2, RIMGUI_HOVERED = 4,
    RIMGUI_BEGIN = 8, RIMGUI_UPDATE = 16, RIMGUI_COMMIT = 32,
    RIMGUI_CANCEL = 64, RIMGUI_VISIBLE = 128, RIMGUI_DEACTIVATED = 256
};
enum rimgui_widget_kind {
    RIMGUI_WINDOW_BEGIN = 1, RIMGUI_WINDOW_END = 2,
    RIMGUI_CHILD_BEGIN = 3, RIMGUI_CHILD_END = 4,
    RIMGUI_ID_PUSH = 5, RIMGUI_ID_POP = 6,
    RIMGUI_TEXT = 7, RIMGUI_BUTTON = 8, RIMGUI_CHECKBOX = 9,
    RIMGUI_INT = 10, RIMGUI_DOUBLE = 11, RIMGUI_COLOR = 12,
    RIMGUI_INPUT_TEXT = 13, RIMGUI_COMBO = 14,
    RIMGUI_TREE_BEGIN = 15, RIMGUI_TREE_END = 16,
    RIMGUI_MENU_BAR_BEGIN = 17, RIMGUI_MENU_BAR_END = 18,
    RIMGUI_MENU_BEGIN = 19, RIMGUI_MENU_END = 20,
    RIMGUI_MENU_ITEM = 21, RIMGUI_SEPARATOR = 22,
    RIMGUI_SAME_LINE = 23, RIMGUI_DISABLED_BEGIN = 24,
    RIMGUI_DISABLED_END = 25, RIMGUI_PROGRESS = 26,
    RIMGUI_IMAGE = 27, RIMGUI_TOOLTIP = 28
};
enum rimgui_widget_options {
    RIMGUI_WINDOW_MENU_BAR = 1, RIMGUI_WINDOW_POSITION = 2,
    RIMGUI_WINDOW_VIEWPORT = 4, /* Route navigation keys to the image consumer. */
    RIMGUI_COLOR_RGB = 8 /* ABI 1.2: COLOR edits bounded RGB without an intensity field. */
};
enum rimgui_key {
    RIMGUI_KEY_A, RIMGUI_KEY_B, RIMGUI_KEY_C, RIMGUI_KEY_D, RIMGUI_KEY_E,
    RIMGUI_KEY_F, RIMGUI_KEY_G, RIMGUI_KEY_H, RIMGUI_KEY_I, RIMGUI_KEY_J,
    RIMGUI_KEY_K, RIMGUI_KEY_L, RIMGUI_KEY_M, RIMGUI_KEY_N, RIMGUI_KEY_O,
    RIMGUI_KEY_P, RIMGUI_KEY_Q, RIMGUI_KEY_R, RIMGUI_KEY_S, RIMGUI_KEY_T,
    RIMGUI_KEY_U, RIMGUI_KEY_V, RIMGUI_KEY_W, RIMGUI_KEY_X, RIMGUI_KEY_Y, RIMGUI_KEY_Z,
    RIMGUI_KEY_0, RIMGUI_KEY_1, RIMGUI_KEY_2, RIMGUI_KEY_3, RIMGUI_KEY_4,
    RIMGUI_KEY_5, RIMGUI_KEY_6, RIMGUI_KEY_7, RIMGUI_KEY_8, RIMGUI_KEY_9,
    RIMGUI_KEY_TAB, RIMGUI_KEY_LEFT, RIMGUI_KEY_RIGHT, RIMGUI_KEY_UP, RIMGUI_KEY_DOWN,
    RIMGUI_KEY_ENTER, RIMGUI_KEY_ESCAPE, RIMGUI_KEY_LEFT_BRACKET, RIMGUI_KEY_RIGHT_BRACKET,
    RIMGUI_KEY_COMMA, RIMGUI_KEY_PERIOD, RIMGUI_KEY_SLASH, RIMGUI_KEY_SPACE,
    RIMGUI_KEY_KEYPAD_ENTER, RIMGUI_KEY_COUNT
};
#define RIMGUI_KEY_BIT(key) (UINT64_C(1) << (key))
enum rimgui_modifiers { RIMGUI_SHIFT=1, RIMGUI_CTRL=2, RIMGUI_ALT=4, RIMGUI_SUPER=8 };
enum rimgui_mouse_button { RIMGUI_MOUSE_LEFT=1, RIMGUI_MOUSE_RIGHT=2, RIMGUI_MOUSE_MIDDLE=4 };
enum rimgui_operation { RIMGUI_TRANSLATE = 1, RIMGUI_ROTATE = 2, RIMGUI_SCALE = 3 };
enum rimgui_mode { RIMGUI_LOCAL = 1, RIMGUI_WORLD = 2 };
enum rimgui_projection { RIMGUI_PERSPECTIVE = 1, RIMGUI_ORTHOGRAPHIC = 2 };
enum rimgui_origin { RIMGUI_TOP_LEFT = 1, RIMGUI_BOTTOM_LEFT = 2 };
enum rimgui_channels { RIMGUI_RGBA = 1, RIMGUI_BGRA = 2 };
enum rimgui_alpha { RIMGUI_OPAQUE = 1, RIMGUI_STRAIGHT = 2 };
enum rimgui_color_space { RIMGUI_DISPLAY_ENCODED = 1 };

typedef struct rimgui_error_v1 {
    uint32_t struct_size;
    int32_t code;
    char *message;
    uint32_t capacity; /* Includes trailing NUL; truncation always terminates. */
} rimgui_error_v1;

typedef struct rimgui_session_desc_v1 {
    uint32_t struct_size;
    uint32_t flags;
    const char *title;
    uint32_t title_length;
    uint32_t width;
    uint32_t height;
} rimgui_session_desc_v1;

typedef struct rimgui_step_v1 {
    uint32_t struct_size;
    uint32_t close_requested;
    uint64_t frame;
    double delta_seconds;
    float width, height;
    float framebuffer_scale_x, framebuffer_scale_y;
    uint32_t focused;
    uint32_t want_keyboard;
} rimgui_step_v1;

typedef struct rimgui_image_v1 {
    uint32_t struct_size;
    uint32_t width, height;
    uint32_t origin, channels, alpha, color_space;
    uint64_t row_stride;
    uint64_t buffer_length;
    uint64_t version; /* Unchanged version means unchanged pixels. */
    const uint8_t *pixels;
} rimgui_image_v1;

typedef struct rimgui_viewport_v1 {
    uint32_t struct_size;
    float x, y, width, height; /* Image bounds in window logical coordinates. */
    float letterbox_x, letterbox_y;
    float mouse_x, mouse_y; /* Relative to the actual image, not its panel. */
    float framebuffer_scale_x, framebuffer_scale_y;
    uint32_t inside, focused, mouse_down, mouse_clicked, mouse_released;
    uint32_t widget_owns_mouse, gizmo_hovered, gizmo_active;
    uint32_t background_click; /* Valid after calling gizmo for this viewport. */
} rimgui_viewport_v1;

/* ABI 1.1: read during a draw callback, after the current window's image/gizmo.
 * Unavailable input has zero key/button masks. Positions come from viewport().
 * pressed is an edge; repeated includes the first press and keyboard repeat.
 */
typedef struct rimgui_input_v1 {
    uint32_t struct_size;
    uint32_t keyboard_available, mouse_available, modifiers;
    uint64_t keys_down, keys_pressed, keys_repeated;
    uint32_t mouse_down, mouse_clicked, mouse_released;
} rimgui_input_v1;

typedef struct rimgui_widget_v1 {
    uint32_t struct_size;
    uint32_t kind;
    uint64_t id; /* Stable application ID; zero allowed for noninteractive layout. */
    const char *label;
    uint32_t label_length;
    uint32_t options;
    int32_t *integers; /* CHECKBOX/INT/COMBO; count 1..4, COMBO uses one. */
    double *numbers; /* DOUBLE: 1..4; COLOR: RGB in [0,1], separate intensity. */
    uint32_t count;
    double minimum, maximum, speed;
    char *text; /* INPUT_TEXT: caller-owned writable UTF-8 with terminating NUL. */
    uint32_t text_capacity;
    const char *const *choices;
    const uint32_t *choice_lengths;
    uint32_t choice_count;
    float width, height, x, y; /* WINDOW_POSITION uses x/y on first appearance. */
    double fraction; /* PROGRESS; or separate COLOR intensity (input/output). */
    rimgui_texture_handle texture;
} rimgui_widget_v1;

typedef struct rimgui_item_v1 {
    uint32_t struct_size;
    uint32_t flags;
} rimgui_item_v1;

typedef struct rimgui_gizmo_v1 {
    uint32_t struct_size;
    uint32_t operation, mode, projection;
    uint64_t id;
    const float *view; /* Borrowed 16 floats, column major, right handed. */
    const float *projection_matrix; /* OpenGL clip Z [-1,1]. */
    float *model; /* 16 floats, column vectors: clip = P * V * M * point. */
    const rimgui_viewport_v1 *viewport;
    uint32_t snap_enabled;
    float snap[3]; /* Translation world units, rotation degrees, scale units. */
} rimgui_gizmo_v1;

typedef struct rimgui_event_v1 {
    uint32_t struct_size;
    uint32_t phase; /* Exactly one of BEGIN/UPDATE/COMMIT/CANCEL. */
    uint64_t application_id;
    uint32_t command;
    uint32_t count;
    double values[4]; /* Copied by event_push; interpreted only by consumer. */
} rimgui_event_v1;

struct rimgui_api_v1;
typedef int32_t (*rimgui_draw_ui_v1)(const struct rimgui_api_v1 *,
    rimgui_session_handle, void *, rimgui_error_v1 *);
typedef struct rimgui_callback_v1 {
    uint32_t struct_size;
    uint32_t abi_major;
    rimgui_draw_ui_v1 draw_ui;
    void *userdata;
} rimgui_callback_v1;

/* Fixed negotiation prefix. Read ONLY this prefix before validating table size.
 * All v1 functions are mandatory. Minor versions may append, never reorder.
 */
typedef struct rimgui_api_header_v1 {
    uint32_t struct_size;
    uint32_t abi_major;
    uint32_t abi_minor;
    uint32_t reserved;
    uint64_t capabilities;
} rimgui_api_header_v1;

typedef struct rimgui_api_v1 {
    rimgui_api_header_v1 header;
    int32_t (*open)(const rimgui_session_desc_v1 *, rimgui_session_handle *, rimgui_error_v1 *);
    int32_t (*close)(rimgui_session_handle *, rimgui_error_v1 *);
    int32_t (*register_callback)(rimgui_session_handle, const rimgui_callback_v1 *, rimgui_callback_handle *, rimgui_error_v1 *);
    int32_t (*unregister_callback)(rimgui_session_handle, rimgui_callback_handle *, rimgui_error_v1 *);
    int32_t (*step)(rimgui_session_handle, rimgui_step_v1 *, rimgui_error_v1 *);
    int32_t (*request_close)(rimgui_session_handle, rimgui_error_v1 *);
    int32_t (*texture_create)(rimgui_session_handle, const rimgui_image_v1 *, rimgui_texture_handle *, rimgui_error_v1 *);
    int32_t (*texture_update)(rimgui_session_handle, rimgui_texture_handle, const rimgui_image_v1 *, rimgui_error_v1 *);
    int32_t (*texture_destroy)(rimgui_session_handle, rimgui_texture_handle *, rimgui_error_v1 *);
    int32_t (*widget)(rimgui_session_handle, rimgui_widget_v1 *, rimgui_item_v1 *, rimgui_error_v1 *);
    int32_t (*gizmo)(rimgui_session_handle, const rimgui_gizmo_v1 *, rimgui_item_v1 *, rimgui_error_v1 *);
    int32_t (*viewport)(rimgui_session_handle, rimgui_viewport_v1 *, rimgui_error_v1 *);
    int32_t (*event_push)(rimgui_session_handle, const rimgui_event_v1 *, rimgui_error_v1 *);
    int32_t (*event_poll)(rimgui_session_handle, rimgui_event_v1 *, rimgui_error_v1 *);
    int32_t (*input)(rimgui_session_handle, rimgui_input_v1 *, rimgui_error_v1 *);
} rimgui_api_v1;

typedef const rimgui_api_v1 *(*rimgui_get_api_v1_fn)(uint32_t major);
#define RIMGUI_MEMBER_SIZE(member) (offsetof(rimgui_api_v1, member) + sizeof(((rimgui_api_v1 *)0)->member))

/* The caller guarantees that candidate points to readable prefix_bytes bytes.
 * A native forged pointer is outside the contract. No functions are called.
 */
static inline int32_t rimgui_validate_api_v1(const void *candidate,
    uint32_t prefix_bytes, uint32_t required_size, uint64_t capabilities) {
    const rimgui_api_header_v1 *header;
    if (!candidate || prefix_bytes < sizeof(rimgui_api_header_v1)) return RIMGUI_ABI;
    header = (const rimgui_api_header_v1 *)candidate;
    if (header->abi_major != 1 || header->struct_size < sizeof(*header) ||
        header->struct_size < required_size ||
        (header->capabilities & capabilities) != capabilities) return RIMGUI_ABI;
    return RIMGUI_OK;
}

#ifdef __cplusplus
}
#endif
#endif
