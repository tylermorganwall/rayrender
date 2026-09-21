#ifndef RAYRENDER_PREVIEW_OBJECT_STATE_H
#define RAYRENDER_PREVIEW_OBJECT_STATE_H
#include <array>
#include <string>
#include <vector>
#include <cstdint>

// UI snapshots only. Renderer pointers and R objects never enter draw callbacks.
// A scalar or short vector plus widget limits. Only the first count values are used.
struct PreviewField {
  std::string name;
  bool color = false;
  bool boolean = false, integer = false, text_input = false;
  std::string text, section, help;
  std::string error, load_error, input_error;
  bool texture_available = false;
  std::vector<std::string> choices;
  // A field may depend on another choice without changing its saved identity.
  std::string condition;
  std::vector<int> visible_choices;
  // Mixed selections display the first value until that property is edited.
  bool mixed = false, changed = false;
  unsigned count = 1;
  std::array<double, 3> values{};
  double minimum = 0, maximum = 1, speed = .01;
};
struct PreviewMaterialTarget {
  size_t row = 0, placement = 0, slot = 0;
};
// Corresponding slots in an instance collection share one inspector panel.
struct PreviewMaterialPanel {
  size_t row = 0, placement = 0, slot = 0;
  bool changed = false;
  std::string label, type;
  std::vector<PreviewField> fields;
  std::vector<PreviewMaterialTarget> targets;
};
// Stable, pointer-free hierarchy nodes. Parent zero denotes the scene root.
struct PreviewHierarchyNode {
  uint64_t id = 0, parent = 0;
  std::string label;
};
// GUI callbacks stage requests here; the renderer consumes them between samples.
// Pending flags retain drafts across frames. Changed inputs request live updates;
// only renderer checkpoints may validate/load textures or replace scene geometry.
struct PreviewObjectState {
  bool enabled = false, selected = false, pick_pending = false, clear_pending = false;
  bool transform_pending = false, transform_active = false, apply_transform = false;
  bool material_pending = false, apply_material = false, revert = false;
  bool projection_valid = false;
  float pick_u = 0, pick_v = 0;
  uint64_t id = 0, select_id = 0;
  bool select_pending = false;
  std::vector<PreviewHierarchyNode> hierarchy;
  std::string label, error;
  int32_t operation = 0, mode = 1, material_slot = 0;
  // Column-major gizmo matrices. model is editable UI state; the renderer keeps
  // its last committed selection transform separately for delta/cancel handling.
  std::array<float, 16> model{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  std::array<float, 16> view{}, projection{};
  uint32_t projection_kind = 1;
  std::array<double, 3> translation{}, rotation{}, scale{1, 1, 1};
  bool numeric_transform = false, cancel_transform = false;
  bool begin_transform = false, end_transform = false, drag_cancelled = false;
  std::vector<PreviewMaterialPanel> materials;
};
#endif
