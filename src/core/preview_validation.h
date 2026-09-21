#ifndef RAYRENDER_PREVIEW_VALIDATION_H
#define RAYRENDER_PREVIEW_VALIDATION_H
#include "preview_object_state.h"
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

// Diagnostics identify fields by their stable material schema names. The GUI
// keeps the invalid draft; exceptions only abort an unpublished scene candidate.
struct PreviewFieldError : std::runtime_error {
  std::vector<std::string> fields;
  PreviewFieldError(std::vector<std::string> names, const std::string& message)
      : std::runtime_error(message), fields(std::move(names)) {
  }
};

inline bool PreviewFieldVisible(const PreviewField& field,
                                const std::vector<PreviewField>& fields) {
  if (field.condition.empty()) {
    return true;
  }
  const auto control = std::find_if(fields.begin(), fields.end(), [&](const auto& f) {
    return f.name == field.condition;
  });
  return control != fields.end() && std::isfinite(control->values[0]) &&
         std::find(field.visible_choices.begin(),
                   field.visible_choices.end(),
                   control->values[0]) != field.visible_choices.end();
}

// Cheap value checks also run inside the draw callback, so red feedback does not
// wait for a rendered sample. Filesystem access and texture decoding run later.
inline bool PreviewValidateMaterialFields(std::vector<PreviewField>& fields) {
  auto find = [&](const char* name) -> PreviewField* {
    auto it = std::find_if(fields.begin(), fields.end(), [&](const auto& f) {
      return f.name == name;
    });
    return it == fields.end() ? nullptr : &*it;
  };
  for (auto& field : fields) {
    field.error.clear();
    if (!PreviewFieldVisible(field, fields)) {
      continue;
    }
    field.error = field.input_error.empty() ? field.load_error : field.input_error;
    if (field.text_input) {
      if (field.text.size() >= 4096 || field.text.find('\0') != std::string::npos) {
        field.error =
            "Texture paths must be shorter than 4096 bytes and contain no NUL.";
      } else if (field.text.empty() && !field.texture_available) {
        field.error = "Choose a texture file.";
      }
      continue;
    }
    for (unsigned i = 0; i < field.count; ++i) {
      const double value = field.values[i];
      if (!std::isfinite(value)) {
        field.error = "Enter a finite number.";
      } else if (value < field.minimum || value > field.maximum) {
        std::ostringstream message;
        message << "Value must be between " << field.minimum << " and " << field.maximum
                << ".";
        field.error = message.str();
      } else if ((field.boolean || field.integer || !field.choices.empty()) &&
                 value != std::floor(value)) {
        field.error = "Enter a whole number.";
      }
    }
  }
  auto range = find("Roughness map range");
  if (range && PreviewFieldVisible(*range, fields) &&
      range->values[0] > range->values[1]) {
    range->error = "Roughness map minimum must not exceed its maximum.";
  }
  auto start = find("Gradient start XYZ"), end = find("Gradient end XYZ");
  if (start && end && PreviewFieldVisible(*start, fields)) {
    double distance = 0;
    for (int i = 0; i < 3; ++i) {
      distance += std::pow(end->values[i] - start->values[i], 2);
    }
    if (distance < 1e-12) {
      start->error = end->error = "Gradient start and end must be different points.";
    }
  }
  auto direction = find("Spotlight direction XYZ");
  if (direction) {
    double length = 0;
    for (double v : direction->values) {
      length += v * v;
    }
    if (length < 1e-12) {
      direction->error = "Spotlight direction must be nonzero.";
    }
  }
  auto width = find("Spotlight width (degrees)"),
       falloff = find("Falloff start (degrees)");
  if (width && falloff && falloff->values[0] > width->values[0]) {
    width->error = falloff->error = "Spotlight falloff must start within its width.";
  }
  return std::none_of(fields.begin(), fields.end(), [](const auto& f) {
    return !f.error.empty();
  });
}

inline void PreviewMaterialChanged(PreviewObjectState& object,
                                   PreviewMaterialPanel& panel, PreviewField& field) {
  // A new draft gets one new validation/load attempt. Failed loads are not
  // repeated every GUI frame while the user is reading the diagnostic.
  for (auto& value : panel.fields) {
    value.load_error.clear();
  }
  field.changed = true;
  field.mixed = false;
  panel.changed = object.material_pending = object.apply_material = true;
  PreviewValidateMaterialFields(panel.fields);
}
#endif
