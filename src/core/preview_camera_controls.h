#ifndef RAYRENDER_PREVIEW_CAMERA_CONTROLS_H
#define RAYRENDER_PREVIEW_CAMERA_CONTROLS_H

#include <Rcpp.h>
#include "../math/float.h"
#include <array>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

// The inspector holds draft numbers, not a live camera. Validation can run while
// workers render; applying a valid draft happens only at a renderer checkpoint.
struct PreviewCameraInputs {
  enum Field { Position, Target, Up, Fov, Aperture, Focus, Ortho, Film, Scale, Count };
  enum Projection { Perspective, Orthographic, Environment, Realistic };
  struct FieldInfo {
    const char* label;
    unsigned count;
    double minimum, maximum, speed;
    std::array<const char*, 3> keys;
  };
  inline static const std::array<FieldInfo, Count> fields{
      {{"Position (X / Y / Z)", 3, -1e12, 1e12, .05, {"x", "y", "z"}},
       {"Target (X / Y / Z)", 3, -1e12, 1e12, .05, {"dx", "dy", "dz"}},
       {"Up direction (X / Y / Z)", 3, -1e6, 1e6, .01, {"upx", "upy", "upz"}},
       {"Field of view (degrees)", 1, 0, 179.9, .25, {"fov", nullptr, nullptr}},
       {"Aperture", 1, 0, 1e12, .01, {"aperture", nullptr, nullptr}},
       {"Focus distance", 1, .001, 1e12, .05, {"focal", nullptr, nullptr}},
       {"Ortho size (width / height)",
        2,
        .001,
        1e12,
        .05,
        {"orthox", "orthoy", nullptr}},
       {"Film diagonal (mm)", 1, .001, 1000, .1, {"film_size", nullptr, nullptr}},
       {"Camera scale", 1, .000001, 1e6, .01, {"camera_scale", nullptr, nullptr}}}};
  std::array<std::array<double, 3>, Count> values{};
  std::array<std::string, Count> errors{}, input_errors{};
  Projection projection = Perspective;
  int32_t model = 0;
  std::vector<std::string> models{
      "Perspective", "Orthographic", "Environment (360 degrees)"};
  std::vector<double> lens_apertures;
  std::string optical_error;

  void Select(int32_t value) {
    const int32_t previous = model;
    model = value;
    projection = model < 3 ? Projection(model) : Realistic;
    values[Fov][0] =
        model == 1   ? 0
        : model == 2 ? 360
        : model >= 3
            ? -1
            : (values[Fov][0] > 0 && values[Fov][0] < 180 ? values[Fov][0] : 45);
    if (model >= 3 && model != previous) {
      values[Aperture][0] = lens_apertures.at(model - 3) / 2;
    } else if (model == 0 && previous >= 3) {
      values[Aperture][0] = 0;
    }
    if ((model == 0 || model >= 3) && values[Focus][0] <= 0) {
      double distance = 0;
      for (size_t axis = 0; axis < 3; ++axis) {
        const double delta = values[Position][axis] - values[Target][axis];
        distance += delta * delta;
      }
      values[Focus][0] = std::max(.001, std::sqrt(distance));
    }
    optical_error.clear();
    pending = true;
  }
  bool available = false, editable = false, pending = false, editing = false;

  bool Visible(size_t field) const {
    if (field == Fov) {
      return projection == Perspective;
    }
    if (field == Aperture || field == Focus) {
      return projection == Perspective || projection == Realistic;
    }
    if (field == Film || field == Scale) {
      return projection == Realistic;
    }
    if (field == Ortho) {
      return projection == Orthographic;
    }
    return true;
  }

  const char* ProjectionName() const {
    switch (projection) {
    case Perspective:
      return "Perspective";
    case Orthographic:
      return "Orthographic";
    case Environment:
      return "Environment (360 degrees)";
    default:
      return "Physical lens camera";
    }
  }

  // Keyframe/export names provide one mapping for navigation, history and edits.
  void Read(const Rcpp::List& state) {
    for (size_t field = 0; field < Count; ++field) {
      for (unsigned axis = 0; axis < fields[field].count; ++axis) {
        const char* key = fields[field].keys[axis];
        if (state.containsElementNamed(key)) {
          values[field][axis] = Rcpp::as<double>(state[key]);
        }
      }
    }
  }

  void Write(Rcpp::List& state) const {
    state["camera_model"] = model;
    state["fov"] = values[Fov][0];
    for (size_t field = 0; field < Count; ++field) {
      if (!Visible(field)) {
        continue;
      }
      for (unsigned axis = 0; axis < fields[field].count; ++axis) {
        state[fields[field].keys[axis]] = values[field][axis];
      }
    }
  }

  Rcpp::NumericMatrix Numbers() const {
    Rcpp::NumericMatrix result(Count, 3);
    for (size_t field = 0; field < Count; ++field) {
      for (int axis = 0; axis < 3; ++axis) {
        result(field, axis) = values[field][axis];
      }
    }
    return result;
  }

  void Restore(const Rcpp::NumericMatrix& saved) {
    for (size_t field = 0; field < Count; ++field) {
      for (int axis = 0; axis < 3; ++axis) {
        values[field][axis] = saved(field, axis);
      }
    }
    input_errors.fill({});
    Validate();
  }

  bool Valid() const {
    return std::all_of(errors.begin(), errors.end(), [](const auto& error) {
      return error.empty();
    });
  }

  bool Validate() {
    errors = input_errors;
    for (size_t field = 0; field < Count; ++field) {
      if (!Visible(field)) {
        errors[field].clear();
        continue;
      }
      const auto& info = fields[field];
      for (unsigned axis = 0; axis < info.count; ++axis) {
        const double value = values[field][axis];
        if (!std::isfinite(value)) {
          errors[field] = "Enter a finite number.";
        } else if (value < info.minimum || value > info.maximum) {
          errors[field] = std::string(info.label) + " must be between " +
                          std::to_string(info.minimum) + " and " +
                          std::to_string(info.maximum) + ".";
        }
      }
    }
    if (projection == Perspective && values[Fov][0] > 0 && values[Fov][0] < .1) {
      errors[Fov] =
          "Use 0 for orthographic, or a field of view from 0.1 to 179.9 degrees.";
    }
    if (projection == Realistic && values[Aperture][0] <= 0) {
      errors[Aperture] = "Lens aperture must be greater than zero (millimeters).";
    }
    if (projection == Realistic && model >= 3 &&
        size_t(model - 3) < lens_apertures.size() &&
        values[Aperture][0] > lens_apertures[model - 3]) {
      errors[Aperture] = "Aperture exceeds the maximum opening of this lens.";
    }
    // Test the pose after Float conversion, as it will be stored by the renderer.
    // Distinct doubles must not collapse to the same point in a float build.
    std::array<double, 3> direction{}, up{};
    double length2 = 0, up2 = 0, cross2 = 0;
    for (int i = 0; i < 3; ++i) {
      direction[i] =
          double(Float(values[Target][i])) - double(Float(values[Position][i]));
      up[i] = double(Float(values[Up][i]));
      length2 += direction[i] * direction[i];
      up2 += up[i] * up[i];
    }
    if (length2 < 1e-12) {
      errors[Position] = errors[Target] =
          "Position and target must be different points.";
    }
    if (up2 < 1e-12) {
      errors[Up] = "Up direction must be nonzero.";
    }
    for (int i = 0; i < 3; ++i) {
      const double cross = direction[(i + 1) % 3] * up[(i + 2) % 3] -
                           direction[(i + 2) % 3] * up[(i + 1) % 3];
      cross2 += cross * cross;
    }
    if (length2 >= 1e-12 && up2 >= 1e-12 && cross2 <= 1e-12 * length2 * up2) {
      errors[Up] = "Up direction must not be parallel to the viewing direction.";
    }
    if (!optical_error.empty()) {
      errors[Focus] = optical_error;
    }
    return Valid();
  }
};

#endif
