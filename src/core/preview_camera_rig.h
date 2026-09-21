#ifndef RAYRENDER_PREVIEW_CAMERA_RIG_H
#define RAYRENDER_PREVIEW_CAMERA_RIG_H

#include "camera.h"
#include <Rcpp.h>
#include <array>
#include <algorithm>
#include <functional>
#include <memory>
#include <string>
#include <vector>

// Own all ordinary projections for the lifetime of a render. Lens cameras are
// prepared on demand; a failed optical configuration never replaces the view.
class PreviewCameraRig {
public:
  struct Lens {
    std::string name, source;
    std::vector<Float> data;
    double aperture = 1;
  };
  std::vector<Lens> lenses;
  int model = 0;
  double film_size = 22, camera_scale = 1;

  explicit PreviewCameraRig(const Rcpp::List& info) {
    auto point = [&](const char* key) {
      Rcpp::NumericVector p = info[key];
      return point3f(p[0], p[1], p[2]);
    };
    const point3f origin = point("lookfrom"), target = point("lookat");
    const vec3f up = convert_to_vec3(point("camera_up"));
    const Float fov = Rcpp::as<Float>(info["fov"]);
    const Float aperture = Rcpp::as<Float>(info["aperture"]);
    const Float focus = Rcpp::as<Float>(info["focal_distance"]);
    width = Rcpp::as<Float>(info["nx"]);
    height = Rcpp::as<Float>(info["ny"]);
    shutter_open = Rcpp::as<Float>(info["shutteropen"]);
    shutter_close = Rcpp::as<Float>(info["shutterclose"]);
    iso = Rcpp::as<Float>(info["iso"]);
    film_size = Rcpp::as<double>(info["film_size"]) * 1000;
    camera_scale = Rcpp::as<double>(info["camera_scale"]);
    Rcpp::NumericVector ortho = info["ortho_dimensions"];
    ordinary[0] = std::make_shared<camera>(origin,
                                           target,
                                           up,
                                           fov > 0 && fov < 180 ? fov : 45,
                                           width / height,
                                           fov < 0 ? 0 : aperture,
                                           focus,
                                           shutter_open,
                                           shutter_close,
                                           iso);
    ordinary[1] = std::make_shared<ortho_camera>(
        origin, target, up, ortho[0], ortho[1], shutter_open, shutter_close, iso);
    ordinary[2] = std::make_shared<environment_camera>(
        origin, target, up, shutter_open, shutter_close, iso);
    start_transform = LookAt(origin, target, up).GetInverseMatrix();

    if (info.containsElementNamed("preview_lenses")) {
      Rcpp::List catalog = info["preview_lenses"];
      for (R_xlen_t i = 0; i < catalog.size(); ++i) {
        Rcpp::List lens = catalog[i];
        AddLens(Rcpp::as<std::string>(lens["name"]),
                Rcpp::as<std::string>(lens["source"]),
                lens["data"]);
      }
    }
    model = ModelForFov(fov);
    if (model == 3) {
      // Retain an arbitrary lens file supplied by the caller alongside presets.
      Rcpp::NumericMatrix original = info["real_camera_info"];
      std::string source = info.containsElementNamed("camera_description_file")
                               ? Rcpp::as<std::string>(info["camera_description_file"])
                               : "";
      auto found = std::find_if(lenses.begin(), lenses.end(), [&](const Lens& lens) {
        return lens.source == source;
      });
      if (found == lenses.end()) {
        AddLens("Original lens", source, original);
        model = int(lenses.size()) + 2;
      } else {
        model = int(found - lenses.begin()) + 3;
      }
      physical =
          MakeLens(model, aperture, focus, film_size, camera_scale, origin, target, up);
      physical_model = model;
    }
    const bool free_rotation =
        info.containsElementNamed("free_rotation") && Rcpp::as<bool>(info["free_rotation"]);
    for (auto &camera : ordinary) {
      camera->set_free_rotation(free_rotation);
    }
    if (physical) {
      physical->set_free_rotation(free_rotation);
    }
  }

  static int ModelForFov(double fov) {
    return fov < 0 ? 3 : fov == 0 ? 1 : fov == 360 ? 2 : 0;
  }

  RayCamera* Active() const {
    return model < 3 ? ordinary[model].get() : physical.get();
  }

  std::string Source() const {
    return model >= 3 ? lenses.at(model - 3).source : "";
  }

  // Allocate/validate optics first. Calling the returned function is the commit
  // point, after renderer workers and other history preparations have finished.
  std::function<RayCamera*()> Prepare(const Rcpp::List& state) {
    const bool free_rotation = Active()->get_free_rotation();
    const Float fov = Rcpp::as<Float>(state["fov"]);
    int next = ModelForFov(fov);
    if (next == 3) {
      next = state.containsElementNamed("camera_model")
                 ? Rcpp::as<int>(state["camera_model"])
                 : model;
      if (next < 3 || size_t(next - 3) >= lenses.size()) {
        throw std::runtime_error(
            "Select a physical lens before using a realistic camera.");
      }
    }
    const double film = state.containsElementNamed("film_size")
                            ? Rcpp::as<double>(state["film_size"])
                            : film_size;
    const double scale = state.containsElementNamed("camera_scale")
                             ? Rcpp::as<double>(state["camera_scale"])
                             : camera_scale;
    const Float aperture = Rcpp::as<Float>(state["aperture"]);
    const Float focus = Rcpp::as<Float>(state["focal"]);
    const point3f origin(Rcpp::as<Float>(state["x"]),
                         Rcpp::as<Float>(state["y"]),
                         Rcpp::as<Float>(state["z"]));
    const point3f target(Rcpp::as<Float>(state["dx"]),
                         Rcpp::as<Float>(state["dy"]),
                         Rcpp::as<Float>(state["dz"]));
    const vec3f up(Rcpp::as<Float>(state["upx"]),
                   Rcpp::as<Float>(state["upy"]),
                   Rcpp::as<Float>(state["upz"]));
    const vec2f ortho(Rcpp::as<Float>(state["orthox"]),
                      Rcpp::as<Float>(state["orthoy"]));
    std::shared_ptr<RealisticCamera> prepared;
    if (next >= 3) {
      if (physical && physical_model == next && film == film_size &&
          scale == camera_scale) {
        // Copy before changing optics so invalid focus/aperture can be rejected
        // without leaving a half-edited live camera or history state.
        prepared = std::make_shared<RealisticCamera>(*physical);
        prepared->update_aperture_absolute(aperture);
        prepared->update_focal_absolute(focus);
      } else {
        prepared = MakeLens(next, aperture, focus, film, scale, origin, target, up);
      }
      prepared->update_pose_absolute(origin, target, up);
    }
    return [this,
            next,
            free_rotation,
            film,
            scale,
            aperture,
            focus,
            fov,
            origin,
            target,
            up,
            ortho,
            prepared]() {
      if (prepared) {
        physical = prepared;
        physical_model = next;
      } else {
        auto& camera = ordinary[next];
        camera->update_focal_absolute(focus);
        camera->update_pose_absolute(origin, target, up);
        camera->update_aperture_absolute(aperture);
        camera->update_fov_absolute(fov);
        camera->update_ortho_absolute(ortho);
      }
      model = next;
      film_size = film;
      camera_scale = scale;
      Active()->set_free_rotation(free_rotation);
      return Active();
    };
  }

private:
  std::array<std::shared_ptr<RayCamera>, 3> ordinary;
  std::shared_ptr<RealisticCamera> physical;
  int physical_model = -1;
  Float width, height, shutter_open, shutter_close, iso;
  Transform start_transform;

  void AddLens(const std::string& name, const std::string& source,
               const Rcpp::NumericMatrix& data) {
    if (data.nrow() == 0 || data.ncol() != 4) {
      throw std::runtime_error(
          "A lens description must have four columns and at least one element.");
    }
    Lens lens;
    lens.name = name;
    lens.source = source;
    for (int row = 0; row < data.nrow(); ++row) {
      for (int col = 0; col < 4; ++col) {
        lens.data.push_back(data(row, col));
      }
      if (data(row, 0) == 0) {
        lens.aperture = data(row, 3);
      }
    }
    lenses.push_back(std::move(lens));
  }

  std::shared_ptr<RealisticCamera> MakeLens(int next, Float aperture, Float focus,
                                            double film, double scale, point3f origin,
                                            point3f target, vec3f up) {
    auto data = lenses.at(next - 3).data;
    AnimatedTransform initial(&start_transform, 0, &start_transform, 0);
    Transform pose = LookAt(origin, target, up).GetInverseMatrix();
    return std::make_shared<RealisticCamera>(initial,
                                             shutter_open,
                                             shutter_close,
                                             aperture,
                                             width,
                                             height,
                                             focus,
                                             false,
                                             data,
                                             film / 1000,
                                             scale,
                                             iso,
                                             up,
                                             pose,
                                             target);
  }
};

#endif
