#ifndef RAYRENDER_PREVIEW_SKY_CONTROLS_H
#define RAYRENDER_PREVIEW_SKY_CONTROLS_H

#include "PreviewDisplay.h"
#include "preview_sky.h"
#include "../hitables/infinite_area_light.h"
#include "../volumes/boundary.h"

// Keep committed sky settings separate from the widgets' draft values. A failed
// model or file load must leave both the active lighting and these settings intact.
struct PreviewSkySettings {
  int model = 0;
  double latitude = 0, longitude = 0, elevation = 0, azimuth = 0;
  std::string datetime;
  bool manual = false, haze = true, altitude = true;
  double base_altitude = 0, meters_per_unit = 1;
};

inline void ConfigureNativeSky(PreviewDisplay& display,
                               std::shared_ptr<InfiniteAreaLight> environment,
                               std::shared_ptr<VolumeScene> volume,
                               TextureCache& textures, const Rcpp::List& controls,
                               const Rcpp::List& original) {
  const R_xlen_t index = Rcpp::as<int>(controls["index"]);
  auto descriptions = std::make_shared<Rcpp::List>(original);
  auto settings = std::make_shared<PreviewSkySettings>();
  settings->model = Rcpp::as<int>(controls["model"]);
  settings->latitude = Rcpp::as<double>(controls["latitude"]);
  settings->longitude = Rcpp::as<double>(controls["longitude"]);
  settings->datetime = Rcpp::as<std::string>(controls["datetime"]);
  settings->elevation = Rcpp::as<double>(controls["elevation"]);
  settings->azimuth = Rcpp::as<double>(controls["azimuth"]);
  settings->manual = controls.containsElementNamed("manual") &&
                     Rcpp::as<bool>(controls["manual"]);
  if (controls.containsElementNamed("haze")) {
    settings->haze = Rcpp::as<bool>(controls["haze"]);
  }
  if (controls.containsElementNamed("altitude")) {
    settings->altitude = Rcpp::as<bool>(controls["altitude"]);
  }
  if (controls.containsElementNamed("base_altitude")) {
    settings->base_altitude = Rcpp::as<double>(controls["base_altitude"]);
  }
  if (controls.containsElementNamed("meters_per_unit")) {
    settings->meters_per_unit = Rcpp::as<double>(controls["meters_per_unit"]);
  }
  display.export_sky = [settings] {
    return Rcpp::List::create(Rcpp::_["model"] = settings->model,
                              Rcpp::_["latitude"] = settings->latitude,
                              Rcpp::_["longitude"] = settings->longitude,
                              Rcpp::_["datetime"] = settings->datetime,
                              Rcpp::_["manual"] = settings->manual,
                              Rcpp::_["elevation"] = settings->elevation,
                              Rcpp::_["azimuth"] = settings->azimuth,
                              Rcpp::_["haze"] = settings->haze,
                              Rcpp::_["altitude"] = settings->altitude,
                              Rcpp::_["base_altitude"] = settings->base_altitude,
                              Rcpp::_["meters_per_unit"] = settings->meters_per_unit);
  };
  Rcpp::List initial = original[index];
  if (Rcpp::as<std::string>(initial["type"]) == "prague") {
    settings->haze = Rcpp::as<bool>(initial["haze"]);
    settings->altitude = Rcpp::as<bool>(initial["query_altitude"]);
    settings->base_altitude = Rcpp::as<double>(initial["altitude"]);
    settings->meters_per_unit = Rcpp::as<double>(initial["meters_per_unit"]);
  }

  // Build the complete replacement before swapping either owner. The native
  // atmosphere pointer belongs to the same light snapshot as emitted radiance.
  auto prepare_publish = [descriptions, environment, volume, &textures](Rcpp::List updated) {
    auto light = BuildInfiniteLights(updated, textures);
    return [descriptions, environment, volume, updated, light] {
      environment->SetLight(light);
      if (volume) volume->atmosphere = light->GetTransportAtmosphere();
      *descriptions = updated;
    };
  };
  auto publish = [prepare_publish](Rcpp::List updated) { prepare_publish(updated)(); };
  auto configure_atmosphere = [&, descriptions, settings, publish, index] {
    Rcpp::List current = (*descriptions)[index];
    const bool native = Rcpp::as<std::string>(current["type"]) == "prague";
    if (native) {
      display.SetAtmosphereControls(
          settings->haze,
          settings->altitude,
          [descriptions, settings, publish, index](bool haze, bool altitude) {
            Rcpp::List updated = Rcpp::clone(*descriptions);
            Rcpp::List sky = updated[index];
            sky["haze"] = haze;
            sky["query_altitude"] = altitude;
            publish(updated);
            settings->haze = haze;
            settings->altitude = altitude;
          });
    } else {
      display.SetAtmosphereControls(false, false, {});
    }
    display.native_gui->has_atmosphere = native;
    display.native_gui->haze = native && settings->haze;
    display.native_gui->altitude = native && settings->altitude;
    display.native_gui->atmosphere_pending = false;
    display.native_gui->base_altitude = settings->base_altitude;
    display.native_gui->meters_per_unit = settings->meters_per_unit;
    display.native_gui->atmosphere_parameters_pending = false;
    display.native_gui->atmosphere_parameters_editing = false;
  };
  Rcpp::Function update = controls["update"];
  auto prepare_restore =
      [settings, update, prepare_publish, configure_atmosphere, index, &display](
          PreviewSkySettings next, bool fast = false) -> std::function<void()> {
    Rcpp::RObject elevation = next.manual ? Rcpp::wrap(next.elevation) : R_NilValue;
    Rcpp::RObject azimuth = next.manual ? Rcpp::wrap(next.azimuth) : R_NilValue;
    Rcpp::List result = update(next.latitude,
                               next.longitude,
                               next.datetime,
                               next.model,
                               elevation,
                               azimuth,
                               Rcpp::_["base_altitude"] = next.base_altitude,
                               Rcpp::_["meters_per_unit"] = next.meters_per_unit,
                               Rcpp::_["fast"] = fast);
    const std::string error = Rcpp::as<std::string>(result["error"]);
    if (!error.empty()) {
      throw std::runtime_error(error);
    }
    // Resolve/clamp the actual sun direction before publishing any lighting.
    next.elevation = Rcpp::as<double>(result["elevation"]);
    next.azimuth = Rcpp::as<double>(result["azimuth"]);
    Rcpp::List updated = result["lights"];
    Rcpp::List sky = updated[index];
    if (Rcpp::as<std::string>(sky["type"]) == "prague") {
      sky["haze"] = next.haze;
      sky["query_altitude"] = next.altitude;
      if (next.manual) {
        updated = PreviewSunDescriptions(updated, index, next.elevation, next.azimuth);
      }
    }
    auto commit_light = prepare_publish(updated);
    return [settings, next, commit_light, configure_atmosphere, &display] {
      commit_light();
      *settings = next;
      display.SetSunPosition(next.elevation, next.azimuth);
      display.native_gui->manual_sun = next.manual;
      configure_atmosphere();
    };
  };
  auto rebuild = [prepare_restore](PreviewSkySettings next,
                                   bool fast = false) -> std::string {
    try {
      prepare_restore(next, fast)();
      return {};
    } catch (const Rcpp::internal::InterruptedException&) {
      throw;
    } catch (const std::exception& error) {
      return error.what();
    }
  };
  // Undo stages the complete prior sky, including automatic/manual sun mode and
  // Prague transport switches, alongside any geometry restoration it needs.
  display.prepare_sky_restore = [prepare_restore](const Rcpp::List& saved) {
    PreviewSkySettings next;
    next.model = Rcpp::as<int>(saved["model"]);
    next.latitude = Rcpp::as<double>(saved["latitude"]);
    next.longitude = Rcpp::as<double>(saved["longitude"]);
    next.datetime = Rcpp::as<std::string>(saved["datetime"]);
    next.manual = Rcpp::as<bool>(saved["manual"]);
    next.elevation = Rcpp::as<double>(saved["elevation"]);
    next.azimuth = Rcpp::as<double>(saved["azimuth"]);
    next.haze = Rcpp::as<bool>(saved["haze"]);
    next.altitude = Rcpp::as<bool>(saved["altitude"]);
    next.base_altitude = Rcpp::as<double>(saved["base_altitude"]);
    next.meters_per_unit = Rcpp::as<double>(saved["meters_per_unit"]);
    return prepare_restore(next);
  };

  display.update_atmosphere_parameters = [settings, rebuild](double base_altitude,
                                                             double meters_per_unit) {
    auto next = *settings;
    next.base_altitude = base_altitude;
    next.meters_per_unit = meters_per_unit;
    return rebuild(next);
  };
  configure_atmosphere();
  display.SetSunControls(
      settings->elevation,
      settings->azimuth,
      [settings, rebuild, &display](double elevation, double azimuth) {
        auto next = *settings;
        next.manual = true;
        next.elevation = elevation;
        next.azimuth = azimuth;
        const auto error = rebuild(next, display.native_gui->sun_editing);
        if (!error.empty()) {
          throw std::runtime_error(error);
        }
      });
  display.SetSkyControls(settings->latitude,
                         settings->longitude,
                         settings->datetime,
                         [settings, rebuild](double latitude,
                                             double longitude,
                                             const std::string& datetime) {
                           auto next = *settings;
                           next.latitude = latitude;
                           next.longitude = longitude;
                           next.datetime = datetime;
                           next.manual = false;
                           return rebuild(next);
                         });
  display.SetSkyModelControls(settings->model, [settings, rebuild](int model) {
    auto next = *settings;
    next.model = model;
    return rebuild(next);
  });
  display.native_gui->has_sun = true;
  display.native_gui->manual_sun = settings->manual;
  // An image-based Prague scene opens with Prague already selected. Upgrade it
  // before the first sample, while the same atomic rebuild handles load errors.
  if (settings->model == 1 && Rcpp::as<std::string>(initial["type"]) != "prague") {
    display.native_gui->sky_error = rebuild(*settings);
  }
}
#endif
