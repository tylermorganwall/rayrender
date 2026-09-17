#ifndef RAYRENDER_PREVIEW_SKY_H
#define RAYRENDER_PREVIEW_SKY_H
#include <Rcpp.h>
#include <cmath>
#include <stdexcept>

// Edit the prepared sky and its unattenuated solar disks together. Atmospheric
// transmission and sampling are rebuilt by the existing BuildInfiniteLights.
inline Rcpp::List PreviewSunDescriptions(const Rcpp::List& descriptions, R_xlen_t index,
                                         double elevation, double azimuth) {
  if (!std::isfinite(elevation) || elevation < -90 || elevation > 90 ||
      !std::isfinite(azimuth) || azimuth < 0 || azimuth > 360) {
    throw std::runtime_error("Invalid preview Sun direction.");
  }

  // Clone before editing: failed sky construction must leave the active R
  // descriptions available for the next frame and for subsequent edits.
  Rcpp::List updated = Rcpp::clone(descriptions);
  Rcpp::List atmosphere = updated[index];
  atmosphere["elevation"] = elevation;
  atmosphere["azimuth"] = azimuth;
  // Convert elevation/azimuth to the sky model's direction convention: Y is up
  // and increasing azimuth turns the horizontal direction toward negative X.
  const double radians = 3.14159265358979323846 / 180;
  Rcpp::NumericVector direction = Rcpp::NumericVector::create(
      -std::sin(azimuth * radians) * std::cos(elevation * radians),
      std::sin(elevation * radians),
      std::cos(azimuth * radians) * std::cos(elevation * radians));
  // Solar disks are separate emitters. Keep their direction and environment
  // rotation in step with the atmosphere while leaving other disks untouched.
  for (R_xlen_t i = 0; i < updated.size(); ++i) {
    Rcpp::List light = updated[i];
    if (Rcpp::as<std::string>(light["type"]) == "disk" &&
        light.containsElementNamed("radiance_spectrum") &&
        Rcpp::as<std::string>(light["radiance_spectrum"]) == "sun") {
      light["direction"] = Rcpp::clone(direction);
      light["rotation"] = atmosphere["rotation"];
    }
  }

  return updated;
}
#endif
