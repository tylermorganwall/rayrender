#ifndef RAYRENDER_PREVIEW_SKY_H
#define RAYRENDER_PREVIEW_SKY_H
#include <Rcpp.h>
#include <cmath>
#include <stdexcept>

// Edit the prepared sky and its unattenuated solar disks together. Atmospheric
// transmission and sampling are rebuilt by the existing BuildInfiniteLights.
inline Rcpp::List PreviewSunDescriptions(const Rcpp::List& descriptions,
                                        R_xlen_t index,double elevation,double azimuth) {
  if(!std::isfinite(elevation)||elevation < -90||elevation > 90||
     !std::isfinite(azimuth)||azimuth < 0||azimuth > 360)
    throw std::runtime_error("Invalid preview Sun direction.");
  Rcpp::List updated=Rcpp::clone(descriptions);
  Rcpp::List atmosphere=updated[index];
  atmosphere["elevation"]=elevation;atmosphere["azimuth"]=azimuth;
  const double radians=3.14159265358979323846/180;
  Rcpp::NumericVector direction=Rcpp::NumericVector::create(
    -std::sin(azimuth*radians)*std::cos(elevation*radians),
    std::sin(elevation*radians),std::cos(azimuth*radians)*std::cos(elevation*radians));
  for(R_xlen_t i=0;i<updated.size();++i) {
    Rcpp::List light=updated[i];
    if(Rcpp::as<std::string>(light["type"])=="disk" &&
       light.containsElementNamed("radiance_spectrum") &&
       Rcpp::as<std::string>(light["radiance_spectrum"])=="sun") {
      light["direction"]=Rcpp::clone(direction);
      light["rotation"]=atmosphere["rotation"];
    }
  }
  return updated;
}
#endif
