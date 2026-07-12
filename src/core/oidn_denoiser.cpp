#include "../core/oidn_denoiser.h"

#ifdef HAS_OIDN

#include <Rcpp.h>

RayOidnDenoiser::RayOidnDenoiser() : prefilter_aux(false), ready(false) {}

oidn::Quality RayOidnDenoiser::OidnQualityValue(RayOidnQuality quality) const {
  switch(quality) {
  case RayOidnQuality::Fast:
    return oidn::Quality::Fast;
  case RayOidnQuality::Balanced:
    return oidn::Quality::Balanced;
  case RayOidnQuality::High:
    return oidn::Quality::High;
  }
  return oidn::Quality::Default;
}

void RayOidnDenoiser::Setup(RayMatrix& color,
                            RayMatrix& albedo,
                            RayMatrix& normal,
                            RayMatrix& output,
                            std::size_t width,
                            std::size_t height,
                            RayOidnQuality quality,
                            bool clean_aux,
                            bool _prefilter_aux) {
  prefilter_aux = _prefilter_aux;
  ready = false;

  device = oidn::newDevice();
  device.commit();
  color_buffer = device.newBuffer(color.begin(), width * height * 3 * sizeof(Float));
  albedo_buffer = device.newBuffer(albedo.begin(), width * height * 3 * sizeof(Float));
  normal_buffer = device.newBuffer(normal.begin(), width * height * 3 * sizeof(Float));
  output_buffer = device.newBuffer(output.begin(), width * height * 3 * sizeof(Float));

  if(prefilter_aux) {
    albedo_filter = device.newFilter("RT");
    albedo_filter.setImage("albedo",
                           albedo_buffer,
                           oidn::Format::Float3,
                           width,
                           height);
    albedo_filter.setImage("output",
                           albedo_buffer,
                           oidn::Format::Float3,
                           width,
                           height);
    albedo_filter.set("quality", OidnQualityValue(quality));
    albedo_filter.commit();

    normal_filter = device.newFilter("RT");
    normal_filter.setImage("normal",
                           normal_buffer,
                           oidn::Format::Float3,
                           width,
                           height);
    normal_filter.setImage("output",
                           normal_buffer,
                           oidn::Format::Float3,
                           width,
                           height);
    normal_filter.set("quality", OidnQualityValue(quality));
    normal_filter.commit();
  } else {
    albedo_filter = oidn::FilterRef();
    normal_filter = oidn::FilterRef();
  }

  beauty_filter = device.newFilter("RT");
  beauty_filter.setImage("color",
                         color_buffer,
                         oidn::Format::Float3,
                         width,
                         height);
  beauty_filter.setImage("albedo",
                         albedo_buffer,
                         oidn::Format::Float3,
                         width,
                         height);
  beauty_filter.setImage("normal",
                         normal_buffer,
                         oidn::Format::Float3,
                         width,
                         height);
  beauty_filter.setImage("output",
                         output_buffer,
                         oidn::Format::Float3,
                         width,
                         height);
  beauty_filter.set("hdr", true);
  beauty_filter.set("cleanAux", clean_aux);
  beauty_filter.set("quality", OidnQualityValue(quality));
  beauty_filter.commit();
  ready = true;
}

void RayOidnDenoiser::Execute() {
  if(!ready) {
    return;
  }
  if(prefilter_aux) {
    albedo_filter.execute();
    normal_filter.execute();
  }
  beauty_filter.execute();
}

bool RayOidnDenoiser::ReportError() {
  if(!ready) {
    return false;
  }
  const char* error_message = nullptr;
  if(device.getError(error_message) != oidn::Error::None) {
    Rcpp::Rcout << "Error: " << error_message << std::endl;
    return true;
  }
  return false;
}

#endif
