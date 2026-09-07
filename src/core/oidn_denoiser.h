#ifndef RAYRENDER_OIDN_DENOISER_H
#define RAYRENDER_OIDN_DENOISER_H

#include <cstddef>

#include "../math/float.h"
#include "../math/RayMatrix.h"

#ifdef HAS_OIDN
#undef None
#include <OpenImageDenoise/oidn.hpp>

enum class RayOidnQuality {
  Fast,
  Balanced,
  High
};

class RayOidnDenoiser {
public:
  RayOidnDenoiser();

  void Setup(RayMatrix& color,
             RayMatrix& albedo,
             RayMatrix& normal,
             RayMatrix& output,
             std::size_t width,
             std::size_t height,
             RayOidnQuality quality,
             bool clean_aux,
             bool prefilter_aux,
             bool use_auxiliary = true);

  void Execute();
  bool ReportError();
  bool Ready() const { return ready; }

private:
  oidn::Quality OidnQualityValue(RayOidnQuality quality) const;

  oidn::DeviceRef device;
  oidn::BufferRef color_buffer;
  oidn::BufferRef albedo_buffer;
  oidn::BufferRef normal_buffer;
  oidn::BufferRef output_buffer;
  oidn::FilterRef beauty_filter;
  oidn::FilterRef albedo_filter;
  oidn::FilterRef normal_filter;
  bool prefilter_aux;
  bool ready;
};

#endif

#endif
