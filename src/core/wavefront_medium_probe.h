#ifndef RAYRENDER_WAVEFRONT_MEDIUM_PROBE_H
#define RAYRENDER_WAVEFRONT_MEDIUM_PROBE_H
#include "wavefront.h"
#if defined(NOT_CRAN) && defined(RAY_HAS_METAL_BVH)
// Development-only probe of successive majorant flights on a fixed world ray.
struct WFMediumProbe {
  WFMedium medium;
  WFVector origin, direction, uniforms[4]; // origin.w is the extinction scale.
};
struct WFMediumResult { WFVector distances[4], majorants[4]; };
std::vector<WFMediumResult> ProbeMetalMedium(const std::vector<WFMediumProbe> &,
                                           const std::vector<float> &);
#endif
#endif
