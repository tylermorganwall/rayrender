#ifndef RAYRENDER_BASE_SAMPLING_H
#define RAYRENDER_BASE_SAMPLING_H

#include "bxdf_flags.h"
#include "color_types.h"

#include "../math/vectypes.h"

#include <optional>

namespace rayrender {
namespace base {

using Vector3f = vec3f;

template <typename SpectrumValue, typename DirectionValue>
struct ScatteringSample {
  SpectrumValue f{};
  DirectionValue wi{};
  Float pdf = 0;
  BxDFFlags flags = BxDFFlags::Unset;
  Float eta = 1;
  bool pdfIsProportional = false;

  bool IsReflection() const {
    return HasFlag(flags, BxDFFlags::Reflection);
  }

  bool IsTransmission() const {
    return HasFlag(flags, BxDFFlags::Transmission);
  }

  bool IsSpecular() const {
    return HasFlag(flags, BxDFFlags::Specular);
  }
};

template <typename SampleType>
using SampleResult = std::optional<SampleType>;

using RGBBSDFSample = ScatteringSample<RGB, Vector3f>;
using RGBBSDFSampleResult = SampleResult<RGBBSDFSample>;

} // namespace base
} // namespace rayrender

#endif
