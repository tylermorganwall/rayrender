#ifndef RAYRENDER_SUN_DIRECTION_H
#define RAYRENDER_SUN_DIRECTION_H

#include <algorithm>

// Keep the solar direction away from the exact zenith in both sky models.
// A tenth of a degree remains distinct from the pole after conversion to Float
// and is visible in the editor's one-decimal input. Keep the R image-generation
// boundary in clamp_sky_sun_elevation() synchronized with this limit.
constexpr double MaxSunElevationDegrees = 89.9;

// Call after validating the physical [-90, 90] range. Night-time elevations
// remain unchanged; only the small cap around the zenith is adjusted.
inline double ClampSunElevation(double elevation) {
  return std::min(elevation, MaxSunElevationDegrees);
}

#endif
