#include "cie.h"
#include "medium.h"
#include <algorithm>
#include <array>

namespace {
point3f integrate_blackbody(double temperature) {
  if (temperature <= 100)
    return point3f(0);
  constexpr double c2 = 0.01438776877;
  double peak = 2.897771955e-3 / temperature;
  double peak_factor = std::expm1(c2 / (peak * temperature));
  double xyz[3] = {0, 0, 0};
  for (int i = 0; i < 471; ++i) {
    double lambda = (360 + i) * 1e-9;
    double exponent = c2 / (lambda * temperature);
    double b = exponent > 700 ? 0 : std::pow(peak / lambda, 5) * peak_factor / std::expm1(exponent);
    xyz[0] += b * CIE_X[i];
    xyz[1] += b * CIE_Y[i];
    xyz[2] += b * CIE_Z[i];
  }
  for (double &v : xyz)
    v /= 106.856895;
  return point3f(std::max(0.0, 3.2404542 * xyz[0] - 1.5371385 * xyz[1] - 0.4985314 * xyz[2]),
                 std::max(0.0, -0.9692660 * xyz[0] + 1.8760108 * xyz[1] + 0.0415560 * xyz[2]),
                 std::max(0.0, 0.0556434 * xyz[0] - 0.2040259 * xyz[1] + 1.0572252 * xyz[2]));
}
} // namespace
point3f BlackbodyRGB(Float kelvin) {
  // One shared table avoids spectral integration at every tracking event.
  static const auto table = [] {
    std::array<point3f, 1991> t;
    for (size_t i = 0; i < t.size(); ++i)
      t[i] = integrate_blackbody(100 + 10 * i);
    return t;
  }();
  if (kelvin <= 100)
    return point3f(0);
  if (kelvin >= 20000)
    return integrate_blackbody(kelvin);
  Float p = (kelvin - 100) / 10;
  size_t i = size_t(p);
  return table[i] * (1 - (p - i)) + table[i + 1] * (p - i);
}
