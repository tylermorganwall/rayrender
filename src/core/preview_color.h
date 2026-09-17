#ifndef RAYRENDER_PREVIEW_COLOR_H
#define RAYRENDER_PREVIEW_COLOR_H

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

// Display conversion for linear sRGB renderer output. Keep the tone curves in
// sync with rayimage::render_tonemap; the GUI receives display-encoded pixels.
class PreviewColorTransform {
public:
  using RGB = std::array<double, 3>;
  enum class ToneMap { Raw, Hbd, Reinhard, Uncharted };

  void SetToneMap(const std::string& name) {
    if (name == "raw") {
      method = ToneMap::Raw;
    } else if (name == "hbd") {
      method = ToneMap::Hbd;
    } else if (name == "reinhard") {
      method = ToneMap::Reinhard;
    } else if (name == "uncharted") {
      method = ToneMap::Uncharted;
    } else {
      throw std::invalid_argument("Unknown preview tone map: " + name);
    }
  }

  bool NeedsLogAverage() const {
    return method == ToneMap::Reinhard || method == ToneMap::Uncharted;
  }

  static double Luminance(const RGB& rgb) {
    // rayimage::CS_SRGB$rgb_to_xyz[2, ], including its D65 precision.
    return 0.2126390058715103 * rgb[0] + 0.7151686787677559 * rgb[1] +
           0.07219231536073371 * rgb[2];
  }

  static double EncodeSRGB(double linear) {
    if (std::isnan(linear)) {
      return 0;
    }
    linear = std::max(0.0, std::min(1.0, linear));
    return linear <= 0.0031308 ? 12.92 * linear
                               : 1.055 * std::pow(linear, 1.0 / 2.4) - 0.055;
  }

  // Normalize luminance for exposure-dependent curves and scale all channels
  // together to preserve chromaticity, then encode pixels for display.
  RGB Apply(RGB rgb) const {
    if (NeedsLogAverage()) {
      const double luminance = Luminance(rgb);
      const double lm = (0.18 / log_average) * luminance;
      const double ld = method == ToneMap::Reinhard
                            ? lm * (1 + lm / (11.2 * 11.2)) / (1 + lm)
                            : Hable(2 * lm) / Hable(2 * 11.2);
      const double scale = ld / std::max(luminance, 1e-6);
      for (double& value : rgb) {
        value *= scale;
      }
    }
    for (double& value : rgb) {
      if (method == ToneMap::Hbd) {
        // The HBD curve already includes the display transfer. rayimage
        // decodes this into linear output, then re-encodes it for display.
        const double x = std::max(value - 0.004, 0.0);
        value = (x * (6.2 * x + 0.5)) / (x * (6.2 * x + 1.7) + 0.06);
      } else {
        value = EncodeSRGB(value);
      }
      if (std::isnan(value)) {
        value = 0;
      }
    }
    return rgb;
  }

  // Prepared from the current image before Apply() is called for individual pixels.
  double log_average = 1;

private:
  static double Hable(double x) {
    return ((x * (0.15 * x + 0.10 * 0.50) + 0.20 * 0.02) /
            (x * (0.15 * x + 0.50) + 0.20 * 0.30)) -
           0.02 / 0.30;
  }

  ToneMap method = ToneMap::Raw;
};

#endif
