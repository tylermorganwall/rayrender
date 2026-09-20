#ifndef RAYRENDER_PREVIEW_SELECTION_OVERLAY_H
#define RAYRENDER_PREVIEW_SELECTION_OVERLAY_H

#include <cstddef>
#include <stdexcept>
#include <cstdint>
#include <vector>

// Selection presentation only: the outline depends on binary coverage,
// never material values, displayed colors, or accumulated radiance.
class PreviewSelectionOverlay {
public:
  // Cache two inward border bands from binary coverage alone. The outer band
  // is dark and the next band is light, giving a continuous contrasting outline
  // without sampling noisy render brightness. Include diagonal neighbors so
  // diagonal boundaries, holes, and thin features cannot leave skipped pixels.
  static std::vector<uint8_t> Outline(const std::vector<uint8_t>& mask, uint32_t width,
                                      uint32_t height) {
    if (!width || !height || mask.size() != size_t(width) * height) {
      throw std::invalid_argument("Selection mask dimensions do not match.");
    }
    std::vector<uint8_t> outline(mask.size(), 0);
    for (uint32_t y = 0; y < height; ++y) {
      for (uint32_t x = 0; x < width; ++x) {
        const size_t index = x + size_t(width) * y;
        if (mask[index] && HasNeighbor(mask, width, height, x, y, 0)) {
          outline[index] = 1;
        }
      }
    }
    // Only look for outer-band neighbors during this pass: adding an inner-band
    // pixel must not grow the border again as the scan moves across the image.
    for (uint32_t y = 0; y < height; ++y) {
      for (uint32_t x = 0; x < width; ++x) {
        const size_t index = x + size_t(width) * y;
        if (mask[index] && !outline[index] &&
            HasNeighbor(outline, width, height, x, y, 1)) {
          outline[index] = 2;
        }
      }
    }
    return outline;
  }

  // Expand the cached border to film resolution. Every non-border pixel stays
  // transparent so selecting an object leaves its rendered interior unchanged.
  static void Compose(uint32_t width, uint32_t height,
                      const std::vector<uint8_t>& outline, uint32_t mask_width,
                      uint32_t mask_height, std::vector<uint8_t>& output) {
    if (!width || !height || !mask_width || !mask_height ||
        outline.size() != size_t(mask_width) * mask_height) {
      throw std::invalid_argument("Selection outline dimensions do not match.");
    }
    output.assign(size_t(width) * height * 4, 0);
    for (uint32_t y = 0; y < height; ++y) {
      const uint32_t my = (2 * uint64_t(y) + 1) * mask_height / (2 * uint64_t(height));
      for (uint32_t x = 0; x < width; ++x) {
        const uint32_t mx = (2 * uint64_t(x) + 1) * mask_width / (2 * uint64_t(width));
        const uint8_t band = outline[mx + size_t(mask_width) * my];
        if (!band) {
          continue;
        }
        const size_t pixel = 4 * (x + size_t(width) * y);
        const uint8_t value = band == 1 ? 10 : 245;
        output[pixel] = output[pixel + 1] = output[pixel + 2] = value;
        output[pixel + 3] = 255;
      }
    }
  }

private:
  static bool HasNeighbor(const std::vector<uint8_t>& values, uint32_t width,
                          uint32_t height, uint32_t x, uint32_t y, uint8_t target) {
    for (int dy = -1; dy <= 1; ++dy) {
      for (int dx = -1; dx <= 1; ++dx) {
        if (dx == 0 && dy == 0) {
          continue;
        }
        const int64_t nx = int64_t(x) + dx, ny = int64_t(y) + dy;
        if (nx < 0 || ny < 0 || nx >= width || ny >= height) {
          if (target == 0) {
            return true; // Outside the image is unselected.
          }
        } else if (values[size_t(nx) + size_t(width) * ny] == target) {
          return true;
        }
      }
    }
    return false;
  }
};

#endif
