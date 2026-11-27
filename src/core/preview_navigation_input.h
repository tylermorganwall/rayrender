#ifndef RAYRENDER_PREVIEW_NAVIGATION_INPUT_H
#define RAYRENDER_PREVIEW_NAVIGATION_INPUT_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>

// Input transitions are collected independently of rendering. A held key may
// accrue frame time only if its current press survived the entire last poll
// interval. A release followed by another press is a new tap, even when both
// events arrive together after a slow render. Each tap gets one 60-Hz step.
class PreviewNavigationInput {
public:
  using Steps = std::array<std::array<double, 6>, 2>; // unshifted/shifted, wasdqz

  void Press(char key, bool shifted) {
    int index = Index(key);
    if(index < 0) return;
    unsigned int bit = 1u << index;
    if(down & bit) return; // OS auto-repeat is not a new physical press.
    down |= bit;
    fresh |= bit;
    taps[shifted ? 1 : 0][index] += 1.0 / 60;
  }

  void Release(char key) {
    int index = Index(key);
    if(index >= 0) down &= ~(1u << index);
  }

  void Reset() {
    down = fresh = 0;
    taps = {};
  }

  Steps Consume(unsigned int physical_keys, bool shifted, double seconds) {
    Steps result = taps;
    // Physical state repairs lost releases. It must never infer a new press:
    // that would resume motion merely because a window regained focus.
    down &= physical_keys;
    double held_seconds = std::isfinite(seconds) ? std::max(0.0, std::min(seconds, .25)) : 0;
    for(int i = 0; i < 6; ++i) {
      if((down & ~fresh) & (1u << i)) result[shifted ? 1 : 0][i] += held_seconds;
    }
    taps = {};
    fresh = 0;
    return result;
  }

private:
  static int Index(char key) {
    if(key == '\0') return -1;
    const char* keys = "wasdqz";
    const char* found = std::strchr(keys, key);
    return found ? static_cast<int>(found - keys) : -1;
  }
  unsigned int down = 0, fresh = 0;
  Steps taps{};
};
#endif
