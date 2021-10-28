#ifndef INTERVALH
#define INTERVALH

#include "mathinline.h"

class Interval {
public:
  // Interval Public Methods
  Interval(Float v) : low(v), high(v) {}
  Interval(Float v0, Float v1)
    : low(std::min(v0, v1)), high(std::max(v0, v1)) {}
  Interval operator+(const Interval &i) const {
    return Interval(low + i.low, high + i.high);
  }
  Interval operator-(const Interval &i) const {
    return Interval(low - i.high, high - i.low);
  }
  Interval operator*(const Interval &i) const {
    return Interval(std::min(std::min(low * i.low, high * i.low),
                             std::min(low * i.high, high * i.high)),
                             std::max(std::max(low * i.low, high * i.low),
                                      std::max(low * i.high, high * i.high)));
  }
  Float low, high;
};

inline Interval Sin(const Interval &i) {
  Float sinLow = std::sin(i.low), sinHigh = std::sin(i.high);
  if (sinLow > sinHigh) std::swap(sinLow, sinHigh);
  if (i.low < M_PI / 2 && i.high > M_PI / 2) sinHigh = 1.;
  if (i.low < (3.f / 2.f) * M_PI && i.high > (3.f / 2.f) * M_PI) sinLow = -1.;
  return Interval(sinLow, sinHigh);
}

inline Interval Cos(const Interval &i) {
  Float cosLow = std::cos(i.low), cosHigh = std::cos(i.high);
  if (cosLow > cosHigh) std::swap(cosLow, cosHigh);
  if (i.low < M_PI && i.high > M_PI) cosLow = -1.;
  return Interval(cosLow, cosHigh);
}

#endif 