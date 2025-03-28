
#include "../math/filter.h"

Float BoxFilter::Evaluate(const vec2f &p) const {
  return 1.0;
}

Float TriangleFilter::Evaluate(const vec2f &p) const {
  return std::max((Float)0, radius.xy.x - std::fabs(p.xy.x)) * 
    std::max((Float)0, radius.xy.y - std::fabs(p.xy.y));
}

Float GaussianFilter::Evaluate(const vec2f &p) const {
  return Gaussian(p.xy.x, expX) * Gaussian(p.xy.y, expY);
}

Float GaussianFilter::Gaussian(Float d, Float expv) const {
  return std::max((Float)0, Float(std::exp(-alpha * d * d) - expv));
}


Float MitchellFilter::Evaluate(const vec2f &p) const {
  return Mitchell1D(p.xy.x * invRadius.xy.x) * Mitchell1D(p.xy.y * invRadius.xy.y);
}

Float MitchellFilter::Mitchell1D(Float x) const {
  x = std::fabs(2 * x);
  if (x > 1) {
    return ((-B - 6*C) * x*x*x + (6*B + 30*C) * x*x +
            (-12*B - 48*C) * x + (8*B + 24*C)) * (1.f/6.f);
  } else {
    return ((12 - 9*B - 6*C) * x*x*x + 
            (-18 + 12*B + 6*C) * x*x +
            (6 - 2*B)) * (1.f/6.f);
  }
}


Float LanczosSincFilter::Evaluate(const vec2f &p) const {
  return WindowedSinc(p.xy.x, radius.xy.x) * WindowedSinc(p.xy.y, radius.xy.y);
}

Float LanczosSincFilter::Sinc(Float x) const {
  x = std::fabs(x);
  if (x < 1e-5)  return 1;
  return std::sin(static_cast<Float>(M_PI)) / (static_cast<Float>(M_PI) * x);
}
Float LanczosSincFilter::WindowedSinc(Float x, Float radius) const {
  x = std::abs(x);
  if (x > radius) return 0;
  Float lanczos = Sinc(x / tau);
  return Sinc(x) * lanczos;
}
