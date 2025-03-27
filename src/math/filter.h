#ifndef FILTERH
#define FILTERH

#include "../math/vec2.h"

class Filter {
public:
  Filter(const vec2f &radius) : radius(radius), invRadius(vec2f(1 / radius.xy.x, 1 / radius.xy.y)) { }
  virtual Float Evaluate(const vec2f &p) const = 0;
  const vec2f radius, invRadius;
    
};

class BoxFilter : public Filter {
public:
  BoxFilter(const vec2f &radius) : Filter(radius) { }
  Float Evaluate(const vec2f &p) const;
};

class TriangleFilter : public Filter {
public:
  TriangleFilter(const vec2f &radius) : Filter(radius) { }
  Float Evaluate(const vec2f &p) const;
};

class GaussianFilter : public Filter {
public:
GaussianFilter(const vec2f &radius, Float alpha)
  : Filter(radius), alpha(alpha),
    expX(std::exp(-alpha * radius.xy.x * radius.xy.x)),
    expY(std::exp(-alpha * radius.xy.y * radius.xy.y)) { }
  Float Evaluate(const vec2f &p) const;
  
private:
  const Float alpha;
  const Float expX, expY;
  Float Gaussian(Float d, Float expv) const;
      
};


class MitchellFilter : public Filter {
  public:
  MitchellFilter(const vec2f &radius, Float B, Float C)
    : Filter(radius), B(B), C(C) {
  }
  Float Evaluate(const vec2f &p) const;
  Float Mitchell1D(Float x) const;
  
  private:
    const Float B, C;
};


class LanczosSincFilter : public Filter {
  public:
  LanczosSincFilter(const vec2f &radius, Float tau)
    : Filter(radius), tau(tau) { }
    Float Evaluate(const vec2f &p) const;
    Float Sinc(Float x) const;
    Float WindowedSinc(Float x, Float radius) const;
    
    private:
      const Float tau;
};


#endif
