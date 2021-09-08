#ifndef EFLOATH
#define EFLOATH

#include "mathinline.h"

// EFloat Declarations
class EFloat {
public:
  // EFloat Public Methods
  EFloat() {}
  EFloat(float v, float err = 0.f) : v(v) {
    if (err == 0.)
      low = high = v;
    else {
      // Compute conservative bounds by rounding the endpoints away
      // from the middle. Note that this will be over-conservative in
      // cases where v-err or v+err are exactly representable in
      // floating-point, but it's probably not worth the trouble of
      // checking this case.
      low = NextFloatDown(v - err);
      high = NextFloatUp(v + err);
    }
  }
  EFloat operator+(EFloat ef) const {
    EFloat r;
    r.v = v + ef.v;
    // Interval arithemetic addition, with the result rounded away from
    // the value r.v in order to be conservative.
    r.low = NextFloatDown(LowerBound() + ef.LowerBound());
    r.high = NextFloatUp(UpperBound() + ef.UpperBound());
    // r.Check();
    return r;
  }
  explicit operator float() const { return v; }
  explicit operator double() const { return v; }
  float GetAbsoluteError() const { return NextFloatUp(std::fmax(std::fabs(high - v),
                                                      std::fabs(v - low))); }
  float UpperBound() const { return high; }
  float LowerBound() const { return low; }
  EFloat operator-(EFloat ef) const {
    EFloat r;
    r.v = v - ef.v;
    r.low = NextFloatDown(LowerBound() - ef.UpperBound());
    r.high = NextFloatUp(UpperBound() - ef.LowerBound());
    // r.Check();
    return r;
  }
  EFloat operator*(EFloat ef) const {
    EFloat r;
    r.v = v * ef.v;
    Float prod[4] = {
      LowerBound() * ef.LowerBound(), UpperBound() * ef.LowerBound(),
      LowerBound() * ef.UpperBound(), UpperBound() * ef.UpperBound()};
    r.low = NextFloatDown(
      std::fmin(std::min(prod[0], prod[1]), std::fmin(prod[2], prod[3])));
    r.high = NextFloatUp(
      std::fmax(std::fmax(prod[0], prod[1]), std::fmax(prod[2], prod[3])));
    // r.Check();
    return r;
  }
  EFloat operator/(EFloat ef) const {
    EFloat r;
    r.v = v / ef.v;
    if (ef.low < 0 && ef.high > 0) {
      // Bah. The interval we're dividing by straddles zero, so just
      // return an interval of everything.
      r.low = -Infinity;
      r.high = Infinity;
    } else {
      Float div[4] = {
        LowerBound() / ef.LowerBound(), UpperBound() / ef.LowerBound(),
        LowerBound() / ef.UpperBound(), UpperBound() / ef.UpperBound()};
      r.low = NextFloatDown(
        std::fmin(std::min(div[0], div[1]), std::fmin(div[2], div[3])));
      r.high = NextFloatUp(
        std::fmax(std::max(div[0], div[1]), std::fmax(div[2], div[3])));
    }
    // r.Check();
    return r;
  }
  EFloat operator-() const {
    EFloat r;
    r.v = -v;
    r.low = -high;
    r.high = -low;
    // r.Check();
    return r;
  }
  inline bool operator==(EFloat fe) const { return v == fe.v; }
  inline void Check() const {
    if (!std::isinf(low) && !std::isnan(low) && !std::isinf(high) &&
        !std::isnan(high)) {}
        
        // CHECK_LE(low, high);
  }
  EFloat(const EFloat &ef) {
    // ef.Check();
    v = ef.v;
    low = ef.low;
    high = ef.high;
  }
  EFloat &operator=(const EFloat &ef) {
    // ef.Check();
    if (&ef != this) {
      v = ef.v;
      low = ef.low;
      high = ef.high;
    }
    return *this;
  }
  bool operator<(const Float &ef) {
    return(high < ef);
  }
  bool operator>(const Float &ef) {
    return(low > ef);
  }
  
  // friend std::ostream &operator<<(std::ostream &os, const EFloat &ef) {
  //   os << StringPrintf("v=%f (%a) - [%f, %f]",
  //                      ef.v, ef.v, ef.low, ef.high);
  //   return os;
  // }
  
private:
  // EFloat Private Data
  float v, low, high;
  friend inline EFloat sqrt(EFloat fe);
  friend inline EFloat abs(EFloat fe);
  friend inline bool Quadratic(EFloat A, EFloat B, EFloat C, EFloat *t0,
                               EFloat *t1);
};

// EFloat Inline Functions
inline EFloat operator*(float f, EFloat fe) { return EFloat(f) * fe; }

inline EFloat operator/(float f, EFloat fe) { return EFloat(f) / fe; }

inline EFloat operator+(float f, EFloat fe) { return EFloat(f) + fe; }

inline EFloat operator-(float f, EFloat fe) { return EFloat(f) - fe; }

inline EFloat sqrt(EFloat fe) {
  EFloat r;
  r.v = std::sqrt(fe.v);
  r.low = NextFloatDown(std::sqrt(fe.low));
  r.high = NextFloatUp(std::sqrt(fe.high));
  // r.Check();
  return r;
}

inline EFloat abs(EFloat fe) {
  if (fe.low >= 0)
    // The entire interval is greater than zero, so we're all set.
    return fe;
  else if (fe.high <= 0) {
    // The entire interval is less than zero.
    EFloat r;
    r.v = -fe.v;
    r.low = -fe.high;
    r.high = -fe.low;
    // r.Check();
    return r;
  } else {
    // The interval straddles zero.
    EFloat r;
    r.v = std::fabs(fe.v);
    r.low = 0;
    r.high = std::fmax(-fe.low, fe.high);
    // r.Check();
    return r;
  }
}

inline bool Quadratic(EFloat A, EFloat B, EFloat C, EFloat *t0, EFloat *t1);
inline bool Quadratic(EFloat A, EFloat B, EFloat C, EFloat *t0, EFloat *t1) {
  // Find quadratic discriminant
  Float discrim = DifferenceOfProducts(B.v, B.v, 4*A.v, C.v);
  // double discrim = (double)B.v * (double)B.v - 4. * (double)A.v * (double)C.v;
  if (discrim < 0.) return false;
  // double rootDiscrim = std::sqrt(discrim);
  Float rootDiscrim = std::sqrt(discrim);
  
  EFloat floatRootDiscrim(rootDiscrim, MachineEpsilon * rootDiscrim);
  
  // Compute quadratic _t_ values
  EFloat q;
  if ((float)B < 0) {
    q = -.5 * (B - floatRootDiscrim);
  } else {
    q = -.5 * (B + floatRootDiscrim);
  }
  *t0 = q / A;
  *t1 = C / q;
  if ((float)*t0 > (float)*t1) {
    std::swap(*t0, *t1);
  }
  return true;
}


#endif