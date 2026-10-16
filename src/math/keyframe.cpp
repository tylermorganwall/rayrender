#include "keyframe.h"
#include "quaternion.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace animation {
namespace {
double Dot(const Value &a, const Value &b) {
  double result = 0;
  for (int k = 0; k < 4; ++k) result += a[k] * b[k];
  return result;
}
Value Unit(Value q) {
  double norm = std::sqrt(Dot(q, q));
  if (!std::isfinite(norm) || norm < 1e-12)
    throw std::runtime_error("Animation quaternion has a zero or unsupported near-zero norm.");
  for (auto &x : q) x /= norm;
  return q;
}
double Cubic(const std::array<double, 4> &c, double u) {
  return ((c[0] * u + c[1]) * u + c[2]) * u + c[3];
}
// A cubic lies in its Bezier control hull. Recursively prove that at least one
// component stays separated from zero on every subinterval. This also rejects
// curves too close to a singularity to certify, before render state is exposed.
bool Separated(const std::array<Value, 4> &b, int depth = 0, double margin = 1e-10) {
  for (int k = 0; k < 4; ++k) {
    double lo = b[0][k], hi = lo;
    for (int i = 1; i < 4; ++i) { lo = std::min(lo, b[i][k]); hi = std::max(hi, b[i][k]); }
    if (lo > margin || hi < -margin) return true;
  }
  if (depth == 24) return false;
  auto midpoint = [](Value a, const Value &c) {
    for (int k = 0; k < 4; ++k) a[k] = (a[k] + c[k]) * .5;
    return a;
  };
  Value a = midpoint(b[0], b[1]), c = midpoint(b[1], b[2]), d = midpoint(b[2], b[3]);
  Value e = midpoint(a, c), f = midpoint(c, d), g = midpoint(e, f);
  return Separated({b[0], a, e, g}, depth + 1, margin) && Separated({g, f, d, b[3]}, depth + 1, margin);
}
} // namespace

ClipClock::ClipClock(double begin, double end, double start, double offset, double speed, Loop loop)
  : begin(begin), end(end), start(start), offset(offset), speed(speed), loop(loop) {
  if (!std::isfinite(begin) || !std::isfinite(end) || end < begin ||
      !std::isfinite(2 * (end - begin)) || !std::isfinite(start) ||
      !std::isfinite(offset) || !std::isfinite(speed))
    throw std::runtime_error("Clip clocks require finite values and a finite ordered interval.");
  if (loop != Loop::Hold && loop != Loop::Repeat && loop != Loop::PingPong)
    throw std::runtime_error("Unknown clip loop mode.");
}

double ClipClock::Phase(double scene_time) const {
  // Multiplication and addition deliberately match the packed evaluator (no
  // fused operation). Do not silently overflow before applying the remainder.
  double phase = speed == 0 ? offset : speed * (scene_time - start) + offset;
  if (!std::isfinite(scene_time) || !std::isfinite(phase))
    throw std::runtime_error("Clip clock exceeds finite time precision.");
  return phase;
}

double ClipClock::At(double scene_time) const {
  double duration = end - begin;
  if (duration == 0) return begin;
  double phase = Phase(scene_time);
  if (loop == Loop::Hold) return begin + std::clamp(phase, 0., duration);
  double period = loop == Loop::Repeat ? duration : 2 * duration;
  // fmod avoids cancellation in x-D*floor(x/D) for many completed cycles.
  double local = std::fmod(phase, period);
  if (local < 0) local += period;
  if (local >= period) local = 0;
  if (loop == Loop::PingPong && local > duration) local = period - local;
  return begin + local;
}

std::vector<std::pair<double, double>> ClipClock::Intervals(double first, double last) const {
  if (!std::isfinite(first) || !std::isfinite(last) || first > last)
    throw std::runtime_error("Clock bounds require finite ordered scene times.");
  double duration = end - begin;
  if (duration == 0 || speed == 0 || first == last) return {{At(first), At(first)}};
  double lo = Phase(first), hi = Phase(last);
  if (lo > hi) std::swap(lo, hi);
  if (loop == Loop::Hold)
    return {{begin + std::clamp(lo, 0., duration), begin + std::clamp(hi, 0., duration)}};
  double period = loop == Loop::Repeat ? duration : 2 * duration;
  // A full cycle has the same bound regardless of cycle count. Never iterate
  // through thousands of seams for a long exposure or a very fast player.
  if (hi - lo >= period) return {{begin, end}};
  double phase = std::fmod(lo, period);
  if (phase < 0) phase += period;
  if (phase >= period) phase = 0;
  double finish = phase + (hi - lo);
  if (loop == Loop::Repeat) {
    if (finish < duration) return {{begin + phase, begin + finish}};
    return {{begin + phase, end}, {begin, begin + (finish - duration)}};
  }
  // Split at each reflection, at most three pieces for less than one cycle.
  std::vector<std::pair<double, double>> result;
  auto reflect = [&](double x) {
    x = std::fmod(x, period);
    return begin + (x <= duration ? x : period - x);
  };
  while (phase < finish) {
    double next = std::min(finish, (std::floor(phase / duration) + 1) * duration);
    double a = reflect(phase), b = reflect(next);
    result.emplace_back(std::min(a, b), std::max(a, b));
    phase = next;
  }
  if (result.empty()) result.emplace_back(At(first), At(first));
  return result;
}

Track::Track(Property property, Interpolation interpolation, std::vector<double> times,
             std::vector<Value> values, std::vector<Value> incoming, std::vector<Value> outgoing)
  : property(property), interpolation(interpolation), times(std::move(times)), values(std::move(values)),
    incoming(std::move(incoming)), outgoing(std::move(outgoing)) {
  const size_t n = this->times.size();
  if (!n || this->values.size() != n) throw std::runtime_error("Tracks require one value per key and at least one key.");
  for (size_t i = 0; i < n; ++i) {
    if (!std::isfinite(this->times[i]) || (i && (!(this->times[i] > this->times[i - 1]) ||
        !std::isfinite(this->times[i] - this->times[i - 1]))))
      throw std::runtime_error("Track times must be finite and strictly increasing.");
    for (double x : this->values[i])
      if (!std::isfinite(x)) throw std::runtime_error("Track values must be finite.");
  }
  if (interpolation == Interpolation::Smooth && property != Property::Rotation)
    throw std::runtime_error("Smooth interpolation requires a rotation track; use cubic for vectors.");
  const bool cubic = interpolation == Interpolation::Cubic || interpolation == Interpolation::Smooth;
  const bool explicit_tangents = !this->incoming.empty() || !this->outgoing.empty();
  if (explicit_tangents && (!cubic || interpolation == Interpolation::Smooth ||
      this->incoming.size() != n || this->outgoing.size() != n))
    throw std::runtime_error("Explicit cubic tangents require incoming and outgoing values at every key.");
  if (property == Property::Rotation) {
    if (interpolation == Interpolation::Cubic && !explicit_tangents)
      throw std::runtime_error("Cubic rotation requires explicit tangents; use smooth for native automatic tangents.");
    for (size_t i = 0; i < n; ++i) {
      double norm = std::sqrt(Dot(this->values[i], this->values[i]));
      if (interpolation == Interpolation::Cubic) {
        // Authored cubic signs and tangent magnitudes are part of the curve.
        if (!std::isfinite(norm) || std::abs(norm - 1) > 1e-5)
          throw std::runtime_error("Explicit cubic quaternion keys must be unit length.");
      } else {
        this->values[i] = Unit(this->values[i]);
        if (i && Dot(this->values[i - 1], this->values[i]) < 0)
          for (auto &x : this->values[i]) x = -x;
      }
    }
  }
  if (cubic && !explicit_tangents) {
    this->incoming.resize(n); this->outgoing.resize(n);
    // Nonuniform centered differences: each neighboring secant is weighted by
    // the opposite interval. Shared derivatives make the Hermite curve C1.
    for (size_t i = 0; i < n; ++i) {
      Value derivative{};
      if (n > 1) {
        size_t left = i ? i - 1 : i, right = i + 1 < n ? i + 1 : i;
        double dl = this->times[i] - this->times[left], dr = this->times[right] - this->times[i];
        for (int k = 0; k < Dimensions(); ++k) {
          double vl = dl ? (this->values[i][k] - this->values[left][k]) / dl : 0;
          double vr = dr ? (this->values[right][k] - this->values[i][k]) / dr : 0;
          derivative[k] = !dl ? vr : !dr ? vl : (dr * vl + dl * vr) / (dl + dr);
        }
        if (property == Property::Rotation) {
          double radial = Dot(derivative, this->values[i]);
          for (int k = 0; k < 4; ++k) derivative[k] -= radial * this->values[i][k];
        }
      }
      this->incoming[i] = this->outgoing[i] = derivative;
    }
  }
  for (const auto *tangents : {&this->incoming, &this->outgoing})
    for (const auto &tangent : *tangents)
      for (double x : tangent)
        if (!std::isfinite(x)) throw std::runtime_error("Track tangents must be finite.");
  if (cubic) for (size_t i = 0; i + 1 < n; ++i) {
    std::array<Value, 4> bezier{};
    for (int k = 0; k < Dimensions(); ++k) {
      auto c = Polynomial(i, k);
      for (double x : c) if (!std::isfinite(x)) throw std::runtime_error("Cubic track coefficients overflow.");
      bezier[0][k] = c[3]; bezier[1][k] = c[3] + c[2] / 3;
      bezier[2][k] = c[3] + 2 * c[2] / 3 + c[1] / 3;
      bezier[3][k] = this->values[i + 1][k];
    }
    if (property == Property::Rotation && !Separated(bezier))
      throw std::runtime_error("Cubic quaternion curve approaches a singularity; add keys or change tangents.");
  }
  if (property == Property::Scale) {
    // STEP can change handedness discontinuously, but continuous segments may
    // not cross zero. Analytic extrema also catch cubic interior crossings.
    for (const auto &v : this->values) for (int k = 0; k < 3; ++k)
      if (std::abs(v[k]) < 1e-8) throw std::runtime_error("Scale trajectory is singular or too close to zero.");
    if (interpolation != Interpolation::Step) for (size_t i = 0; i + 1 < n; ++i) {
      auto range = Range(this->times[i], this->times[i + 1]);
      for (int k = 0; k < 3; ++k)
        if (range[0][k] <= 1e-8 && range[1][k] >= -1e-8)
          throw std::runtime_error("Scale trajectory crosses or approaches zero between keys.");
    }
  }
  if (property != Property::Rotation) {
    auto range = Range(this->times.front(), this->times.back());
    for (int k = 0; k < Dimensions(); ++k)
      if (!std::isfinite(range[0][k]) || !std::isfinite(range[1][k]) ||
          std::abs(range[0][k]) > 1e12 || std::abs(range[1][k]) > 1e12)
        throw std::runtime_error("Track trajectory exceeds the supported magnitude of 1e12.");
  }
}

std::array<double, 4> Track::Polynomial(size_t i, int k) const {
  double a = values[i][k], b = values[i + 1][k], dt = times[i + 1] - times[i];
  if (interpolation == Interpolation::Step) return {0, 0, 0, a};
  if (interpolation == Interpolation::Linear) return {0, 0, b - a, a};
  double m0 = dt * outgoing[i][k], m1 = dt * incoming[i + 1][k];
  return {2 * a - 2 * b + m0 + m1, -3 * a + 3 * b - 2 * m0 - m1, m0, a};
}

void Track::ValidateFloatPacking() const {
  if (property != Property::Scale && property != Property::Rotation) return;
  for (size_t i = 0; i + 1 < times.size(); ++i) {
    if (interpolation == Interpolation::Step) continue;
    if (property == Property::Rotation && interpolation == Interpolation::Linear) continue;
    std::array<Value, 4> bezier{};
    double magnitude = 0;
    for (int k = 0; k < Dimensions(); ++k) {
      auto c = Polynomial(i, k);
      for (double x : c) magnitude = std::max(magnitude, std::abs(x));
      bezier[0][k] = c[3]; bezier[1][k] = c[3] + c[2] / 3;
      bezier[2][k] = c[3] + 2 * c[2] / 3 + c[1] / 3;
      bezier[3][k] = values[i + 1][k];
    }
    // Leave room for both float polynomial arithmetic and relative key-time
    // rounding (export allows at most 1e-5 of an adjacent interval).
    double margin = std::max(1e-12, magnitude * 1e-4);
    if (property == Property::Rotation && (magnitude > 1e12 || !Separated(bezier, 0, margin)))
      throw std::runtime_error("Device quaternion curve cannot be certified away from a float singularity; use CPU.");
    if (property == Property::Scale) {
      auto range = Range(times[i], times[i + 1]);
      for (int k = 0; k < 3; ++k)
        if (range[0][k] <= margin && range[1][k] >= -margin)
          throw std::runtime_error("Device scale curve approaches float singularity; use CPU.");
    }
  }
}

Value Track::At(double time) const {
  if (!std::isfinite(time)) throw std::runtime_error("Evaluation time must be finite.");
  if (time <= times.front()) return property == Property::Rotation ? Unit(values.front()) : values.front();
  if (time >= times.back()) return property == Property::Rotation ? Unit(values.back()) : values.back();
  size_t i = std::upper_bound(times.begin(), times.end(), time) - times.begin() - 1;
  if (interpolation == Interpolation::Step || time == times[i])
    return property == Property::Rotation ? Unit(values[i]) : values[i];
  double u = (time - times[i]) / (times[i + 1] - times[i]);
  Value result{};
  if (property == Property::Rotation && interpolation == Interpolation::Linear) {
    double cosine = std::clamp(Dot(values[i], values[i + 1]), -1., 1.);
    double a = 1 - u, b = u;
    if (cosine < .9995) {
      double theta = std::acos(cosine), sine = std::sin(theta);
      a = std::sin((1 - u) * theta) / sine; b = std::sin(u * theta) / sine;
    }
    for (int k = 0; k < 4; ++k) result[k] = a * values[i][k] + b * values[i + 1][k];
  } else {
    for (int k = 0; k < Dimensions(); ++k) result[k] = Cubic(Polynomial(i, k), u);
  }
  return property == Property::Rotation ? Unit(result) : result;
}

std::array<Value, 2> Track::Range(double begin, double end) const {
  if (!std::isfinite(begin) || !std::isfinite(end) || begin > end)
    throw std::runtime_error("Bounds require finite ordered times.");
  if (property == Property::Rotation)
    throw std::runtime_error("Quaternion bounds use a rotation-invariant sphere.");
  Value first = At(begin), last = At(end);
  std::array<Value, 2> range{};
  for (int k = 0; k < Dimensions(); ++k) {
    range[0][k] = std::min(first[k], last[k]); range[1][k] = std::max(first[k], last[k]);
  }
  for (size_t i = 0; i + 1 < times.size(); ++i) {
    if (times[i] > end || times[i + 1] < begin) continue;
    double dt = times[i + 1] - times[i];
    double lo = std::clamp((begin - times[i]) / dt, 0., 1.);
    double hi = std::clamp((end - times[i]) / dt, 0., 1.);
    for (int k = 0; k < Dimensions(); ++k) {
      auto c = Polynomial(i, k);
      auto include = [&](double u) {
        if (u < lo || u > hi) return;
        double v = Cubic(c, u);
        range[0][k] = std::min(range[0][k], v); range[1][k] = std::max(range[1][k], v);
      };
      include(lo); include(hi);
      double coefficient_scale = std::max({std::abs(c[0]), std::abs(c[1]), std::abs(c[2]), 1.});
      double a = 3 * (c[0] / coefficient_scale), b = 2 * (c[1] / coefficient_scale), d = c[2] / coefficient_scale;
      if (a == 0) { if (b != 0) include(-d / b); }
      else {
        double disc = b * b - 4 * a * d;
        if (disc >= 0) {
          double q = -.5 * (b + std::copysign(std::sqrt(disc), b));
          include(q / a);
          if (q != 0) include(d / q);
        }
      }
    }
  }
  return range;
}

Transform TransformTrack::At(double time) const {
  return AtLocal(clock ? clock->At(time) : time);
}

Transform TransformTrack::AtLocal(double time) const {
  Value p = position ? position->At(time) : Value{}, s = scale ? scale->At(time) : Value{1, 1, 1, 0};
  Value q = rotation ? rotation->At(time) : Value{0, 0, 0, 1};
  Quaternion orientation; orientation.v = vec3f(q[0], q[1], q[2]); orientation.w = q[3];
  return Translate(vec3f(p[0] + pivot[0], p[1] + pivot[1], p[2] + pivot[2])) *
    orientation.ToTransform() * Scale(s[0], s[1], s[2]) * Translate(vec3f(-pivot[0], -pivot[1], -pivot[2]));
}

aabb TransformTrack::Bounds(const aabb &box, double begin, double end) const {
  if (!std::isfinite(begin) || !std::isfinite(end) || begin > end)
    throw std::runtime_error("Bounds require finite ordered times.");
  if (begin == end) return At(begin)(box);
  if (!clock) return LocalBounds(box, begin, end);
  aabb result;
  for (const auto &interval : clock->Intervals(begin, end))
    result = surrounding_box(result, LocalBounds(box, interval.first, interval.second));
  return result;
}

aabb TransformTrack::LocalBounds(const aabb &box, double begin, double end) const {
  if (begin == end) return AtLocal(begin)(box);
  std::array<Value, 2> translation{}, scaling{Value{1, 1, 1, 0}, Value{1, 1, 1, 0}};
  if (position) translation = position->Range(begin, end);
  if (scale) scaling = scale->Range(begin, end);
  // Every scaled point is in a sphere around the pivot with this radius.
  // Rotation preserves that sphere, independently of asynchronous key times.
  // Translation's analytic component range then encloses its entire sweep.
  double radius2 = 0;
  for (int k = 0; k < 3; ++k) {
    double extent = std::max(std::abs(box.min()[k] - pivot[k]), std::abs(box.max()[k] - pivot[k]));
    double scale_max = std::max(std::abs(scaling[0][k]), std::abs(scaling[1][k]));
    radius2 += extent * extent * scale_max * scale_max;
  }
  double radius = std::sqrt(radius2);
  point3f lo, hi;
  for (int k = 0; k < 3; ++k) {
    double lower = translation[0][k] + pivot[k] - radius, upper = translation[1][k] + pivot[k] + radius;
    double pad = 64 * std::numeric_limits<Float>::epsilon() * std::max({1., std::abs(lower), std::abs(upper)});
    lo[k] = lower - pad; hi[k] = upper + pad;
  }
  return aabb(lo, hi);
}
} // namespace animation
