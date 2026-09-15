#include "haze.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace {
double next_knot(const DensityIndexRay &ray, double t, double end) {
  for (int a = 0; a < 3; ++a) {
    double d = ray.direction[a];
    if (d == 0) continue;
    double x = std::fma(t, d, ray.origin[a]);
    double plane = d > 0 ? std::floor(x) + 1 : std::ceil(x) - 1;
    double next = (plane - ray.origin[a]) / d;
    if (next <= t) next = (plane + (d > 0 ? 1 : -1) - ray.origin[a]) / d;
    if (next > t) end = std::min(end, next);
  }
  return end;
}
}

MediumHazeIterator::MediumHazeIterator(const Medium *m, const Ray &local, double end,
                                       bool enabled, const std::atomic<bool> *cancel)
    : medium(m), cancel(cancel), end(end), enabled(enabled) {
  if (!enabled || !m) return;
  if (!m->haze) { this->enabled = false; return; }
  if (m->haze_density_threshold == 0) return;
  if (m->density_scale == 0) return;
  threshold = std::max(std::numeric_limits<double>::denorm_min(),
                        m->haze_density_threshold / m->density_scale);
  if (m->IsHomogeneous()) { this->enabled = 1 < threshold; return; }
  uniform = false;
  index_ray = m->DensityRay(local);
  majorants.emplace(m->SampleRay(local, end));
}


// Along one interpolation cell, density is a cubic polynomial in normalized
// ray distance u. Corner values construct that polynomial directly, avoiding
// noisy point probes and missed thin intervals between null events.
void MediumHazeIterator::SplitCell(double a, double b) {
  std::array<double, 3> cell, start, delta;
  for (int axis = 0; axis < 3; ++axis) {
    cell[axis] = std::floor(std::fma(a + (b - a) / 2, index_ray.direction[axis], index_ray.origin[axis]));
    start[axis] = std::fma(a, index_ray.direction[axis], index_ray.origin[axis]) - cell[axis];
    delta[axis] = (b - a) * index_ray.direction[axis];
  }
  auto values = medium->DensityCorners(cell);
  auto bounds = std::minmax_element(values.begin(), values.end());
  piece = 0;
  if (*bounds.second < threshold || *bounds.first >= threshold) {
    pieces[0] = {a, b, *bounds.second < threshold};
    piece_count = 1;
    return;
  }
  std::array<double, 4> polynomial{};
  for (int corner = 0; corner < 8; ++corner) {
    std::array<double, 4> term{values[corner], 0, 0, 0};
    for (int axis = 0; axis < 3; ++axis) {
      bool upper = (corner >> axis) & 1;
      double c = upper ? start[axis] : 1 - start[axis];
      double d = upper ? delta[axis] : -delta[axis];
      for (int degree = axis + 1; degree >= 0; --degree)
        term[degree] = term[degree] * c + (degree ? term[degree - 1] * d : 0);
    }
    for (int degree = 0; degree < 4; ++degree) polynomial[degree] += term[degree];
  }
  polynomial[0] -= threshold;
  double magnitude = 0;
  for (double c : polynomial) magnitude = std::max(magnitude, std::abs(c));
  if (magnitude > 0) for (double &c : polynomial) c /= magnitude;
  auto value = [&](double u) {
    return ((polynomial[3] * u + polynomial[2]) * u + polynomial[1]) * u + polynomial[0];
  };


  // Split at derivative roots first: each resulting interval is monotone and
  // contains at most one crossing. This also handles three crossings in one
  // voxel and tangent contacts without assuming density is monotone along rays.
  // Unused slots remain at the upper endpoint, after every interior cut. Sorting
  // the fixed-size array also lets compilers verify the sort's access bounds.
  std::array<double, 4> cuts{0, 1, 1, 1};
  size_t count = 2;
  auto add_cut = [&](double u) { if (u > 0 && u < 1) cuts[count++] = u; };
  double qa = 3 * polynomial[3], qb = 2 * polynomial[2], qc = polynomial[1];
  if (std::abs(qa) < 1e-14) {
    if (std::abs(qb) >= 1e-14) add_cut(-qc / qb);
  } else {
    double discriminant = qb * qb - 4 * qa * qc;
    if (discriminant >= 0) {
      double q = -.5 * (qb + std::copysign(std::sqrt(discriminant), qb));
      if (q != 0) { add_cut(q / qa); add_cut(qc / q); }
      else add_cut(-qb / (2 * qa));
    }
  }
  std::sort(cuts.begin(), cuts.end());
  std::array<double, 5> edges{0};
  size_t edge_count = 1;
  auto add_edge = [&](double u) {
    if (u > edges[edge_count - 1] + 1e-12 && u < 1 - 1e-12) edges[edge_count++] = u;
  };
  for (size_t i = 1; i < count; ++i) {
    double lo = cuts[i - 1], hi = cuts[i], left = value(lo), right = value(hi);
    if (left == 0) add_edge(lo);
    if ((left < 0 && right > 0) || (left > 0 && right < 0)) {
      for (int step = 0; step < 42; ++step) {
        double mid = lo + (hi - lo) / 2;
        if ((value(mid) < 0) == (left < 0)) lo = mid;
        else hi = mid;
      }
      add_edge(lo + (hi - lo) / 2);
    }
    if (right == 0) add_edge(cuts[i]);
  }
  edges[edge_count++] = 1;
  piece_count = edge_count - 1;
  for (size_t i = 0; i < piece_count; ++i)
    pieces[i] = {a + (b - a) * edges[i], a + (b - a) * edges[i + 1],
                 value((edges[i] + edges[i + 1]) / 2) < 0};
}


std::optional<MediumHazeSegment> MediumHazeIterator::Piece() {
  if (cancel && cancel->load(std::memory_order_relaxed)) return {};
  if (piece < piece_count) return pieces[piece++];
  if (!(cursor < end)) return {};
  if (uniform) { double start = cursor; cursor = end; return MediumHazeSegment{start, end, enabled}; }
  if (!majorant || cursor >= majorant->t_max) majorant = majorants->Next();
  if (!majorant || cursor < majorant->t_min) {
    double start = cursor;
    cursor = majorant ? majorant->t_min : end;
    return MediumHazeSegment{start, cursor, true};
  }
  double start = cursor;
  if (majorant->density_max < threshold) {
    cursor = majorant->t_max;
    return MediumHazeSegment{start, cursor, true};
  }
  cursor = next_knot(index_ray, cursor, majorant->t_max);
  SplitCell(start, cursor);
  return pieces[piece++];
}

std::optional<MediumHazeSegment> MediumHazeIterator::Next(bool merge) {
  auto result = lookahead ? lookahead : Piece();
  lookahead.reset();
  if (!result) return {};
  if (!merge) return result;
  while (auto next = Piece()) {
    if (result->integrate != next->integrate) { lookahead = next; break; }
    result->t_max = next->t_max;
  }
  return result;
}
