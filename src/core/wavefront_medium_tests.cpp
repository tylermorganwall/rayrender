#if defined(NOT_CRAN) && defined(RAY_HAS_METAL_BVH)
#include "wavefront_medium_probe.h"
#include "hlbvh.h"
#include <algorithm>
#include <array>
#include <random>
#include <testthat.h>

namespace {
// Independent reference: enumerate every grid plane, sort the clipped intervals,
// and integrate their optical depths in double precision. No incremental DDA.
std::vector<std::array<double, 3>> medium_intervals(const WFMediumProbe &p,
                                                   const std::vector<float> &values) {
  double origin[3], direction[3], enter = 0, leave = INFINITY;
  const auto &m = p.medium;
  unsigned n = m.majorant.nx;
  std::vector<double> planes;
  for (unsigned a = 0; a < 3; ++a) {
    const auto &r = m.inverse[a];
    double extent = (&m.hi.x)[a] - (&m.lo.x)[a];
    origin[a] = (r.x * double(p.origin.x) + r.y * double(p.origin.y) +
                 r.z * double(p.origin.z) + r.w - (&m.lo.x)[a]) / extent;
    direction[a] = (r.x * double(p.direction.x) + r.y * double(p.direction.y) +
                    r.z * double(p.direction.z)) / extent;
    if (direction[a] == 0) {
      if (origin[a] < 0 || origin[a] > 1) return {};
    } else {
      double start = -origin[a] / direction[a], end = (1 - origin[a]) / direction[a];
      enter = std::max(enter, std::min(start, end));
      leave = std::min(leave, std::max(start, end));
      for (unsigned j = 0; j <= n; ++j)
        planes.push_back((double(j) / n - origin[a]) / direction[a]);
    }
  }
  if (!(enter < leave) || p.origin.w == 0) return {};
  planes.push_back(enter);
  planes.push_back(leave);
  std::sort(planes.begin(), planes.end());
  std::vector<std::array<double, 3>> result;
  for (size_t j = 1; j < planes.size(); ++j) {
    double start = std::max(enter, planes[j-1]), end = std::min(leave, planes[j]);
    if (!(start < end)) continue;
    unsigned cell[3];
    for (unsigned a = 0; a < 3; ++a)
      cell[a] = std::clamp(int(std::floor((origin[a] + (start+end)*.5*direction[a])*n)), 0, int(n)-1);
    double rate = p.origin.w * double(values[m.majorant.offset + cell[0] + n*(cell[1]+n*cell[2])]);
    result.push_back({start, end, rate});
  }
  return result;
}
}
context("Metal retained majorant flights") {
  test_that("successive flights invert the world-distance optical-depth CDF") {
    if (!MetalBVHAvailable()) return;
    std::vector<float> values(2*8*8*8);
    for (size_t j = 0; j < 512; ++j) {
      values[j] = 1;
      values[512+j] = j%5 == 0 ? 0 : float((j%7)+1)*.2f;
    }
    std::mt19937 rng(513);
    std::uniform_real_distribution<float> uniform(0, 1);
    std::vector<WFMediumProbe> probes;
    for (unsigned transform = 0; transform < 2; ++transform)
      for (unsigned field = 0; field < 2; ++field)
        for (unsigned ray = 0; ray < 160; ++ray) {
          WFMediumProbe p = {};
          p.medium.inverse[0] = transform ? WFVector{.5f, .125f, 0, -.1f} : WFVector{1,0,0,0};
          p.medium.inverse[1] = {0, transform ? 1.5f : 1.f, 0, 0};
          p.medium.inverse[2] = {0, 0, transform ? .75f : 1.f, 0};
          p.medium.hi = {1, 1, 1, 0};
          p.medium.majorant = {field*512, 8, 8, 8};
          p.origin = {uniform(rng), uniform(rng), uniform(rng), ray%9 == 0 ? 0.f : 20.f};
          p.direction = {2*uniform(rng)-1, 2*uniform(rng)-1, 2*uniform(rng)-1, 0};
          if (ray < 8) { // exact grid planes, axis-parallel and outward rays
            p.origin = {float(ray)/8, .5f, .5f, 20};
            p.direction = {ray%2 ? -1.f : 1.f, 0, 0, 0};
          } else if (ray < 16) { // tied two/three-axis crossings, both signs
            p.origin = {.5f, .5f, .5f, 20};
            p.direction = {ray%2 ? -1.f : 1.f, ray%2 ? -1.f : 1.f, ray%3 ? 1.f : 0.f, 0};
          } else if (ray < 24) { // enter from outside, or miss a parallel slab
            p.origin.x = -1;
            p.direction = {ray%2 ? 1.f : 0.f, 0, 0, 0};
          } else if (ray < 32) { // many very short flights in the same cell
            p.origin.w = 100000;
          }
          for (auto &u : p.uniforms) u = {uniform(rng), uniform(rng), uniform(rng), uniform(rng)};
          p.uniforms[0].x = 0;
          p.uniforms[0].y = 1e-6f;
          p.uniforms[3].w = .999999f;
          probes.push_back(p);
        }
    auto actual = ProbeMetalMedium(probes, values);
    size_t invalid = 0, wrong_rate = 0, finite = 0;
    double maximum_error = 0;
    for (size_t i = 0; i < probes.size(); ++i) {
      auto intervals = medium_intervals(probes[i], values);
      double optical = 0;
      for (unsigned j = 0; j < 16; ++j) {
        optical -= std::log1p(-double((&probes[i].uniforms[j/4].x)[j%4]));
        double residual = optical, expected = INFINITY, rate = 0;
        for (auto segment : intervals) {
          double capacity = segment[2] * (segment[1]-segment[0]);
          if (segment[2] > 0 && residual < capacity) {
            expected = segment[0]+residual/segment[2]; rate = segment[2]; break;
          }
          residual -= capacity;
        }
        double distance = (&actual[i].distances[j/4].x)[j%4];
        double majorant = (&actual[i].majorants[j/4].x)[j%4];
        invalid += std::isfinite(distance) != std::isfinite(expected) || std::isnan(distance);
        if (std::isfinite(expected) && std::isfinite(distance)) {
          ++finite;
          maximum_error = std::max(maximum_error, std::abs(distance-expected)/std::max(1., std::abs(expected)));
          wrong_rate += std::abs(majorant-rate) > 1e-5*std::max(1., rate);
        } else wrong_rate += majorant != 0;
      }
    }
    expect_true(finite > 1000);
    expect_true(invalid == 0);
    expect_true(wrong_rate == 0);
    expect_true(maximum_error < 3e-6);
  }
}
#endif
