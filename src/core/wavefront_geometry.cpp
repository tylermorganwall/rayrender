#include "../hitables/curve.h"
#include "../hitables/cylinder.h"
#include "../hitables/disk.h"
#include "wavefront.h"
#include <array>
#include <cmath>
#include <stdexcept>

namespace {
// Regular lathe construction, as in rayvertex's cylindrical primitives,
// generated directly into the GPU snapshot to avoid an R round trip per shape.
// 128 azimuth segments bound radial chord error by 1-cos(pi/128) (~0.03%).
constexpr int radial_segments = 128;
struct Vertex {
  point3f p;
  normal3f n;
  point2f uv;
};
template <class T> WFVector Pack(const T &v) { return {float(v[0]), float(v[1]), float(v[2]), 0}; }

class Tessellator {
public:
  std::vector<WFTriangle> triangles;
  Transform transform;
  bool flip;
  Tessellator(const hitable &shape, const Transform &placement)
      : transform(placement * *shape.ObjectToWorld), flip(shape.reverseOrientation ^ shape.transformSwapsHandedness) {}
  void Triangle(const Vertex &a, const Vertex &b, const Vertex &c, bool smooth = true) {
    auto g = convert_to_normal3(cross(b.p - a.p, c.p - a.p));
    if (g.squared_length() == 0)
      return;
    if (dot(g, a.n + b.n + c.n) < 0)
      g = -g;
    WFTriangle out;
    out.geometric = Pack(unit_vector(transform(g)) * (flip ? -1 : 1));
    out.flags = smooth ? 1 : 0;
    std::array<Vertex, 3> v{a, b, c};
    for (int k = 0; k < 3; ++k) {
      out.p[k] = Pack(transform(v[k].p));
      out.n[k] = Pack(unit_vector(transform(v[k].n)) * (flip ? -1 : 1));
    }
    out.uv01 = {float(a.uv[0]), float(a.uv[1]), float(b.uv[0]), float(b.uv[1])};
    out.uv2 = {float(c.uv[0]), float(c.uv[1]), 0, 0};
    triangles.push_back(out);
  }
  void Quad(const Vertex &a, const Vertex &b, const Vertex &c, const Vertex &d, bool smooth = true) {
    Triangle(a, b, c, smooth);
    Triangle(a, c, d, smooth);
  }
};

void Disk(Tessellator &out, Float radius, Float inner, Float y, Float normal_sign) {
  auto vertex = [&](int i, Float r) {
    double phi = 2 * M_PI * (i % radial_segments) / radial_segments;
    point3f p(r * std::cos(phi), y, r * std::sin(phi));
    return Vertex{p, normal3f(0, normal_sign, 0), {Float(.5) - p[0] / (2 * radius), Float(.5) + p[2] / (2 * radius)}};
  };
  for (int i = 0; i < radial_segments; ++i)
    out.Quad(vertex(i, inner), vertex(i, radius), vertex(i + 1, radius), vertex(i + 1, inner), false);
}
} // namespace

// Export cached native subcurves without rebuilding a tessellated tube. Geometry
// remains object-local, so nonuniform/mirrored/animated transforms act on the ray.
WFCurve ExportWavefrontCurve(const curve &shape) {
  const auto &common = *shape.common;
  WFCurve out;
  for (int i = 0; i < 4; ++i)
    out.cp[i] = Pack(shape.segmentPoints[i]);
  if (common.type == CurveType::Ribbon)
    for (int i = 0; i < 2; ++i)
      out.normal[i] = Pack(common.n[i]);
  out.normal[0].w = common.normalAngle;
  out.normal[1].w = common.invSinNormalAngle;
  out.width = {float(common.width[0]), float(common.width[1]), float(shape.uMin), float(shape.uMax)};
  out.fallback = Pack(common.cpObj[3] - common.cpObj[0]);
  out.type = uint32_t(common.type);
  return out;
}

std::vector<WFTriangle> TriangulateWavefrontPrimitive(const hitable &shape, const Transform &placement) {
  if (!shape.ObjectToWorld)
    return {};
  Tessellator out(shape, placement);
  if (auto *s = dynamic_cast<const disk *>(&shape)) {
    Disk(out, s->radius, s->inner_radius, 0, 1);
  } else if (auto *s = dynamic_cast<const cylinder *>(&shape)) {
    // Split exactly at the UV seam phi=pi as well as clipped angular endpoints.
    std::vector<double> angles{s->phi_min, s->phi_max};
    for (int i = 0; i <= radial_segments; ++i) {
      double phi = 2 * M_PI * i / radial_segments;
      if (phi > s->phi_min && phi < s->phi_max)
        angles.push_back(phi);
    }
    std::sort(angles.begin(), angles.end());
    for (size_t i = 1; i < angles.size(); ++i) {
      bool second_half = .5 * (angles[i - 1] + angles[i]) > M_PI;
      auto vertex = [&](double angle, Float v) {
        double phi = angle >= 2 * M_PI - 1e-7 ? 0 : angle;
        normal3f n(std::cos(phi), 0, std::sin(phi));
        point2f uv(1 - (angle - (second_half ? 2 * M_PI : 0) + M_PI) / (2 * M_PI), v);
        return Vertex{point3f(s->radius * n[0], (v - Float(.5)) * s->length, s->radius * n[2]), n, uv};
      };
      out.Quad(vertex(angles[i - 1], 0), vertex(angles[i], 0), vertex(angles[i], 1), vertex(angles[i - 1], 1));
    }
    // Native cylinder caps are full disks even for a clipped cylindrical side.
    if (s->has_caps) {
      Disk(out, s->radius, 0, -s->length / 2, -1);
      Disk(out, s->radius, 0, s->length / 2, 1);
    }
  }
  return std::move(out.triangles);
}
