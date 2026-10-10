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
  Float curve_radius = 0;
};
template <class T> WFVector Pack(const T &v) { return {float(v[0]), float(v[1]), float(v[2]), 0}; }

class Tessellator {
public:
  std::vector<WFTriangle> triangles;
  Transform transform;
  bool flip;
  Tessellator(const hitable &shape, const Transform &placement)
      : transform(placement * *shape.ObjectToWorld),
        flip(shape.reverseOrientation ^ shape.transformSwapsHandedness) {}
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
      out.n[k].w = v[k].curve_radius * transform(convert_to_vec3(v[k].n)).length();
    }
    out.uv01 = {float(a.uv[0]), float(a.uv[1]), float(b.uv[0]), float(b.uv[1])};
    out.uv2 = {float(c.uv[0]), float(c.uv[1]), 0, 0};
    triangles.push_back(out);
  }
  void Quad(const Vertex &a, const Vertex &b, const Vertex &c, const Vertex &d,
            bool smooth = true) {
    Triangle(a, b, c, smooth);
    Triangle(a, c, d, smooth);
  }
};

void Disk(Tessellator &out, Float radius, Float inner, Float y, Float normal_sign) {
  auto vertex = [&](int i, Float r) {
    double phi = 2 * M_PI * (i % radial_segments) / radial_segments;
    point3f p(r * std::cos(phi), y, r * std::sin(phi));
    return Vertex{p,
                  normal3f(0, normal_sign, 0),
                  {Float(.5) - p[0] / (2 * radius), Float(.5) + p[2] / (2 * radius)}};
  };
  for (int i = 0; i < radial_segments; ++i)
    out.Quad(vertex(i, inner), vertex(i, radius), vertex(i + 1, radius), vertex(i + 1, inner),
             false);
}
} // namespace

// Curve data remain encapsulated: only the cold export can read control points.
class WavefrontCurveTessellator {
public:
  static void Generate(const curve &shape, Tessellator &out) {
    const auto &c = *shape.common;
    if (c.type == CurveType::Flat)
      throw std::runtime_error("ray-facing flat curves (use cylindrical or ribbon curves)");
    const int steps = std::max(2, int(std::ceil(64 * (shape.uMax - shape.uMin))));
    const int sides = c.type == CurveType::Ribbon ? 1 : 16;
    auto tangent_at = [&](Float u) {
      Float a = 1 - u;
      vec3f t = 3 * a * a * (c.cpObj[1] - c.cpObj[0]) + 6 * a * u * (c.cpObj[2] - c.cpObj[1]) +
                3 * u * u * (c.cpObj[3] - c.cpObj[2]);
      if (t.squared_length() == 0)
        t = c.cpObj[3] - c.cpObj[0];
      if (t.squared_length() == 0)
        throw std::runtime_error("degenerate Bezier curve tangent");
      return unit_vector(t);
    };
    auto transport = [](vec3f right, const vec3f &from, const vec3f &to) {
      vec3f axis = cross(from, to);
      Float cosine = std::clamp(dot(from, to), Float(-1), Float(1));
      // Minimal rotation (Rodrigues). At a U-turn choose the previous right
      // axis for the otherwise ambiguous half turn, then re-orthogonalize.
      if (cosine > Float(-.999999))
        right += cross(axis, right) + cross(axis, cross(axis, right)) / (1 + cosine);
      right -= dot(right, to) * to;
      return unit_vector(right);
    };
    // Construct a common frame field over the WHOLE Bezier curve, so adjacent
    // PBRT subcurves meet at identical ring vertices. Changing a reference axis
    // independently at each ring would abruptly twist the triangulated tube.
    std::array<vec3f, 129> tangents, right;
    if (c.type == CurveType::Cylinder) {
      tangents[0] = tangent_at(0);
      vec3f axis = std::abs(tangents[0][0]) > .9 ? vec3f(0, 1, 0) : vec3f(1, 0, 0);
      right[0] = unit_vector(cross(axis, tangents[0]));
      for (int j = 1; j <= 128; ++j) {
        tangents[j] = tangent_at(Float(j) / 128);
        right[j] = transport(right[j - 1], tangents[j - 1], tangents[j]);
      }
    }
    auto vertex = [&](int step, int side) {
      Float u = step == 0 ? shape.uMin
                          : (step == steps
                                 ? shape.uMax
                                 : shape.uMin + (shape.uMax - shape.uMin) * Float(step) / steps);
      Float v = Float(side) / sides, a = 1 - u;
      point3f p = a * a * a * c.cpObj[0] + 3 * a * a * u * c.cpObj[1] + 3 * a * u * u * c.cpObj[2] +
                  u * u * u * c.cpObj[3];
      vec3f tangent = tangent_at(u);
      Float radius = .5 * ((1 - u) * c.width[0] + u * c.width[1]);
      normal3f n;
      if (c.type == CurveType::Ribbon) {
        n = c.normalAngle > 1e-5
                ? Float(std::sin((1 - u) * c.normalAngle) * c.invSinNormalAngle) * c.n[0] +
                      Float(std::sin(u * c.normalAngle) * c.invSinNormalAngle) * c.n[1]
                : c.n[0];
        vec3f width = cross(convert_to_vec3(n), tangent);
        if (width.squared_length() == 0)
          throw std::runtime_error("ribbon normal parallel to tangent");
        width.make_unit_vector();
        p += radius * (2 * v - 1) * width;
      } else {
        int anchor = std::min(128, int(u * 128));
        vec3f x = transport(right[anchor], tangents[anchor], tangent), up = cross(tangent, x);
        double angle = 2 * M_PI * (side % sides) / sides;
        n = convert_to_normal3(Float(std::cos(angle)) * x + Float(std::sin(angle)) * up);
        p += radius * convert_to_vec3(n);
      }
      return Vertex{p, unit_vector(n), {u, v}, radius};
    };
    for (int i = 0; i < steps; ++i)
      for (int j = 0; j < sides; ++j)
        out.Quad(vertex(i, j), vertex(i + 1, j), vertex(i + 1, j + 1), vertex(i, j + 1));
    // Hair integrates scattering THROUGH a fiber. Preserve its ray-dependent
    // impact parameter and native width-sized spawn error, so a tessellated
    // back face does not become an additional hair scattering event.
    for (auto &triangle : out.triangles) {
      triangle.flags |= c.type == CurveType::Cylinder ? 2 : 4;
      triangle.geometric.w = c.width[0];
      triangle.uv2.z = c.width[1];
      for (int axis = 0; axis < 3; ++axis)
        for (int j = 0; j < 3; ++j)
          triangle.p[axis].w += std::abs(out.transform.GetMatrix().m[axis][j]);
    }
  }
};

std::vector<WFTriangle> TriangulateWavefrontPrimitive(const hitable &shape,
                                                      const Transform &placement) {
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
        return Vertex{point3f(s->radius * n[0], (v - Float(.5)) * s->length, s->radius * n[2]), n,
                      uv};
      };
      out.Quad(vertex(angles[i - 1], 0), vertex(angles[i], 0), vertex(angles[i], 1),
               vertex(angles[i - 1], 1));
    }
    // Native cylinder caps are full disks even for a clipped cylindrical side.
    if (s->has_caps) {
      Disk(out, s->radius, 0, -s->length / 2, -1);
      Disk(out, s->radius, 0, s->length / 2, 1);
    }
  } else if (auto *s = dynamic_cast<const curve *>(&shape)) {
    WavefrontCurveTessellator::Generate(*s, out);
  }
  return std::move(out.triangles);
}
