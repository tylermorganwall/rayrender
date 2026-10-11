#include "../hitables/csg.h"
#include "wavefront.h"
#include <stdexcept>

namespace {
template <class T> WFVector pack(const T &p, float w = 0) { return {float(p[0]), float(p[1]), float(p[2]), w}; }
// Postfix distance expressions use a value stack. Coordinate modifiers bracket
// a child program with save/restore instructions; sibling inputs always see the
// same point. Limits are checked before upload, never by unchecked GPU indexing.
class CsgCompiler {
public:
  std::vector<WFCsgInstruction> &code;
  explicit CsgCompiler(std::vector<WFCsgInstruction> &code) : code(code) {}
  void emit(const ImplicitShape &shape, unsigned depth = 0) {
    if (depth >= 32)
      throw std::runtime_error("CSG expression exceeds the Metal nesting limit (32)");
    WFCsgInstruction n;
    auto child = [&](const std::shared_ptr<ImplicitShape> &s) { emit(*s, depth + 1); };
    auto binary = [&](const auto *s, uint32_t op, float parameter = 0) {
      child(s->shape1);
      child(s->shape2);
      n.operation = op;
      n.a.w = parameter;
    };
    if (auto s = dynamic_cast<const Union *>(&shape))
      binary(s, 30);
    else if (auto s = dynamic_cast<const Subtract *>(&shape))
      binary(s, 31);
    else if (auto s = dynamic_cast<const Intersect *>(&shape))
      binary(s, 32);
    else if (auto s = dynamic_cast<const Blend *>(&shape))
      binary(s, 33, s->op.k);
    else if (auto s = dynamic_cast<const BlendMinus *>(&shape))
      binary(s, 34, s->op.k);
    else if (auto s = dynamic_cast<const Mix *>(&shape))
      binary(s, 35, s->op.t);
    else if (auto s = dynamic_cast<const csg_list *>(&shape)) {
      if (s->shapes.empty())
        throw std::runtime_error("empty CSG group");
      child(s->shapes[0]);
      for (size_t i = 1; i < s->shapes.size(); ++i) {
        child(s->shapes[i]);
        n.operation = 30;
        code.push_back(n);
      }
      return;
    } else if (auto s = dynamic_cast<const csg_sphere *>(&shape)) {
      n.operation = 2;
      n.a = pack(s->center, s->radius);
    } else if (auto s = dynamic_cast<const csg_plane *>(&shape)) {
      n.operation = 3;
      n.a = pack(s->pointOnPlane);
      n.b = pack(s->axis.u(), s->width_x);
      n.c = pack(s->axis.v());
      n.d = pack(s->axis.w(), s->width_z);
    } else if (auto s = dynamic_cast<const csg_box *>(&shape)) {
      n.operation = 4;
      n.a = pack(s->center);
      n.b = pack(s->width);
    } else if (auto s = dynamic_cast<const csg_rounded_box *>(&shape)) {
      n.operation = 4;
      n.a = pack(s->center, s->radius);
      n.b = pack(s->width);
    } else if (auto s = dynamic_cast<const csg_torus *>(&shape)) {
      n.operation = 7;
      n.a = pack(s->center, s->ring_radius);
      n.b.w = s->cross_radius;
    } else if (auto s = dynamic_cast<const csg_capsule *>(&shape)) {
      n.operation = 8;
      n.a = pack(s->start, s->radius);
      n.b = pack(s->end);
    } else if (auto s = dynamic_cast<const csg_cylinder *>(&shape)) {
      n.operation = 9;
      n.a = pack(s->start, s->radius);
      n.b = pack(s->end, s->corner_radius);
    } else if (auto s = dynamic_cast<const csg_ellipsoid *>(&shape)) {
      n.operation = 10;
      n.a = pack(s->center);
      n.b = pack(s->axes);
    } else if (auto s = dynamic_cast<const csg_rounded_cone *>(&shape)) {
      n.operation = 11;
      n.a = pack(s->start, s->r1);
      n.b = pack(s->end, s->r2);
    } else if (auto s = dynamic_cast<const csg_cone *>(&shape)) {
      n.operation = 12;
      n.a = pack(s->start, s->radius);
      n.b = pack(s->axis.u(), s->height);
      n.c = pack(s->axis.v());
      n.d = pack(s->axis.w());
    } else if (auto s = dynamic_cast<const csg_pyramid *>(&shape)) {
      n.operation = 13;
      n.a = pack(s->center_bottom, s->h);
      n.b = pack(s->base_inv);
    } else if (auto s = dynamic_cast<const csg_triangle *>(&shape)) {
      n.operation = 14;
      n.a = pack(s->a);
      n.b = pack(s->b);
      n.c = pack(s->c);
    } else if (auto s = dynamic_cast<const csg_round *>(&shape)) {
      child(s->shape);
      n.operation = 17;
      n.a.w = s->r;
    } else if (auto s = dynamic_cast<const csg_onion *>(&shape)) {
      child(s->shape);
      n.operation = 18;
      n.a.w = s->thickness;
    } else {
      std::shared_ptr<ImplicitShape> operand;
      if (auto s = dynamic_cast<const csg_scale *>(&shape)) {
        if (!(s->scale > 0))
          throw std::runtime_error("CSG scale must be positive on Metal");
        n.operation = 119;
        n.a.w = s->scale;
        operand = s->shape;
      } else if (auto s = dynamic_cast<const csg_rotate *>(&shape)) {
        n.operation = 120;
        n.a = pack(s->pivot_point);
        n.b = pack(s->axis.u());
        n.c = pack(s->axis.v());
        n.d = pack(s->axis.w());
        operand = s->shape;
      } else if (auto s = dynamic_cast<const csg_translate *>(&shape)) {
        n.operation = 121;
        n.a = pack(s->translate);
        operand = s->shape;
      } else if (auto s = dynamic_cast<const csg_elongate *>(&shape)) {
        n.operation = 115;
        n.a = pack(s->center);
        n.b = pack(s->elongate);
        operand = s->shape;
      } else if (auto s = dynamic_cast<const csg_elongate_robust *>(&shape)) {
        n.operation = 116;
        n.a = pack(s->center);
        n.b = pack(s->elongate);
        operand = s->shape;
      } else
        throw std::runtime_error("unknown CSG distance operation");
      code.push_back(n);
      child(operand);
      n.operation += 100;
    }
    for (const auto &v : {n.a, n.b, n.c, n.d})
      if (!std::isfinite(v.x) || !std::isfinite(v.y) || !std::isfinite(v.z) || !std::isfinite(v.w))
        throw std::runtime_error("non-finite CSG parameter");
    code.push_back(n);
  }
};
} // namespace
void ExportWavefrontCSG(const ImplicitShape &shape, std::vector<WFCsgInstruction> &code) {
  CsgCompiler(code).emit(shape);
}
void WavefrontCSGBounds(const WFQuadric &q, WFVector &lo, WFVector &hi) {
  // Affine interval arithmetic encloses all eight local corners, including
  // mirrored/nonuniform scales. Marching tolerance must also fit inside the AS.
  for (int a = 0; a < 3; ++a) {
    const auto &r = q.forward[a];
    double lower = r.w, upper = r.w;
    for (int b = 0; b < 3; ++b) {
      double x = (&r.x)[b] * double((&q.csg_lo.x)[b]), y = (&r.x)[b] * double((&q.csg_hi.x)[b]);
      lower += std::min(x, y);
      upper += std::max(x, y);
    }
    double pad = 8 * std::numeric_limits<float>::epsilon() * std::max({1., std::abs(lower), std::abs(upper)});
    (&lo.x)[a] = std::nextafter(float(lower - pad), -INFINITY);
    (&hi.x)[a] = std::nextafter(float(upper + pad), INFINITY);
  }
}
