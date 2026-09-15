#include "lights.h"
#include "../hitables/box.h"
#include "../hitables/ellipsoid.h"
#include "../hitables/infinite_area_light.h"
#include "../hitables/instance.h"
#include "../hitables/mesh3d.h"
#include "../hitables/plymesh.h"
#include "../hitables/raymesh.h"
#include "../hitables/rectangle.h"
#include "../hitables/sphere.h"
#include "../hitables/triangle.h"
#include "../hitables/trimesh.h"
#include "boundary.h"
#include "intersections.h"
#include <algorithm>
#include <array>
#include <unordered_set>

namespace {
using Adapter = std::unique_ptr<VolumeLightAdapter>;
Adapter make_adapter(std::shared_ptr<hitable>);
class Group final : public VolumeLightAdapter {
  std::vector<Adapter> children;

public:
  explicit Group(const hitable_list &list) : Group(list.objects) {}
  explicit Group(std::span<const std::shared_ptr<hitable>> objects) {
    for (const auto &object : objects)
      children.push_back(make_adapter(object));
  }
  vec3f Sample(const point3f &p, Sampler *s, Float t) const override {
    if (children.empty())
      return vec3f(0);
    size_t i = std::min(children.size() - 1, size_t(s->Get1D() * children.size()));
    return children[i]->Sample(p, s, t);
  }
  Float Pdf(const point3f &p, const vec3f &wi, random_gen &rng, Float t) const override {
    double sum = 0;
    for (const auto &light : children)
      sum += light->Pdf(p, wi, rng, t);
    return children.empty() ? 0 : Float(sum / children.size());
  }
};
class ExistingLight final : public VolumeLightAdapter {
  std::shared_ptr<hitable> shape;

public:
  explicit ExistingLight(std::shared_ptr<hitable> s) : shape(std::move(s)) {}
  vec3f Sample(const point3f &p, Sampler *s, Float t) const override {
    return shape->random(p, s, t);
  }
  Float Pdf(const point3f &p, const vec3f &w, random_gen &r, Float t) const override {
    return shape->pdf_value(p, w, r, t);
  }
};
// A parallelogram or triangle sampled uniformly in surface area.
class PlanarLight final : public VolumeLightAdapter {
  point3f origin;
  vec3f u, v, n;
  double uu, uv, vv, gram, area;
  bool triangle;

public:
  PlanarLight(const point3f &p, const vec3f &u, const vec3f &v, bool tri = false)
      : origin(p), u(u), v(v), n(cross(u, v)), triangle(tri) {
    uu = dot(u, u);
    uv = dot(u, v);
    vv = dot(v, v);
    gram = uu * vv - uv * uv;
    area = n.length() * (tri ? .5 : 1);
    n = unit_vector(n);
  }
  vec3f Sample(const point3f &p, Sampler *s, Float) const override {
    vec2f z = s->Get2D();
    Float a = z.xy.x, b = z.xy.y;
    if (triangle) {
      a = std::sqrt(a);
      b *= a;
      a = 1 - a;
    }
    return origin + u * a + v * b - p;
  }
  Float Pdf(const point3f &p, const vec3f &wi, random_gen &, Float) const override {
    double cosine = dot(n, wi), distance = dot(origin - p, n) / cosine;
    if (!(distance > 0) || !(area > 0) || !std::isfinite(distance))
      return 0;
    vec3f q = p + wi * Float(distance) - origin;
    double qu = dot(q, u), qv = dot(q, v);
    double a = (qu * vv - qv * uv) / gram, b = (qv * uu - qu * uv) / gram;
    if (a < 0 || b < 0 || a > 1 || b > 1 || (triangle && a + b > 1))
      return 0;
    return Float(distance * distance * wi.squared_length() * wi.length() /
                 (std::abs(cosine) * area));
  }
};
class SphereLight final : public VolumeLightAdapter {
  Float radius;
  static double SquaredDistance(const point3f &p) {
    double result = 0;
    for (int i = 0; i < 3; ++i)
      result += double(p[i]) * p[i];
    return result;
  }

public:
  explicit SphereLight(Float radius) : radius(std::abs(radius)) {}
  vec3f Sample(const point3f &p, Sampler *s, Float) const override {
    vec2f u = s->Get2D();
    double d2 = SquaredDistance(p);
    if (d2 <= double(radius) * radius) {
      Float z = 1 - 2 * u.xy.x, r = std::sqrt(std::max(Float(0), 1 - z * z)),
            phi = 2 * M_PI * u.xy.y;
      return vec3f(r * std::cos(phi), r * std::sin(phi), z);
    }
    double sin2 = double(radius) * radius / d2;
    double one_minus_cos = sin2 / (1 + std::sqrt(std::max(0.0, 1 - sin2)));
    double cos = 1 - u.xy.x * one_minus_cos;
    double sin = std::sqrt(std::max(0.0, (1 - cos) * (1 + cos)));
    onb frame;
    frame.build_from_w(-convert_to_vec3(p));
    return frame.local(sin * std::cos(2 * M_PI * u.xy.y), sin * std::sin(2 * M_PI * u.xy.y), cos);
  }
  Float Pdf(const point3f &p, const vec3f &wi, random_gen &, Float) const override {
    double d2 = SquaredDistance(p);
    if (d2 <= double(radius) * radius)
      return 1 / (4 * M_PI);
    double sin2 = double(radius) * radius / d2;
    double one_minus_cos = sin2 / (1 + std::sqrt(std::max(0.0, 1 - sin2)));
    // A Float cosine test loses precision at the edge of a distant light's
    // narrow cone. Use the same roots as the actual NEE sphere intersection;
    // otherwise valid light hits can incorrectly receive a zero sampling PDF.
    Float t_near, t_far;
    if (!VolumeSphereRoots(Ray(p, wi), radius, t_near, t_far) || t_far <= 0)
      return 0;
    return Float(1 / (2 * M_PI * one_minus_cos));
  }
};
class TransformedLight final : public VolumeLightAdapter {
  Adapter child;
  Transform fixed;
  std::shared_ptr<AnimatedHitable> animated;
  Transform At(Float t) const {
    Transform a = fixed;
    if (animated)
      animated->PrimitiveToWorld.Interpolate(t, &a);
    return a;
  }

public:
  TransformedLight(Adapter child, const Transform &matrix)
      : child(std::move(child)), fixed(matrix) {}
  TransformedLight(Adapter child, std::shared_ptr<AnimatedHitable> a)
      : child(std::move(child)), animated(std::move(a)) {}
  vec3f Sample(const point3f &p, Sampler *s, Float time) const override {
    Transform a = At(time);
    return a(child->Sample(Inverse(a)(p), s, time));
  }
  Float Pdf(const point3f &p, const vec3f &wi, random_gen &rng, Float time) const override {
    Transform a = At(time), inv = Inverse(a);
    // Preserve the ray direction used for intersection. Normalizing again can
    // round it across the silhouette of a small emitter.
    vec3f local = inv(wi);
    vec3f x = a(vec3f(1, 0, 0)), y = a(vec3f(0, 1, 0)), z = a(vec3f(0, 0, 1));
    double determinant = std::abs(dot(x, cross(y, z)));
    double jacobian = determinant * std::pow(double(local.length()) / wi.length(), 3);
    return jacobian > 0 ? Float(child->Pdf(inv(p), local, rng, time) / jacobian) : 0;
  }
};
Adapter make_adapter(std::shared_ptr<hitable> shape) {
  if (auto b = std::dynamic_pointer_cast<MediumBoundary>(shape))
    return make_adapter(b->geometry);
  if (auto a = std::dynamic_pointer_cast<AnimatedHitable>(shape))
    return std::make_unique<TransformedLight>(make_adapter(a->primitive), a);
  if (auto *i = dynamic_cast<instance *>(shape.get()))
    return std::make_unique<TransformedLight>(
        std::make_unique<Group>(*i->importance_sampled_objects), *i->ObjectToWorld);
  if (auto *b = dynamic_cast<box *>(shape.get()))
    return std::make_unique<Group>(b->list);
  if (auto *m = dynamic_cast<trimesh *>(shape.get()))
    return std::make_unique<Group>(m->tri_mesh_bvh->Primitives());
  if (auto *m = dynamic_cast<mesh3d *>(shape.get()))
    return std::make_unique<Group>(m->mesh_bvh->Primitives());
  if (auto *m = dynamic_cast<plymesh *>(shape.get()))
    return std::make_unique<Group>(m->ply_mesh_bvh->Primitives());
  if (auto *m = dynamic_cast<raymesh *>(shape.get()))
    return std::make_unique<Group>(m->tri_mesh_bvh->Primitives());
  if (auto *t = dynamic_cast<triangle *>(shape.get())) {
    const point3f &a = t->mesh->p[t->v[0]], &b = t->mesh->p[t->v[1]], &c = t->mesh->p[t->v[2]];
    return std::make_unique<PlanarLight>(a, b - a, c - a, true);
  }
  if (auto *r = dynamic_cast<xy_rect *>(shape.get())) {
    const Transform &a = *r->ObjectToWorld;
    return std::make_unique<PlanarLight>(a(point3f(r->x0, r->y0, r->k)),
                                         a(vec3f(r->x1 - r->x0, 0, 0)),
                                         a(vec3f(0, r->y1 - r->y0, 0)));
  }
  if (auto *r = dynamic_cast<xz_rect *>(shape.get())) {
    const Transform &a = *r->ObjectToWorld;
    return std::make_unique<PlanarLight>(a(point3f(r->x0, r->k, r->z0)),
                                         a(vec3f(r->x1 - r->x0, 0, 0)),
                                         a(vec3f(0, 0, r->z1 - r->z0)));
  }
  if (auto *r = dynamic_cast<yz_rect *>(shape.get())) {
    const Transform &a = *r->ObjectToWorld;
    return std::make_unique<PlanarLight>(a(point3f(r->k, r->y0, r->z0)),
                                         a(vec3f(0, r->y1 - r->y0, 0)),
                                         a(vec3f(0, 0, r->z1 - r->z0)));
  }
  if (auto *s = dynamic_cast<sphere *>(shape.get()))
    return std::make_unique<TransformedLight>(std::make_unique<SphereLight>(s->radius),
                                              *s->ObjectToWorld);
  if (auto *s = dynamic_cast<ellipsoid *>(shape.get()))
    return std::make_unique<TransformedLight>(std::make_unique<SphereLight>(1),
                                              *s->ObjectToWorld *
                                                  Scale(s->axes[0], s->axes[1], s->axes[2]));
  if (dynamic_cast<InfiniteAreaLight *>(shape.get()))
    return std::make_unique<ExistingLight>(std::move(shape));
  // For remaining finite emitters use a conservative enclosing sphere. Sampling
  // and PDF then agree even when the legacy shape sampler ignores affine scale.
  // Tracing the sampled direction resolves the actual first emitter endpoint.
  aabb bounds;
  if (shape->bounding_box(0, 1, bounds)) {
    point3f center = (bounds.min() + bounds.max()) / 2;
    Float radius = (bounds.max() - center).length();
    if (radius > 0 && std::isfinite(radius))
      return std::make_unique<TransformedLight>(std::make_unique<SphereLight>(radius),
                                                Translate(convert_to_vec3(center)));
  }
  return std::make_unique<ExistingLight>(std::move(shape));
}


// These bounds estimate importance, not radiance or visibility. The light tree
// always mixes in the original distribution, so a dark texture probe, loose
// motion bound or imperfect directional estimate can only affect variance.
struct LightBounds {
  aabb box;
  point3f center{0};
  vec3f axis{0, 1, 0};
  double radius2 = 0, power = 0;
  double cos_orientation = -1, sin_orientation = 0, cos_emission = 0;

  void SetBox(const aabb &b) {
    box = b;
    center = b.min() * Float(.5) + b.max() * Float(.5);
    radius2 = 0;
    for (int i = 0; i < 3; ++i) {
      double radius = .5 * (double(b.max()[i]) - b.min()[i]);
      radius2 += radius * radius;
    }
  }

  static LightBounds Union(const LightBounds &a, const LightBounds &b) {
    LightBounds result;
    result.SetBox(surrounding_box(a.box, b.box));
    result.power = a.power + b.power;
    result.cos_emission = std::min(a.cos_emission, b.cos_emission);
    // Retain one axis and enlarge its cone to enclose the other cone. This
    // conservative merge also handles opposing faces and full-sphere emitters.
    const LightBounds &wide = a.cos_orientation < b.cos_orientation ? a : b;
    const LightBounds &other = &wide == &a ? b : a;
    result.axis = wide.axis;
    double separation = std::acos(std::clamp(double(dot(wide.axis, other.axis)), -1.0, 1.0));
    double angle = std::max(std::acos(wide.cos_orientation),
                           separation + std::acos(other.cos_orientation));
    if (angle < M_PI) {
      result.cos_orientation = std::cos(angle);
      result.sin_orientation = std::sin(angle);
    }
    return result;
  }

  // Maximum cosine after allowing an angular deviation. Evaluating this form
  // needs no inverse trigonometry during stochastic tree traversal.
  static double MaxCos(double cosine, double cos_angle, double sin_angle) {
    if (cosine >= cos_angle) return 1;
    return cosine * cos_angle + std::sqrt(std::max(0.0, 1 - cosine * cosine)) * sin_angle;
  }

  double Importance(const VolumeLightSampler::Context &ctx) const {
    double v[3], distance2 = 0;
    for (int i = 0; i < 3; ++i) {
      v[i] = double(ctx.p[i]) - center[i];
      distance2 += v[i] * v[i];
    }
    if (!(distance2 > radius2)) return power / std::max(radius2, 1e-30);
    double inv_distance = 1 / std::sqrt(distance2);
    double sin_angle = std::sqrt(radius2 / distance2);
    double cos_angle = std::sqrt(std::max(0.0, 1 - sin_angle * sin_angle));
    double emission = 1;
    if (cos_orientation > -1) {
      double cosine = 0;
      for (int i = 0; i < 3; ++i) cosine += axis[i] * v[i] * inv_distance;
      cosine = MaxCos(std::clamp(cosine, -1.0, 1.0), cos_orientation, sin_orientation);
      emission = MaxCos(std::clamp(cosine, -1.0, 1.0), cos_angle, sin_angle);
      if (emission <= cos_emission) return 0;
    }
    double receiver = 1;
    if (ctx.n.squared_length() > 0) {
      double cosine = 0;
      for (int i = 0; i < 3; ++i) cosine += ctx.n[i] * v[i] * inv_distance;
      receiver = MaxCos(std::min(1.0, std::abs(cosine)), cos_angle, sin_angle);
    }
    return power * std::max(0.0, emission) * receiver / distance2;
  }
};

bool finite_bounds(const aabb &box) {
  for (int i = 0; i < 3; ++i)
    if (!std::isfinite(box.min()[i]) || !std::isfinite(box.max()[i]) ||
        box.min()[i] > box.max()[i]) return false;
  return true;
}

// A few deterministic probes are sufficient for a power estimate. In
// particular, these do not replace the emitter's actual texture evaluation.
double emission_estimate(const material *mat, const point3f &p) {
  const texture *emit = nullptr;
  double intensity = 1;
  if (auto *light = dynamic_cast<const diffuse_light *>(mat)) {
    emit = light->emit.get(); intensity = light->intensity;
  } else if (auto *light = dynamic_cast<const spot_light *>(mat)) {
    emit = light->emit.get(); intensity = light->intensity;
  }
  if (!emit) return 1;
  double sum = 0;
  for (int y = 0; y < 4; ++y) for (int x = 0; x < 4; ++x) {
    point3f rgb = emit->value((x + .5f) / 4, (y + .5f) / 4, p);
    for (int c = 0; c < 3; ++c)
      if (std::isfinite(rgb[c])) sum += std::max(0.0, double(rgb[c]));
  }
  return std::max(0.0, intensity) * sum / 48;
}
} // namespace
// The directional mixture above remains available for atmospheric source terms.
// Finite emission instead lives on (emitter, direction): a selected lamp cannot
// collect another lamp's emission when that lamp blocks the connection.
struct VolumeLightSampler::Emitters {
  struct Key {
    const hitable *shape;
    uint64_t placement;
    bool operator==(const Key &other) const {
      return shape == other.shape && placement == other.placement;
    }
  };
  struct Hash {
    size_t operator()(const Key &key) const {
      return std::hash<const hitable *>{}(key.shape) ^
             (std::hash<uint64_t>{}(key.placement) + 0x9e3779b9);
    }
  };
  struct Placement {
    std::shared_ptr<hitable> owner;
    LightPlacementMap *ids;
    Transform At(Float time) const {
      if (auto *a = dynamic_cast<AnimatedHitable *>(owner.get())) {
        Transform result;
        a->PrimitiveToWorld.Interpolate(time, &result);
        return result;
      }
      return *owner->ObjectToWorld;
    }
  };
  struct Emitter {
    std::shared_ptr<hitable> shape;
    Adapter proposal;
    std::vector<Placement> placements;
    Key key;
    double pmf;
    size_t node = size_t(-1);
  };
  struct Node {
    LightBounds bounds;
    double fixed_mass = 0;
    size_t parent = size_t(-1), left = size_t(-1), right = size_t(-1), emitter = size_t(-1);
  };
  std::vector<Emitter> entries;
  std::vector<double> cdf;
  std::vector<Node> nodes;
  std::vector<size_t> unbounded;
  std::vector<double> unbounded_cdf;
  double finite_mass = 0;
  static constexpr double fallback_probability = .05;
  std::unordered_map<Key, size_t, Hash> lookup;
  std::unordered_set<LightPlacementMap *> initialized;
  uint64_t next_placement = 1;

  void AddGroup(const hitable_list &group, double pmf, std::vector<Placement> &path) {
    AddGroup(group.objects, pmf, path);
  }
  void AddGroup(std::span<const std::shared_ptr<hitable>> objects, double pmf,
                std::vector<Placement> &path) {
    if (objects.empty()) return;
    for (const auto &shape : objects) Add(shape, pmf / objects.size(), path);
  }
  void Add(std::shared_ptr<hitable> shape, double pmf, std::vector<Placement> &path) {
    if (auto b = std::dynamic_pointer_cast<MediumBoundary>(shape)) {
      Add(b->geometry, pmf, path);
      return;
    }
    auto *animated = dynamic_cast<AnimatedHitable *>(shape.get());
    auto *placed = dynamic_cast<instance *>(shape.get());
    if (animated || placed) {
      auto *ids = animated ? &animated->light_placements : &placed->light_placements;
      if (initialized.insert(ids).second) ids->Reset();
      path.push_back({shape, ids});
      if (animated) Add(animated->primitive, pmf, path);
      else if (placed->importance_sampled_objects)
        AddGroup(*placed->importance_sampled_objects, pmf, path);
      path.pop_back();
      return;
    }
    if (auto *b = dynamic_cast<box *>(shape.get())) { AddGroup(b->list, pmf, path); return; }
    if (auto *m = dynamic_cast<trimesh *>(shape.get())) { AddGroup(m->tri_mesh_bvh->Primitives(), pmf, path); return; }
    if (auto *m = dynamic_cast<mesh3d *>(shape.get())) { AddGroup(m->mesh_bvh->Primitives(), pmf, path); return; }
    if (auto *m = dynamic_cast<plymesh *>(shape.get())) { AddGroup(m->ply_mesh_bvh->Primitives(), pmf, path); return; }
    if (auto *m = dynamic_cast<raymesh *>(shape.get())) { AddGroup(m->tri_mesh_bvh->Primitives(), pmf, path); return; }

    uint64_t placement = 0;
    Adapter proposal = make_adapter(shape);
    for (auto i = path.rbegin(); i != path.rend(); ++i) {
      placement = i->ids->Register(placement, next_placement);
      if (auto a = std::dynamic_pointer_cast<AnimatedHitable>(i->owner))
        proposal = std::make_unique<TransformedLight>(std::move(proposal), a);
      else proposal = std::make_unique<TransformedLight>(std::move(proposal), i->At(0));
    }
    Key key{shape.get(), placement};
    auto found = lookup.find(key);
    if (found != lookup.end()) {
      // Repeated references to the same physical emitter are one MIS outcome.
      entries[found->second].pmf += pmf;
      return;
    }
    lookup.emplace(key, entries.size());
    entries.push_back({std::move(shape), std::move(proposal), path, key, pmf});
  }


  bool Bounds(const Emitter &entry, LightBounds &bounds) const {
    // Distant lights have no useful spatial attenuation bound. Preserve their
    // original probability, including the balance between sky, sun and moon.
    if (dynamic_cast<InfiniteAreaLight *>(entry.shape.get())) return false;
    aabb box;
    if (!entry.shape->bounding_box(0, 1, box) || !finite_bounds(box)) return false;
    Transform placement;
    bool animated = false;
    for (auto p = entry.placements.rbegin(); p != entry.placements.rend(); ++p) {
      if (auto *a = dynamic_cast<AnimatedHitable *>(p->owner.get())) {
        box = a->PrimitiveToWorld.MotionBounds(box);
        animated = true;
      } else box = p->At(0)(box);
      placement = p->At(0) * placement;
    }
    if (!finite_bounds(box)) return false;
    bounds.SetBox(box);

    const material *mat = entry.shape->mat_ptr.get();
    vec3f u(0), v(0), normal(0);
    double area = 0;
    if (auto *t = dynamic_cast<triangle *>(entry.shape.get())) {
      const point3f &a = t->mesh->p[t->v[0]], &b = t->mesh->p[t->v[1]], &c = t->mesh->p[t->v[2]];
      u = placement(b - a); v = placement(c - a);
      area = .5 * cross(u, v).length();
      mat = t->mesh->mesh_materials[t->mesh->face_material_id[t->face_number]].get();
      // Vertex normals, alpha masks and consistent-normal interpolation can
      // change which side emits. Use an unrestricted cone in those cases.
      if (!t->mesh->has_normals && !t->mesh->alpha_textures[t->mesh->face_material_id[t->face_number]])
        normal = convert_to_vec3(placement(convert_to_normal3(cross(b - a, c - a))));
    } else if (dynamic_cast<xy_rect *>(entry.shape.get()) ||
               dynamic_cast<xz_rect *>(entry.shape.get()) ||
               dynamic_cast<yz_rect *>(entry.shape.get()) ||
               dynamic_cast<sphere *>(entry.shape.get()) ||
               dynamic_cast<ellipsoid *>(entry.shape.get())) {
      Transform a = placement * *entry.shape->ObjectToWorld;
      if (auto *r = dynamic_cast<xy_rect *>(entry.shape.get())) {
        u = a(vec3f(r->x1 - r->x0, 0, 0)); v = a(vec3f(0, r->y1 - r->y0, 0));
        if (!r->alpha_mask) normal = convert_to_vec3(a(normal3f(0, 0, 1)));
      } else if (auto *r = dynamic_cast<xz_rect *>(entry.shape.get())) {
        u = a(vec3f(r->x1 - r->x0, 0, 0)); v = a(vec3f(0, 0, r->z1 - r->z0));
        if (!r->alpha_mask) normal = convert_to_vec3(a(normal3f(0, 1, 0)));
      } else if (auto *r = dynamic_cast<yz_rect *>(entry.shape.get())) {
        u = a(vec3f(0, r->y1 - r->y0, 0)); v = a(vec3f(0, 0, r->z1 - r->z0));
        if (!r->alpha_mask) normal = convert_to_vec3(a(normal3f(1, 0, 0)));
      }
      // Masked rectangles face their shading normal toward the incoming ray
      // and can emit on both sides, just like masked triangles above.
      area = cross(u, v).length();
      if (entry.shape->reverseOrientation) normal = -normal;
      vec3f axes(0);
      if (auto *s = dynamic_cast<sphere *>(entry.shape.get())) axes = vec3f(std::abs(s->radius));
      else if (auto *s = dynamic_cast<ellipsoid *>(entry.shape.get())) axes = s->axes;
      if (axes.squared_length() > 0) {
        vec3f x = a(vec3f(axes[0], 0, 0)), y = a(vec3f(0, axes[1], 0)), z = a(vec3f(0, 0, axes[2]));
        area = (4 * M_PI / 3) * (cross(x, y).length() + cross(x, z).length() + cross(y, z).length());
      }
    }
    if (!(area > 0) || !std::isfinite(area)) area = 4 * M_PI * bounds.radius2;
    bounds.power = area * emission_estimate(mat, bounds.center);
    if (!std::isfinite(bounds.power)) bounds.power = 0;
    if (!animated && normal.squared_length() > 0) {
      bounds.axis = unit_vector(normal);
      bounds.cos_orientation = 1;
    }
    // Spotlight directions live in the same coordinates as emitted(), rather
    // than in the primitive's object space. Placement transforms are applied
    // later to the hit record, so they do not rotate the material's direction.
    if (auto *spot = dynamic_cast<const spot_light *>(mat)) {
      bounds.axis = spot->spot_direction;
      bounds.cos_orientation = 1;
      bounds.cos_emission = std::max(0.0, double(spot->cosTotalWidth));
    }
    return true;
  }

  size_t Build(std::vector<size_t> &order, size_t begin, size_t end,
               const std::vector<LightBounds> &bounds, size_t parent = size_t(-1)) {
    size_t index = nodes.size();
    nodes.emplace_back();
    nodes[index].parent = parent;
    if (end - begin == 1) {
      size_t emitter = order[begin];
      nodes[index].emitter = emitter;
      nodes[index].bounds = bounds[emitter];
      nodes[index].fixed_mass = entries[emitter].pmf;
      entries[emitter].node = index;
      return index;
    }
    aabb centers;
    for (size_t i = begin; i < end; ++i) centers = surrounding_box(centers, bounds[order[i]].center);
    int axis = centers.MaxDimension();
    size_t middle = begin + (end - begin) / 2;
    std::nth_element(order.begin() + begin, order.begin() + middle, order.begin() + end,
      [&](size_t a, size_t b) {
        Float ca = bounds[a].center[axis], cb = bounds[b].center[axis];
        return ca == cb ? a < b : ca < cb;
      });
    size_t left = Build(order, begin, middle, bounds, index);
    size_t right = Build(order, middle, end, bounds, index);
    nodes[index].left = left; nodes[index].right = right;
    nodes[index].bounds = LightBounds::Union(nodes[left].bounds, nodes[right].bounds);
    nodes[index].fixed_mass = nodes[left].fixed_mass + nodes[right].fixed_mass;
    return index;
  }

  double LeftProbability(size_t index, const Context &ctx) const {
    const Node &node = nodes[index], &left = nodes[node.left], &right = nodes[node.right];
    double a = left.bounds.Importance(ctx), b = right.bounds.Importance(ctx);
    if (a + b > 0 && std::isfinite(a + b)) return a / (a + b);
    return left.fixed_mass / node.fixed_mass;
  }

  double TreePmf(const Context &ctx, size_t emitter) const {
    double pmf = finite_mass;
    for (size_t child = entries[emitter].node; child != 0;) {
      size_t parent = nodes[child].parent;
      double p = LeftProbability(parent, ctx);
      pmf *= nodes[parent].left == child ? p : 1 - p;
      child = parent;
    }
    return pmf;
  }

  double Pmf(const Context &ctx, size_t emitter) const {
    const Emitter &entry = entries[emitter];
    if (nodes.empty() || entry.node == size_t(-1)) return entry.pmf;
    return fallback_probability * entry.pmf + (1 - fallback_probability) * TreePmf(ctx, emitter);
  }

  size_t Select(const Context &ctx, double u, double &pmf) const {
    if (nodes.empty() || u < fallback_probability) {
      if (!nodes.empty()) u /= fallback_probability;
      size_t index = std::upper_bound(cdf.begin(), cdf.end(), u) - cdf.begin();
      pmf = index < entries.size() ? Pmf(ctx, index) : 0;
      return index;
    }
    u = (u - fallback_probability) / (1 - fallback_probability);
    if (u >= finite_mass) {
      size_t i = std::upper_bound(unbounded_cdf.begin(), unbounded_cdf.end(), u) - unbounded_cdf.begin();
      pmf = i < unbounded.size() ? entries[unbounded[i]].pmf : 0;
      return i < unbounded.size() ? unbounded[i] : entries.size();
    }
    u /= finite_mass;
    size_t index = 0;
    // Median splits bound the depth by the number of bits in size_t. Reuse
    // these branch probabilities instead of evaluating importance a second
    // time, multiplying in the same leaf-to-root order as emitter-hit MIS.
    std::array<double, std::numeric_limits<size_t>::digits> probabilities;
    size_t depth = 0;
    while (nodes[index].emitter == size_t(-1)) {
      double p = LeftProbability(index, ctx);
      if (u < p) {
        probabilities[depth++] = p;
        u /= p; index = nodes[index].left;
      } else {
        probabilities[depth++] = 1 - p;
        u = p < 1 ? (u - p) / (1 - p) : 0; index = nodes[index].right;
      }
      u = std::min(u, std::nextafter(1.0, 0.0));
    }
    pmf = finite_mass;
    while (depth) pmf *= probabilities[--depth];
    size_t emitter = nodes[index].emitter;
    pmf = fallback_probability * entries[emitter].pmf + (1 - fallback_probability) * pmf;
    return emitter;
  }


  explicit Emitters(const hitable_list &list, SelectionMethod method) {
    std::vector<Placement> path;
    AddGroup(list, 1, path);
    double sum = 0;
    for (const auto &entry : entries) cdf.push_back(sum += entry.pmf);
    // Empty sampled subgroups retain their original null probability. They
    // must not silently increase the probabilities of the remaining emitters.
    if (method == SelectionMethod::Fixed) return;
    std::vector<LightBounds> bounds(entries.size());
    std::vector<size_t> order;
    for (size_t i = 0; i < entries.size(); ++i) {
      if (Bounds(entries[i], bounds[i])) {
        order.push_back(i);
        finite_mass += entries[i].pmf;
      } else unbounded.push_back(i);
    }
    if (order.size() < 2) return;
    nodes.reserve(2 * order.size() - 1);
    Build(order, 0, order.size(), bounds);
    sum = finite_mass;
    for (size_t i : unbounded) unbounded_cdf.push_back(sum += entries[i].pmf);
  }
};

VolumeLightSampler::Selection VolumeLightSampler::SampleEmitter(
    const Context &ctx, Sampler *sampler, random_gen &rng) const {
  Selection sample;
  double selection_pmf = 0;
  sample.index = emitters->Select(ctx, sampler->Get1D(), selection_pmf);
  if (sample.index >= emitters->entries.size()) return sample;
  const auto &entry = emitters->entries[sample.index];
  sample.wi = entry.proposal->Sample(ctx.p, sampler, ctx.time);
  if (!(sample.wi.squared_length() > 0)) return sample;
  sample.wi = unit_vector(sample.wi);
  sample.pdf = Float(selection_pmf * entry.proposal->Pdf(ctx.p, sample.wi, rng, ctx.time));
  return sample;
}

double VolumeLightSampler::SelectionPmf(const Context &ctx, size_t index) const {
  return index < emitters->entries.size() ? emitters->Pmf(ctx, index) : 0;
}

bool VolumeLightSampler::Endpoint(const Selection &sample, const Ray &ray,
                                  hit_record &endpoint, random_gen &rng) const {
  const auto &entry = emitters->entries[sample.index];
  Ray local = ray;
  for (const auto &placement : entry.placements) local = Inverse(placement.At(ray.time()))(local);
  if (!entry.shape->hit(local, 0, MaxT, endpoint, rng)) return false;
  for (auto p = entry.placements.rbegin(); p != entry.placements.rend(); ++p)
    endpoint = p->At(ray.time())(endpoint);
  endpoint.light_placement = entry.key.placement;
  return true;
}

bool VolumeLightSampler::Matches(const Selection &sample, const hit_record &hit) const {
  const auto &key = emitters->entries[sample.index].key;
  return key.shape == hit.shape && key.placement == hit.light_placement;
}

Float VolumeLightSampler::EmitterPdf(const Context &ctx, const vec3f &wi,
                                     const hit_record &hit, random_gen &rng) const {
  auto found = emitters->lookup.find({hit.shape, hit.light_placement});
  if (found == emitters->lookup.end()) return 0;
  const auto &entry = emitters->entries[found->second];
  return Float(emitters->Pmf(ctx, found->second) * entry.proposal->Pdf(ctx.p, wi, rng, ctx.time));
}

VolumeLightSampler::VolumeLightSampler(const hitable_list &list, SelectionMethod method)
    : emitters(std::make_unique<Emitters>(list, method)), lights(std::make_unique<Group>(list)) {}
VolumeLightSampler::~VolumeLightSampler() = default;
vec3f VolumeLightSampler::Sample(const point3f &p, Sampler *sampler, Float t) const {
  return lights->Sample(p, sampler, t);
}
Float VolumeLightSampler::Pdf(const point3f &p, const vec3f &wi, random_gen &rng, Float t) const {
  return lights->Pdf(p, wi, rng, t);
}
