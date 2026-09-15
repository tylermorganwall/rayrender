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
    Float near, far;
    if (!VolumeSphereRoots(Ray(p, wi), radius, near, far) || far <= 0)
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
  };
  std::vector<Emitter> entries;
  std::vector<double> cdf;
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
  explicit Emitters(const hitable_list &list) {
    std::vector<Placement> path;
    AddGroup(list, 1, path);
    double sum = 0;
    for (const auto &entry : entries) cdf.push_back(sum += entry.pmf);
  }
};

VolumeLightSampler::Selection VolumeLightSampler::SampleEmitter(
    const point3f &p, Sampler *sampler, random_gen &rng, Float time) const {
  Selection sample;
  double u = sampler->Get1D();
  auto selected = std::upper_bound(emitters->cdf.begin(), emitters->cdf.end(), u);
  if (selected == emitters->cdf.end()) return sample;
  sample.index = selected - emitters->cdf.begin();
  const auto &entry = emitters->entries[sample.index];
  sample.wi = entry.proposal->Sample(p, sampler, time);
  if (!(sample.wi.squared_length() > 0)) return sample;
  sample.wi = unit_vector(sample.wi);
  sample.pdf = Float(entry.pmf * entry.proposal->Pdf(p, sample.wi, rng, time));
  return sample;
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

Float VolumeLightSampler::EmitterPdf(const point3f &p, const vec3f &wi,
                                     const hit_record &hit, random_gen &rng, Float time) const {
  auto found = emitters->lookup.find({hit.shape, hit.light_placement});
  if (found == emitters->lookup.end()) return 0;
  const auto &entry = emitters->entries[found->second];
  return Float(entry.pmf * entry.proposal->Pdf(p, wi, rng, time));
}

VolumeLightSampler::VolumeLightSampler(const hitable_list &list)
    : emitters(std::make_unique<Emitters>(list)), lights(std::make_unique<Group>(list)) {}
VolumeLightSampler::~VolumeLightSampler() = default;
vec3f VolumeLightSampler::Sample(const point3f &p, Sampler *sampler, Float t) const {
  return lights->Sample(p, sampler, t);
}
Float VolumeLightSampler::Pdf(const point3f &p, const vec3f &wi, random_gen &rng, Float t) const {
  return lights->Pdf(p, wi, rng, t);
}
