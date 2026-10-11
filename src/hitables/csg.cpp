#include "../hitables/csg.h"
#include "../utils/raylog.h"
#include "../math/vectypes.h"

namespace {
// Object-space slab test against the cached CSG bounds. CSG distance fields are
// comparatively expensive to sample, so rays that never enter the bounds can
// skip marching entirely, and rays that do enter only march until the bbox exit.
inline bool ray_box_interval(const Ray& r, const aabb& box, Float t_min, Float t_max,
                             Float& t_enter, Float& t_exit) {
  t_enter = t_min;
  t_exit = t_max;

  for(int axis = 0; axis < 3; axis++) {
    const Float origin = r.o.e[axis];
    const Float direction = r.d.e[axis];
    const Float min_bound = box.min().e[axis];
    const Float max_bound = box.max().e[axis];

    if(direction == 0) {
      if(origin < min_bound || origin > max_bound) {
        return(false);
      }
      continue;
    }

    Float t0 = (min_bound - origin) * r.inv_dir_pad.e[axis];
    Float t1 = (max_bound - origin) * r.inv_dir_pad.e[axis];
    if(t0 > t1) {
      std::swap(t0,t1);
    }

    t_enter = t0 > t_enter ? t0 : t_enter;
    t_exit = t1 < t_exit ? t1 : t_exit;
    if(t_enter > t_exit) {
      return(false);
    }
  }

  return(true);
}

inline void populate_hit_record(const csg& object, const Ray& object_ray,
                                const point3f& hit_point, Float march_t,
                                Float ray_length, Float threshold,
                                Float delta, hit_record& rec) {
  rec.normal = normal3f(
    object.shapes->getDistance(hit_point + vec3f(delta, 0, 0)) -
      object.shapes->getDistance(hit_point + vec3f(-delta, 0, 0)),
    object.shapes->getDistance(hit_point + vec3f(0, delta, 0)) -
      object.shapes->getDistance(hit_point + vec3f(0, -delta, 0)),
    object.shapes->getDistance(hit_point + vec3f(0, 0, delta)) -
      object.shapes->getDistance(hit_point + vec3f(0, 0, -delta))
  );

  // Degenerate gradients can occur on exact CSG seams.
  if(rec.normal.xyz.x == 0 && rec.normal.xyz.y == 0 && rec.normal.xyz.z == 0) {
    rec.normal = convert_to_normal3(-object_ray.direction());
  }

  rec.p = hit_point;
  rec.normal.make_unit_vector();
  rec.bump_normal = rec.normal;
  rec.t = march_t / ray_length;
  rec.u = 0.5;
  rec.v = 0.5;
  rec.dpdu = 0.5;
  rec.dpdv = 0.5;
  // This implicit surface currently has no differentiable UV chart.
  rec.has_differentials = false;
  rec.dpdx = rec.dpdy = vec3f(0);
  rec.dudx = rec.dvdx = rec.dudy = rec.dvdy = 0;
  rec.dndu = rec.dndv = normal3f(0);
  rec.mat_ptr = object.mat_ptr.get();
  rec.has_bump = false;
  rec.pError = vec3f(threshold,threshold,threshold);
  rec.texture_object_p = rec.p;
  rec.texture_object_normal = rec.normal * (object.reverseOrientation ? -1 : 1);
  rec = (*object.ObjectToWorld)(rec);
  rec.normal *= object.reverseOrientation ? -1 : 1;
  rec.bump_normal *= object.reverseOrientation ? -1 : 1;
  rec.geometric_normal = rec.normal;
  rec.physical_shading_normal = normal3f(0);
  rec.base_shading_normal = normal3f(0);
  rec.shape = &object;
  rec.alpha_miss = false;
}

inline bool march_csg(const csg &object, const Ray &r, Float t_min, Float t_max, hit_record *rec) {
  Ray r2 = (*object.WorldToObject)(r);
  const Float threshold =
      std::max(Float(1e-7), Float(1e-6) * (object.object_bounds.max() - object.object_bounds.min()).length());
  const Float delta = threshold;
  const Float ray_length = r2.direction().length();

  if (ray_length == 0) {
    return (false);
  }

  Float t_enter;
  Float t_exit;
  if (!ray_box_interval(r2, Expand(object.object_bounds, 2 * threshold), t_min, t_max, t_enter, t_exit)) {
    return (false);
  }

  const vec3f dir = r2.direction() / ray_length;
  Float t = t_enter * ray_length;
  const Float max_t = t_exit * ray_length;
  Float cut_lo = INFINITY, cut_hi = -INFINITY;
  // Convex cutters give a whole interval known to contain no resulting solid.
  // This particularly avoids tiny steps parallel to a subtracted box face.
  if (const auto *difference = dynamic_cast<const Subtract *>(object.shapes.get())) {
    Ray normalized(r2.o, dir);
    if (const auto *box = dynamic_cast<const csg_box *>(difference->shape2.get())) {
      vec3f half = box->width / Float(2) - vec3f(2 * threshold);
      if (MinComponent(half) > 0) {
        aabb interior(box->center - half, box->center + half);
        if (!ray_box_interval(normalized, interior, -INFINITY, INFINITY, cut_lo, cut_hi))
          cut_lo = INFINITY;
      }
    } else if (const auto *sphere = dynamic_cast<const csg_sphere *>(difference->shape2.get())) {
      vec3f v = r2.o - convert_to_point3(sphere->center);
      Float b = dot(v, dir), radius = sphere->radius - 2 * threshold;
      vec3f closest = v - b * dir;
      Float discriminant = radius * radius - closest.squared_length();
      if (radius > 0 && discriminant > 0) {
        Float root = std::sqrt(discriminant);
        cut_lo = -b - root;
        cut_hi = -b + root;
      }
    }
  }
  // Start at the clipped interval, not the original origin. This is essential
  // for containment replay and diffusion chords that request later crossings.
  bool departing = t_enter == t_min && std::abs(object.shapes->getDistance(r2.o + t * dir)) <= threshold;
  while (t <= max_t) {
    if (t >= cut_lo && t < cut_hi) {
      t = cut_hi;
      if (t > max_t)
        return false;
    }
    point3f from = r2.origin() + t * dir;
    Float distance = std::abs(object.shapes->getDistance(from));
    if (!std::isfinite(distance))
      return false;
    if (distance > threshold)
      departing = false;
    if (distance <= threshold && !departing) {
      if (rec)
        populate_hit_record(object, r2, from, t, ray_length, 2 * threshold, delta, *rec);
      return true;
    }
    if (departing)
      distance = std::max(distance, 4 * threshold);
    Float next = t + distance;
    if (!(next > t))
      return false;
    t = next;
  }

  return (false);
}
} // namespace

const bool csg::hit(const Ray &r, Float t_min, Float t_max, hit_record &rec, random_gen &rng) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("CSG");
  return (march_csg(*this, r, t_min, t_max, &rec));
}

const bool csg::hit(const Ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("CSG");
  return(march_csg(*this, r, t_min, t_max, &rec));
}


bool csg::HitP(const Ray& r, Float t_min, Float t_max, random_gen& rng) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("CSG");
  return(march_csg(*this, r, t_min, t_max, nullptr));
}


bool csg::HitP(const Ray& r, Float t_min, Float t_max, Sampler* sampler) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("CSG");
  return(march_csg(*this, r, t_min, t_max, nullptr));
}

bool csg::bounding_box(Float t0, Float t1, aabb& box) const {
  box = world_bounds;
  return(true);
}

namespace {
// Uniform solid-angle sampling of a conservative enclosing sphere. Samples in
// holes legitimately miss; resampling them would bias the unknown projected area.
Float csg_cone_width(const csg &shape, const point3f &origin, vec3f &axis) {
  point3f center = (shape.world_bounds.min() + shape.world_bounds.max()) / Float(2);
  axis = center - origin;
  Float radius2 = (shape.world_bounds.max() - shape.world_bounds.min()).squared_length() / 4;
  Float distance2 = axis.squared_length();
  if (distance2 == 0)
    axis = vec3f(0, 0, 1);
  else
    axis /= std::sqrt(distance2);
  if (distance2 <= radius2)
    return 2;
  Float ratio = radius2 / distance2;
  return ratio / (1 + std::sqrt(1 - ratio));
}
vec3f sample_csg_light(const csg &shape, const point3f &o, vec2f u) {
  vec3f axis;
  Float dc = u[0] * csg_cone_width(shape, o, axis);
  onb basis;
  basis.build_from_w(axis);
  Float r = std::sqrt(dc * (2 - dc)), phi = 2 * M_PI * u[1];
  return basis.local(vec3f(r * std::cos(phi), r * std::sin(phi), 1 - dc));
}
Float csg_light_pdf(const csg &shape, const point3f &o, const vec3f &v) {
  vec3f axis;
  Float width = csg_cone_width(shape, o, axis);
  if (v.squared_length() == 0 || dot(unit_vector(v), axis) < 1 - width)
    return 0;
  return 1 / (2 * M_PI * width);
}
} // namespace
Float csg::pdf_value(const point3f &o, const vec3f &v, random_gen &, Float) { return csg_light_pdf(*this, o, v); }
Float csg::pdf_value(const point3f &o, const vec3f &v, Sampler *, Float) { return csg_light_pdf(*this, o, v); }
vec3f csg::random(const point3f &o, random_gen &r, Float) {
  return sample_csg_light(*this, o, vec2f(r.unif_rand(), r.unif_rand()));
}
vec3f csg::random(const point3f &o, Sampler *r, Float) { return sample_csg_light(*this, o, r->Get2D()); }
