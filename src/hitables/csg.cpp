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
  rec.mat_ptr = object.mat_ptr.get();
  rec.has_bump = false;
  rec.pError = vec3f(threshold,threshold,threshold);
  rec = (*object.ObjectToWorld)(rec);
  rec.normal *= object.reverseOrientation ? -1 : 1;
  rec.bump_normal *= object.reverseOrientation ? -1 : 1;
  rec.shape = &object;
  rec.alpha_miss = false;
}

inline bool march_csg(const csg& object, const Ray& r, Float t_min, Float t_max,
                      hit_record* rec) {
  Ray r2 = (*object.WorldToObject)(r);
  const Float threshold = 0.001;
  const Float delta = 10e-5 * object.max_dist/100;
  const Float ray_length = r2.direction().length();

  if(ray_length == 0) {
    return(false);
  }

  Float t_enter;
  Float t_exit;
  if(!ray_box_interval(r2, object.object_bounds, t_min, t_max, t_enter, t_exit)) {
    return(false);
  }

  vec3f dir = r2.direction() / ray_length;
  // Keep the original march start at the object-space ray origin. The cached
  // interval is used for rejection and as a finite miss limit; hit acceptance
  // still happens through t_min/t_max after converting march distance back to a
  // ray parameter.
  Float t = 0;
  const Float max_t = t_exit * ray_length;
  bool first = true;

  while(t < max_t) {
    point3f from = r2.origin() + t * dir;

    // Use the interior distance as well so dielectric exits can be found.
    Float minDistance = ffabs(object.shapes->getDistance(from));

    // Refraction often starts just inside a surface; skip the self-hit.
    if(first && minDistance < threshold) {
      t += 0.01;
      first = false;
      continue;
    }

    first = false;
    if(minDistance <= threshold) {
      Float tval = t / ray_length;
      if(tval > t_min && tval < t_max) {
        if(rec) {
          populate_hit_record(object, r2, from, t, ray_length, threshold, delta, *rec);
        }
        return(true);
      }
      return(false);
    }

    if(!std::isfinite(minDistance)) {
      return(false);
    }
    t += minDistance;
  }

  return(false);
}
}

const bool csg::hit(const Ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("CSG");
  return(march_csg(*this, r, t_min, t_max, &rec));
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
