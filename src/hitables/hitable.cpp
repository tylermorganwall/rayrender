#include "../hitables/hitable.h"
#include "../utils/raylog.h"
#include <cmath>

void SetEllipsoidDerivatives(hit_record &h, const vec3f &axes) {
  // The existing ellipsoid UV chart maps its unit normal, not its radial
  // position. Parameterize p=A^2 n / |A n| so tangents match that chart.
  vec3f q(h.p[0] / (axes[0] * axes[0]), h.p[1] / (axes[1] * axes[1]), h.p[2] / (axes[2] * axes[2]));
  const auto n = unit_vector(q);
  const Float horizontal = std::sqrt(n[0] * n[0] + n[2] * n[2]);
  const Float cosine = horizontal > 0 ? n[0] / horizontal : 1,
              sine = horizontal > 0 ? n[2] / horizontal : 0;
  const vec3f nu = 2 * Float(M_PI) * vec3f(n[2], 0, -n[0]);
  const vec3f nv = Float(M_PI) * vec3f(-n[1] * cosine, horizontal, -n[1] * sine);
  const vec3f squared(axes[0] * axes[0], axes[1] * axes[1], axes[2] * axes[2]);
  const auto an = squared * n;
  const Float denominator = std::sqrt(dot(n, an));
  h.dpdu = (squared * nu - an * (dot(an, nu) / (denominator * denominator))) / denominator;
  h.dpdv = (squared * nv - an * (dot(an, nv) / (denominator * denominator))) / denominator;
  h.dndu = convert_to_normal3(nu);
  h.dndv = convert_to_normal3(nv);
}

void hit_record::ComputeDifferentials(const Ray &ray) {
  has_differentials = false;
  dpdx = dpdy = vec3f(0);
  dudx = dvdx = dudy = dvdy = 0;
  if (!ray.has_differentials)
    return;
  const normal3f n = geometric_normal.squared_length() > 0 ? geometric_normal : normal;
  const double dx = dot(n, ray.rx_direction), dy = dot(n, ray.ry_direction);
  if (std::abs(dx) < 1e-10 || std::abs(dy) < 1e-10)
    return;
  const double tx = dot(n, p - ray.rx_origin) / dx, ty = dot(n, p - ray.ry_origin) / dy;
  if (!std::isfinite(tx) || !std::isfinite(ty))
    return;
  const auto px = ray.rx_origin + ray.rx_direction * Float(tx) - p;
  const auto py = ray.ry_origin + ray.ry_direction * Float(ty) - p;
  // Solve in the two coordinates least parallel to the geometric normal.
  const int drop = MaxDimension(Abs(convert_to_vec3(n)));
  const int a = (drop + 1) % 3, b = (drop + 2) % 3;
  const double det = double(dpdu[a]) * dpdv[b] - double(dpdu[b]) * dpdv[a];
  const double scale =
      std::max(std::abs(double(dpdu[a]) * dpdv[b]), std::abs(double(dpdu[b]) * dpdv[a]));
  if (!(std::abs(det) > 1e-12 * scale))
    return;
  auto solve = [&](const vec3f &d, Float &du, Float &dv) {
    const double u = (double(d[a]) * dpdv[b] - double(d[b]) * dpdv[a]) / det;
    const double v = (double(dpdu[a]) * d[b] - double(dpdu[b]) * d[a]) / det;
    if (!std::isfinite(u) || !std::isfinite(v) || std::abs(u) > 1e8 || std::abs(v) > 1e8)
      return false;
    du = Float(u);
    dv = Float(v);
    return true;
  };
  if (!solve(px, dudx, dvdx) || !solve(py, dudy, dvdy)) {
    dudx = dvdx = dudy = dvdy = 0;
    return;
  }
  dpdx = px;
  dpdy = py;
  has_differentials = true;
}

void PropagateRayDifferentials(const Ray &incoming, const hit_record &h, Ray &outgoing,
                               bool transmission, Float eta, bool passthrough) {
  outgoing.has_differentials = false;
  if (!incoming.has_differentials)
    return;
  if (passthrough) {
    // Keep the same two geometric lines through invisible/priority interfaces.
    // Their origins need not move with the central origin's numerical offset.
    outgoing.rx_origin = incoming.rx_origin;
    outgoing.ry_origin = incoming.ry_origin;
    outgoing.rx_direction = incoming.rx_direction;
    outgoing.ry_direction = incoming.ry_direction;
    outgoing.has_differentials = true;
    return;
  }
  if (!h.has_differentials || !(eta > 0))
    return;
  // eta is transmitted / incident IOR. Redirect adjacent rays using the
  // first-order change in the surface normal, then anchor their differences
  // to the actual sampled central direction. As in PBRT, bump second
  // derivatives are not included in these normal derivatives. Neighboring
  // transmission that becomes TIR invalidates this local approximation.
  normal3f n = h.has_bump                                       ? h.bump_normal
               : h.physical_shading_normal.squared_length() > 0 ? h.physical_shading_normal
                                                                : h.normal;
  n = unit_vector(n);
  normal3f nx = unit_vector(n + h.dndu * h.dudx + h.dndv * h.dvdx);
  normal3f ny = unit_vector(n + h.dndu * h.dudy + h.dndv * h.dvdy);
  auto redirect = [&](vec3f d, normal3f normal, vec3f &result) {
    d = unit_vector(d);
    if (dot(d, normal) > 0)
      normal = -normal;
    if (transmission)
      return refract(d, normal, 1 / eta, result);
    result = d - 2 * dot(d, normal) * convert_to_vec3(normal);
    return true;
  };
  vec3f center, x, y;
  if (!redirect(incoming.d, n, center) || !redirect(incoming.rx_direction, nx, x) ||
      !redirect(incoming.ry_direction, ny, y))
    return;
  outgoing.rx_origin = outgoing.o + h.dpdx;
  outgoing.ry_origin = outgoing.o + h.dpdy;
  outgoing.rx_direction = outgoing.d + (x - center);
  outgoing.ry_direction = outgoing.d + (y - center);
  for (int i = 0; i < 3; ++i)
    if (!std::isfinite(outgoing.rx_direction[i]) || !std::isfinite(outgoing.ry_direction[i]))
      return;
  outgoing.has_differentials = true;
}

// Translate implementation

void get_sphere_uv(const vec3f &p, Float &u, Float &v) {
  // phi returns the angle from the +x axis to the to the point (x, z),
  Float phi = atan2(p.xyz.z,p.xyz.x);
  Float theta = asin(p.xyz.y);
  u = 1 - (phi * M_1_PI + 1) * 0.5;
  v = (theta * M_1_PI + 0.5);
}

void get_sphere_uv_z(const vec3f &p, Float &u, Float &v) {
  // phi returns the angle from the +z axis to the to the point (x, z),
  Float phi = atan2(-p.xyz.x, -p.xyz.z);
  Float theta = asin(p.xyz.y);
  u = 1 - (phi * M_1_PI + 1) * 0.5; //Converts to 0..1
  v = (theta * M_1_PI + 0.5);
}

void get_sphere_uv_z(const normal3f &p, Float &u, Float &v) {
  Float phi = atan2(-p.xyz.x, -p.xyz.z);
  Float theta = asin(p.xyz.y);
  u = 1 - (phi * M_1_PI + 1) * 0.5; // Converts to 0..1
  v = (theta * M_1_PI + 0.5);
}

void get_sphere_uv(const normal3f& p, Float& u, Float& v) {
  Float phi = atan2(p.xyz.z,p.xyz.x);
  Float theta = asin(p.xyz.y);
  u = 1 - (phi * M_1_PI + 1) * 0.5;
  v = (theta * M_1_PI + 0.5);
}

Float AnimatedHitable::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) {
  Transform InterpolatedPrimToWorld;
  PrimitiveToWorld.Interpolate(time, &InterpolatedPrimToWorld);
  return(primitive->pdf_value(Inverse(InterpolatedPrimToWorld)(o),
                              Inverse(InterpolatedPrimToWorld)(v),
                              rng, time));
}
Float AnimatedHitable::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) {
  Transform InterpolatedPrimToWorld;
  PrimitiveToWorld.Interpolate(time, &InterpolatedPrimToWorld);
  return(primitive->pdf_value(Inverse(InterpolatedPrimToWorld)(o),
                              Inverse(InterpolatedPrimToWorld)(v),
                              sampler, time));
}
vec3f AnimatedHitable::random(const point3f& o, random_gen& rng, Float time) {
  Transform InterpolatedPrimToWorld;
  PrimitiveToWorld.Interpolate(time, &InterpolatedPrimToWorld);
  return(InterpolatedPrimToWorld(primitive->random(Inverse(InterpolatedPrimToWorld)(o), rng, time)));
}
vec3f AnimatedHitable::random(const point3f& o, Sampler* sampler, Float time) {
  Transform InterpolatedPrimToWorld;
  PrimitiveToWorld.Interpolate(time, &InterpolatedPrimToWorld);
  
  return(InterpolatedPrimToWorld(primitive->random(Inverse(InterpolatedPrimToWorld)(o), sampler, time)));
}
std::string AnimatedHitable::GetName() const {
  return(std::string("AnimatedHitable"));
}



bool AnimatedHitable::bounding_box(Float t0, Float t1, aabb& box) const {
  primitive->bounding_box(t0, t1, box);
  box = PrimitiveToWorld.MotionBounds(box);
  return(true);
}

const bool AnimatedHitable::hit(const Ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("Animation");
  
  Transform InterpolatedPrimToWorld;
  PrimitiveToWorld.Interpolate(r.time(), &InterpolatedPrimToWorld);
  Ray ray_interp = Inverse(InterpolatedPrimToWorld)(r);
  
  if (!primitive->hit(ray_interp, t_min, t_max, rec, rng)) {
    return false;
  }
  r.tMax = ray_interp.tMax;
  rec.light_placement = light_placements.Resolve(rec.light_placement);
  if (!InterpolatedPrimToWorld.IsIdentity()) {
    rec = InterpolatedPrimToWorld(rec);
  }
  return true;
}

const bool AnimatedHitable::hit(const Ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("Animation");
  
  Transform InterpolatedPrimToWorld;
  PrimitiveToWorld.Interpolate(r.time(), &InterpolatedPrimToWorld);
  Ray ray_interp = Inverse(InterpolatedPrimToWorld)(r);
  
  if (!primitive->hit(ray_interp, t_min, t_max, rec, sampler)) {
    return false;
  }
  r.tMax = ray_interp.tMax;
  rec.light_placement = light_placements.Resolve(rec.light_placement);
  if (!InterpolatedPrimToWorld.IsIdentity()) {
    rec = InterpolatedPrimToWorld(rec);
  }
  return true;
}

#ifdef NOT_CRAN
#include <testthat.h>
context("Ray differential transport") {
  test_that("[surface footprints solve actual UV units and reject degenerate charts]") {
    Ray ray(point3f(0, 0, -2), vec3f(0, 0, 1));
    ray.has_differentials = true;
    ray.rx_origin = ray.ry_origin = ray.o;
    ray.rx_direction = vec3f(.01, 0, 1);
    ray.ry_direction = vec3f(0, .02, 1);
    hit_record hit;
    hit.p = point3f(0);
    hit.geometric_normal = normal3f(0, 0, -1);
    hit.dpdu = vec3f(2, 0, 0);
    hit.dpdv = vec3f(0, 4, 0);
    hit.ComputeDifferentials(ray);
    expect_true(hit.has_differentials);
    expect_true(hit.dudx == Approx(.01));
    expect_true(hit.dvdy == Approx(.01));
    expect_true(hit.dvdx == 0);
    expect_true(hit.dudy == 0);
    hit.dpdv = hit.dpdu;
    hit.ComputeDifferentials(ray);
    expect_true(!hit.has_differentials);
    expect_true(hit.dudx == 0);
  }
  test_that("[reflection refraction and priority pass-through preserve adjacent rays]") {
    Ray ray(point3f(0, 0, -2), vec3f(0, 0, 1));
    ray.has_differentials = true;
    ray.rx_origin = ray.ry_origin = ray.o;
    ray.rx_direction = unit_vector(vec3f(.01, 0, 1));
    ray.ry_direction = unit_vector(vec3f(0, .02, 1));
    hit_record hit;
    hit.p = point3f(0);
    hit.normal = hit.geometric_normal = normal3f(0, 0, -1);
    hit.dpdu = vec3f(1, 0, 0);
    hit.dpdv = vec3f(0, 1, 0);
    hit.has_bump = false;
    hit.ComputeDifferentials(ray);
    Ray reflected(hit.p, vec3f(0, 0, -1));
    PropagateRayDifferentials(ray, hit, reflected, false, 1);
    expect_true(reflected.has_differentials);
    expect_true(reflected.rx_direction[0] == Approx(ray.rx_direction[0]));
    expect_true(reflected.rx_direction[2] == Approx(-ray.rx_direction[2]));
    Ray transmitted(hit.p, vec3f(0, 0, 1));
    PropagateRayDifferentials(ray, hit, transmitted, true, 1.5);
    expect_true(transmitted.has_differentials);
    expect_true(transmitted.rx_direction[0] == Approx(ray.rx_direction[0] / 1.5));
    expect_true(transmitted.rx_origin[0] == Approx(.02));
    Ray skipped(point3f(0, 0, .001), ray.d);
    PropagateRayDifferentials(ray, hit, skipped, true, 1, true);
    expect_true(skipped.has_differentials);
    expect_true((skipped.rx_origin - ray.rx_origin).squared_length() == 0);
    expect_true((skipped.rx_direction - ray.rx_direction).squared_length() == 0);
    expect_true(!Ray(hit.p, vec3f(1, 0, 1)).has_differentials);
  }
  test_that("[nonuniform transforms preserve differential ray lines]") {
    Ray ray(point3f(1, 2, 3), vec3f(.2, .3, 1));
    ray.has_differentials = true;
    ray.rx_origin = ray.o + vec3f(.01, 0, 0);
    ray.ry_origin = ray.o + vec3f(0, .01, 0);
    ray.rx_direction = ray.d + vec3f(.01, 0, 0);
    ray.ry_direction = ray.d + vec3f(0, .01, 0);
    Transform transform = Translate(vec3f(3, 4, 5)) * Scale(2, 3, 4);
    Ray moved = transform(ray);
    expect_true(moved.has_differentials);
    expect_true((moved.rx_origin - transform(ray.rx_origin)).squared_length() == 0);
    expect_true((moved.ry_direction - transform(ray.ry_direction)).squared_length() == 0);
    Ray restored = Inverse(transform)(moved);
    expect_true((restored.rx_origin - ray.rx_origin).squared_length() < 1e-10);
    ray.ScaleDifferentials(.5);
    expect_true((ray.rx_origin - ray.o)[0] == Approx(.005));
  }
}
#endif

#ifdef NOT_CRAN
#include "sphere.h"
#include "ellipsoid.h"
#include "cylinder.h"
#include "rectangle.h"
#include "instance.h"
context("Primitive differential geometry") {
  test_that("[analytic UV footprints agree with neighboring intersections]") {
    Transform identity, scaled = Translate(vec3f(.1, .2, .3)) * Scale(1.3, .7, 1.2);
    auto mat = std::make_shared<diffuse_material>(std::make_shared<constant_texture>(point3f(.5)));
    sphere ball(1, mat, nullptr, nullptr, &identity, &identity, false);
    ellipsoid oval(point3f(0), 1, vec3f(2, 1, 3), mat, nullptr, nullptr, &identity, &identity,
                   false);
    cylinder tube(1, 2, 0, 2 * M_PI, true, mat, nullptr, nullptr, &identity, &identity, false);
    xy_rect plane(-2, 2, -3, 3, 0, mat, nullptr, nullptr, &identity, &identity, false);
    instance placed(&ball, &scaled, nullptr, 0);
    random_gen rng(22);
    for (hitable *shape : {static_cast<hitable *>(&ball), static_cast<hitable *>(&oval),
                           static_cast<hitable *>(&tube), static_cast<hitable *>(&plane),
                           static_cast<hitable *>(&placed)}) {
      Ray ray(point3f(.3, .2, -5), vec3f(0, 0, 1));
      ray.has_differentials = true;
      ray.rx_origin = ray.o + vec3f(.0001, 0, 0);
      ray.ry_origin = ray.o + vec3f(0, .0001, 0);
      ray.rx_direction = ray.ry_direction = ray.d;
      hit_record center, x, y;
      expect_true(shape->hit(ray, 0, 100, center, rng));
      expect_true(shape->hit(Ray(ray.rx_origin, ray.d), 0, 100, x, rng));
      expect_true(shape->hit(Ray(ray.ry_origin, ray.d), 0, 100, y, rng));
      expect_true(center.has_differentials);
      expect_true(std::abs(center.dudx - (x.u - center.u)) < 2e-7);
      expect_true(std::abs(center.dvdx - (x.v - center.v)) < 2e-7);
      expect_true(std::abs(center.dudy - (y.u - center.u)) < 2e-7);
      expect_true(std::abs(center.dvdy - (y.v - center.v)) < 2e-7);
    }
  }
}
#endif
