#ifdef NOT_CRAN
#include "subsurface.h"
#include "boundary.h"
#include "../math/rng.h"
#include "../materials/material.h"
#include "../hitables/sphere.h"
#include "../hitables/rectangle.h"
#include "../core/bvh.h"
#include "intersections.h"
#include <testthat.h>
#include <cmath>

namespace {
class OrderedTriangleProbe : public hitable {
public:
  OrderedTriangleProbe(point3f a, point3f b, point3f c) : a(a), b(b), c(c) {}
  const bool hit(const Ray &r, Float lo, Float hi, hit_record &h, random_gen &) const override {
    return Intersect(r, lo, hi, h);
  }
  const bool hit(const Ray &r, Float lo, Float hi, hit_record &h, Sampler *) const override {
    return Intersect(r, lo, hi, h);
  }
  bool bounding_box(Float, Float, aabb &bounds) const override {
    bounds = surrounding_box(aabb(a, b), c);
    return true;
  }
  std::string GetName() const override { return "Ordered triangle probe"; }
  size_t GetSize() override { return sizeof(*this); }
  void hitable_info_bounds(Float, Float) const override {}
private:
  bool Intersect(const Ray &r, Float lo, Float hi, hit_record &h) const {
    Float b0, b1, b2;
    if (!VolumeTriangleIntersection(r, a, b, c, lo, hi, h.t, b0, b1, b2, &h.precise_t)) return false;
    h.p = r.o + r.d * h.t;
    h.geometric_normal = unit_vector(convert_to_normal3(cross(b - a, c - a)));
    h.normal = h.bump_normal = h.geometric_normal;
    h.pError = h.dpdu = h.dpdv = vec3f(0);
    h.u = h.v = 0;
    h.mat_ptr = nullptr;
    h.shape = this;
    return true;
  }
  point3f a, b, c;
};
}

context("Subsurface proposal and boundary") {
  test_that("volume triangle coverage includes edges and rejects points just outside") {
    struct Probe { Float x, y; bool inside; };
    const Float step = std::numeric_limits<Float>::epsilon();
    const Probe probes[] = {
      {0, 0, true}, {-Float(0), Float(0), true}, {1, 0, true}, {0, 1, true},
      {.5f, 0, true}, {0, .5f, true}, {.5f, .5f, true}, {.25f, .25f, true},
      {-step, .5f, false}, {.5f, -step, false}, {.5f, .5f + step, false},
      {1, 1, false}, {-1, -1, false}
    };
    // Rotate the plane through all three dominant ray axes and test either
    // winding from both sides. Expected coverage comes from x>=0,y>=0,x+y<=1.
    for (int axis = 0; axis < 3; ++axis) {
      auto point = [axis](Float x, Float y, Float z) {
        point3f p(0);
        p[axis] = z;
        p[(axis + 1) % 3] = x;
        p[(axis + 2) % 3] = y;
        return p;
      };
      for (bool reverse : {false, true}) {
        point3f a = point(0, 0, 0), b = point(1, 0, 0), c = point(0, 1, 0);
        if (reverse) std::swap(b, c);
        for (Float side : {-1.f, 1.f}) {
          vec3f direction(0);
          direction[axis] = -side;
          for (const auto &probe : probes) {
            Ray ray(point(probe.x, probe.y, side), direction);
            Float t = -1, b0 = -1, b1 = -1, b2 = -1;
            double precise = -1;
            bool hit = VolumeTriangleIntersection(ray, a, b, c, 1, 1, t, b0, b1, b2, &precise);
            expect_true(hit == probe.inside);
            if (hit) {
              expect_true(t == 1);
              expect_true(precise == 1);
              expect_true(b0 == Float(1) - probe.x - probe.y);
              expect_true(b1 == (reverse ? probe.y : probe.x));
              expect_true(b2 == (reverse ? probe.x : probe.y));
              // A double lower bound can exclude the hit even when its Float
              // representation rounds back to exactly the same distance.
              ray.medium_t_min = std::nextafter(1.0, INFINITY);
              expect_false(VolumeTriangleIntersection(ray, a, b, c, 1, 1, t, b0, b1, b2));
            }
          }
        }
      }
    }
  }

  test_that("volume triangle origin contacts and degenerate planes retain their rules") {
    const point3f a(0, 0, 0), b(1, 0, 0), c(0, 1, 0);
    Float t, b0, b1, b2;
    double precise;
    for (Float zero : {Float(0), -Float(0)}) {
      for (Float direction : {-1.f, 1.f}) {
        Ray ray(point3f(.25f, .25f, zero), vec3f(0, 0, direction));
        expect_true(VolumeTriangleIntersection(ray, a, b, c, 0, 0, t, b0, b1, b2, &precise));
        expect_true(t == 0);
        expect_true(precise == 0);
        ray.medium_t_min = std::nextafter(0.0, INFINITY);
        expect_false(VolumeTriangleIntersection(ray, a, b, c, 0, 1, t, b0, b1, b2));
      }
    }
    Ray normal(point3f(.25f, .25f, 1), vec3f(0, 0, -1));
    expect_false(VolumeTriangleIntersection(normal, a, a, c, 0, 2, t, b0, b1, b2));
    expect_false(VolumeTriangleIntersection(normal, a, b, point3f(2, 0, 0), 0, 2, t, b0, b1, b2));
    Ray parallel(point3f(.25f, .25f, 1), vec3f(1, 0, 0));
    expect_false(VolumeTriangleIntersection(parallel, a, b, c, 0, 2, t, b0, b1, b2));
  }

  test_that("bounded flights retain all surfaces before the collision and near its endpoint") {
    // Two surfaces stand for a nested object and the enclosing body's exit.
    // The query must not assume that the active medium owns the first surface.
    std::vector<std::shared_ptr<hitable>> faces;
    for (Float z : {.25f, 1.f})
      faces.push_back(std::make_shared<OrderedTriangleProbe>(
          point3f(-10, -10, z), point3f(10, -10, z), point3f(0, 10, z)));
    BVHAggregate world(faces, 0, 1, 1, true);
    Ray ray(point3f(0), vec3f(0, 0, 1));
    ray.segment_absorption = true;
    random_gen rng(824);
    hit_record full, bounded;
    expect_true(world.hit(ray, 0, MaxT, full, rng));
    for (double t : {0.0, .001, .1, .249, .24999999, .25, .25000001, .8, 2.0}) {
      Float limit = SubsurfaceFlightLimit(ray, t);
      expect_true(double(limit) > t);
      bool hit = world.hit(ray, 0, limit, bounded, rng);
      expect_true((t < full.t) == (!hit || t < bounded.t));
      if (hit) expect_true(bounded.t == full.t);
    }
    expect_true(world.hit(ray, 0, SubsurfaceFlightLimit(ray, .24999999), bounded, rng));
    expect_true(SubsurfaceFlightLimit(ray, INFINITY) == MaxT);
    expect_true(SubsurfaceFlightLimit(ray, double(MaxT)) == MaxT);
    // A large coordinate increases the guard rather than rounding below t.
    ray.o = point3f(1000000, -1000000, 1000000);
    expect_true(SubsurfaceFlightLimit(ray, .001) > .001);
  }

  test_that("containment probes do not jump across pointed mesh corners") {
    Transform identity;
    Rcpp::NumericMatrix transform(4, 4);
    for (int i = 0; i < 4; ++i) transform(i, i) = 1;
    auto medium = std::make_shared<Medium>(Rcpp::List::create(
        Rcpp::Named("sigma_a") = Rcpp::NumericVector::create(.1, .1, .1),
        Rcpp::Named("sigma_s") = Rcpp::NumericVector::create(0, 0, 0),
        Rcpp::Named("density_scale") = 1, Rcpp::Named("g") = 0,
        Rcpp::Named("emission") = Rcpp::NumericVector::create(0, 0, 0),
        Rcpp::Named("temperature") = R_NilValue, Rcpp::Named("emission_scale") = 1,
        Rcpp::Named("temperature_scale") = 1, Rcpp::Named("temperature_offset") = 0,
        Rcpp::Named("medium_transform") = transform));
    // The pointed end of the star extrusion. Offsetting the first probe hit
    // along its normal used to jump past the adjacent face's exit.
    point3f vertices[] = {
      point3f(.18834585807153054, 0, -.57966894668189606),
      point3f(.93036954353118950, 0, -.67595304013634394),
      point3f(.6095, 0, 0),
      point3f(.18834585807153054, .65, -.57966894668189606),
      point3f(.93036954353118950, .65, -.67595304013634394),
      point3f(.6095, .65, 0)};
    int faces[][3] = {{0,1,2}, {3,5,4}, {0,3,4}, {0,4,1},
                      {1,4,5}, {1,5,2}, {2,5,3}, {2,3,0}};
    auto geometry = std::make_shared<hitable_list>();
    for (const auto &f : faces) geometry->add(std::make_shared<OrderedTriangleProbe>(
        vertices[f[0]], vertices[f[1]], vertices[f[2]]));
    VolumeScene scene;
    scene.boundaries.add(std::make_shared<MediumBoundary>(geometry, medium, identity, false,
                                                         scene.NextBoundaryId()));
    scene.Finish(0, 1);
    Ray ray(point3f(.93036812543869019, .36028149724006653, -.67595314979553223),
            vec3f(-.71635448932647705, -.4611370861530304, -.52363055944442749));
    expect_true(scene.InitialState(ray, nullptr).media.empty());
    ray.o[2] = -.6759526;
    expect_true(scene.InitialState(ray, nullptr).media.size() == 1);
    ray.segment_absorption = true;
    ray.medium_t_min = 1.234567890123;
    vec3f a(0), b(0), c(0), d(0);
    expect_true(identity(ray).medium_t_min == ray.medium_t_min);
    expect_true(identity(ray, &a, &b).medium_t_min == ray.medium_t_min);
    expect_true(identity(ray, a, b, &c, &d).medium_t_min == ray.medium_t_min);
  }
  test_that("a grazing mesh entry precedes its exit even when Float distances tie") {
    // Exact primary ray from the 256-sample holed extrusion render. It enters
    // the cap only 4.96e-8 ray units before leaving the wall, at t ~= 10.2.
    Ray ray(point3f(4, 4.5, 7), vec3f(-.49616354703903198, -.37251961231231689, -.78425180912017822));
    ray.segment_absorption = true;
    auto cap = std::make_shared<OrderedTriangleProbe>(point3f(-2, .7, -1), point3f(-2, .7, 1), point3f(2, .7, -1));
    auto wall = std::make_shared<OrderedTriangleProbe>(point3f(-2, 0, -1), point3f(-2, .7, -1), point3f(2, .7, -1));
    random_gen rng(19);
    RandomSampler sampler(rng);
    hit_record entry, exit;
    expect_true(cap->hit(ray, 0, MaxT, entry, rng));
    expect_true(wall->hit(ray, 0, MaxT, exit, rng));
    expect_true(entry.OrderedDistance() < exit.OrderedDistance());
    if (sizeof(Float) == sizeof(float)) expect_true(entry.t == exit.t);
    for (bool reversed : {false, true}) {
      std::vector<std::shared_ptr<hitable>> faces = reversed
          ? std::vector<std::shared_ptr<hitable>>{wall, cap}
          : std::vector<std::shared_ptr<hitable>>{cap, wall};
      hitable_list list;
      list.objects = faces;
      BVHAggregate bvh(faces, 0, 1, 1, true);
      for (bool sampled : {false, true}) {
        hit_record h;
        expect_true((sampled ? list.hit(ray, 0, MaxT, h, &sampler) : list.hit(ray, 0, MaxT, h, rng)));
        expect_true(h.shape == cap.get());
        expect_true((sampled ? bvh.hit(ray, 0, MaxT, h, &sampler) : bvh.hit(ray, 0, MaxT, h, rng)));
        expect_true(h.shape == cap.get());
        expect_true(dot(ray.d, h.geometric_normal) < 0);
        // Instances retain ray parameters when transforming hit records.
        Transform transform = Translate(vec3f(3, 5, 7));
        expect_true(transform(h).OrderedDistance() == h.OrderedDistance());
        const hit_record constant = h;
        expect_true(transform(constant).OrderedDistance() == h.OrderedDistance());
      }
    }
    hit_record rounded_down;
    rounded_down.t = Float(10.2);
    rounded_down.precise_t = double(rounded_down.t) + .25 *
        (double(std::nextafter(rounded_down.t, Float(INFINITY))) - rounded_down.t);
    expect_true(double(rounded_down.DistanceUpperBound()) >= rounded_down.precise_t);
  }
  test_that("volume triangle crossings include exact surface origins consistently") {
    Float height = 2.715;
    point3f a(-1, height, -1), b(1, height, -1), c(0, height, 1);
    point3f origin(.2, height, .1);
    Float t, b0, b1, b2;
    for (Float sign : {-1.f, 1.f}) {
      Ray ray(origin, unit_vector(vec3f(.317, sign, .129)));
      expect_true(VolumeTriangleIntersection(ray, a, b, c, 0, MaxT, t, b0, b1, b2));
      expect_true(t == 0);
      expect_true(std::abs(b0 + b1 + b2 - 1) < 1e-6);
      expect_false(VolumeTriangleIntersection(ray, a, b, c, 1e-8, MaxT, t, b0, b1, b2));
    }
    origin[1] = std::nextafter(height, Float(-INFINITY));
    expect_true(VolumeTriangleIntersection(Ray(origin, vec3f(0,1,0)), a, b, c,
                                           0, MaxT, t, b0, b1, b2));
    expect_true(t > 0);
    expect_false(VolumeTriangleIntersection(Ray(origin, vec3f(0,-1,0)), a, b, c,
                                            0, MaxT, t, b0, b1, b2));
  }
  test_that("glass priority suppresses SSS and selects the actual contact IOR") {
    Rcpp::NumericMatrix transform(4, 4);
    for (int i = 0; i < 4; ++i) transform(i, i) = 1;
    auto medium = std::make_shared<Medium>(Rcpp::List::create(
        Rcpp::Named("sigma_a") = Rcpp::NumericVector::create(.01, .01, .01),
        Rcpp::Named("sigma_s") = Rcpp::NumericVector::create(8, 8, 8),
        Rcpp::Named("density_scale") = 1, Rcpp::Named("g") = 0,
        Rcpp::Named("emission") = Rcpp::NumericVector::create(0, 0, 0),
        Rcpp::Named("temperature") = R_NilValue, Rcpp::Named("emission_scale") = 1,
        Rcpp::Named("temperature_scale") = 1, Rcpp::Named("temperature_offset") = 0,
        Rcpp::Named("medium_transform") = transform));
    medium->subsurface = true;
    Transform identity;
    auto geometry = std::make_shared<sphere>(1, nullptr, nullptr, nullptr,
                                             &identity, &identity, false);
    MediumBoundary boundary(geometry, medium, identity, true, 1);
    dielectric glass(point3f(1), 1.5, point3f(0), 0);
    dielectric milk(point3f(1), 1.333, point3f(0), 1);
    hit_record wall, liquid;
    wall.mat_ptr = &glass;
    wall.geometric_normal = wall.normal = normal3f(1, 0, 0);
    liquid = wall;
    liquid.mat_ptr = &milk;
    liquid.medium_boundary = &boundary;
    liquid.boundary_id = 1;
    vec3f enter(-1, 0, 0), leave(1, 0, 0);
    VolumePathState state;
    state.CrossDielectric(wall, enter);
    expect_true(state.Interface(liquid, enter).Hidden());
    // At a grazing hidden wedge, a rounded origin can skip the losing solid's
    // entry. An observed exit still leaves glass as the sole physical medium.
    state.CrossDielectric(liquid, leave);
    expect_true(state.glass.size() == 1);
    expect_true(state.media.empty());
    state.CrossDielectric(liquid, enter);
    state.CrossDielectric(liquid, enter);
    expect_true(state.glass.size() == 2);
    expect_true(state.media.size() == 1);
    expect_true(state.ContainsSubsurface());
    expect_true(state.Active() == nullptr);
    auto contact = state.Interface(wall, leave);
    expect_false(contact.Hidden());
    expect_true(std::abs(contact.Eta() - 1.333 / 1.5) < 1e-6);
    // Merely evaluating reflection must not mutate either membership list.
    expect_true(state.glass.size() == 2);
    expect_true(state.media.size() == 1);
    state.CrossDielectric(wall, leave);
    expect_true((state.Active() && state.Active()->boundary_id == 1));
    expect_false(state.Active()->guide_valid);
    state.CrossDielectric(liquid, leave);
    expect_true((state.glass.empty() && state.media.empty()));
    bool rejected = false;
    try { state.CrossDielectric(liquid, leave); }
    catch (const std::runtime_error &) { rejected = true; }
    expect_true(rejected);
    // Traverse the same overlap in reverse: the submerged milk boundary exits
    // while glass still wins, so it must be a null interface, not milk-to-air.
    state.CrossDielectric(liquid, enter);
    state.CrossDielectric(wall, enter);
    expect_true(state.Interface(liquid, leave).Hidden());
    state.CrossDielectric(liquid, leave);
    expect_true(state.Active() == nullptr);
    expect_true(state.ActiveDielectric() == &glass);
    state.CrossDielectric(wall, leave);
    expect_true((state.glass.empty() && state.media.empty()));
    // Equal priorities select the most recently entered surface.
    milk.priority = 0;
    state.CrossDielectric(wall, enter);
    expect_false(state.Interface(liquid, enter).Hidden());
    state.CrossDielectric(liquid, enter);
    expect_true((state.Active() && state.Active()->boundary_id == 1));
    expect_true(state.Interface(wall, leave).Hidden());
    state.CrossDielectric(wall, leave);
    state.CrossDielectric(liquid, leave);
    expect_true((state.glass.empty() && state.media.empty()));
    // A hidden SSS shell must not mask an independently authored enclosing
    // explicit medium. Its existing transport through glass remains intact.
    auto enclosing_medium = std::make_shared<Medium>(*medium);
    enclosing_medium->subsurface = false;
    MediumBoundary enclosing(geometry, enclosing_medium, identity, false, 2);
    hit_record fog = liquid;
    fog.medium_boundary = &enclosing;
    fog.boundary_id = 2;
    milk.priority = 1;
    state.Cross(fog, enter);
    state.CrossDielectric(wall, enter);
    expect_true(state.Active()->boundary_id == 2);
    state.CrossDielectric(liquid, enter);
    expect_true((state.Active() && state.Active()->boundary_id == 2));
    state.CrossDielectric(liquid, leave);
    expect_true(state.Active()->boundary_id == 2);
    state.CrossDielectric(wall, leave);
    state.Cross(fog, leave);
    expect_true((state.glass.empty() && state.media.empty()));
  }
  test_that("a collision inside endpoint error bounds stays on its incoming side") {
    // Rounded position from the high-sample thick-sphere regression, outside
    // the unit sphere by 3.9e-9 despite a sampled distance below the hit's t.
    point3f p(.087437309324741364, .74083441495895386, .66597229242324829);
    hit_record h;
    h.p = p;
    h.p *= 1 / h.p.length();
    h.pError = convert_to_vec3(gamma(5) * Abs(h.p));
    h.normal = h.geometric_normal = unit_vector(convert_to_normal3(h.p));
    h.t = 1;
    Ray ray(point3f(0), unit_vector(convert_to_vec3(p)));
    bool adjusted = false;
    expect_true(Float(.9999999) < h.t);
    auto inside = SubsurfaceCollisionPoint(p, ray, .9999999, h, adjusted);
    expect_true(adjusted);
    double radius2 = 0;
    for (int a=0; a<3; ++a) radius2 += double(inside[a])*inside[a];
    expect_true(radius2 < 1);
    auto unchanged = SubsurfaceCollisionPoint(point3f(0), ray, .5, h, adjusted);
    expect_false(adjusted);
    expect_true(unchanged.squared_length() == 0);
  }
  test_that("grazing collisions stay inside a plane without jumping to its distant endpoint") {
    const Float plane = Float(2.715);
    for (Float side : {Float(-1), Float(1)}) {
      const Float start = std::nextafter(plane, side < 0 ? Float(-INFINITY) : Float(INFINITY));
      Ray ray(point3f(0, start, 0), vec3f(1, plane - start, 0));
      hit_record h;
      h.p = point3f(1, plane, 0);
      h.pError = gamma(7) * convert_to_vec3(Abs(h.p));
      h.normal = h.geometric_normal = normal3f(0, -side, 0);
      h.t = 1;
      const double t = .75;
      point3f collision(Float(t), Float(double(start) + double(ray.d[1]) * t), 0);
      expect_true(collision[1] == plane);
      bool adjusted = false;
      auto inside = SubsurfaceCollisionPoint(collision, ray, t, h, adjusted);
      expect_true(adjusted);
      expect_true(side * (inside[1] - plane) > 0);
      expect_true(inside[0] == collision[0]);
      expect_true(inside[2] == collision[2]);
    }
  }
  test_that("coordinate-zero boundary planes get a resolvable scale-aware offset") {
    hit_record h;
    h.p = point3f(0, .3, .2);
    h.pError = vec3f(0, 1e-7, 1e-7);
    h.geometric_normal = h.normal = normal3f(1, 0, 0);
    auto outside = OffsetMediumOrigin(h, vec3f(1, 0, 0));
    auto inside = OffsetMediumOrigin(h, vec3f(-1, 0, 0));
    expect_true(outside[0] > 1e-8);
    expect_true(inside[0] < -1e-8);
    h.p *= .01;
    h.pError *= .01;
    auto scaled = OffsetMediumOrigin(h, vec3f(1, 0, 0));
    expect_true(std::abs(scaled[0] / outside[0] - .01) < 1e-8);
  }
  test_that("extreme optical scales and physical phase parameters remain finite") {
    SubsurfaceProposal p;
    p.guided = true;
    p.extinction = {1e-18, 1, 1e18};
    p.pole.fill(SubsurfaceProposal::Pole(.01));
    for (int c = 0; c < 3; ++c) {
      double b = 1 / p.extinction[c];
      expect_true(std::isfinite(p.LogDistancePdf(c, vec3f(0,0,1), b, true)));
      expect_true(std::isfinite(p.LogDistancePdf(c, vec3f(0,0,1), b, false)));
      expect_true(std::isfinite(p.SampleDistance(c, vec3f(0,0,1), .3, .999999)));
    }
    expect_true(SubsurfaceProposal::Pole(0) == 0);
    expect_true(SubsurfaceProposal::Pole(1) == 0);
    for (Float g : {-.9999f, .9999f}) {
      HGPhaseFunction phase(g);
      auto s = phase.Sample(vec3f(0,0,-1), .4, .9);
      expect_true(std::isfinite(s.p));
      expect_true(s.p > 0);
    }
  }
  test_that("direction and joint flight densities normalize including endpoint survival") {
    SubsurfaceProposal p;
    p.guided = true;
    p.extinction = {0, 1, 8};
    p.pole = {0, SubsurfaceProposal::Pole(.8), SubsurfaceProposal::Pole(.99)};
    p.wo = vec3f(0, 0, -1);
    p.axis = vec3f(0, 0, 1);
    const int n = 12000;
    double norm[3] = {}, joint[3] = {};
    for (int i = 0; i < n; ++i) {
      double mu = 2 * (i + .5) / n - 1;
      vec3f wi(std::sqrt(1 - mu * mu), 0, mu);
      for (int c = 0; c < 3; ++c) {
        double d = p.DirectionPdf(c, wi);
        norm[c] += d * 4 * M_PI / n;
        double exit = std::exp(p.LogDistancePdf(c, wi, .7, false));
        // Integrate continuous collision mass independently with midpoint quadrature.
        double collisions = 0;
        for (int k = 0; k < 60; ++k)
          collisions += std::exp(p.LogDistancePdf(c, wi, .7 * (k + .5) / 60, true)) * .7 / 60;
        joint[c] += d * (exit + collisions) * 4 * M_PI / n;
      }
    }
    for (int c = 0; c < 3; ++c) {
      expect_true(std::abs(norm[c] - 1) < 2e-6);
      expect_true(std::abs(joint[c] - 1) < .0006);
    }
    expect_true(p.SampleDistance(0, vec3f(0,0,1), .2, .5) == INFINITY);
  }
  test_that("sampled joint moments agree with independently integrated densities") {
    SubsurfaceProposal p;
    p.guided = true;
    p.extinction = {1, 1, 1};
    p.pole.fill(SubsurfaceProposal::Pole(.7));
    p.wo = vec3f(0, 0, -1);
    p.axis = vec3f(0, 0, 1);
    random_gen rng(729);
    double expected_exit = 0, expected_mu_exit = 0;
    for (int i = 0; i < 10000; ++i) {
      double mu = 2 * (i + .5) / 10000 - 1;
      vec3f wi(std::sqrt(1-mu*mu), 0, mu);
      double q = p.DirectionPdf(0,wi) * std::exp(p.LogDistancePdf(0,wi,.8,false));
      expected_exit += q * 4*M_PI/10000;
      expected_mu_exit += mu * q * 4*M_PI/10000;
    }
    double exits=0, mu_exit=0;
    for (int i=0; i<100000; ++i) {
      double strategy=rng.unif_rand(), u=rng.unif_rand(), v=rng.unif_rand();
      auto wi=p.SampleDirection(0,strategy,u,v);
      strategy=rng.unif_rand(); u=rng.unif_rand();
      double t=p.SampleDistance(0,wi,strategy,u);
      if(t>=.8) { ++exits; mu_exit+=wi[2]; }
    }
    expect_true(std::abs(exits/100000-expected_exit)<.006);
    expect_true(std::abs(mu_exit/100000-expected_mu_exit)<.006);
  }
  test_that("HG retains the physical anisotropy sign under guidance") {
    for (Float g : {-.9f, 0.f, .9f}) {
      HGPhaseFunction phase(g);
      double mean=0;
      for(int i=0;i<50000;++i) mean+=phase.Sample(vec3f(0,0,-1), (i+.5)/50000, .4).wi[2];
      expect_true(std::abs(mean/50000-g)<1e-4);
    }
  }
  test_that("smooth dielectric has exact normal Fresnel, TIR and radiance eta weights") {
    SubsurfaceBoundaryBSDF b(vec3f(0,0,-1),normal3f(0,0,1),1.5,0);
    auto reflection=b.Sample(.01,0,0), transmission=b.Sample(.5,0,0);
    expect_true(std::abs(reflection.pdf-.04)<1e-6);
    expect_true(reflection.weight==1);
    expect_true(std::abs(transmission.weight-1/2.25)<1e-6);
    expect_true(transmission.transmission);
    SubsurfaceBoundaryBSDF exit(vec3f(.9,0,std::sqrt(.19)),normal3f(0,0,1),1/1.5,0);
    expect_false(exit.Sample(.99,0,0).transmission);
    expect_true(exit.Sample(.99,0,0).weight==1);
    SubsurfaceBoundaryBSDF matched(vec3f(0,0,-1),normal3f(0,0,1),1,.8);
    expect_true(matched.Sample(.5,.4,.5).weight==1);
  }
  test_that("rough dielectric samples match PDFs and adjoint reciprocity") {
    random_gen rng(188);
    for (double roughness : {.15, .5, 1.}) {
      vec3f incoming=unit_vector(vec3f(.3,0,-1));
      SubsurfaceBoundaryBSDF b(incoming,normal3f(0,0,1),1.4,roughness);
      double energy=0;
      int transmitted=0;
      double pdf_error=0, weight_error=0, reciprocity_error=0;
      bool finite=true;
      for(int i=0;i<30000;++i) {
        double branch=rng.unif_rand(),u=rng.unif_rand(),v=rng.unif_rand();
        auto s=b.Sample(branch,u,v);
        if(s.weight==0) continue;
        finite = finite && std::isfinite(s.weight);
        pdf_error = std::max(pdf_error, std::abs(s.pdf-b.Pdf(s.wi)));
        weight_error = std::max(weight_error, std::abs(s.weight*s.pdf-b.Evaluate(s.wi)));
        energy+=s.weight*(s.transmission?1.4*1.4:1);
        if(s.transmission && i%100==0) {
          ++transmitted;
          SubsurfaceBoundaryBSDF reverse(-s.wi,normal3f(0,0,1),1/1.4,roughness);
          double forward=b.Evaluate(s.wi)/std::abs(s.wi[2]);
          double backward=reverse.Evaluate(-incoming)/std::abs(incoming[2]);
          reciprocity_error = std::max(reciprocity_error, std::abs(forward*1.4*1.4-backward)/std::max(1.,backward));
        }
      }
      expect_true(finite);
      expect_true(pdf_error<1e-8);
      expect_true(weight_error<1e-8);
      expect_true(reciprocity_error<.004);
      expect_true(transmitted>20);
      expect_true((energy/30000>.65 && energy/30000<1.01));
    }
  }
}

context("Normalized diffusion") {
  test_that("coplanar projection rays do not produce NaN rectangle hits") {
    Transform identity;
    std::vector<std::shared_ptr<hitable>> faces{
      std::make_shared<yz_rect>(-1, 1, -1, 1, 0, nullptr, nullptr, nullptr, &identity, &identity, false),
      std::make_shared<xz_rect>(-1, 1, -1, 1, 0, nullptr, nullptr, nullptr, &identity, &identity, false),
      std::make_shared<xy_rect>(-1, 1, -1, 1, 0, nullptr, nullptr, nullptr, &identity, &identity, false)};
    random_gen rng(72);
    RandomSampler sampler(rng);
    for (int axis = 0; axis < 3; ++axis) {
      for (bool volume_ray : {false, true}) {
        vec3f direction(0);
        direction[(axis + 1) % 3] = 1;
        for (Float plane_distance : {Float(0), Float(.5)}) {
          point3f origin(0);
          origin[axis] = plane_distance;
          Ray ray(origin, direction);
          ray.segment_absorption = volume_ray;
          hit_record hit;
          expect_false(faces[axis]->hit(ray, 0, MaxT, hit, rng));
          expect_false(faces[axis]->hit(ray, 0, MaxT, hit, &sampler));
          expect_false(faces[axis]->HitP(ray, 0, MaxT, rng));
          expect_false(faces[axis]->HitP(ray, 0, MaxT, &sampler));
        }
      }
    }
  }

  test_that("normalized radial profiles integrate to one and sample the full tail") {
    NormalizedDiffusionProfile profile{{.2, .7, 1}, {.01, .2, 3}};
    random_gen rng(927);
    for (int c = 0; c < 3; ++c) {
      double integral = 0, mean = 0;
      int below = 0;
      constexpr int count = 100000;
      for (int i = 0; i < count; ++i) {
        double r = profile.radius[c] * 90 * (i + .5) / count;
        integral += 2 * M_PI * r * std::exp(profile.LogAreaPdf(c, r)) *
                    profile.radius[c] * 90 / count;
        double component = rng.unif_rand(), uniform = rng.unif_rand();
        double sample = profile.SampleRadius(c, component, uniform);
        mean += sample / count;
        below += sample <= profile.radius[c];
      }
      expect_true(std::abs(integral - 1) < 1e-6);
      expect_true(std::abs(mean / profile.radius[c] - 2.5) < .035);
      expect_true(std::abs(double(below) / count - profile.Cdf(c, profile.radius[c])) < .006);
      expect_true(profile.SampleRadius(c, .5, 1 - 1e-12) > 80 * profile.radius[c]);
    }
  }
  test_that("exit lobes normalize Fresnel for air and higher-index neighboring glass") {
    for (double ratio : {1.0, 1.333, 1.333 / 1.5, 1.5, .5, 2.5}) {
      DiffusionExitBSDF lobe(normal3f(0, 0, 1), ratio, 1);
      double integral = 0;
      constexpr int count = 100000;
      for (int i = 0; i < count; ++i) {
        double cosine = (i + .5) / count;
        vec3f wi(Float(std::sqrt(1 - cosine * cosine)), 0, Float(cosine));
        integral += 2 * M_PI * lobe.Evaluate(wi) / count / lobe.EtaSquared();
      }
      expect_true(std::abs(integral - 1) < 2e-5);
      expect_true(lobe.Evaluate(vec3f(0, 0, -1)) == 0);
      expect_true(lobe.Pdf(vec3f(0, 0, -1)) == 0);
      vec3f sample = lobe.Sample(.3, .7);
      expect_true(std::abs(sample.length() - 1) < 1e-6);
      expect_true(lobe.Pdf(sample) > 0);
    }
  }
  test_that("Fresnel exit sampling matches its PDF and removes angular weight variance") {
    for (double ratio : {.25, .5, 1.333 / 1.5, .9999, 1.0, 1.0001, 1.333, 2.5, 10.0}) {
      DiffusionExitBSDF lobe(normal3f(0, 0, 1), ratio, 1);
      constexpr int count = 100000;
      double integral = 0, pdf_moment = 0, sample_moment = 0;
      double weight_sum = 0, weight_square_sum = 0;
      double old_sum = 0, old_square_sum = 0;
      for (int i = 0; i < count; ++i) {
        double u = (i + .5) / count;
        vec3f direction(Float(std::sqrt(1 - u * u)), 0, Float(u));
        double pdf = lobe.Pdf(direction);
        integral += 2 * M_PI * pdf / count;
        pdf_moment += 2 * M_PI * pdf * u / count;
        vec3f sample = lobe.Sample(u, .37);
        double sample_pdf = lobe.Pdf(sample);
        double weight = sample_pdf > 0 ? lobe.Evaluate(sample) / sample_pdf / lobe.EtaSquared() : 0;
        sample_moment += sample[2] / count;
        weight_sum += weight / count;
        weight_square_sum += weight * weight / count;
        // The old cosine estimator evaluates exactly the same physical lobe.
        double cosine = std::sqrt(1 - u);
        vec3f old_sample(Float(std::sqrt(u)), 0, Float(cosine));
        double old_weight = lobe.Evaluate(old_sample) * M_PI / cosine / lobe.EtaSquared();
        old_sum += old_weight / count;
        old_square_sum += old_weight * old_weight / count;
      }
      expect_true(std::abs(integral - 1) < 3e-4);
      expect_true(std::abs(pdf_moment - sample_moment) < 3e-4);
      expect_true(std::abs(weight_sum - 1) < 3e-5);
      double variance = weight_square_sum - weight_sum * weight_sum;
      double old_variance = old_square_sum - old_sum * old_sum;
      expect_true(variance < 1e-5);
      if (ratio < .9 || ratio > 1.2) expect_true(variance < old_variance * .001);
      for (double endpoint : {0.0, 1.0}) {
        vec3f sample = lobe.Sample(endpoint, endpoint);
        expect_true(std::isfinite(lobe.Pdf(sample)));
        expect_true(lobe.Pdf(sample) > 0);
        expect_true(std::abs(sample.length() - 1) < 1e-6);
      }
    }
    // Frame transforms and critical-angle rounding must not create NaN weights.
    DiffusionExitBSDF rotated(unit_vector(normal3f(1, 2, 3)), 1.333, 1.5);
    for (double u : {0.0, .1, .9, .999, 1.0}) {
      vec3f sample = rotated.Sample(u, .6);
      double pdf = rotated.Pdf(sample);
      expect_true((pdf >= 0 && std::isfinite(pdf)));
      expect_true(std::isfinite(rotated.Evaluate(sample)));
    }
  }
  test_that("diffusion samples effective glass boundaries and excludes hidden liquid faces") {
    Rcpp::NumericMatrix transform(4, 4);
    for (int i = 0; i < 4; ++i) transform(i, i) = 1;
    auto medium = std::make_shared<Medium>(Rcpp::List::create(
        Rcpp::Named("sigma_a") = Rcpp::NumericVector::create(0, 0, 0),
        Rcpp::Named("sigma_s") = Rcpp::NumericVector::create(0, 0, 0),
        Rcpp::Named("density_scale") = 1, Rcpp::Named("g") = 0,
        Rcpp::Named("emission") = Rcpp::NumericVector::create(0, 0, 0),
        Rcpp::Named("temperature") = R_NilValue, Rcpp::Named("emission_scale") = 1,
        Rcpp::Named("temperature_scale") = 1, Rcpp::Named("temperature_offset") = 0,
        Rcpp::Named("medium_transform") = transform,
        Rcpp::Named("subsurface") = Rcpp::List::create(
          Rcpp::Named("method") = "diffusion", Rcpp::Named("refraction") = 1.333,
          Rcpp::Named("roughness") = 0,
          Rcpp::Named("color") = Rcpp::NumericVector::create(1, 1, 1),
          Rcpp::Named("radius") = Rcpp::NumericVector::create(.3, .3, .3))));
    Transform identity, glass_transform = Translate(vec3f(.6, 0, 0));
    Transform glass_inverse = Inverse(glass_transform);
    auto milk = std::make_shared<dielectric>(point3f(1), 1.333, point3f(0), 1);
    auto glass = std::make_shared<dielectric>(point3f(1), 1.5, point3f(0), 0);
    auto milk_shape = std::make_shared<sphere>(1, milk, nullptr, nullptr, &identity, &identity, false);
    auto glass_shape = std::make_shared<sphere>(.8, glass, nullptr, nullptr,
                                              &glass_transform, &glass_inverse, false);
    VolumeScene scene;
    auto milk_boundary = std::make_shared<MediumBoundary>(milk_shape, medium, identity, true,
                                                         scene.NextBoundaryId());
    auto glass_boundary = std::make_shared<MediumBoundary>(glass_shape, nullptr, glass_transform,
                                                          true, scene.NextBoundaryId());
    scene.boundaries.add(milk_boundary);
    scene.boundaries.add(glass_boundary);
    scene.Finish(0, 1);
    random_gen rng(1927);
    RandomSampler sampler(rng);
    auto state = scene.InitialState(Ray(point3f(-.99, 0, 0), vec3f(1, 0, 0)), nullptr);
    expect_true(state.Active() != nullptr);
    MediumEntry body = *state.Active();
    int own_hits = 0, glass_hits = 0;
    for (int i = 0; i < 3000; ++i) {
      auto sample = SampleDiffusionSurface(scene, body, point3f(-1, 0, 0),
                                           normal3f(-1, 0, 0), 0, &sampler, rng, nullptr);
      if (!sample) continue;
      expect_true((!sample->outside.Active() || sample->outside.Active()->boundary_id != body.boundary_id));
      bool on_glass = sample->hit.medium_boundary == glass_boundary.get();
      glass_hits += on_glass;
      own_hits += !on_glass;
      if (on_glass) {
        expect_true(sample->outside.ActiveDielectric() == glass.get());
        expect_true(dot(sample->outward, sample->hit.geometric_normal) < -.99);
        expect_true((sample->hit.p - point3f(0)).length() < 1.00001);
      } else {
        expect_true((sample->hit.p - point3f(.6, 0, 0)).length() >= .79999);
      }
      for (double w : sample->weight) expect_true((std::isfinite(w) && w >= 0));
    }
    expect_true(own_hits > 100);
    expect_true(glass_hits > 100);
    // Reversing priority restores the uncut milk surface, including the overlap.
    glass->priority = 2;
    int buried = 0;
    for (int i = 0; i < 500; ++i) {
      auto sample = SampleDiffusionSurface(scene, body, point3f(-1, 0, 0),
                                           normal3f(-1, 0, 0), 0, &sampler, rng, nullptr);
      if (!sample) continue;
      expect_true(sample->hit.medium_boundary == milk_boundary.get());
      buried += (sample->hit.p - point3f(.6, 0, 0)).length() < .8;
    }
    expect_true(buried > 5);
    // A white sphere in a white furnace has response CDF(2R), approximately
    // one for this short profile. Check both its mean and spatial-estimator
    // second moment: uniform far-side selection formerly gave about four.
    medium->diffusion_radius = {.01, .01, .01};
    double mean = 0, second_moment = 0;
    constexpr int samples = 30000;
    for (int i = 0; i < samples; ++i) {
      auto sample = SampleDiffusionSurface(scene, body, point3f(-1, 0, 0),
                                           normal3f(-1, 0, 0), 0, &sampler, rng, nullptr);
      double w = sample ? sample->weight[0] : 0;
      mean += w / samples;
      second_moment += w * w / samples;
    }
    expect_true(std::abs(mean - 1) < .025);
    expect_true(second_moment < 2.6);
    std::atomic<bool> cancelled{true};
    expect_false(SampleDiffusionSurface(scene, body, point3f(-1, 0, 0), normal3f(-1, 0, 0),
                                        0, &sampler, rng, &cancelled).has_value());
  }
}
#endif
