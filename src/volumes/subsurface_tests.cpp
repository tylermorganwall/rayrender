#ifdef NOT_CRAN
#include "subsurface.h"
#include "boundary.h"
#include "../math/rng.h"
#include "../materials/material.h"
#include "../hitables/sphere.h"
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
#endif
