#ifdef NOT_CRAN
#include "../core/bvh.h"
#include "../hitables/box.h"
#include "../hitables/instance.h"
#include "../hitables/mesh3d.h"
#include "../hitables/infinite_area_light.h"
#include "../hitables/rectangle.h"
#include "../hitables/sphere.h"
#include "../hitables/triangle.h"
#include "../materials/material.h"
#include "../materials/texture.h"
#include "boundary.h"
#include "medium.h"
#include "haze.h"
#include "picking.h"
#include "volpath.h"
#include <chrono>
#include <set>
#include <stdexcept>
#include <testthat.h>
#include <thread>

namespace {
Rcpp::List medium_description(Float scattering = 1) {
  Rcpp::NumericMatrix transform(4, 4);
  for (int i = 0; i < 4; ++i)
    transform(i, i) = 1;
  return Rcpp::List::create(
      Rcpp::Named("sigma_a") = Rcpp::NumericVector::create(0, 0, 0),
      Rcpp::Named("sigma_s") = Rcpp::NumericVector::create(scattering, scattering, scattering),
      Rcpp::Named("density_scale") = 1, Rcpp::Named("g") = 0,
      Rcpp::Named("emission") = Rcpp::NumericVector::create(0, 0, 0),
      Rcpp::Named("temperature") = R_NilValue, Rcpp::Named("emission_scale") = 1,
      Rcpp::Named("temperature_scale") = 1, Rcpp::Named("temperature_offset") = 0,
      Rcpp::Named("medium_transform") = transform);
}
struct PickingScene {
  Transform identity;
  VolumeScene scene;
  hitable_list world;
  std::shared_ptr<material> mat =
      std::make_shared<lambertian>(std::make_shared<constant_texture>(point3f(1)));
  void Add(const std::shared_ptr<const Medium> &medium, Float radius, bool surface = false) {
    auto geometry = std::make_shared<box>(vec3f(-radius), vec3f(radius), mat, nullptr, nullptr,
                                         &identity, &identity, false);
    auto boundary = std::make_shared<MediumBoundary>(geometry, medium, identity, surface,
                                                    scene.NextBoundaryId());
    scene.boundaries.add(boundary);
    world.add(boundary);
    scene.Finish(0, 1);
  }
};
Rcpp::List grid_description(const Rcpp::NumericVector &values, int nx, int ny, int nz) {
  Rcpp::List description = medium_description();
  Rcpp::NumericVector density = Rcpp::clone(values);
  density.attr("dim") = Rcpp::IntegerVector::create(nx, ny, nz);
  Rcpp::NumericMatrix bounds(2, 3);
  for (int a = 0; a < 3; ++a) {
    bounds(0, a) = -1;
    bounds(1, a) = 1;
  }
  description["density"] = density;
  description["bounds"] = bounds;
  return description;
}
} // namespace
context("Hit record medium placement") {
  test_that("ordinary records defer matrices and transformed boundaries retain their placement") {
    hit_record ordinary;
    expect_false(ordinary.medium_to_world.has_value());
    expect_true(ordinary.MediumToWorld().IsIdentity());
    expect_false(ordinary.medium_to_world.has_value());

    PickingScene scene;
    scene.Add(std::make_shared<Medium>(medium_description()), 1);
    Ray ray(point3f(-3, 0, 0), vec3f(1, 0, 0));
    ray.segment_absorption = true;
    random_gen rng(19);
    hit_record boundary;
    expect_true(scene.world.hit(ray, 0, 10, boundary, rng));
    expect_true(boundary.medium_to_world.has_value());
    Transform placement = Translate(vec3f(2, 3, 4)) * Scale(2, 3, 4);
    const hit_record transformed = placement(boundary);
    expect_true(transformed.MediumToWorld() == placement);
    hit_record copied = transformed;
    expect_true(copied.MediumToWorld() == placement);
    VolumePathState state;
    state.Cross(copied, vec3f(1, 0, 0));
    expect_true((state.Active() && state.Active()->medium_to_world == placement));
    copied.medium_to_world.reset();
    expect_true(copied.MediumToWorld().IsIdentity());
  }
}
context("Grid medium emission") {
  test_that("zero RGB fields preserve density and return zero emission") {
    for (bool array_field : {false, true}) {
      auto description = grid_description(Rcpp::NumericVector(8, .35), 2, 2, 2);
      if (array_field) {
        Rcpp::NumericVector emission(24, 0.0);
        emission.attr("dim") = Rcpp::IntegerVector::create(2, 2, 2, 3);
        description["emission"] = emission;
      }
      GridMedium medium(description);
      for (point3f p : {point3f(0), point3f(.75), point3f(2)}) {
        auto properties = medium.SamplePoint(p);
        for (int c = 0; c < 3; ++c) {
          expect_true(properties.Le[c] == 0);
          expect_true(medium.Emission(p)[c] == 0);
          expect_true(properties.sigma_s[c] == medium.Density(p));
        }
      }
    }
  }
  test_that("a zero first voxel and zero absorption do not disable nonzero RGB emission") {
    auto description = grid_description(Rcpp::NumericVector(8, 1.0), 2, 2, 2);
    Rcpp::NumericVector emission = Rcpp::NumericVector::create(0, 2, 0, 4, 0, 6);
    emission.attr("dim") = Rcpp::IntegerVector::create(2, 1, 1, 3);
    description["emission"] = emission;
    description["emission_scale"] = 2;
    GridMedium medium(description);
    expect_false(medium.IsEmissive()); // sigma_a is zero; Le itself is still nonzero.
    for (int c = 0; c < 3; ++c) {
      expect_true(medium.SamplePoint(point3f(0)).Le[c] == 2 * (c + 1));
      expect_true(medium.Emission(point3f(-.5, 0, 0))[c] == 0);
      expect_true(medium.Emission(point3f(.5, 0, 0))[c] == 4 * (c + 1));
    }
  }
  test_that("temperature emission survives zero RGB and an emission scale of zero disables it") {
    auto description = grid_description(Rcpp::NumericVector(8, 1.0), 2, 2, 2);
    Rcpp::NumericVector temperature = Rcpp::NumericVector::create(3000, 6000);
    temperature.attr("dim") = Rcpp::IntegerVector::create(2, 1, 1);
    description["temperature"] = temperature;
    description["temperature_offset"] = 500;
    description["temperature_scale"] = 2;
    for (double scale : {1.5, 0.}) {
      description["emission_scale"] = scale;
      GridMedium medium(description);
      point3f expected = Float(scale) * BlackbodyRGB(8000);
      for (int c = 0; c < 3; ++c) {
        expect_true(medium.Emission(point3f(0))[c] == expected[c]);
        expect_true(medium.SamplePoint(point3f(0)).Le[c] == expected[c]);
      }
    }
  }
}
context("Density-selected atmospheric haze") {
  test_that("density metadata preserves fractional RGB tracking bounds") {
    auto description = grid_description(Rcpp::NumericVector(8, .35), 2, 2, 2);
    description["sigma_s"] = Rcpp::NumericVector::create(.2, 1.1, 2.2);
    GridMedium medium(description);
    auto iterator = medium.SampleRay(Ray(point3f(0, 0, -1), vec3f(0, 0, 1)), 2);
    auto segment = iterator.Next();
    for (int c = 0; c < 3; ++c) {
      double expected = double(medium.sigma_s[c]) * Float(.35);
      expect_true(segment->sigma_maj[c] >= expected);
      expect_true(std::abs(segment->sigma_maj[c] - expected) < 1e-6);
    }
    expect_true(std::abs(segment->density_max - .35) < 1e-6);
  }
  test_that("threshold crossings preserve empty space and scaled grid ramps without collisions") {
    auto description = grid_description(Rcpp::NumericVector::create(0, 2), 1, 1, 2);
    description["sigma_s"] = Rcpp::NumericVector::create(0, 0, 0);
    description["haze"] = true;
    description["haze_density_threshold"] = .5;
    for (double scale : {1., 2.}) {
      description["density_scale"] = scale;
      GridMedium medium(description);
      for (double speed : {1., 2.}) {
        MediumHazeIterator iterator(&medium, Ray(point3f(0, 0, -3), vec3f(0, 0, speed)),
                                      5 / speed, true);
        auto air = iterator.Next(), dense = iterator.Next(), tail = iterator.Next();
        expect_true((bool(air) && bool(dense) && bool(tail)));
        expect_true((air->integrate && !dense->integrate && tail->integrate));
        expect_true(air->t_min == 0);
        expect_true(std::abs(air->t_max - (2.5 + .25 / scale) / speed) < 1e-10);
        expect_true((dense->t_min == air->t_max && dense->t_max == 4 / speed));
        expect_true((tail->t_min == dense->t_max && tail->t_max == 5 / speed));
        expect_false(bool(iterator.Next()));
      }
      MediumHazeIterator reverse(&medium, Ray(point3f(0, 0, 3), vec3f(0, 0, -1)), 6, true);
      auto air = reverse.Next(), dense = reverse.Next(), tail = reverse.Next();
      expect_true(air->t_max == 2);
      expect_true(std::abs(dense->t_max - (3.5 - .25 / scale)) < 1e-10);
      expect_true((tail->integrate && tail->t_max == 6));
      MediumHazeIterator disabled(&medium, Ray(point3f(0), vec3f(0, 0, 1)), 2, false);
      auto off = disabled.Next();
      expect_true((!off->integrate && off->t_max == 2));
      expect_false(bool(disabled.Next()));
    }
  }
  test_that("cubic density crosses a threshold three times inside one cell") {
    class CubicMedium : public Medium {
    public:
      CubicMedium() : Medium(medium_description(0)) { haze_density_threshold = 1; }
      bool IsHomogeneous() const override { return false; }
      DensityIndexRay DensityRay(const Ray &) const override { return {{0, 0, 0}, {1, 1, 1}}; }
      std::array<double, 8> DensityCorners(const std::array<double, 3> &) const override {
        // rho(u) = 1 + (u - .2)(u - .5)(u - .8), in Bernstein form.
        return {.92, 1.14, 1.14, .86, 1.14, .86, .86, 1.08};
      }
    } medium;
    MediumHazeIterator iterator(&medium, Ray(point3f(0), vec3f(1)), 1, true);
    double endpoints[4] = {.2, .5, .8, 1};
    for (int i = 0; i < 4; ++i) {
      auto segment = iterator.Next();
      expect_true(bool(segment));
      expect_true(segment->integrate == (i % 2 == 0));
      expect_true(std::abs(segment->t_max - endpoints[i]) < 1e-10);
    }
    expect_false(bool(iterator.Next()));
    std::atomic<bool> cancel(true);
    MediumHazeIterator cancelled(&medium, Ray(point3f(0), vec3f(1)), 1, true, &cancel);
    expect_false(bool(cancelled.Next()));
  }
  test_that("constant equality, tag overrides, and sparse NanoVDB transitions are deterministic") {
    auto description = medium_description();
    description["haze"] = true;
    description["haze_density_threshold"] = 1.;
    for (double scale : {0., .5, 1., 2.}) {
      description["density_scale"] = scale;
      Medium medium(description);
      MediumHazeIterator iterator(&medium, Ray(point3f(0), vec3f(0, 0, 1)), 2, true);
      auto interval = iterator.Next();
      expect_true(interval->integrate == (scale < 1));
      expect_true(interval->t_max == 2);
      medium.haze = false;
      MediumHazeIterator tagged(&medium, Ray(point3f(0), vec3f(0, 0, 1)), 2, true);
      expect_false(tagged.Next()->integrate);
    }
    Rcpp::Function test_path = Rcpp::Environment::namespace_env("testthat")["test_path"];
    description["density_scale"] = 1.;
    description["haze_density_threshold"] = .5;
    description["filename"] = test_path("fixtures", "volumes", "tiles.nvdb");
    description["density_grid"] = "density";
    description["temperature_grid"] = R_NilValue;
    NanoVDBMedium medium(description);
    MediumHazeIterator iterator(&medium, Ray(point3f(3, 3, -3), vec3f(0, 0, 1)), 14, true);
    auto air = iterator.Next(), dense = iterator.Next(), tail = iterator.Next();
    expect_true((air->integrate && !dense->integrate && tail->integrate));
    expect_true(std::abs(air->t_max - 2.5) < 1e-10);
    expect_true(std::abs(dense->t_max - 10.5) < 1e-10);
    expect_true(tail->t_max == 14);
  }
  test_that("disabled haze and an absent cutoff bypass heterogeneous density traversal") {
    class CountingMedium : public Medium {
    public:
      CountingMedium() : Medium(medium_description(0)) {}
      bool IsHomogeneous() const override { return false; }
      DensityIndexRay DensityRay(const Ray &) const override { ++queries; return {}; }
      mutable int queries = 0;
    } medium;
    medium.haze = false;
    medium.haze_density_threshold = .05;
    Ray ray(point3f(0), vec3f(0, 0, 1));
    MediumHazeIterator off(&medium, ray, 2, true);
    auto interval = off.Next();
    expect_true((!interval->integrate && interval->t_min == 0 && interval->t_max == 2));
    expect_false(bool(off.Next()));
    expect_true(medium.queries == 0);

    medium.haze = true;
    MediumHazeIterator global_off(&medium, ray, 2, false);
    expect_false(global_off.Next()->integrate);
    expect_true(medium.queries == 0);

    medium.haze_density_threshold = 0;
    MediumHazeIterator full(&medium, ray, 2, true);
    expect_true(full.Next()->integrate);
    expect_false(bool(full.Next()));
    expect_true(medium.queries == 0);
  }
}
context("Deterministic volume picking") {
  test_that("thin media reach the background while visible matter takes priority") {
    PickingScene s;
    s.Add(std::make_shared<Medium>(medium_description(.01)), 1);
    auto texture = std::make_shared<constant_texture>(point3f(1));
    auto light = std::make_shared<diffuse_light>(texture, 1, false);
    s.world.add(std::make_shared<InfiniteAreaLight>(16, 8, 100, point3f(0), texture,
                                                   light, &s.identity, &s.identity, false));
    Ray ray(point3f(0, 0, -3), vec3f(0, 0, 1));
    auto target = PickRay(ray, &s.world, &s.scene);
    expect_true(bool(target));
    expect_true(target->background);
    expect_false(target->volume);
    s.Add(std::make_shared<Medium>(medium_description(1)), .5);
    target = PickRay(ray, &s.world, &s.scene);
    expect_true(target->volume);
    expect_false(target->background);
    s.world.add(std::make_shared<sphere>(.75, s.mat, nullptr, nullptr,
                                        &s.identity, &s.identity, false));
    target = PickRay(ray, &s.world, &s.scene);
    expect_false(target->background);
    expect_false(target->volume);
    expect_true(target->p[2] == Approx(-.75));
    expect_false(bool(PickRay(ray, &s.world, &s.scene, .15, [] { return true; })));
  }
  test_that("homogeneous RGB opacity uses world distance and repeats exactly") {
    PickingScene s;
    auto description = medium_description();
    description["sigma_s"] = Rcpp::NumericVector::create(1, 0, 0);
    s.Add(std::make_shared<Medium>(description), 1);
    Ray ray(point3f(0, 0, -3), vec3f(0, 0, 7));
    auto pick = PickRay(ray, &s.world, &s.scene);
    expect_true(bool(pick));
    expect_true(pick->volume);
    expect_true(std::abs(pick->p[2] - (-1 - std::log(1 - 3 * .15))) < 1e-5);
    for (int i = 0; i < 8; ++i) {
      auto repeated = PickRay(ray, &s.world, &s.scene);
      expect_true(repeated->p[2] == pick->p[2]);
    }
    auto inside = PickRay(Ray(point3f(0), vec3f(0, 0, 1)), &s.world, &s.scene);
    expect_true(std::abs(inside->p[2] + std::log(1 - 3 * .15)) < 1e-5);
    // A single extinguished channel cannot reach 50% mean RGB opacity.
    expect_false(bool(PickRay(ray, &s.world, &s.scene, .5)));
  }
  test_that("thin media fall through to surfaces and vacuum cavities restore their enclosing medium") {
    PickingScene thin;
    thin.Add(std::make_shared<Medium>(medium_description(.01)), 1);
    Ray ray(point3f(0, 0, -3), vec3f(0, 0, 1));
    expect_false(bool(PickRay(ray, &thin.world, &thin.scene)));
    thin.world.add(std::make_shared<sphere>(.25, thin.mat, nullptr, nullptr,
                                          &thin.identity, &thin.identity, false));
    auto surface = PickRay(ray, &thin.world, &thin.scene);
    expect_true(bool(surface));
    expect_false(surface->volume);
    expect_true(surface->p[2] == Approx(-.25));
    PickingScene nested;
    nested.Add(std::make_shared<Medium>(medium_description(.1)), 2);
    nested.Add(std::make_shared<Medium>(medium_description(0)), 1);
    auto pick = PickRay(ray, &nested.world, &nested.scene);
    expect_true(pick->volume);
    expect_true(std::abs(pick->p[2] - (1 + (-std::log(.85) - .1) / .1)) < 1e-5);
    PickingScene visible;
    visible.Add(std::make_shared<Medium>(medium_description()), 1, true);
    auto container = PickRay(ray, &visible.world, &visible.scene);
    expect_false(container->volume);
    expect_true(container->p[2] == Approx(-1));
  }
  test_that("grid ramps and isolated voxels are integrated between interpolation knots") {
    PickingScene ramp;
    ramp.Add(std::make_shared<GridMedium>(
                 grid_description(Rcpp::NumericVector::create(0, 2), 1, 1, 2)), 1);
    Ray ray(point3f(0, 0, -3), vec3f(0, 0, 1));
    auto pick = PickRay(ray, &ramp.world, &ramp.scene);
    expect_true(std::abs(pick->p[2] - (std::sqrt(-std::log(.85)) - .5)) < 1e-5);
    PickingScene spike;
    Rcpp::NumericVector density(128);
    density[64] = 100;
    spike.Add(std::make_shared<GridMedium>(grid_description(density, 1, 1, 128)), 1);
    pick = PickRay(ray, &spike.world, &spike.scene);
    double start = -.0078125, width = 2.0 / 128;
    expect_true(std::abs(pick->p[2] - (start + std::sqrt(-2 * std::log(.85) * width / 100))) < 1e-5);
    PickingScene constant;
    constant.Add(std::make_shared<GridMedium>(
                     grid_description(Rcpp::NumericVector(8, 1.0), 2, 2, 2)), 1);
    pick = PickRay(ray, &constant.world, &constant.scene);
    expect_true(std::abs(pick->p[2] - (-1 - std::log(.85))) < 1e-5);
    PickingScene empty;
    empty.Add(std::make_shared<GridMedium>(
                  grid_description(Rcpp::NumericVector(8), 2, 2, 2)), 1);
    expect_false(bool(PickRay(ray, &empty.world, &empty.scene)));
  }
  test_that("diagonal trilinear fields use cubic density integrals") {
    PickingScene s;
    Rcpp::NumericVector density(8);
    density[7] = 1;
    s.Add(std::make_shared<GridMedium>(grid_description(density, 2, 2, 2)), 1);
    auto pick = PickRay(Ray(point3f(-2), vec3f(1)), &s.world, &s.scene);
    // On the diagonal inside [-.5,.5], density = (x+.5)^3.
    double expected = std::pow(-4 * std::log(.85) / std::sqrt(3.0), .25) - .5;
    expect_true(std::abs(pick->p[0] - expected) < 1e-5);
  }
  test_that("NanoVDB picking respects sparse tiles and interpolation boundaries") {
    Rcpp::Function test_path = Rcpp::Environment::namespace_env("testthat")["test_path"];
    Rcpp::List description = medium_description();
    description["filename"] = test_path("fixtures", "volumes", "tiles.nvdb");
    description["density_grid"] = "density";
    description["temperature_grid"] = R_NilValue;
    PickingScene s;
    s.Add(std::make_shared<NanoVDBMedium>(description), 10);
    auto pick = PickRay(Ray(point3f(3, 3, -12), vec3f(0, 0, 1)), &s.world, &s.scene);
    expect_true(bool(pick));
    expect_true(std::abs(pick->p[2] - (-1 + std::sqrt(-2 * std::log(.85)))) < 1e-5);
    pick = PickRay(Ray(point3f(3), vec3f(0, 0, 1)), &s.world, &s.scene);
    expect_true(std::abs(pick->p[2] - (3 - std::log(.85))) < 1e-5);
    expect_false(bool(PickRay(Ray(point3f(-3, -3, -12), vec3f(0, 0, 1)), &s.world, &s.scene)));
  }
  test_that("majorant partition and magnitude do not move the opacity target") {
    class RepartitionedMedium : public Medium {
      GridMedium field;
      MajorantGrid grid;
    public:
      RepartitionedMedium(int resolution, Float bound)
          : Medium(medium_description()),
            field(grid_description(Rcpp::NumericVector::create(0, 2), 1, 1, 2)) {
        grid.lo = point3f(-1);
        grid.hi = point3f(1);
        grid.resolution = resolution;
        grid.density.assign(resolution * resolution * resolution, bound);
      }
      bool IsHomogeneous() const override { return false; }
      Float Density(const point3f &p) const override { return field.Density(p); }
      DensityIndexRay DensityRay(const Ray &r) const override { return field.DensityRay(r); }
      RayMajorantIterator SampleRay(const Ray &r, double end) const override {
        return RayMajorantIterator(r, end, point3f(1), &grid);
      }
    };
    for (int resolution : {1, 5, 32}) {
      PickingScene s;
      s.Add(std::make_shared<RepartitionedMedium>(resolution, resolution * 10), 1);
      auto pick = PickRay(Ray(point3f(0, 0, -3), vec3f(0, 0, 1)), &s.world, &s.scene);
      expect_true(std::abs(pick->p[2] - (std::sqrt(-std::log(.85)) - .5)) < 1e-5);
    }
  }
  test_that("absorption inside glass adds to explicit medium extinction") {
    PickingScene s;
    s.mat = std::make_shared<dielectric>(point3f(1), 1.5, point3f(.5), 0);
    s.Add(std::make_shared<Medium>(medium_description(.5)), 1, true);
    auto pick = PickRay(Ray(point3f(0), vec3f(0, 0, 1)), &s.world, &s.scene);
    expect_true(pick->volume);
    expect_true(std::abs(pick->p[2] + std::log(.85)) < 1e-5);
    pick = PickRay(Ray(point3f(0, 0, -3), vec3f(0, 0, 1)), &s.world, &s.scene);
    expect_false(pick->volume);
    expect_true(pick->p[2] == Approx(-1));
  }
  test_that("picking preserves animated instance transforms and can cancel voxel traversal") {
    PickingScene s;
    auto medium = std::make_shared<GridMedium>(
        grid_description(Rcpp::NumericVector::create(0, 2), 1, 1, 2));
    s.Add(medium, 1);
    Transform transform = Translate(vec3f(4, 0, 0)) * Scale(1, 1, 3);
    Transform inverse = Inverse(transform), end = Translate(vec3f(0, 2, 0));
    AnimatedTransform motion(&s.identity, 0, &end, 1);
    hitable_list lights, world;
    std::shared_ptr<hitable> placed = std::make_shared<instance>(
        &s.world, &transform, &inverse, &lights, 0);
    auto animated = std::make_shared<AnimatedHitable>(placed, motion);
    world.add(animated);
    VolumeScene scene;
    scene.boundaries.add(animated);
    scene.Finish(0, 1);
    for (Float time : {0.f, .5f, 1.f}) {
      auto pick = PickRay(Ray(point3f(4, 2 * time, -6), vec3f(0, 0, 1), time), &world, &scene);
      expect_true(std::abs(pick->p[2] - (std::sqrt(-3 * std::log(.85)) - 1.5)) < 1e-5);
    }
    int polls = 0;
    auto cancel = [&] { return ++polls >= 12; };
    auto pick = PickRay(Ray(point3f(0, 0, -3), vec3f(0, 0, 1)), &s.world, &s.scene, .99, cancel);
    expect_false(bool(pick));
    expect_true(polls >= 12);
    polls = 12;
    expect_false(bool(PickRay(Ray(point3f(0), vec3f(0, 0, 1)), &s.world, &s.scene, .15, cancel)));
  }
}
context("Participating media geometry and sampling") {
  test_that("NEE scalar and vector requests use distinct sampler coordinates") {
    random_gen rng1(7), rng2(7);
    SobolBlueNoiseSampler a(rng1), b(rng2);
    a.independent_dimensions = b.independent_dimensions = true;
    a.StartPixel(127, 127);
    b.StartPixel(127, 127);
    bool same = true, finite = true;
    for (int sample = 0; sample < 32; ++sample) {
      for (int dim = 0; dim < 300; ++dim) {
        Float x = a.Get1D(), y = a.Get1D();
        vec2f xy = b.Get2D();
        same &= x == xy[0] && y == xy[1];
        finite &= x >= 0 && x < 1 && y >= 0 && y < 1;
      }
      a.StartNextSample();
      b.StartNextSample();
    }
    expect_true(same);
    expect_true(finite);
  }
  test_that("cancelled containment traversal returns promptly") {
    VolumeScene scene;
    Transform identity;
    auto mat = std::make_shared<lambertian>(std::make_shared<constant_texture>(point3f(1)));
    auto geometry = std::make_shared<sphere>(1, mat, nullptr, nullptr, &identity, &identity, false);
    auto medium = std::make_shared<Medium>(medium_description());
    scene.boundaries.add(std::make_shared<MediumBoundary>(geometry, medium, identity, false,
                                                         scene.NextBoundaryId()));
    scene.Finish(0, 1);
    expect_true(scene.InitialState(Ray(point3f(0), vec3f(0, 0, 1)), nullptr).media.size() == 1);
    std::atomic<bool> cancel(true);
    auto state = scene.InitialState(Ray(point3f(0), vec3f(0, 0, 1)), &cancel);
    expect_true(state.media.empty());
  }
  test_that("boundary ID ranges reject exhaustion without wrapping or consuming IDs") {
    VolumeScene scene;
    const uint64_t maximum = std::numeric_limits<uint64_t>::max();
    expect_true(scene.NextBoundaryId() == 1);
    expect_true(scene.ReserveBoundaryIds(maximum - 2) == 1);
    bool rejected = false;
    try {
      scene.ReserveBoundaryIds(2);
    } catch (const std::overflow_error &) {
      rejected = true;
    }
    expect_true(rejected);
    expect_true(scene.NextBoundaryId() == maximum);
    expect_true(scene.ReserveBoundaryIds(0) == maximum);
    rejected = false;
    try {
      scene.NextBoundaryId();
    } catch (const std::overflow_error &) {
      rejected = true;
    }
    expect_true(rejected);
    expect_true(scene.BoundaryCount() == maximum);
  }
  test_that("shared nested boundaries have distinct IDs in render and containment traversal") {
    Transform identity, larger = Scale(2, 2, 2), larger_inverse = Inverse(larger);
    Transform side = Translate(vec3f(0, 8, 0)), side_inverse = Inverse(side);
    auto mat = std::make_shared<lambertian>(std::make_shared<constant_texture>(point3f(.5)));
    auto medium = std::make_shared<Medium>(medium_description());
    VolumeScene leaf, branch, scene;
    auto geometry = std::make_shared<sphere>(1, mat, nullptr, nullptr, &identity, &identity, false);
    auto boundary = std::make_shared<MediumBoundary>(geometry, medium, identity, false,
                                                    leaf.NextBoundaryId());
    leaf.boundaries.add(boundary);
    leaf.Finish(0, 1);
    hitable_list leaf_world, branch_world, world, lights;
    leaf_world.add(boundary);
    auto side_geometry = std::make_shared<sphere>(1, mat, nullptr, nullptr, &side, &side_inverse, false);
    auto side_boundary = std::make_shared<MediumBoundary>(side_geometry, medium, side, false,
                                                         branch.NextBoundaryId());
    branch.boundaries.add(side_boundary);
    branch_world.add(side_boundary);
    for (int i = 0; i < 2; ++i) {
      Transform *transform = i == 0 ? &larger : &identity;
      Transform *inverse = i == 0 ? &larger_inverse : &identity;
      uint64_t offset = branch.ReserveBoundaryIds(leaf.BoundaryCount());
      branch_world.add(std::make_shared<instance>(&leaf_world, transform, inverse, &lights, offset));
      branch.boundaries.add(std::make_shared<instance>(leaf.boundary_bvh.get(), transform,
                                                       inverse, &lights, offset));
    }
    branch.Finish(0, 1);
    Transform far_transform = Translate(vec3f(30, 0, 0)), far_inverse = Inverse(far_transform);
    auto far_geometry = std::make_shared<sphere>(1, mat, nullptr, nullptr, &far_transform, &far_inverse, false);
    auto far_boundary = std::make_shared<MediumBoundary>(far_geometry, medium, far_transform, false,
                                                        scene.NextBoundaryId());
    scene.boundaries.add(far_boundary);
    world.add(far_boundary);
    Transform placement[2] = {Translate(vec3f(-5, 0, 0)), Translate(vec3f(5, 0, 0))};
    Transform inverse[2] = {Inverse(placement[0]), Inverse(placement[1])};
    Transform motion_end = Translate(vec3f(0, 2, 0));
    AnimatedTransform motion(&identity, 0, &motion_end, 1);
    for (int i = 0; i < 2; ++i) {
      uint64_t offset = scene.ReserveBoundaryIds(branch.BoundaryCount());
      std::shared_ptr<hitable> render_instance = std::make_shared<instance>(
          &branch_world, &placement[i], &inverse[i], &lights, offset);
      std::shared_ptr<hitable> probe_instance = std::make_shared<instance>(
          branch.boundary_bvh.get(), &placement[i], &inverse[i], &lights, offset);
      world.add(std::make_shared<AnimatedHitable>(render_instance, motion));
      scene.boundaries.add(std::make_shared<AnimatedHitable>(probe_instance, motion));
    }
    scene.Finish(0, 1);
    expect_true(scene.BoundaryCount() == 7);
    random_gen rng(17);
    RandomSampler sampler(rng);
    for (Float time : {0.f, .5f, 1.f}) {
      std::set<uint64_t> ids;
      auto direct = scene.InitialState(Ray(point3f(30, 0, 0), vec3f(0, 0, 1), time), nullptr);
      expect_true(direct.media.size() == 1);
      ids.insert(direct.media.at(0).boundary_id);
      for (int i = 0; i < 2; ++i) {
        point3f center(i == 0 ? -5 : 5, 2 * time, 0);
        Ray ray(center, vec3f(0, 0, 1), time);
        auto state = scene.InitialState(ray, nullptr);
        expect_true(state.media.size() == 2);
        expect_true(state.media.at(0).boundary == boundary.get());
        expect_true(state.media.at(1).boundary == boundary.get());
        for (const auto &entry : state.media)
          ids.insert(entry.boundary_id);
        while (!state.media.empty()) {
          state.SetRay(ray);
          hit_record render_hit, probe_hit, sampled_hit;
          expect_true(world.hit(ray, 0, MaxT, render_hit, rng));
          expect_true(world.hit(ray, 0, MaxT, sampled_hit, &sampler));
          expect_true(scene.boundary_bvh->hit(ray, 0, MaxT, probe_hit, rng));
          expect_true(render_hit.boundary_id == state.media.back().boundary_id);
          expect_true(render_hit.boundary_id == probe_hit.boundary_id);
          expect_true(render_hit.boundary_id == sampled_hit.boundary_id);
          state.Cross(render_hit, ray.d);
          ray = Ray(OffsetMediumOrigin(render_hit, ray.d), ray.d, time);
        }
        // Direct boundaries inside an instanced scene use the same ID range.
        auto side_state = scene.InitialState(Ray(center + vec3f(0, 8, 0), vec3f(0, 0, 1), time), nullptr);
        expect_true(side_state.media.size() == 1);
        ids.insert(side_state.media.at(0).boundary_id);
      }
      expect_true(ids.size() == 7);
      expect_true(*ids.begin() == 1);
      expect_true(*ids.rbegin() == 7);
    }
    // An ordinary surface must keep the zero sentinel even with an ID offset.
    instance ordinary(geometry.get(), &identity, &identity, &lights, 100);
    Ray ray(point3f(0, 0, -3), vec3f(0, 0, 1));
    ray.segment_absorption = true;
    hit_record h;
    expect_true(ordinary.hit(ray, 0, MaxT, h, rng));
    expect_true(h.boundary_id == 0);
    expect_true(ordinary.hit(ray, 0, MaxT, h, &sampler));
    expect_true(h.boundary_id == 0);
  }
  test_that("cancellation stops a long null-event tracking loop") {
    class NullOnlyMedium final : public Medium {
    public:
      mutable std::atomic<size_t> calls{0};
      NullOnlyMedium() : Medium(medium_description(1e9)) {}
      MediumProperties SamplePoint(const point3f &) const override {
        ++calls;
        MediumProperties properties;
        return properties;
      }
    };
    auto medium = std::make_shared<NullOnlyMedium>();
    Transform identity;
    auto mat = std::make_shared<lambertian>(std::make_shared<constant_texture>(point3f(1)));
    auto geometry = std::make_shared<sphere>(1, mat, nullptr, nullptr, &identity, &identity, false);
    auto scene = std::make_shared<VolumeScene>();
    auto boundary = std::make_shared<MediumBoundary>(geometry, medium, identity, false,
                                                    scene->NextBoundaryId());
    scene->boundaries.add(boundary);
    scene->Finish(0, 1);
    // Skip the separate alpha walk to exercise radiance tracking directly.
    hitable_list world, lights;
    world.add(boundary);
    lights.volume_scene = scene;
    random_gen rng(19);
    RandomSampler sampler(rng);
    Float transparency;
    point3f radiance, albedo;
    normal3f normal;
    std::atomic<bool> cancel(false);
    auto start = std::chrono::steady_clock::now();
    std::thread stop([&] {
      while (medium->calls < 1000 &&
             std::chrono::steady_clock::now() - start < std::chrono::seconds(1))
        std::this_thread::yield();
      cancel = true;
    });
    try {
      color_volume(Ray(point3f(0), vec3f(0, 0, 1)), &world, &lights, 10, 5, rng, &sampler,
                   transparency, radiance, normal, albedo, &cancel);
    } catch (...) {
      cancel = true;
      stop.join();
      throw;
    }
    stop.join();
    expect_true(medium->calls >= 1000);
    expect_true(std::chrono::steady_clock::now() - start < std::chrono::seconds(2));
  }
  test_that("HG samples preserve the requested forward mean cosine") {
    for (Float g : {-0.7f, 0.f, 0.7f}) {
      HGPhaseFunction phase(g);
      double mean = 0;
      bool valid = true;
      constexpr int n = 32768;
      for (int i = 0; i < n; ++i) {
        auto sample = phase.Sample(vec3f(0, 0, -1), (i + .5f) / n, 0.371f);
        mean += sample.wi[2];
        valid &= std::isfinite(sample.pdf) && sample.pdf > 0;
      }
      expect_true(valid);
      expect_true(std::abs(mean / n - g) < 0.0002);
    }
  }
  test_that("DDA covers a transformed world-distance interval once") {
    MajorantGrid grid;
    grid.lo = point3f(-1);
    grid.hi = point3f(1);
    grid.resolution = 16;
    grid.density.assign(4096, 1);
    RayMajorantIterator it(Ray(point3f(-2, 0, 0), vec3f(.5, 0, 0)), 10, point3f(2), &grid);
    double length = 0, previous = 2;
    int count = 0;
    while (auto s = it.Next()) {
      expect_true(std::abs(s->t_min - previous) < 1e-8);
      expect_true(s->t_max > s->t_min);
      length += s->t_max - s->t_min;
      previous = s->t_max;
      ++count;
    }
    expect_true(length == Approx(4));
    expect_true(count == 16);
  }
  test_that("cell-centered fields interpolate channels in R array order") {
    SampledField f;
    f.dims = {2, 2, 2};
    f.channels = 3;
    for (int channel = 0; channel < 3; ++channel)
      for (int z = 0; z < 2; ++z)
        for (int y = 0; y < 2; ++y)
          for (int x = 0; x < 2; ++x)
            f.values.push_back(channel * 10 + x + 2 * y + 4 * z);
    expect_true(f.Lookup(point3f(.5), 2) == Approx(23.5));
    expect_true(f.Lookup(point3f(0), 1) == Approx(10));
    expect_true(f.Lookup(point3f(1), 0) == Approx(7));
  }
  test_that("box edge contacts are not medium crossings") {
    Transform identity;
    auto material = std::make_shared<lambertian>(std::make_shared<constant_texture>(point3f(.5)));
    auto geometry = std::make_shared<box>(vec3f(-1), vec3f(1), material, nullptr, nullptr,
                                          &identity, &identity, false);
    auto medium = std::make_shared<Medium>(medium_description());
    VolumeScene scene;
    MediumBoundary boundary(geometry, medium, identity, false, scene.NextBoundaryId());
    Ray contact(point3f(-2, 0, 0), vec3f(1, 0, 1));
    contact.segment_absorption = true;
    random_gen rng(4);
    hit_record h;
    expect_false(boundary.hit(contact, 0, 100, h, rng));
    Ray crossing(point3f(-2, 0, -2), vec3f(1, 0, 1));
    crossing.segment_absorption = true;
    expect_true(boundary.hit(crossing, 0, 100, h, rng));
    expect_true(dot(crossing.d, h.geometric_normal) < 0);
  }
  test_that("invisible boxes retain grazing entries before their exits") {
    Transform transform = Scale(-1, 1, -1) * Translate(vec3f(5.889312909, 11.514999986, 0));
    Transform inverse = Inverse(transform);
    auto material = std::make_shared<lambertian>(std::make_shared<constant_texture>(point3f(.5)));
    vec3f half_size(16.956099017, 5, 23.314636148);
    auto geometry = std::make_shared<box>(-half_size, half_size, material, nullptr, nullptr,
                                         &transform, &inverse, false);
    auto medium = std::make_shared<Medium>(medium_description(0));
    VolumeScene scene;
    auto boundary = std::make_shared<MediumBoundary>(geometry, medium, transform, false,
                                                    scene.NextBoundaryId());
    scene.boundaries.add(boundary);
    scene.Finish(0, 1);
    // Primary ray from the geographic cloud example, crossing an x/z edge.
    Ray ray(point3f(48.439773559570312, 53.806549072265625, -76.688796997070312),
            vec3f(-0.49673172831535339, -0.5, 0.70940649509429932));
    auto state = scene.InitialState(ray, nullptr);
    expect_true(state.media.empty());
    state.SetRay(ray);
    random_gen rng(4);
    RandomSampler sampler(rng);
    hit_record random_hit, sampled_hit;
    expect_true(boundary->hit(ray, 0, MaxT, random_hit, rng));
    expect_true(boundary->hit(ray, 0, MaxT, sampled_hit, &sampler));
    expect_true(dot(ray.d, random_hit.geometric_normal) < 0);
    expect_true(random_hit.t == Approx(75.23777));
    expect_true(random_hit.t == sampled_hit.t);
    state.Cross(random_hit, ray.d);
    Ray continued(OffsetMediumOrigin(random_hit, ray.d), ray.d);
    state.SetRay(continued);
    expect_true(boundary->hit(continued, 0, MaxT, random_hit, rng));
    expect_true(dot(continued.d, random_hit.geometric_normal) > 0);
    state.Cross(random_hit, continued.d);
    expect_true(state.media.empty());
  }
  test_that("camera containment resolves points close to a curved boundary") {
    Transform transform = Translate(vec3f(278, 220, 250)), inverse = Inverse(transform);
    auto material = std::make_shared<lambertian>(std::make_shared<constant_texture>(point3f(.5)));
    auto geometry =
        std::make_shared<sphere>(110, material, nullptr, nullptr, &transform, &inverse, false);
    auto medium = std::make_shared<Medium>(medium_description());
    VolumeScene scene;
    scene.boundaries.add(std::make_shared<MediumBoundary>(geometry, medium, transform, false,
                                                         scene.NextBoundaryId()));
    scene.Finish(0, 1);
    auto outside = scene.InitialState(
        Ray(point3f(329.332824707, 137.041244507, 199.178741455), vec3f(0, 0, 1)), nullptr);
    auto inside = scene.InitialState(Ray(point3f(278, 220, 250), vec3f(0, 0, 1)), nullptr);
    expect_true(outside.media.empty());
    expect_true(inside.media.size() == 1);
  }
  test_that("rays starting on a box face count the initial crossing once") {
    Transform transform = Translate(vec3f(0, 10.5, 0)), inverse = Inverse(transform);
    auto material = std::make_shared<lambertian>(std::make_shared<constant_texture>(point3f(.5)));
    auto geometry = std::make_shared<box>(vec3f(-50, -6.5, -38), vec3f(50, 6.5, 38),
                                         material, nullptr, nullptr, &transform, &inverse, false);
    auto medium = std::make_shared<Medium>(medium_description());
    VolumeScene scene;
    auto boundary = std::make_shared<MediumBoundary>(geometry, medium, transform, false,
                                                    scene.NextBoundaryId());
    scene.boundaries.add(boundary);
    scene.Finish(0, 1);
    random_gen rng(4);
    for (int face : {-1, 1}) {
      for (int direction : {-1, 1}) {
        Ray ray(point3f(0, face > 0 ? 17 : 4, 0), vec3f(0, direction, 0));
        auto state = scene.InitialState(ray, nullptr);
        bool outgoing = direction == face;
        expect_true(state.media.size() == (outgoing ? 1 : 0));
        state.SetRay(ray);
        hit_record h;
        expect_true(boundary->hit(ray, 0, MaxT, h, rng));
        expect_true(h.t == 0);
        state.Cross(h, ray.d);
        expect_true(state.media.size() == (outgoing ? 0 : 1));
        if (!outgoing) {
          Ray continued(OffsetMediumOrigin(h, ray.d), ray.d);
          state.SetRay(continued);
          expect_true(boundary->hit(continued, 0, MaxT, h, rng));
          expect_true(h.t > 0);
          state.Cross(h, continued.d);
          expect_true(state.media.empty());
        }
      }
      auto tangent = scene.InitialState(
          Ray(point3f(0, face > 0 ? 17 : 4, 0), vec3f(0, 0, 1)), nullptr);
      expect_true(tangent.media.empty());
    }
  }
  test_that("distant sphere light PDFs agree with their actual intersections") {
    for (Transform transform : {Translate(vec3f(1000, 2000, 3000)),
                                Translate(vec3f(1000, 2000, 3000)) * Scale(.7, 1.8, 1.2)}) {
      Transform inverse = Inverse(transform);
      auto material = std::make_shared<lambertian>(std::make_shared<constant_texture>(point3f(.5)));
      auto geometry = std::make_shared<sphere>(50, material, nullptr, nullptr,
                                              &transform, &inverse, false);
      hitable_list lights;
      lights.add(geometry);
      VolumeLightSampler light_sampler(lights);
      random_gen rng(31);
      RandomSampler sampler(rng);
      point3f p(22, 17, 7);
      int missing_pdf = 0, hits = 0;
      for (int i = 0; i < 100000; ++i) {
        vec3f wi = unit_vector(light_sampler.Sample(p, &sampler, 0));
        Ray ray(p, wi);
        ray.segment_absorption = true;
        hit_record h;
        if (geometry->hit(ray, 0, MaxT, h, rng)) {
          ++hits;
          if (!(light_sampler.Pdf(p, wi, rng, 0) > 0))
            ++missing_pdf;
        }
      }
      expect_true(hits > 99000);
      expect_true(missing_pdf == 0);
    }
  }
}

context("Conservative bounds around medium boundaries") {
  test_that("slab rounding preserves clipped hits and parallel rays") {
    aabb bounds(point3f(-1), point3f(1));
    random_gen rng(2026);
    RandomSampler sampler(rng);
    auto check = [&](const Ray &ray, Float lo, Float hi, bool expected) {
      expect_true(bounds.hit(ray, lo, hi, rng) == expected);
      expect_true(bounds.hit(ray, lo, hi, &sampler) == expected);
#ifdef RAYSIMD
      BBox4 packed(bounds, bounds, bounds, bounds);
      IVec4 hits;
      FVec4 entries;
      rayBBoxIntersect4(RayBBox4(ray), packed, lo, hi, hits, entries);
      for (int lane = 0; lane < 4; ++lane) {
        expect_true((hits.xyzw[lane] != 0) == expected);
        if (expected) {
          expect_true(entries.xyzw[lane] >= lo);
          expect_true(entries.xyzw[lane] <= hi);
        }
      }
#endif
    };
    for (int axis = 0; axis < 3; ++axis) {
      for (Float sign : {Float(-1), Float(1)}) {
        point3f origin(0);
        vec3f direction(0);
        origin.e[axis] = -2 * sign;
        direction.e[axis] = sign;
        check(Ray(origin, direction), 0, 1, true);
        check(Ray(origin, -direction), -1, -1, true);
        check(Ray(origin, -direction), 0, MaxT, false);

        // Parallel rays on either face generate 0 * infinity in one slab.
        // Signed zero must work as well as ordinary zero, and rays just
        // outside the parallel slab must still miss.
        for (Float face : {Float(-1), Float(1)}) {
          origin.e[(axis + 1) % 3] = face;
          direction.e[(axis + 1) % 3] = std::copysign(Float(0), face);
          check(Ray(origin, direction), 0, 1, true);
          origin.e[(axis + 1) % 3] = face * 1.01;
          check(Ray(origin, direction), 0, MaxT, false);
        }
      }
    }
  }

  test_that("a surface just behind a cloud face cannot hide its entry") {
    Transform transform = Translate(vec3f(0, 13.4, 63)), inverse = Inverse(transform), identity;
    auto mat = std::make_shared<lambertian>(std::make_shared<constant_texture>(point3f(.5)));
    auto geometry = std::make_shared<box>(vec3f(-45, -.9, -12), vec3f(45, .9, 12),
                                         mat, nullptr, nullptr, &transform, &inverse, false);
    auto medium = std::make_shared<Medium>(medium_description(0));
    VolumeScene scene;
    auto boundary = std::make_shared<MediumBoundary>(geometry, medium, transform, false,
                                                    scene.NextBoundaryId());
    for (int side : {-1, 1}) {
      const Float face = side < 0 ? Float(12.5) : Float(13.4) + Float(.9);
      auto surface = std::make_shared<xz_rect>(
          -1, 1, 62, 64, std::nextafter(face, Float(13.4)), mat, nullptr, nullptr,
          &identity, &identity, false);
      hitable_list reference;
      reference.add(boundary);
      reference.add(surface);
      // Separate leaves force traversal to compare the cloud's entry bound
      // with the surface hit. A padded reciprocal used to reverse their order.
      BVHAggregate world(reference.objects, 0, 1, 1, true);
      Ray ray(point3f(0, side < 0 ? 0 : 26.8, 63), vec3f(0, -side, 0));
      ray.segment_absorption = true;
      random_gen rng(2026);
      RandomSampler sampler(rng);
      hit_record expected, actual;
      expect_true(reference.hit(ray, 0, MaxT, expected, rng));
      expect_true(expected.medium_boundary == boundary.get());
      expect_true(world.hit(ray, 0, MaxT, actual, rng));
      expect_true(actual.medium_boundary == boundary.get());
      expect_true(actual.t == expected.t);
      expect_true(world.hit(ray, 0, MaxT, actual, &sampler));
      expect_true(actual.medium_boundary == boundary.get());
      expect_true(actual.t == expected.t);
      expect_true(world.hit(ray, 0, expected.t, actual, rng));
      expect_true(actual.medium_boundary == boundary.get());
      expect_true(world.HitP(ray, 0, expected.t, rng));
      expect_true(world.HitP(ray, 0, expected.t, &sampler));
    }
  }
}


namespace {
class OrderedShadowWorld : public hitable_list {
public:
  OpaqueShadowType ShadowType() const { return OpaqueShadowType::Unsupported; }
};
}

context("Restricted opaque shadow connections") {
  test_that("visibility preserves radiance and both random streams, including near-light ties") {
    Transform identity;
    auto matte = std::make_shared<lambertian>(std::make_shared<constant_texture>(point3f(.6)));
    auto emission = std::make_shared<diffuse_light>(
        std::make_shared<constant_texture>(point3f(1, .8, .4)), 3, false);
    auto ground = std::make_shared<xz_rect>(-20, 20, -20, 20, 0, matte, nullptr, nullptr,
                                          &identity, &identity, false);
    auto light = std::make_shared<xz_rect>(-3, 3, -3, 3, 5, emission, nullptr, nullptr,
                                         &identity, &identity, true);
    for (Float height : {Float(2), Float(5), std::nextafter(Float(5), Float(0)),
                         std::nextafter(Float(5), Float(10))}) {
      auto blocker = std::make_shared<xz_rect>(-.7, .7, -.7, .7, height, matte, nullptr,
                                              nullptr, &identity, &identity, false);
      std::vector<std::shared_ptr<hitable>> objects{ground, blocker, light};
      auto bvh = std::make_shared<BVHAggregate>(objects, 0, 1, 1, true);
      hitable_list fast, lights;
      OrderedShadowWorld ordered;
      fast.add(bvh); ordered.add(bvh); lights.add(light);
      lights.volume_scene = std::make_shared<VolumeScene>();
      lights.volume_scene->light_sampler = std::make_shared<VolumeLightSampler>(lights);
      expect_true(fast.ShadowType() == OpaqueShadowType::Mixed);
      bool same = true;
      for (int seed = 0; seed < 512; ++seed) {
        random_gen a(seed), b(seed);
        RandomSampler sa(a), sb(b);
        sa.independent_dimensions = sb.independent_dimensions = true;
        Ray ray(point3f(8, 3, -8), unit_vector(vec3f(-8 + (seed % 17) * .15,
                                                     -3, 8 + (seed % 13) * .2)));
        Float ta, tb;
        point3f ra, rb, aa, ab;
        normal3f na, nb;
        color_volume(ray, &fast, &lights, 8, 5, a, &sa, ta, ra, na, aa, nullptr);
        color_volume(ray, &ordered, &lights, 8, 5, b, &sb, tb, rb, nb, ab, nullptr);
        same &= ta == tb && a.unif_rand() == b.unif_rand() && sa.Get1D() == sb.Get1D();
        for (int c = 0; c < 3; ++c) same &= ra[c] == rb[c] && aa[c] == ab[c] && na[c] == nb[c];
      }
      expect_true(same);
    }
  }

  test_that("unsupported materials, masks, and boundaries opt out") {
    Transform identity;
    auto texture = std::make_shared<constant_texture>(point3f(.5));
    auto matte = std::make_shared<lambertian>(texture);
    auto glass = std::make_shared<dielectric>(point3f(1), 1.5, point3f(0), 0);
    auto emitter = std::make_shared<diffuse_light>(texture, 1, false);
    auto invisible = std::make_shared<diffuse_light>(texture, 1, true);
    unsigned char pixels[4] = {255, 255, 255, 128};
    auto alpha = std::make_shared<alpha_texture>(pixels, 1, 1, 4);
    auto shape = [&](std::shared_ptr<material> m, std::shared_ptr<alpha_texture> a = nullptr) {
      return std::make_shared<sphere>(1, m, a, nullptr, &identity, &identity, false);
    };
    expect_true(shape(matte)->ShadowType() == OpaqueShadowType::Opaque);
    expect_true(shape(glass)->ShadowType() == OpaqueShadowType::Unsupported);
    expect_true(shape(matte, alpha)->ShadowType() == OpaqueShadowType::Unsupported);
    expect_true(shape(invisible)->ShadowType() == OpaqueShadowType::Unsupported);
    auto boundary = std::make_shared<MediumBoundary>(shape(matte), nullptr, identity, false, 1);
    expect_true(boundary->ShadowType() == OpaqueShadowType::Unsupported);
    BVHAggregate many({shape(emitter), shape(emitter)}, 0, 1, 1, true);
    expect_true(many.ShadowType() == OpaqueShadowType::Mixed);
    BVHAggregate mixed({shape(matte), shape(matte, alpha)}, 0, 1, 1, true);
    expect_true(mixed.ShadowType() == OpaqueShadowType::Unsupported);
  }

  test_that("triangle predicates match NEE hit distances at edges and tiny offsets") {
    Transform identity;
    auto matte = std::make_shared<lambertian>(std::make_shared<constant_texture>(point3f(.5)));
    float vertices[9] = {-1, -1, 0, 1, -1, 0, 0, 1, 0};
    int indices[3] = {0, 1, 2};
    TriangleMesh mesh(vertices, indices, nullptr, nullptr, 3, 3, nullptr, nullptr,
                      matte, &identity, &identity, false);
    triangle tri(&mesh, mesh.vertexIndices.data(), mesh.normalIndices.data(),
                 mesh.texIndices.data(), 0, &identity, &identity, false);
    expect_true(tri.ShadowType() == OpaqueShadowType::Opaque);
    bool same = true;
    for (Float z : {Float(1), Float(1e-9), Float(1e-20)}) {
      for (int i = -15; i <= 15; ++i) {
        for (int j = -15; j <= 15; ++j) {
          Ray ray(point3f(i / 10.f, j / 10.f, -z), vec3f(0, 0, 1));
          ray.segment_absorption = true;
          for (Float end : {std::nextafter(z, Float(0)), z, std::nextafter(z, Float(2))}) {
            random_gen a(4), b(4);
            hit_record rec;
            same &= tri.hit(ray, 0, end, rec, a) == tri.OpaqueHit(ray, 0, end, b);
            same &= a.unif_rand() == b.unif_rand();
          }
        }
      }
    }
    expect_true(same);
    // A collapsed triangle must not become a blocker in the reduced predicate.
    mesh.p[2] = mesh.p[1];
    Ray ray(point3f(0, -1, -1), vec3f(0, 0, 1)); ray.segment_absorption = true;
    random_gen rng(1);
    expect_false(tri.OpaqueHit(ray, 0, MaxT, rng));
  }
}

context("Explicit emitter sampling") {
  test_that("light-tree probabilities normalize and match sampling at different vertices") {
    Transform identity;
    std::vector<Transform> transforms, inverses;
    transforms.reserve(8); inverses.reserve(8);
    hitable_list lights;
    for (int i = 0; i < 8; ++i) {
      transforms.push_back(Translate(vec3f((i - 3.5f) * 3, i % 2 ? 5 : 1, 0)));
      inverses.push_back(Inverse(transforms.back()));
      auto emission = std::make_shared<diffuse_light>(
          std::make_shared<constant_texture>(point3f(1)), i == 3 ? 20 : 1, false);
      lights.add(std::make_shared<sphere>(.25, emission, nullptr, nullptr,
                                         &transforms.back(), &inverses.back(), false));
    }
    VolumeLightSampler sampler(lights), fixed(lights, VolumeLightSampler::SelectionMethod::Fixed);
    using Context = VolumeLightSampler::Context;
    Context left{point3f(-10, 0, -2), normal3f(0), 0};
    Context right{point3f(10, 0, -2), normal3f(0), 0};
    expect_true(sampler.SelectionPmf(left, 0) > sampler.SelectionPmf(right, 0) * 5);
    expect_true(sampler.SelectionPmf(right, 7) > sampler.SelectionPmf(left, 7) * 5);
    bool correct = true;
    for (const Context &ctx : {left, right, Context{point3f(0), normal3f(0, 1, 0), .5f},
                              Context{point3f(0), normal3f(0, -1, 0), .5f}}) {
      double total = 0;
      std::vector<int> counts(8);
      random_gen rng(845); RandomSampler random(rng);
      const int count = 32768;
      for (int i = 0; i < count; ++i) {
        auto sample = sampler.SampleEmitter(ctx, &random, rng);
        correct &= sample.index < 8;
        if (sample.index < 8) ++counts[sample.index];
        hit_record hit;
        hit.shape = lights.objects[sample.index].get();
        correct &= sample.pdf == sampler.EmitterPdf(ctx, sample.wi, hit, rng);
      }
      for (size_t i = 0; i < 8; ++i) {
        double p = sampler.SelectionPmf(ctx, i);
        total += p;
        correct &= p >= .05 / 8 && p <= 1;
        correct &= std::abs(double(counts[i]) / count - p) < 6 * std::sqrt(p * (1 - p) / count) + 5.0 / count;
        correct &= fixed.SelectionPmf(ctx, i) == 1.0 / 8;
      }
      correct &= std::abs(total - 1) < 1e-12;
      correct &= sampler.SelectionPmf(ctx, 8) == 0;
    }
    expect_true(correct);
  }

  test_that("power, facing and receiver estimates retain support for missed texture features") {
    Transform identity;
    class SmallBrightPatch final : public texture {
    public:
      point3f value(Float u, Float v, const point3f &) const override {
        return point3f(u < .01f && v < .01f ? 10000 : 0);
      }
    };
    auto unit = std::make_shared<constant_texture>(point3f(1));
    auto dim = std::make_shared<diffuse_light>(unit, 1, false);
    auto bright = std::make_shared<diffuse_light>(unit, 10, false);
    auto patch = std::make_shared<diffuse_light>(std::make_shared<SmallBrightPatch>(), 1, false);
    hitable_list lights;
    for (auto mat : {dim, bright, patch})
      lights.add(std::make_shared<xz_rect>(-1, 1, -1, 1, 3, mat, nullptr, nullptr,
                                          &identity, &identity, true));
    lights.add(std::make_shared<xz_rect>(-1, 1, -1, 1, 3, bright, nullptr, nullptr,
                                        &identity, &identity, false));
    VolumeLightSampler sampler(lights);
    VolumeLightSampler::Context ctx{point3f(0), normal3f(0, 1, 0), 0};
    expect_true(sampler.SelectionPmf(ctx, 1) > sampler.SelectionPmf(ctx, 0) * 5);
    expect_true(sampler.SelectionPmf(ctx, 2) >= .05 / 4);
    expect_true(sampler.SelectionPmf(ctx, 3) >= .05 / 4);
    hitable_list facing;
    facing.add(lights.objects[1]); facing.add(lights.objects[3]);
    VolumeLightSampler facing_sampler(facing);
    expect_true(facing_sampler.SelectionPmf(ctx, 0) > facing_sampler.SelectionPmf(ctx, 1) * 5);
    auto opposite = ctx; opposite.n = -ctx.n;
    for (size_t i = 0; i < 4; ++i)
      expect_true(sampler.SelectionPmf(ctx, i) == sampler.SelectionPmf(opposite, i));

    unsigned char pixel[4] = {255, 255, 255, 255};
    auto alpha = std::make_shared<alpha_texture>(pixel, 1, 1, 4);
    hitable_list masked;
    for (bool flipped : {false, true})
      masked.add(std::make_shared<xz_rect>(-1, 1, -1, 1, 3, bright, alpha, nullptr,
                                          &identity, &identity, flipped));
    VolumeLightSampler masked_sampler(masked);
    expect_true(masked_sampler.SelectionPmf(ctx, 0) == .5);
    ctx.p = point3f(0, 6, 0);
    expect_true(masked_sampler.SelectionPmf(ctx, 0) == .5);
  }

  test_that("distant lights and empty groups retain their original probability mass") {
    Transform identity;
    auto texture = std::make_shared<constant_texture>(point3f(1));
    auto material = std::make_shared<diffuse_light>(texture, 1, false);
    hitable_list lights, empty;
    lights.add(std::make_shared<sphere>(1, material, nullptr, nullptr, &identity, &identity, false));
    lights.add(std::make_shared<sphere>(2, material, nullptr, nullptr, &identity, &identity, false));
    lights.add(std::make_shared<InfiniteAreaLight>(4, 2, 100, point3f(0), texture,
                                                  material, &identity, &identity, false));
    lights.add(std::make_shared<instance>(&empty, &identity, &identity, &empty));
    VolumeLightSampler sampler(lights);
    VolumeLightSampler::Context ctx{point3f(0), normal3f(0), 0};
    double sum = 0;
    for (size_t i = 0; i < 3; ++i) sum += sampler.SelectionPmf(ctx, i);
    expect_true(std::abs(sum - .75) < 1e-12);
    expect_true(sampler.SelectionPmf(ctx, 2) == .25);
    random_gen rng(423); RandomSampler random(rng);
    int distant = 0, null = 0;
    const int count = 16384;
    for (int i = 0; i < count; ++i) {
      auto sample = sampler.SampleEmitter(ctx, &random, rng);
      distant += sample.index == 2;
      null += sample.index == 3 && sample.pdf == 0;
    }
    expect_true(std::abs(double(distant) / count - .25) < .02);
    expect_true(std::abs(double(null) / count - .25) < .02);
  }

  test_that("shared emitters retain distinct nested, mirrored and animated placement IDs") {
    Transform identity, left = Translate(vec3f(-2, 0, 0)), right = Translate(vec3f(2, 0, 0));
    Transform left_inv = Inverse(left), right_inv = Inverse(right);
    auto emission = std::make_shared<diffuse_light>(
        std::make_shared<constant_texture>(point3f(1)), 2, false);
    auto lamp = std::make_shared<sphere>(.25, emission, nullptr, nullptr, &identity, &identity, false);
    hitable_list leaf(lamp), inner;
    inner.add(std::make_shared<instance>(&leaf, &left, &left_inv, &leaf));
    inner.add(std::make_shared<instance>(&leaf, &right, &right_inv, &leaf));
    Transform upper = Translate(vec3f(0, 3, 0)) * Scale(-1, 1.4, 1);
    Transform lower = Translate(vec3f(0, -3, 0)), upper_inv = Inverse(upper), lower_inv = Inverse(lower);
    hitable_list lights;
    lights.add(std::make_shared<instance>(&inner, &upper, &upper_inv, &inner));
    std::shared_ptr<hitable> moving = std::make_shared<instance>(&inner, &lower, &lower_inv, &inner);
    Transform end = Translate(vec3f(.4, 0, 0));
    lights.add(std::make_shared<AnimatedHitable>(moving, AnimatedTransform(&identity, 0, &end, 1)));
    // Rebuilding the sampler must rebuild the intern table consistently as well.
    VolumeLightSampler first(lights);
    VolumeLightSampler sampler(lights);
    random_gen rng(2026);
    RandomSampler random(rng);
    std::set<uint64_t> ids;
    bool correct = true;
    int hits = 0;
    point3f p(0, 0, -15);
    for (int i = 0; i < 2048; ++i) {
      Float time = Float(i % 17) / 16;
      auto selected = sampler.SampleEmitter({p, normal3f(0), time}, &random, rng);
      if (!(selected.pdf > 0)) continue;
      Ray ray(p, selected.wi, time); ray.segment_absorption = true;
      hit_record endpoint, reached;
      if (!sampler.Endpoint(selected, ray, endpoint, rng) || !lights.hit(ray, 0, MaxT, reached, rng)) {
        correct = false;
        continue;
      }
      ++hits;
      ids.insert(reached.light_placement);
      correct &= reached.shape == lamp.get() && sampler.Matches(selected, reached);
      correct &= endpoint.light_placement == reached.light_placement;
      Float pdf = sampler.EmitterPdf({p, normal3f(0), time}, selected.wi, reached, rng);
      correct &= pdf == selected.pdf;
      reached.light_placement = 0;
      correct &= !sampler.Matches(selected, reached);
      correct &= sampler.EmitterPdf({p, normal3f(0), time}, selected.wi, reached, rng) == 0;
    }
    expect_true(correct);
    expect_true(hits > 2000);
    expect_true(ids.size() == 4);
  }

  test_that("mesh lights sample the BVH-owned triangles after releasing the build list") {
    Transform identity;
    auto material = std::make_shared<diffuse_light>(
        std::make_shared<constant_texture>(point3f(1)), 3, false);
    float vertices[12] = {-1, 3, -1, 1, 3, -1, 1, 3, 1, -1, 3, 1};
    int indices[6] = {0, 1, 2, 0, 2, 3};
    auto mesh = std::make_shared<mesh3d>();
    mesh->mesh = std::make_unique<TriangleMesh>(vertices, indices, nullptr, nullptr, 4, 6,
                                              nullptr, nullptr, material, &identity, &identity, false);
    for (int i = 0; i < 6; i += 3)
      mesh->triangles.add(std::make_shared<triangle>(mesh->mesh.get(),
          &mesh->mesh->vertexIndices[i], &mesh->mesh->normalIndices[i],
          &mesh->mesh->texIndices[i], i / 3, &identity, &identity, false));
    mesh->mesh_bvh = std::make_shared<BVHAggregate>(mesh->triangles.objects, 0, 1, 1, true);
    mesh->triangles.objects.clear();
    hitable_list lights(mesh);
    VolumeLightSampler sampler(lights);
    random_gen rng(42);
    RandomSampler random(rng);
    bool correct = true;
    for (int i = 0; i < 1024; ++i) {
      auto sample = sampler.SampleEmitter({point3f(0), normal3f(0), 0}, &random, rng);
      Ray ray(point3f(0), sample.wi); ray.segment_absorption = true;
      hit_record hit;
      correct &= sample.pdf > 0 && sampler.Pdf(ray.o, ray.d, rng, 0) > 0;
      if (!lights.hit(ray, 0, MaxT, hit, rng)) { correct = false; continue; }
      correct &= sampler.Matches(sample, hit);
      correct &= sampler.EmitterPdf({ray.o, normal3f(0), 0}, ray.d, hit, rng) == sample.pdf;
    }
    expect_true(correct);
  }

  test_that("overlapping emitter proposals preserve mean radiance with matching MIS") {
    Transform identity;
    auto matte = std::make_shared<lambertian>(std::make_shared<constant_texture>(point3f(.6)));
    auto warm = std::make_shared<diffuse_light>(
        std::make_shared<constant_texture>(point3f(1, .4, .2)), 3, false);
    auto cool = std::make_shared<diffuse_light>(
        std::make_shared<constant_texture>(point3f(.2, .4, 1)), 5, false);
    auto floor = std::make_shared<xz_rect>(-20, 20, -20, 20, 0, matte, nullptr, nullptr,
                                         &identity, &identity, false);
    auto a = std::make_shared<xz_rect>(-1, 1, -1, 1, 3, warm, nullptr, nullptr,
                                     &identity, &identity, true);
    auto b = std::make_shared<xz_rect>(-3, 3, -3, 3, 5, cool, nullptr, nullptr,
                                     &identity, &identity, true);
    // The second configuration adds real scattering, null candidates and
    // boundary crossings to the same overlapping-emitter comparison.
    class TestMedium final : public Medium {
    public:
      TestMedium() : Medium(medium_description(.4)) {}
      bool IsHomogeneous() const override { return false; }
      MediumProperties SamplePoint(const point3f &p) const override {
        auto properties = Medium::SamplePoint(p);
        properties.sigma_s *= Float(.4 + .2 * std::sin(p[0]));
        return properties;
      }
    };
    for (bool fog : {false, true}) {
      hitable_list objects, selected, mixture;
      objects.add(floor); objects.add(a); objects.add(b);
      selected.add(a); selected.add(b); mixture.add(a); mixture.add(b);
      selected.volume_scene = std::make_shared<VolumeScene>();
      selected.volume_scene->light_sampler = std::make_shared<VolumeLightSampler>(selected);
      if (fog) {
        auto container = std::make_shared<sphere>(10, matte, nullptr, nullptr,
                                                 &identity, &identity, false);
        auto boundary = std::make_shared<MediumBoundary>(container, std::make_shared<TestMedium>(),
                                                         identity, false, 1);
        objects.add(boundary);
        mixture.volume_scene = std::make_shared<VolumeScene>();
        for (auto scene : {selected.volume_scene, mixture.volume_scene}) {
          scene->boundaries.add(boundary);
          scene->has_media = true;
          scene->Finish(0, 1);
        }
      }
      BVHAggregate world(objects.objects, 0, 1, 1, true);
      double sum[3] = {}, squares[3] = {}, reference[3] = {};
      const int count = 40000;
      for (int i = 0; i < count; ++i) {
        random_gen ra(i), rb(i);
        RandomSampler sa(ra), sb(rb);
        sa.independent_dimensions = sb.independent_dimensions = true;
        Ray ray(point3f(0, 1, -2), unit_vector(vec3f(0, -1, 2)));
        Float ta, tb;
        point3f ca, cb, albedo;
        normal3f normal;
        color_volume(ray, &world, &selected, 1, 5, ra, &sa, ta, ca, normal, albedo, nullptr);
        color_volume(ray, &world, &mixture, 1, 5, rb, &sb, tb, cb, normal, albedo, nullptr);
        for (int c = 0; c < 3; ++c) {
          double difference = double(ca[c]) - cb[c];
          sum[c] += difference;
          squares[c] += difference * difference;
          reference[c] += cb[c];
        }
      }
      for (int c = 0; c < 3; ++c) {
        double mean = sum[c] / count;
        double se = std::sqrt((squares[c] / count - mean * mean) / count);
        expect_true(std::abs(mean) < 6 * se + 1e-6);
        expect_true(std::abs(sum[c] / reference[c]) < .02);
      }
    }
  }
}

#endif
