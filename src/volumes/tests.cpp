#ifdef NOT_CRAN
#include "../hitables/box.h"
#include "../hitables/instance.h"
#include "../hitables/sphere.h"
#include "../materials/material.h"
#include "../materials/texture.h"
#include "boundary.h"
#include "medium.h"
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
} // namespace
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
    Transform far = Translate(vec3f(30, 0, 0)), far_inverse = Inverse(far);
    auto far_geometry = std::make_shared<sphere>(1, mat, nullptr, nullptr, &far, &far_inverse, false);
    auto far_boundary = std::make_shared<MediumBoundary>(far_geometry, medium, far, false,
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
#endif
