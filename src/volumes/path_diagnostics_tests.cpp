#ifdef NOT_CRAN
#include "path_diagnostics.h"
#include "boundary.h"
#include "volpath.h"
#include "../materials/material.h"
#include <testthat.h>
#include <thread>

namespace {
class EmittingPassthrough : public material {
public:
  const std::string GetName() override { return "EmittingPassthrough"; }
  size_t GetSize() override { return sizeof(*this); }
  point3f emitted(const Ray &, const hit_record &, Float, Float, const point3f &, bool &) override {
    return point3f(1, 2, 3);
  }
  bool scatter(const Ray &r, const hit_record &h, scatter_record &s, Sampler *) override {
    s.is_specular = true;
    s.is_passthrough = true;
    s.attenuation = point3f(1);
    s.specular_ray = Ray(h.p, r.d, r.time());
    return true;
  }
};
class FailingPathWorld : public hitable {
public:
  explicit FailingPathWorld(int mode) : mode(mode) {}
  const bool hit(const Ray &r, Float, Float, hit_record &h, random_gen &) const override {
    if (calls++) {
      if (mode == 1) throw PathFailure(PathFailureKind::RepeatedEntry, "synthetic repeated entry");
      if (mode == 2) throw std::runtime_error("unclassified failure");
      return false;
    }
    h.t = 1;
    h.p = r(1);
    h.normal = h.geometric_normal = normal3f(0, 0, -1);
    h.pError = vec3f(0);
    h.u = h.v = 0;
    h.mat_ptr = mat.get();
    h.shape = this;
    return true;
  }
  const bool hit(const Ray &r, Float lo, Float hi, hit_record &h, Sampler *) const override {
    random_gen rng(0);
    return hit(r, lo, hi, h, rng);
  }
  bool bounding_box(Float, Float, aabb &) const override { return false; }
  std::string GetName() const override { return "FailingPathWorld"; }
  size_t GetSize() override { return sizeof(*this); }
  void hitable_info_bounds(Float, Float) const override {}
private:
  int mode;
  mutable int calls = 0;
  std::shared_ptr<material> mat = std::make_shared<EmittingPassthrough>();
};
}
context("Recoverable transport failures") {
  test_that("a failed continuation retains earlier emission and later paths run") {
    hitable_list lights;
    lights.volume_scene = std::make_shared<VolumeScene>();
    random_gen rng(123);
    RandomSampler sampler(rng);
    Ray input(point3f(0), vec3f(0, 0, 1));
    Float alpha;
    point3f radiance, albedo;
    normal3f normal;
    for (int mode : {1, 0}) {
      FailingPathWorld world(mode);
      color_volume(input, &world, &lights, 5, 5, rng, &sampler, alpha,
                   radiance, normal, albedo, nullptr);
      expect_true((radiance[0] == 1 && radiance[1] == 2 && radiance[2] == 3));
      expect_true(std::isfinite(alpha));
      auto record = lights.volume_scene->path_diagnostics.Take();
      expect_true(Rcpp::as<double>(record["terminated_paths"]) == mode);
    }
    FailingPathWorld fatal(2);
    bool propagated = false;
    try {
      color_volume(input, &fatal, &lights, 5, 5, rng, &sampler, alpha,
                   radiance, normal, albedo, nullptr);
    } catch (const std::runtime_error &e) {
      propagated = std::string(e.what()) == "unclassified failure";
    }
    expect_true(propagated);
    expect_true(Rcpp::as<double>(lights.volume_scene->path_diagnostics.Take()["terminated_paths"]) == 0);
  }
  test_that("concurrent failures are counted exactly and example memory is bounded") {
    PathDiagnostics diagnostics;
    std::vector<std::thread> workers;
    for (int t = 0; t < 8; ++t) workers.emplace_back([&]() {
      Ray ray(point3f(1, 2, 3), vec3f(0, 1, 0));
      for (int i = 0; i < 1000; ++i)
        diagnostics.Record(PathFailure(PathFailureKind::InvalidExit, "invalid exit"),
                           ray, ray, 2, 17, 4, "radiance");
    });
    for (auto &worker : workers) worker.join();
    auto record = diagnostics.Take();
    expect_true(Rcpp::as<double>(record["terminated_paths"]) == 8000);
    Rcpp::NumericVector counts = record["counts"];
    expect_true(double(counts["invalid_exit"]) == 8000);
    Rcpp::List examples = record["examples"];
    Rcpp::CharacterVector exits = examples["invalid_exit"];
    expect_true(exits.size() == 8);
    expect_true(Rcpp::as<std::string>(exits[0]).find("active_boundary=4") != std::string::npos);
    expect_true(Rcpp::as<double>(diagnostics.Take()["terminated_paths"]) == 0);
  }
}
#endif
