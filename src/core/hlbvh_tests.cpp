#ifdef NOT_CRAN
#include "../hitables/hitablelist.h"
#include "bvh.h"
#include <random>
#include <set>
#include <thread>
#include <testthat.h>

namespace {
bool check_tree(const HLBVHTree &tree, std::span<const HLBVHBounds> boxes) {
  if (boxes.empty())
    return tree.order.empty() && !tree.nodes;
  std::vector<unsigned> seen(boxes.size(), 0);
  std::vector<uint32_t> pending{tree.root};
  size_t visited = 0;
  while (!pending.empty()) {
    uint32_t index = pending.back();
    pending.pop_back();
    if (index >= 2 * boxes.size() + 4096 || ++visited > 2 * boxes.size())
      return false;
    const auto &n = tree.nodes[index];
    auto contains = [&](const HLBVHBounds &b) {
      for (int d = 0; d < 3; ++d)
        if (n.bounds.lo[d] > b.lo[d] || n.bounds.hi[d] < b.hi[d])
          return false;
      return true;
    };
    if (n.count) {
      if (n.count > 4 || size_t(n.first) + n.count > boxes.size())
        return false;
      for (uint32_t j = n.first; j < n.first + n.count; ++j) {
        uint32_t p = tree.order[j].index;
        if (p >= boxes.size() || ++seen[p] != 1 || !contains(boxes[p]))
          return false;
      }
    } else {
      if (n.left >= 2 * boxes.size() + 4096 || n.right >= 2 * boxes.size() + 4096)
        return false;
      if (!contains(tree.nodes[n.left].bounds) || !contains(tree.nodes[n.right].bounds))
        return false;
      pending.push_back(n.left);
      pending.push_back(n.right);
    }
  }
  return std::all_of(seen.begin(), seen.end(), [](unsigned n) { return n == 1; });
}
bool same_tree(const HLBVHTree &a, const HLBVHTree &b) {
  if (a.order.size() != b.order.size() || a.root != b.root)
    return false;
  for (size_t i = 0; i < a.order.size(); ++i)
    if (a.order[i].code != b.order[i].code || a.order[i].index != b.order[i].index)
      return false;
  if (a.order.empty())
    return true;
  std::vector<uint32_t> pending{a.root};
  while (!pending.empty()) {
    uint32_t i = pending.back();
    pending.pop_back();
    const auto &x = a.nodes[i];
    const auto &y = b.nodes[i];
    if (x.count != y.count || x.axis != y.axis)
      return false;
    for (int d = 0; d < 3; ++d)
      if (x.bounds.lo[d] != y.bounds.lo[d] || x.bounds.hi[d] != y.bounds.hi[d])
        return false;
    if (x.count) {
      if (x.first != y.first)
        return false;
    } else {
      if (x.left != y.left || x.right != y.right)
        return false;
      pending.push_back(x.left);
      pending.push_back(x.right);
    }
  }
  return true;
}

// Moving analytic boxes exercise shutter-interval bounds and inside/outside
// rays against a brute-force oracle, independently of triangle intersection.
class MovingBoxProbe : public hitable {
public:
  point3f center;
  OpaqueShadowType shadow_type;
  explicit MovingBoxProbe(point3f c, OpaqueShadowType type = OpaqueShadowType::Opaque)
      : center(c), shadow_type(type) {}
  bool bounding_box(Float t0, Float t1, aabb &b) const override {
    b = aabb(center + vec3f(t0 - .4, -.4, -.4), center + vec3f(t1 + .4, .4, .4));
    return true;
  }
  const bool hit(const Ray &r, Float lo, Float hi, hit_record &rec, random_gen &) const override {
    return intersect(r, lo, hi, rec);
  }
  const bool hit(const Ray &r, Float lo, Float hi, hit_record &rec, Sampler *) const override {
    return intersect(r, lo, hi, rec);
  }
  OpaqueShadowType ShadowType() const override { return shadow_type; }
  bool OpaqueHit(const Ray &r, Float lo, Float hi, random_gen &) const override {
    hit_record rec;
    return intersect(r, lo, hi, rec);
  }
  std::string GetName() const override { return "Moving HLBVH test box"; }
  size_t GetSize() override { return sizeof(*this); }
  void hitable_info_bounds(Float, Float) const override {}

private:
  bool intersect(const Ray &r, Float lo, Float hi, hit_record &rec) const {
    double entry = -INFINITY, exit = INFINITY;
    for (int d = 0; d < 3; ++d) {
      double c = center[d] + (d == 0 ? r.time() : 0);
      if (r.d[d] == 0) {
        if (r.o[d] < c - .4 || r.o[d] > c + .4)
          return false;
        continue;
      }
      double a = (c - .4 - r.o[d]) / r.d[d], b = (c + .4 - r.o[d]) / r.d[d];
      entry = std::max(entry, std::min(a, b));
      exit = std::min(exit, std::max(a, b));
    }
    double t = entry >= lo ? entry : exit;
    if (entry > exit || t < lo || t > hi)
      return false;
    rec.t = Float(t);
    rec.precise_t = t;
    return true;
  }
};
struct PreparationAudit {
  const std::thread::id caller = std::this_thread::get_id();
  std::atomic<unsigned> bounds_calls{0}, shadow_calls{0}, worker_calls{0};
};

class PreparationProbe final : public MovingBoxProbe {
public:
  PreparationAudit &audit;
  bool has_bounds = true;
  PreparationProbe(point3f center, PreparationAudit &audit)
      : MovingBoxProbe(center), audit(audit) {}
  bool bounding_box(Float t0, Float t1, aabb &b) const override {
    ++audit.bounds_calls;
    if (std::this_thread::get_id() != audit.caller)
      ++audit.worker_calls;
    return has_bounds && MovingBoxProbe::bounding_box(t0, t1, b);
  }
  OpaqueShadowType ShadowType() const override {
    ++audit.shadow_calls;
    if (std::this_thread::get_id() != audit.caller)
      ++audit.worker_calls;
    return MovingBoxProbe::ShadowType();
  }
};
} // namespace

context("Parallel HLBVH construction") {
  test_that("primitive query parallelism requires opt-in and preserves moving bounds and hits") {
    bool correct = true;
    for (int n : {32767, 32768, 32769}) {
      PreparationAudit audit;
      hitable_list brute;
      for (int i = 0; i < n; ++i)
        brute.add(std::make_shared<PreparationProbe>(
            point3f(2 * (i % 300), 0, 2 * (i / 300)), audit));
      aabb expected_bounds;
      brute.bounding_box(-.5, 1.5, expected_bounds);
      for (auto option : {BVHBuildOptions{BVHBuildMethod::HLBVH, 6, false},
                          BVHBuildOptions{BVHBuildMethod::HLBVH, 1, true},
                          BVHBuildOptions{BVHBuildMethod::HLBVH, 6, true}}) {
        audit.bounds_calls = audit.shadow_calls = audit.worker_calls = 0;
        BVHAggregate bvh(brute.objects, -.5, 1.5, 1, true, option);
        correct &= audit.bounds_calls == unsigned(n) && audit.shadow_calls == unsigned(n);
        const bool parallel = option.parallel_primitive_queries && option.threads > 1 && n >= 32768;
        correct &= parallel ? audit.worker_calls > 0 : audit.worker_calls == 0;
        aabb actual_bounds;
        correct &= bvh.bounding_box(-.5, 1.5, actual_bounds);
        for (int d = 0; d < 3; ++d)
          correct &= actual_bounds.min()[d] == expected_bounds.min()[d] &&
                     actual_bounds.max()[d] == expected_bounds.max()[d];
        random_gen rng(1);
        for (int i : {0, 4095, 4096, n - 1}) {
          Ray ray(point3f(2 * (i % 300) + .5f, -4, 2 * (i / 300)), vec3f(0, 1, 0), .5f);
          hit_record expected, actual;
          bool hit = brute.hit(ray, 0, 100, expected, rng);
          correct &= bvh.hit(ray, 0, 100, actual, rng) == hit;
          if (hit)
            correct &= actual.t == expected.t;
        }
      }
    }
    expect_true(correct);
  }
  test_that("parallel block reduction preserves shadow types and primitive ordering") {
    bool correct = true;
    std::vector<std::shared_ptr<hitable>> objects(32769);
    auto opaque = std::make_shared<MovingBoxProbe>(point3f(0, 0, 0));
    for (auto type : {OpaqueShadowType::Light, OpaqueShadowType::Mixed, OpaqueShadowType::Unsupported})
      for (int position : {0, 4096, 32768}) {
        std::fill(objects.begin(), objects.end(), opaque);
        objects[position] = std::make_shared<MovingBoxProbe>(point3f(2, 0, 0), type);
        BVHAggregate serial(objects, 0, 1, 1, true, {BVHBuildMethod::HLBVH, 6, false});
        BVHAggregate parallel(objects, 0, 1, 1, true, {BVHBuildMethod::HLBVH, 6, true});
        correct &= serial.ShadowType() == parallel.ShadowType();
        correct &= std::equal(serial.Primitives().begin(), serial.Primitives().end(),
                              parallel.Primitives().begin());
      }
    expect_true(correct);
  }
  test_that("parallel primitive query errors are propagated without leaving workers running") {
    bool correct = true;
    PreparationAudit audit;
    auto good = std::make_shared<PreparationProbe>(point3f(0, 0, 0), audit);
    auto bad = std::make_shared<PreparationProbe>(point3f(1, 0, 0), audit);
    bad->has_bounds = false;
    std::vector<std::shared_ptr<hitable>> objects(32769, good);
    for (int position : {0, 4096, 32768}) {
      objects[position] = bad;
      bool caught = false;
      try {
        BVHAggregate bvh(objects, 0, 1, 1, true, {BVHBuildMethod::HLBVH, 6, true});
      } catch (const std::runtime_error &error) {
        caught = std::string(error.what()) == "BVH primitive has no bounding box.";
      }
      correct &= caught;
      objects[position] = good;
    }
    BVHAggregate recovered(objects, 0, 1, 1, true, {BVHBuildMethod::HLBVH, 6, true});
    correct &= recovered.ShadowType() == OpaqueShadowType::Opaque;
    expect_true(correct);
  }
  test_that("CPU treelets are deterministic and Metal matches CPU topology and bounds") {
    std::mt19937 rng(41);
    std::uniform_real_distribution<float> random(-10, 10);
    bool valid = true, equal = true;
    const bool metal = MetalBVHAvailable();
    for (unsigned leaf_size : {1u, 4u})
      for (int n : {0, 1, 2, 3, 257, 10000, 70000})
        for (int shape = 0; shape < 3; ++shape) {
          std::vector<HLBVHBounds> boxes(n);
          for (auto &b : boxes)
            for (int d = 0; d < 3; ++d) {
              float c = shape == 0 ? 0 : shape == 1 ? random(rng) : random(rng) * 1e36f;
              float radius = shape == 2 ? 1e35f : .01f;
              b.lo[d] = c - radius;
              b.hi[d] = c + radius;
            }
          auto serial = BuildHLBVH(boxes, leaf_size, {BVHBuildMethod::HLBVH, 1});
          auto parallel = BuildHLBVH(boxes, leaf_size, {BVHBuildMethod::HLBVH, 6});
          valid &= check_tree(serial, boxes) && check_tree(parallel, boxes);
          equal &= same_tree(serial, parallel);
          if (metal) {
            auto gpu = BuildHLBVH(boxes, leaf_size, {BVHBuildMethod::Metal, 6});
            valid &= check_tree(gpu, boxes);
            equal &= same_tree(serial, gpu);
          }
        }
    expect_true(valid);
    expect_true(equal);
    Rprintf("HLBVH topology backend validation: CPU%s\n", metal ? " + Metal" : " (Metal unavailable)");
  }
  test_that("all builders preserve moving closest hits, shadows, and containment filtering") {
    std::mt19937 random(79);
    std::uniform_real_distribution<float> uniform(-8, 8);
    hitable_list brute;
    for (int i = 0; i < 300; ++i)
      brute.add(std::make_shared<MovingBoxProbe>(point3f(uniform(random), uniform(random), uniform(random))));
    std::vector<BVHBuildMethod> methods{BVHBuildMethod::SAH, BVHBuildMethod::Equal, BVHBuildMethod::HLBVH};
    if (MetalBVHAvailable())
      methods.push_back(BVHBuildMethod::Metal);
    bool correct = true;
    for (auto method : methods) {
      BVHAggregate bvh(brute.objects, 0, 1, 4, true, {method, 6});
      random_gen rng(1);
      RandomSampler sampler(rng);
      for (int i = 0; i < 4000; ++i) {
        Ray ray(point3f(uniform(random), uniform(random), uniform(random)),
                unit_vector(vec3f(uniform(random), uniform(random), uniform(random))), Float(i % 11) / 10);
        hit_record expected, actual;
        bool found = brute.hit(ray, 0, 100, expected, rng);
        for (bool sampled : {false, true}) {
          bool hit = sampled ? bvh.hit(ray, 0, 100, actual, &sampler) : bvh.hit(ray, 0, 100, actual, rng);
          correct &= hit == found && (!hit || std::abs(actual.t - expected.t) < 1e-5);
          correct &= (sampled ? bvh.HitP(ray, 0, 100, &sampler) : bvh.HitP(ray, 0, 100, rng)) == found;
        }
        correct &= bvh.OpaqueHit(ray, 0, 100, rng) == found;
        hitable_list containing;
        for (const auto &object : brute.objects) {
          aabb b;
          object->bounding_box(ray.time(), ray.time(), b);
          bool inside = true;
          for (int d = 0; d < 3; ++d)
            inside &= ray.o[d] >= b.min()[d] && ray.o[d] <= b.max()[d];
          if (inside)
            containing.add(object);
        }
        found = containing.hit(ray, 0, 100, expected, rng);
        bool hit = bvh.HitContainingObjects(ray, 0, 100, actual, rng);
        correct &= hit == found && (!hit || std::abs(actual.t - expected.t) < 1e-5);
      }
    }
    expect_true(correct);
  }
  test_that("invalid bounds and unavailable Metal report an actionable error") {
    std::vector<HLBVHBounds> boxes(1);
    boxes[0].lo[0] = NAN;
    bool rejected = false;
    try {
      BuildHLBVH(boxes, 1, {BVHBuildMethod::HLBVH, 1});
    } catch (const std::runtime_error &) {
      rejected = true;
    }
    expect_true(rejected);
    boxes[0] = {};
    if (!MetalBVHAvailable()) {
      rejected = false;
      try {
        BuildHLBVH(boxes, 1, {BVHBuildMethod::Metal, 1});
      } catch (const std::runtime_error &) {
        rejected = true;
      }
      expect_true(rejected);
    }
  }
  test_that("parallel centroid validation rejects invalid blocks and recovers") {
    // Exercise both sides of the parallel threshold and the partial last block.
    // Invalid inputs must fail before a Metal device or any partial result is used.
    bool correct = true;
    for (size_t count : {32767u, 32768u, 32769u}) {
      std::vector<HLBVHBounds> boxes(count);
      for (size_t position : {size_t(0), size_t(4095), size_t(4096), count - 1}) {
        for (int kind = 0; kind < 3; ++kind) {
          boxes[position] = {};
          if (kind == 0)
            boxes[position].lo[0] = NAN;
          else if (kind == 1)
            boxes[position].hi[1] = INFINITY;
          else
            boxes[position].lo[2] = 1;
          for (auto method : {BVHBuildMethod::HLBVH, BVHBuildMethod::Metal}) {
            bool rejected = false;
            try {
              BuildHLBVH(boxes, 1, {method, 6});
            } catch (const std::runtime_error &error) {
              rejected = std::string(error.what()) == "HLBVH requires finite, ordered primitive bounds.";
            }
            correct &= rejected;
          }
        }
        boxes[position] = {};
      }
      auto serial = BuildHLBVH(boxes, 1, {BVHBuildMethod::HLBVH, 1});
      auto recovered = BuildHLBVH(boxes, 1, {BVHBuildMethod::HLBVH, 6});
      correct &= check_tree(recovered, boxes) && same_tree(serial, recovered);
    }
    expect_true(correct);
  }
  test_that("packed small trees and parallel frontiers preserve axis-aligned hits") {
    std::vector<BVHBuildMethod> methods{BVHBuildMethod::HLBVH};
    if (MetalBVHAvailable())
      methods.push_back(BVHBuildMethod::Metal);
    bool correct = true;
    // Exercise root leaves, partial four-wide nodes, and both sides of the
    // parallel packing threshold. Zero direction components check masked lanes.
    for (int n : {0, 1, 2, 3, 17, 32767, 32768, 70000}) {
      hitable_list brute;
      for (int i = 0; i < n; ++i)
        brute.add(std::make_shared<MovingBoxProbe>(point3f(2 * (i % 300), 0, 2 * (i / 300))));
      for (auto method : methods)
        for (unsigned threads : {1u, 6u})
          for (int leaf_size : {1, 4}) {
            BVHAggregate bvh(brute.objects, 0, 1, leaf_size, true, {method, threads});
            random_gen rng(1);
            for (int j = 0; j < 16; ++j) {
              int index = n ? (int64_t(j) * (n - 1)) / 15 : 0;
              // Alternate a direct hit with a ray passing between columns.
              Ray ray(point3f(2 * (index % 300) + (j % 2 ? .5f : 1.5f), -4, 2 * (index / 300)),
                      vec3f(0, 1, 0), .5f);
              hit_record expected, actual;
              bool found = brute.hit(ray, 0, 100, expected, rng);
              bool hit = bvh.hit(ray, 0, 100, actual, rng);
              correct &= hit == found && (!hit || std::abs(actual.t - expected.t) < 1e-5);
              correct &= bvh.OpaqueHit(ray, 0, 100, rng) == found;
            }
          }
    }
    expect_true(correct);
  }
  test_that("fused primitive preparation preserves aggregate shadow classification") {
    std::vector<BVHBuildMethod> methods{BVHBuildMethod::SAH, BVHBuildMethod::HLBVH};
    if (MetalBVHAvailable())
      methods.push_back(BVHBuildMethod::Metal);
    bool correct = true;
    for (auto type : {OpaqueShadowType::Opaque, OpaqueShadowType::Light,
                      OpaqueShadowType::Mixed, OpaqueShadowType::Unsupported})
      for (int position = 0; position < 3; ++position) {
        hitable_list objects;
        for (int i = 0; i < 3; ++i)
          objects.add(std::make_shared<MovingBoxProbe>(point3f(i, 0, 0),
                      i == position ? type : OpaqueShadowType::Opaque));
        auto expected = type == OpaqueShadowType::Light ? OpaqueShadowType::Mixed : type;
        for (auto method : methods) {
          BVHAggregate bvh(objects.objects, 0, 1, 1, true, {method, 6});
          correct &= bvh.ShadowType() == expected;
        }
      }
    expect_true(correct);
  }
  test_that("primitive permutation transfers each shared ownership slot once") {
    std::vector<BVHBuildMethod> methods{BVHBuildMethod::HLBVH};
    if (MetalBVHAvailable())
      methods.push_back(BVHBuildMethod::Metal);
    bool correct = true;
    for (auto method : methods) {
      auto object = std::make_shared<MovingBoxProbe>(point3f(0, 0, 0));
      // Distinct input slots may intentionally share one primitive/control block.
      std::vector<std::shared_ptr<hitable>> objects{object, object};
      const auto before = object.use_count();
      {
        BVHAggregate bvh(objects, 0, 1, 1, true, {method, 6});
        correct &= object.use_count() == before + 2;
        random_gen rng(1);
        hit_record rec;
        correct &= bvh.hit(Ray(point3f(0, 0, -2), vec3f(0, 0, 1)), 0, 10, rec, rng);
        correct &= std::abs(rec.t - 1.6f) < 1e-5;
      }
      correct &= object.use_count() == before;
    }
    expect_true(correct);
  }
  test_that("all 4096 Morton cells can be occupied without losing boundary cells") {
    std::vector<HLBVHBounds> boxes(4096);
    for (unsigned i = 0; i < boxes.size(); ++i)
      for (unsigned d = 0; d < 3; ++d)
        boxes[i].lo[d] = boxes[i].hi[d] = float((i >> (4 * d)) & 15);
    bool correct = true;
    for (unsigned leaf_size : {1u, 4u}) {
      auto cpu = BuildHLBVH(boxes, leaf_size, {BVHBuildMethod::HLBVH, 6});
      correct &= cpu.treelets.size() == 4096 && check_tree(cpu, boxes);
      if (MetalBVHAvailable()) {
        auto gpu = BuildHLBVH(boxes, leaf_size, {BVHBuildMethod::Metal, 6});
        correct &= gpu.treelets.size() == 4096 && check_tree(gpu, boxes) && same_tree(cpu, gpu);
      }
    }
    expect_true(correct);
  }
}
#endif
