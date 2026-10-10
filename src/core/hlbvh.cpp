#include "hlbvh.h"
#include <RcppThread.h>
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace {
HLBVHBounds empty_bounds() {
  HLBVHBounds b{};
  for (int d = 0; d < 3; ++d) {
    b.lo[d] = INFINITY;
    b.hi[d] = -INFINITY;
  }
  return b;
}
HLBVHBounds unite(HLBVHBounds a, const HLBVHBounds &b) {
  for (int d = 0; d < 3; ++d) {
    a.lo[d] = std::min(a.lo[d], b.lo[d]);
    a.hi[d] = std::max(a.hi[d], b.hi[d]);
  }
  return a;
}
double area(const HLBVHBounds &b) {
  double x = std::max(0., double(b.hi[0]) - b.lo[0]);
  double y = std::max(0., double(b.hi[1]) - b.lo[1]);
  double z = std::max(0., double(b.hi[2]) - b.lo[2]);
  return 2 * (x * y + x * z + y * z);
}
float centroid(const HLBVHBounds &b, int d) { return .5f * b.lo[d] + .5f * b.hi[d]; }
uint32_t spread(uint32_t v) {
  v = (v | (v << 16)) & 0x030000FF;
  v = (v | (v << 8)) & 0x0300F00F;
  v = (v | (v << 4)) & 0x030C30C3;
  return (v | (v << 2)) & 0x09249249;
}

// Each worker owns a disjoint node arena and sorted primitive interval. The
// stable ordering removes scheduling-dependent primitive permutations.
uint32_t emit_treelet(HLBVHTree &tree, std::span<const HLBVHBounds> bounds, uint32_t first, uint32_t count, int bit,
                      unsigned max_leaf, uint32_t &next) {
  const uint32_t index = next++;
  auto &node = tree.nodes[index];
  if (count <= max_leaf) {
    // Nonempty leaves start with a real bound. The usual singleton leaf needs
    // only a copy, avoiding an infinity fill and three redundant min/max pairs.
    node.bounds = bounds[tree.order[first].index];
    for (uint32_t i = first + 1; i < first + count; ++i)
      node.bounds = unite(node.bounds, bounds[tree.order[i].index]);
    node.first = first;
    node.count = count;
    node.axis = 0;
    return index;
  }
  // Sorted endpoints share every bit above their highest differing bit.
  // countl_zero(0) is 32, yielding -1 for coincident Morton codes.
  const uint32_t different = tree.order[first].code ^ tree.order[first + count - 1].code;
  bit = std::min(bit, 31 - int(std::countl_zero(different)));
  uint32_t split = first + count / 2;
  if (bit >= 0) {
    uint32_t lo = first, hi = first + count;
    while (lo < hi) {
      uint32_t mid = lo + (hi - lo) / 2;
      if (tree.order[mid].code & (1u << bit))
        hi = mid;
      else
        lo = mid + 1;
    }
    split = lo;
  }
  // Unlike PBRT's unlimited coincident-code leaf, balanced index splits keep
  // scalar leaf counts representable even for >65535 coincident primitives.
  node.count = 0;
  node.axis = bit >= 0 ? unsigned(bit % 3) : 0;
  node.left = emit_treelet(tree, bounds, first, split - first, bit - 1, max_leaf, next);
  node.right = emit_treelet(tree, bounds, split, first + count - split, bit - 1, max_leaf, next);
  node.bounds = unite(tree.nodes[node.left].bounds, tree.nodes[node.right].bounds);
  return index;
}

uint32_t stitch(HLBVHTree &tree, std::span<uint32_t> roots, uint32_t &next) {
  if (roots.size() == 1)
    return roots[0];
  HLBVHBounds centers = empty_bounds();
  for (uint32_t root : roots)
    for (int d = 0; d < 3; ++d) {
      float c = centroid(tree.nodes[root].bounds, d);
      centers.lo[d] = std::min(centers.lo[d], c);
      centers.hi[d] = std::max(centers.hi[d], c);
    }
  int axis = 0;
  for (int d = 1; d < 3; ++d)
    if (double(centers.hi[d]) - centers.lo[d] > double(centers.hi[axis]) - centers.lo[axis])
      axis = d;
  size_t middle = roots.size() / 2;
  const double extent = double(centers.hi[axis]) - centers.lo[axis];
  if (extent > 0) {
    HLBVHBounds boxes[12];
    unsigned counts[12]{};
    for (auto &box : boxes)
      box = empty_bounds();
    auto bucket = [&](uint32_t root) {
      return std::clamp(int(12 * (double(centroid(tree.nodes[root].bounds, axis)) - centers.lo[axis]) / extent), 0, 11);
    };
    for (uint32_t root : roots) {
      int b = bucket(root);
      ++counts[b];
      boxes[b] = unite(boxes[b], tree.nodes[root].bounds);
    }
    double best = INFINITY;
    int split = 0;
    for (int s = 0; s < 11; ++s) {
      HLBVHBounds left = empty_bounds(), right = empty_bounds();
      unsigned nl = 0, nr = 0;
      for (int b = 0; b < 12; ++b) {
        if (b <= s) {
          left = unite(left, boxes[b]);
          nl += counts[b];
        } else {
          right = unite(right, boxes[b]);
          nr += counts[b];
        }
      }
      double cost = nl * area(left) + nr * area(right);
      if (nl && nr && cost < best) {
        best = cost;
        split = s;
      }
    }
    auto mid = std::partition(roots.begin(), roots.end(), [&](uint32_t root) { return bucket(root) <= split; });
    if (mid != roots.begin() && mid != roots.end())
      middle = mid - roots.begin();
  }
  const uint32_t index = next++;
  auto &node = tree.nodes[index];
  node.count = 0;
  node.axis = axis;
  node.left = stitch(tree, roots.first(middle), next);
  node.right = stitch(tree, roots.subspan(middle), next);
  node.bounds = unite(tree.nodes[node.left].bounds, tree.nodes[node.right].bounds);
  return index;
}
} // namespace

HLBVHTree BuildHLBVH(std::span<const HLBVHBounds> bounds, unsigned max_leaf, BVHBuildOptions options) {
  HLBVHTree tree;
  if (bounds.empty())
    return tree;
  if (bounds.size() > (size_t(std::numeric_limits<int>::max()) - 4096) / 2)
    throw std::runtime_error("HLBVH exceeds the renderer's node index range.");
  max_leaf = std::clamp(max_leaf, 1u, 255u);
  auto scan_centroids = [&](size_t first, size_t end) {
    HLBVHBounds centers = empty_bounds();
    for (size_t i = first; i < end; ++i)
      for (int d = 0; d < 3; ++d) {
        const auto &b = bounds[i];
        if (!std::isfinite(b.lo[d]) || !std::isfinite(b.hi[d]) || b.lo[d] > b.hi[d])
          throw std::runtime_error("HLBVH requires finite, ordered primitive bounds.");
        float c = centroid(b, d);
        centers.lo[d] = std::min(centers.lo[d], c);
        centers.hi[d] = std::max(centers.hi[d], c);
      }
    return centers;
  };
  HLBVHBounds centers;
  if (options.threads > 1 && bounds.size() >= 32768) {
    // Compact bounds are immutable native records, independent of primitive
    // query safety. Reduce contiguous blocks before either CPU or Metal work;
    // input-order merging preserves ties independently of worker scheduling.
    constexpr size_t block_size = 4096;
    const size_t blocks = (bounds.size() + block_size - 1) / block_size;
    std::vector<HLBVHBounds> partial(blocks);
    RcppThread::ThreadPool pool(options.threads);
    pool.parallelFor(size_t(0), blocks, [&](size_t block) {
      partial[block] = scan_centroids(block * block_size, std::min((block + 1) * block_size, bounds.size()));
    });
    pool.join(); // Invalid bounds must fail before consuming partial results.
    centers = empty_bounds();
    for (const auto &p : partial)
      centers = unite(centers, p);
  } else {
    centers = scan_centroids(0, bounds.size());
  }
  tree.order.resize(bounds.size());

  if (options.method == BVHBuildMethod::Metal) {
    BuildMetalTreelets(bounds, centers, max_leaf, tree);
  } else {
    tree.nodes.reset(new HLBVHNode[2 * bounds.size() + 4096]);
    RcppThread::ThreadPool pool(std::max(1u, options.threads));
    constexpr size_t block_size = 4096;
    const size_t blocks = (bounds.size() + block_size - 1) / block_size;
    // Work stealing reserves an iteration with an atomic operation. Reserve
    // blocks, not individual primitives, so this cheap loop can run locally.
    pool.parallelFor(size_t(0), blocks, [&](size_t block) {
      const size_t end = std::min(bounds.size(), (block + 1) * block_size);
      for (size_t i = block * block_size; i < end; ++i) {
        uint32_t code = 0;
        for (int d = 0; d < 3; ++d) {
          float extent = .5f * centers.hi[d] - .5f * centers.lo[d];
          float offset = extent > 0 ? (.5f * centroid(bounds[i], d) - .5f * centers.lo[d]) / extent : 0;
          code |= spread(uint32_t(std::clamp(offset * 1024.f, 0.f, 1023.f))) << d;
        }
        tree.order[i] = {code, uint32_t(i)};
      }
    });
    pool.wait();
    // Stable parallel LSD radix sort. Per-block histograms and disjoint
    // scatter intervals avoid atomics and preserve coincident-code order.
    std::vector<std::array<size_t, 64>> offsets(blocks);
    std::vector<MortonPrimitive> scratch(bounds.size());
    for (unsigned shift = 0; shift < 30; shift += 6) {
      pool.parallelFor(size_t(0), blocks, [&](size_t block) {
        offsets[block].fill(0);
        size_t end = std::min(bounds.size(), (block + 1) * block_size);
        for (size_t i = block * block_size; i < end; ++i)
          ++offsets[block][(tree.order[i].code >> shift) & 63];
      });
      pool.wait();
      size_t base = 0;
      for (int bucket = 0; bucket < 64; ++bucket)
        for (size_t block = 0; block < blocks; ++block) {
          size_t count = offsets[block][bucket];
          offsets[block][bucket] = base;
          base += count;
        }
      pool.parallelFor(size_t(0), blocks, [&](size_t block) {
        auto cursor = offsets[block];
        size_t end = std::min(bounds.size(), (block + 1) * block_size);
        for (size_t i = block * block_size; i < end; ++i) {
          auto p = tree.order[i];
          scratch[cursor[(p.code >> shift) & 63]++] = p;
        }
      });
      pool.wait();
      tree.order.swap(scratch);
    }
    for (size_t first = 0, end = 1; end <= bounds.size(); ++end) {
      if (end == bounds.size() || (tree.order[first].code >> 18) != (tree.order[end].code >> 18)) {
        tree.treelets.push_back({uint32_t(first), uint32_t(end - first)});
        first = end;
      }
    }
    pool.parallelFor(size_t(0), tree.treelets.size(), [&](size_t i) {
      const auto t = tree.treelets[i];
      uint32_t next = 2 * t.first;
      emit_treelet(tree, bounds, t.first, t.count, 17, max_leaf, next);
    });
    pool.join();
  }
  std::vector<uint32_t> roots;
  for (const auto &t : tree.treelets)
    roots.push_back(2 * t.first);
  uint32_t next = uint32_t(2 * bounds.size());
  tree.root = stitch(tree, roots, next);
  return tree;
}

#ifndef RAY_HAS_METAL_BVH
bool MetalBVHAvailable() { return false; }
void BuildMetalTreelets(std::span<const HLBVHBounds>, const HLBVHBounds &, unsigned, HLBVHTree &) {
  throw std::runtime_error("Metal BVH construction is unavailable in this build; use bvh_type = 'hlbvh'.");
}
#endif
