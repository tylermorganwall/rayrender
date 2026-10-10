#include "../core/bvh.h"
#include "bvh_timing.h"
#include "../utils/assert.h"
#include "../math/mathinline.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <new>
#include <cstring>
#include "../utils/raylog.h"
#include "../math/aabb.h"
#include "../math/simd.h"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <functional>
#include <queue>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>
#include <RcppThread.h>

#ifdef NOT_CRAN
#include <testthat.h>
#endif

namespace {

constexpr int kBVH4Empty = -1;

BVH4NodeStorage allocateBVH4NodeStorage(size_t count) {
  if (!count)
    return {};
  // allocator starts the array lifetime and supplies its required alignment;
  // each element's lifetime begins when its packing worker constructs it.
  return BVH4NodeStorage(std::allocator<LinearBVHNode4>{}.allocate(count), BVH4NodeDeleter{count});
}

bool isBVH4Leaf(int reference) { return reference < kBVH4Empty; }
int bvh4LeafReference(int index) { return -2 - index; }
int bvh4LeafIndex(int reference) { return -(reference + 2); }

struct BVHNodeEntry {
    int nodeIndex;
    float tEnter;

    bool operator<(const BVHNodeEntry& other) const {
        return tEnter > other.tEnter;
    }
};

static_assert(sizeof(BVHNodeEntry) == 8,
              "BVHNodeEntry should remain compact: int + float");

template <size_t MaxSize> class StaticPriorityQueue {
public:
  StaticPriorityQueue() : size_(0) {}

  bool try_push(const BVHNodeEntry &entry) {
    if (size_ >= MaxSize) {
      return false;
    }

    int i = static_cast<int>(size_) - 1;

    while (i >= 0 && data_[i].tEnter < entry.tEnter) {
      data_[i + 1] = data_[i];
      --i;
    }

    data_[i + 1] = entry;
    ++size_;
    return true;
  }

  // Merge up to four children in one pass instead of shifting the same pending
  // nodes once per child. Equal-distance children retain insertion order here,
  // so popping still visits the newest equal-distance entry first.
  bool try_push_batch(BVHNodeEntry* entries, size_t count) {
    if (count > MaxSize - size_) return false;
    if (count == 1) return try_push(entries[0]);
    for (size_t i = 1; i < count; ++i) {
      BVHNodeEntry entry = entries[i];
      size_t j = i;
      while (j > 0 && entries[j - 1].tEnter < entry.tEnter) {
        entries[j] = entries[j - 1];
        --j;
      }
      entries[j] = entry;
    }
    size_t old = size_, added = count, dest = size_ + count;
    while (added > 0) {
      if (old > 0 && data_[old - 1].tEnter < entries[added - 1].tEnter) {
        data_[--dest] = data_[--old];
      } else {
        data_[--dest] = entries[--added];
      }
    }
    size_ += count;
    return true;
  }

  void push(const BVHNodeEntry &entry) {
    const bool ok = try_push(entry);
    if (!ok) {
      ASSERT(false && "Priority queue overflow");
    }
  }

  BVHNodeEntry pop() {
    ASSERT(size_ > 0 && "Priority queue underflow");
    return data_[--size_];
  }

  const BVHNodeEntry &top() const {
    ASSERT(size_ > 0 && "Priority queue is empty");
    return data_[size_ - 1];
  }

  bool empty() const { return size_ == 0; }

  size_t size() const { return size_; }

private:
  BVHNodeEntry data_[MaxSize];
  size_t size_;
};

} // namespace

BVHAggregate::BVHAggregate(std::vector<std::shared_ptr<hitable>> prims, float t_min, float t_max, int maxPrimsInNode,
                           bool sah, Transform *ObjectToWorld, Transform *WorldToObject, bool reverseOrientation,
                           BVHBuildOptions options)
    : hitable(ObjectToWorld, WorldToObject, nullptr, reverseOrientation),
      maxPrimsInNode(std::clamp(maxPrimsInNode, 1, 255)), primitives(std::move(prims)) {
  build(t_min, t_max, options);
}

BVHAggregate::BVHAggregate(std::vector<std::shared_ptr<hitable>> prims, float t_min, float t_max, int maxPrimsInNode,
                           bool sah, BVHBuildOptions options)
    : maxPrimsInNode(std::clamp(maxPrimsInNode, 1, 255)), primitives(std::move(prims)) {
  build(t_min, t_max, options);
}

void BVHAggregate::build(float t_min, float t_max, BVHBuildOptions options) {
  SCOPED_CONTEXT("Initialization");
  SCOPED_TIMER_COUNTER("BVH Build");
  ScopedBVHBuildTiming build_timing;
  build_method = options.method;
  // Opt-in development timings separate parallel construction from the serial
  // primitive interface and final CPU traversal-layout preparation.
  const bool profile = std::getenv("RAYRENDER_HLBVH_PROFILE") != nullptr &&
                       (build_method == BVHBuildMethod::HLBVH || build_method == BVHBuildMethod::Metal);
  using Clock = std::chrono::steady_clock;
  auto checkpoint = profile ? Clock::now() : Clock::time_point{};
  auto stage = [&](const char *name) {
    if (!profile)
      return;
    const auto now = Clock::now();
    std::fprintf(stderr, "HLBVHHost n=%zu method=%u workers=%u stage=%s wall=%.6f\n", primitives.size(),
                 unsigned(build_method), options.threads, name,
                 std::chrono::duration<double>(now - checkpoint).count());
    checkpoint = now;
  };

  if (primitives.empty())
    return;
  std::vector<std::shared_ptr<hitable>> orderedPrims(primitives.size());
  if (build_method == BVHBuildMethod::HLBVH || build_method == BVHBuildMethod::Metal) {
    std::vector<HLBVHBounds> bounds(primitives.size());
    struct Preparation {
      aabb bounds;
      OpaqueShadowType shadow = OpaqueShadowType::Opaque;
    };
    auto combine_shadow = [](OpaqueShadowType a, OpaqueShadowType b) {
      if (a == OpaqueShadowType::Unsupported || b == OpaqueShadowType::Unsupported)
        return OpaqueShadowType::Unsupported;
      return a == OpaqueShadowType::Opaque && b == OpaqueShadowType::Opaque
                 ? OpaqueShadowType::Opaque : OpaqueShadowType::Mixed;
    };
    // Each block owns its compact bounds range and its reduction state. Keep
    // native-precision scene bounds separate from conservatively rounded GPU
    // bounds, and classify shadows while the primitive is already in cache.
    auto prepare = [&](size_t first, size_t end) {
      Preparation result;
      for (size_t i = first; i < end; ++i) {
        aabb box;
        if (!primitives[i]->bounding_box(t_min, t_max, box))
          throw std::runtime_error("BVH primitive has no bounding box.");
        result.bounds = surrounding_box(result.bounds, box);
        for (int d = 0; d < 3; ++d) {
          Float lo = box.min()[d], hi = box.max()[d];
          bounds[i].lo[d] = static_cast<float>(lo);
          bounds[i].hi[d] = static_cast<float>(hi);
          if (bounds[i].lo[d] > lo)
            bounds[i].lo[d] = std::nextafter(bounds[i].lo[d], -INFINITY);
          if (bounds[i].hi[d] < hi)
            bounds[i].hi[d] = std::nextafter(bounds[i].hi[d], INFINITY);
        }
        if (result.shadow != OpaqueShadowType::Unsupported)
          result.shadow = combine_shadow(result.shadow, primitives[i]->ShadowType());
      }
      return result;
    };
    Preparation prepared;
    if (options.parallel_primitive_queries && options.threads > 1 && primitives.size() >= 32768) {
      // Scheduling blocks amortizes pool overhead. Only explicitly audited
      // callers opt in; arbitrary/lazy hitables keep their original serial calls.
      constexpr size_t block_size = 4096;
      const size_t nblocks = (primitives.size() + block_size - 1) / block_size;
      std::vector<Preparation> blocks(nblocks);
      RcppThread::ThreadPool pool(options.threads);
      pool.parallelFor(0, nblocks, [&](size_t b) {
        blocks[b] = prepare(b * block_size, std::min((b + 1) * block_size, primitives.size()));
      });
      pool.join(); // Propagate missing-bounds errors before consuming any block.
      // Reduce in input order, independent of worker scheduling. Finite bounds
      // and the conservative shadow classification agree with the serial walk.
      for (const auto &block : blocks) {
        prepared.bounds = surrounding_box(prepared.bounds, block.bounds);
        prepared.shadow = combine_shadow(prepared.shadow, block.shadow);
      }
    } else {
      prepared = prepare(0, primitives.size());
    }
    scene_bounds = prepared.bounds;
    shadow_type = prepared.shadow;
    stage("primitive_preparation");
    HLBVHTree tree = BuildHLBVH(bounds, maxPrimsInNode, options);
    stage("tree_construction");
    // The sorted order is a permutation, so every slot transfers its ownership
    // exactly once. Moving avoids an atomic retain/release pair per primitive.
    for (size_t i = 0; i < primitives.size(); ++i)
      orderedPrims[i] = std::move(primitives[tree.order[i].index]);
    primitives.swap(orderedPrims);
    stage("primitive_reordering");
    packHLBVH(tree, options.threads);
    stage("layout_packing");
  } else {
    std::vector<BVHPrimitive> bvhPrimitives(primitives.size());
    for (size_t i = 0; i < primitives.size(); ++i) {
      aabb box;
      if (!primitives[i]->bounding_box(t_min, t_max, box))
        throw std::runtime_error("BVH primitive has no bounding box.");
      bvhPrimitives[i] = BVHPrimitive(i, box);
      scene_bounds = surrounding_box(scene_bounds, box);
    }
    std::atomic<int> totalNodes{0}, orderedPrimsOffset{0};
    BVHBuildNode *root = buildRecursive(bvhPrimitives, &totalNodes, &orderedPrimsOffset, orderedPrims);
    primitives.swap(orderedPrims);
    std::vector<BVHPrimitive>().swap(bvhPrimitives);
    std::vector<std::shared_ptr<hitable>>().swap(orderedPrims);
#ifndef RAYSIMD
    nodes.reset(new LinearBVHNode[totalNodes]);
    n_nodes = totalNodes;
    int offset = 0;
    flattenBVH(root, &offset);
    delete root;
#else
    buildBVH4(root);
#endif
  }
  stage("temporary_release");
  if (build_method != BVHBuildMethod::HLBVH && build_method != BVHBuildMethod::Metal)
    classifyOpaqueShadow();
}

// Pack indexed construction nodes directly, avoiding two pointer-tree copies
// and millions of individual allocations. Collapse the same two binary levels
// as the existing BVH4 conversion; all traversal paths keep their current ABI.
void BVHAggregate::packHLBVH(const HLBVHTree &tree, unsigned threads) {
  // Count only reachable nodes: the construction arena also has gaps between
  // treelets. Keep this pass free of bounds work and type-erased recursion.
#ifndef RAYSIMD
  auto count = [&](auto &&self, uint32_t i) -> void {
    ++n_nodes;
    const auto &n = tree.nodes[i];
    if (!n.count) {
      self(self, n.left);
      self(self, n.right);
    }
  };
  count(count, tree.root);
  nodes.reset(new LinearBVHNode[n_nodes]);
  n_nodes = 0;
  auto pack = [&](auto &&self, uint32_t i) -> int {
    const int index = n_nodes++;
    const auto &n = tree.nodes[i];
    auto &out = nodes[index];
    out.bounds.bounds[0] = point3f(n.bounds.lo[0], n.bounds.lo[1], n.bounds.lo[2]);
    out.bounds.bounds[1] = point3f(n.bounds.hi[0], n.bounds.hi[1], n.bounds.hi[2]);
    out.nPrimitives = n.count;
    out.axis = n.axis;
    if (n.count) {
      out.primitivesOffset = n.first;
    } else {
      self(self, n.left);
      out.secondChildOffset = self(self, n.right);
    }
    return index;
  };
  pack(pack, tree.root);
#else
  // Collapse two binary levels into one four-wide node, preserving the old
  // depth-first child order (including its tie-breaking during traversal).
  auto children = [&](uint32_t i, uint32_t (&indices)[4]) {
    int count = 0;
    const auto &n = tree.nodes[i];
    for (uint32_t child : {n.left, n.right}) {
      const auto &c = tree.nodes[child];
      if (c.count) {
        indices[count++] = child;
      } else {
        indices[count++] = c.left;
        indices[count++] = c.right;
      }
    }
    return count;
  };
  const bool profile = std::getenv("RAYRENDER_HLBVH_PROFILE") != nullptr;
  using Clock = std::chrono::steady_clock;
  auto checkpoint = profile ? Clock::now() : Clock::time_point{};
  auto stage = [&](const char *name) {
    if (!profile)
      return;
    const auto now = Clock::now();
    std::fprintf(stderr, "HLBVHPack n=%zu workers=%u stage=%s wall=%.6f\n", tree.order.size(), threads, name,
                 std::chrono::duration<double>(now - checkpoint).count());
    checkpoint = now;
  };
  struct Counts {
    int nodes = 0, leaves = 0;
  };
  struct PackTask {
    uint32_t root;
    Counts size, offset;
  };
  std::vector<PackTask> tasks;
  // Four BVH4 levels give at most 256 independent subtrees. Count them in
  // parallel, then reserve their exact depth-first spans on the main thread.
  // Workers write disjoint ranges without atomics or schedule-dependent order.
  // Small trees stay serial to avoid paying for another thread pool.
  constexpr unsigned frontier_depth = 4;
  const bool parallel = threads > 1 && tree.order.size() >= 32768;
  auto count = [&](auto &&self, uint32_t i, Counts &size, unsigned depth) -> void {
    if (tree.nodes[i].count) {
      ++size.leaves;
      return;
    }
    if (parallel && depth == frontier_depth) {
      tasks.push_back({i, {}, {}});
      return;
    }
    ++size.nodes;
    uint32_t indices[4];
    const int nchildren = children(i, indices);
    for (int j = 0; j < nchildren; ++j)
      self(self, indices[j], size, depth + 1);
  };
  Counts size;
  count(count, tree.root, size, 0);
  std::unique_ptr<RcppThread::ThreadPool> pool;
  if (!tasks.empty()) {
    pool = std::make_unique<RcppThread::ThreadPool>(threads);
    pool->parallelFor(0, tasks.size(), [&](size_t i) {
      // Start below the frontier to count each task completely.
      count(count, tasks[i].root, tasks[i].size, frontier_depth + 1);
    });
    pool->wait();
    for (const auto &task : tasks) {
      size.nodes += task.size.nodes;
      size.leaves += task.size.leaves;
    }
  }
  totalNodes4 = size.nodes;
  totalLeaves4 = size.leaves;
  stage("count");
  nodes4 = allocateBVH4NodeStorage(totalNodes4);
  leaves4.reset(new LinearBVHLeaf4[totalLeaves4]);
#ifdef NOT_CRAN
  const bool poison_storage = std::getenv("RAYRENDER_HLBVH_POISON_PACKING") != nullptr;
  if (poison_storage) {
    // Development-only check: stale bounds become NaNs and stale references
    // become -1. Every field used by traversal must be overwritten below.
    if (totalNodes4)
      std::memset(static_cast<void*>(nodes4.get()), 0xff, sizeof(LinearBVHNode4) * totalNodes4);
    std::memset(leaves4.get(), 0xff, sizeof(LinearBVHLeaf4) * totalLeaves4);
  }
#endif
  stage("allocate");
  size_t next_task = 0;
  auto pack = [&](auto &&self, uint32_t i, Counts &next, unsigned depth) -> int {
    const auto &n = tree.nodes[i];
    if (n.count) {
      const int leaf = next.leaves++;
      leaves4[leaf] = {int(n.first), int(n.count)};
      return bvh4LeafReference(leaf);
    }
    if (parallel && depth == frontier_depth) {
      auto &task = tasks[next_task++];
      task.offset = next;
      next.nodes += task.size.nodes;
      next.leaves += task.size.leaves;
      return task.offset.nodes;
    }
    const int index = next.nodes++;
    // Every reserved slot has exactly one writer. Start its object lifetime on
    // that worker so serial allocation does not first-touch the entire array.
    ::new (static_cast<void*>(nodes4.get() + index)) LinearBVHNode4;
    auto &out = nodes4[index];
    uint32_t indices[4];
    out.nChildren = children(i, indices);
    for (int j = 0; j < out.nChildren; ++j) {
      out.childOffsets[j] = self(self, indices[j], next, depth + 1);
      const auto &bounds = tree.nodes[indices[j]].bounds;
      // Construction already guarantees ordered, conservative float bounds.
      // Write the traversal lanes directly instead of converting through
      // four aabbs and reordering their min/max coordinates a second time.
      for (int d = 0; d < 3; ++d) {
        out.bbox4.corners[2 * d].xyzw[j] = bounds.lo[d];
        out.bbox4.corners[2 * d + 1].xyzw[j] = bounds.hi[d];
      }
    }
    for (int j = out.nChildren; j < 4; ++j) {
      out.childOffsets[j] = kBVH4Empty;
      // Initialize inactive lanes for SIMD loads; their absent-child references
      // mask out the intersection results, including any parallel-ray NaNs.
      for (auto &corner : out.bbox4.corners)
        corner.xyzw[j] = 0;
    }
    return index;
  };
  Counts next;
  root4 = pack(pack, tree.root, next, 0);
  if (pool) {
    pool->parallelFor(0, tasks.size(), [&](size_t i) {
      Counts cursor = tasks[i].offset;
      pack(pack, tasks[i].root, cursor, frontier_depth + 1);
    });
    pool->join();
  }
  stage("write");
#ifdef NOT_CRAN
  if (poison_storage) {
    validateBVH4();
    for (int i = 0; i < totalNodes4; ++i)
      for (const auto &corner : nodes4[i].bbox4.corners)
        for (int lane = 0; lane < 4; ++lane)
          if (!std::isfinite(corner.xyzw[lane]) ||
              (lane >= nodes4[i].nChildren && corner.xyzw[lane] != 0))
            throw std::runtime_error("HLBVH packing left an invalid bounds lane.");
  }
#endif
#ifndef NDEBUG
  validateBVH4();
#endif
#endif
}

bool BVHAggregate::bounding_box(Float t0, Float t1, aabb &box) const {
  if (primitives.empty()) {
    return (false);
  }
  box = scene_bounds;
  return (true);
}

void BVHAggregate::classifyOpaqueShadow() {
    shadow_type = OpaqueShadowType::Opaque;
    for (const auto& primitive : primitives) {
        const auto type = primitive->ShadowType();
        if (type == OpaqueShadowType::Unsupported) {
            shadow_type = type;
            return;
        }
        if (type != OpaqueShadowType::Opaque) shadow_type = OpaqueShadowType::Mixed;
    }
}

// Opaque visibility needs any blocker, so order only the current node's children
// and keep pending work on a small LIFO stack. Spill safely for unusually deep
// trees; ordinary closest-hit and stochastic HitP traversal remain unchanged.
bool BVHAggregate::OpaqueHit(const Ray& r, Float t_min, Float t_max, random_gen& rng) const {
    int pending[64];
    size_t count = 0;
    std::vector<int> overflow;
    auto push = [&](int index) {
        if (count < 64 && overflow.empty()) pending[count++] = index;
        else overflow.push_back(index);
    };
    auto pop = [&]() {
        if (overflow.empty()) return pending[--count];
        int index = overflow.back();
        overflow.pop_back();
        return index;
    };
#ifdef RAYSIMD
    if (root4 == kBVH4Empty) return false;
    const RayBBox4 rbox(r);
    push(root4);
    while (count || !overflow.empty()) {
        const int index = pop();
        if (isBVH4Leaf(index)) {
            const auto& leaf = leaves4[bvh4LeafIndex(index)];
            for (int i = 0; i < leaf.nPrimitives; ++i)
                if (primitives[leaf.primitivesOffset + i]->OpaqueHit(r, t_min, t_max, rng))
                    return true;
        } else {
            const auto& node = nodes4[index];
            IVec4 hits;
            FVec4 enters;
            rayBBoxIntersect4(rbox, node.bbox4, t_min, t_max, hits, enters);
            int mask = simd_extract_hitmask(simd_and(hits, simd_not_equals_minus_one(node.childOffsets)));
            if (!mask) continue;
            float distances[4];
            simd_extract_fvec4(enters, distances);
            BVHNodeEntry children[4];
            int n = 0;
            for (int lane = 0; lane < 4; ++lane) {
                if (!(mask & (1 << lane))) continue;
                BVHNodeEntry child{node.childOffsets[lane], distances[lane]};
                int j = n++;
                while (j && children[j - 1].tEnter < child.tEnter) {
                    children[j] = children[j - 1];
                    --j;
                }
                children[j] = child;
            }
            for (int i = 0; i < n; ++i) push(children[i].nodeIndex);
        }
    }
#else
    if (!nodes) return false;
    push(0);
    while (count || !overflow.empty()) {
        const int index = pop();
        const auto& node = nodes[index];
        if (!node.bounds.hit(r, t_min, t_max, rng)) continue;
        if (node.nPrimitives) {
            for (int i = 0; i < node.nPrimitives; ++i)
                if (primitives[node.primitivesOffset + i]->OpaqueHit(r, t_min, t_max, rng))
                    return true;
        } else if (r.d[node.axis] < 0) {
            push(index + 1);
            push(node.secondChildOffset);
        } else {
            push(node.secondChildOffset);
            push(index + 1);
        }
    }
#endif
    return false;
}

BVHBuildNode *BVHAggregate::buildRecursive(std::span<BVHPrimitive> bvhPrimitives,
                                           std::atomic<int> *totalNodes,
                                           std::atomic<int> *orderedPrimsOffset,
                                           std::vector<std::shared_ptr<hitable> > &orderedPrims) {
    // DCHECK_NE(bvhPrimitives.size(), 0);
    // Allocator alloc = threadAllocators.Get();
    // BVHBuildNode *node = alloc.new_object<BVHBuildNode>();
    // std::unique_ptr<BVHBuildNode> node = std::make_unique<BVHBuildNode>();
    BVHBuildNode *node = new BVHBuildNode();

    // Initialize _BVHBuildNode_ for primitive range
    bool isRoot = (*totalNodes == 0);
    (*totalNodes)++;
    // Compute bounds of all primitives in BVH node
    aabb bounds;
    for (const auto &prim : bvhPrimitives) {
        bounds = surrounding_box(bounds, prim.bounds);
    }

    if ((bounds.surface_area() == 0 || bvhPrimitives.size() <= static_cast<size_t>(maxPrimsInNode)) && !isRoot) {
        // Create leaf _BVHBuildNode_
        int firstPrimOffset = orderedPrimsOffset->fetch_add(bvhPrimitives.size());
        for (size_t i = 0; i < bvhPrimitives.size(); ++i) {
            int index = bvhPrimitives[i].primitiveIndex;
            orderedPrims[firstPrimOffset + i] = primitives[index];
        }
        node->InitLeaf(firstPrimOffset, bvhPrimitives.size(), bounds);
        return node;
    } else {
        // Compute bound of primitive centroids and choose split dimension _dim_
        aabb centroidBounds;
        for (const auto &prim : bvhPrimitives)
            centroidBounds = surrounding_box(centroidBounds, prim.Centroid());
        int dim = centroidBounds.MaxDimension();

        // Partition primitives into two sets and build children
        if (centroidBounds.max()[dim] == centroidBounds.min()[dim]) {
            // Create leaf _BVHBuildNode_
            int firstPrimOffset = orderedPrimsOffset->fetch_add(bvhPrimitives.size());
            for (size_t i = 0; i < bvhPrimitives.size(); ++i) {
                int index = bvhPrimitives[i].primitiveIndex;
                orderedPrims[firstPrimOffset + i] = primitives[index];
            }
            node->InitLeaf(firstPrimOffset, bvhPrimitives.size(), bounds);
            return node;

        } else {
            int mid = bvhPrimitives.size() / 2;
            // Partition primitives using approximate SAH
            if (build_method == BVHBuildMethod::Equal || bvhPrimitives.size() <= 2) {
                // Partition primitives into equally sized subsets
                mid = bvhPrimitives.size() / 2;
                std::nth_element(bvhPrimitives.begin(), bvhPrimitives.begin() + mid,
                                    bvhPrimitives.end(),
                                    [dim](const BVHPrimitive &a, const BVHPrimitive &b) {
                                        return a.Centroid()[dim] < b.Centroid()[dim];
                                    });

            } else {
                // Allocate _BVHSplitBucket_ for SAH partition buckets
                constexpr int nBuckets = 12;
                BVHSplitBucket buckets[nBuckets];

                // Initialize _BVHSplitBucket_ for SAH partition buckets
                for (const auto &prim : bvhPrimitives) {
                    int b = nBuckets * centroidBounds.offset(prim.Centroid())[dim];
                    if (b == nBuckets) {
                        b = nBuckets - 1;
                    }
                    // DCHECK_GE(b, 0);
                    // DCHECK_LT(b, nBuckets);
                    buckets[b].count++;
                    buckets[b].bounds = surrounding_box(buckets[b].bounds, prim.bounds);
                }

                // Compute costs for splitting after each bucket
                constexpr int nSplits = nBuckets - 1;
                Float costs[nSplits] = {};
                // Partially initialize _costs_ using a forward scan over splits
                int countBelow = 0;
                aabb boundBelow;
                for (int i = 0; i < nSplits; ++i) {
                    boundBelow = surrounding_box(boundBelow, buckets[i].bounds);
                    countBelow += buckets[i].count;
                    costs[i] += countBelow * boundBelow.surface_area();
                }

                // Finish initializing _costs_ using a backward scan over splits
                int countAbove = 0;
                aabb boundAbove;
                for (int i = nSplits; i >= 1; --i) {
                    boundAbove = surrounding_box(boundAbove, buckets[i].bounds);
                    countAbove += buckets[i].count;
                    costs[i - 1] += countAbove * boundAbove.surface_area();
                }

                // Find bucket to split at that minimizes SAH metric
                int minCostSplitBucket = -1;
                Float minCost = INFINITY;
                for (int i = 0; i < nSplits; ++i) {
                    // Compute cost for candidate split and update minimum if
                    // necessary
                    if (costs[i] < minCost) {
                        minCost = costs[i];
                        minCostSplitBucket = i;
                    }
                }
                // Compute leaf cost and SAH split cost for chosen split
                Float leafCost = bvhPrimitives.size();
                minCost = 1.f / 2.f + minCost / bounds.surface_area();

                // Either create leaf or split primitives at selected SAH bucket
                if (bvhPrimitives.size() > static_cast<size_t>(maxPrimsInNode) || minCost < leafCost) {
                    auto midIter = std::partition(
                        bvhPrimitives.begin(), bvhPrimitives.end(),
                        [=](const BVHPrimitive &bp) {
                            int b = nBuckets * centroidBounds.offset(bp.Centroid())[dim];
                            if (b == nBuckets) {
                                b = nBuckets - 1;
                            }
                            return b <= minCostSplitBucket;
                        });
                    mid = static_cast<int>(midIter - bvhPrimitives.begin());
                    if (mid == 0 ||
                        mid == static_cast<int>(bvhPrimitives.size())) {
                      mid = static_cast<int>(bvhPrimitives.size() / 2);
                      std::nth_element(
                          bvhPrimitives.begin(), bvhPrimitives.begin() + mid,
                          bvhPrimitives.end(),
                          [dim](const BVHPrimitive &a, const BVHPrimitive &b) {
                            return a.Centroid()[dim] < b.Centroid()[dim];
                          });
                    }
                } else {
                    // Create leaf _BVHBuildNode_
                    int firstPrimOffset =
                        orderedPrimsOffset->fetch_add(bvhPrimitives.size());
                    for (size_t i = 0; i < bvhPrimitives.size(); ++i) {
                        int index = bvhPrimitives[i].primitiveIndex;
                        orderedPrims[firstPrimOffset + i] = primitives[index];
                    }
                    node->InitLeaf(firstPrimOffset, bvhPrimitives.size(), bounds);
                    return node;
                }
            }
            BVHBuildNode *children[2];

            // Recursively build child BVHs sequentially
            children[0] =
                buildRecursive(bvhPrimitives.subspan(0, mid),
                               totalNodes, orderedPrimsOffset, orderedPrims);
            children[1] =
                buildRecursive(bvhPrimitives.subspan(mid),
                                totalNodes, orderedPrimsOffset, orderedPrims);

            node->InitInterior(dim, children[0], children[1]);
        }
    }
    return node;
}

int BVHAggregate::flattenBVH(BVHBuildNode *node, int *offset) {
    LinearBVHNode *linearNode = &nodes[*offset];
    linearNode->bounds = node->bounds;
    int nodeOffset = (*offset)++;
    if (node->nPrimitives > 0) {
        // CHECK_LT(node->nPrimitives, 65536);
        linearNode->primitivesOffset = node->firstPrimOffset;
        linearNode->nPrimitives = node->nPrimitives;
    } else {
        // Create interior flattened BVH node
        linearNode->axis = node->splitAxis;
        linearNode->nPrimitives = 0;
        flattenBVH(node->children[0], offset);
        linearNode->secondChildOffset = flattenBVH(node->children[1], offset);
    }
    return nodeOffset;
}


#ifndef RAYSIMD
namespace {
inline bool rayInvDirIsNeg(const Ray& r, int axis) {
    return r.inv_dir_is_neg[axis] != 0;
}

inline bool rayBoundsHitTEnter(
    const Ray& r,
    const aabb& bounds,
    Float t_min,
    Float t_max,
    float& tEnter
) {
    const int sx = r.inv_dir_is_neg[0] ? 1 : 0;
    const int sy = r.inv_dir_is_neg[1] ? 1 : 0;
    const int sz = r.inv_dir_is_neg[2] ? 1 : 0;

    const Float txNear =
        (bounds.bounds[sx].e[0] - r.o.e[0]) * r.inv_dir_pad.e[0];
    const Float txFar =
        (bounds.bounds[1 - sx].e[0] - r.o.e[0]) * r.inv_dir_pad.e[0];

    const Float tyNear =
        (bounds.bounds[sy].e[1] - r.o.e[1]) * r.inv_dir_pad.e[1];
    const Float tyFar =
        (bounds.bounds[1 - sy].e[1] - r.o.e[1]) * r.inv_dir_pad.e[1];

    const Float tzNear =
        (bounds.bounds[sz].e[2] - r.o.e[2]) * r.inv_dir_pad.e[2];
    const Float tzFar =
        (bounds.bounds[1 - sz].e[2] - r.o.e[2]) * r.inv_dir_pad.e[2];

    const Float t_near = std::fmax(std::fmax(txNear, tyNear), tzNear);
    const Float t_far = std::fmin(std::fmin(txFar, tyFar), tzFar);
    t_min = std::fmax(t_min, bounds_entry_lower(t_near));
    t_max = std::fmin(t_max, bounds_exit_upper(t_far));

    tEnter = static_cast<float>(t_min);
    return t_min <= t_max;
}

template <typename HitPrimitive>
bool traverseAnyPriorityBVH2(
    const LinearBVHNode* nodes,
    const Ray& r,
    Float t_min,
    Float t_max,
    HitPrimitive&& hitPrimitive
) {
    if (!nodes) {
        return false;
    }

    float rootTEnter;
    if (!rayBoundsHitTEnter(r, nodes[0].bounds, t_min, t_max, rootTEnter)) {
        return false;
    }

    std::priority_queue<BVHNodeEntry> frontier;
    frontier.push({0, rootTEnter});

    while (!frontier.empty()) {
        const BVHNodeEntry entry = frontier.top();
        frontier.pop();

        const LinearBVHNode* node = &nodes[entry.nodeIndex];
        if (node->nPrimitives > 0) {
            for (int i = 0; i < node->nPrimitives; ++i) {
                const int primIndex = node->primitivesOffset + i;

                if (hitPrimitive(primIndex, t_min, t_max)) {
                    return true;
                }
            }

            continue;
        }

        const int firstChildOffset = entry.nodeIndex + 1;
        const int secondChildOffset = node->secondChildOffset;
        float firstTEnter;
        float secondTEnter;

        if (rayBoundsHitTEnter(
                r, nodes[firstChildOffset].bounds, t_min, t_max, firstTEnter
            )) {
            frontier.push({firstChildOffset, firstTEnter});
        }

        if (rayBoundsHitTEnter(
                r, nodes[secondChildOffset].bounds, t_min, t_max, secondTEnter
            )) {
            frontier.push({secondChildOffset, secondTEnter});
        }
    }

    return false;
}
} // namespace

const bool BVHAggregate::hit(const Ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const {
    // SCOPED_CONTEXT("Hit");
    // SCOPED_TIMER_COUNTER("BVH Serial");
    if (!nodes) {
        return false;
    }
    // Follow ray through BVH nodes to find primitive intersections
    int toVisitOffset = 0, currentNodeIndex = 0;
    int nodesToVisit[64];
    // int nodesVisited = 0;
    bool any_hit = false;
    while (true) {
        // ++nodesVisited;
        const LinearBVHNode *node = &nodes[currentNodeIndex];
        // Check ray against BVH node
        if (node->bounds.hit(r, t_min, t_max, rng)) {
            if (node->nPrimitives > 0) {
                // Intersect ray with primitives in leaf BVH node
                for (int i = 0; i < node->nPrimitives; ++i) {
                    // Check for intersection with primitive in BVH node
                    hit_record hrec_temp;
                    bool prim_hrec = primitives[node->primitivesOffset + i]->hit(r, t_min, t_max, hrec_temp, rng);
                    if (prim_hrec && (!r.segment_absorption || !any_hit ||
                                      hrec_temp.OrderedDistance() <= rec.OrderedDistance())) {
                        any_hit = true;
                        rec = hrec_temp;
                        t_max = r.segment_absorption ? rec.DistanceUpperBound() : rec.t;
                    }
                }
                if (toVisitOffset == 0) {
                    break;
                }
                currentNodeIndex = nodesToVisit[--toVisitOffset];
            } else {
                // Put far BVH node on _nodesToVisit_ stack, advance to near node
                if (rayInvDirIsNeg(r, node->axis)) {
                    nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
                    currentNodeIndex = node->secondChildOffset;
                } else {
                    nodesToVisit[toVisitOffset++] = node->secondChildOffset;
                    currentNodeIndex = currentNodeIndex + 1;
                }
            }
        } else {
            if (toVisitOffset == 0) {
                break;
            }
            currentNodeIndex = nodesToVisit[--toVisitOffset];
        }
    }

    // bvhNodesVisited += nodesVisited;
    return any_hit;
}

const bool BVHAggregate::hit(const Ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const {
    // SCOPED_CONTEXT("Hit");
    // SCOPED_TIMER_COUNTER("BVH Serial");
    if (!nodes) {
        return false;
    }
    // Follow ray through BVH nodes to find primitive intersections
    int toVisitOffset = 0, currentNodeIndex = 0;
    int nodesToVisit[64];
    // int nodesVisited = 0;
    bool any_hit = false;
    while (true) {
        // ++nodesVisited;
        const LinearBVHNode *node = &nodes[currentNodeIndex];
        // Check ray against BVH node
        if (node->bounds.hit(r, t_min, t_max, sampler)) {
            if (node->nPrimitives > 0) {
                // Intersect ray with primitives in leaf BVH node
                for (int i = 0; i < node->nPrimitives; ++i) {
                    // Check for intersection with primitive in BVH node
                    hit_record hrec_temp;
                    bool prim_hrec = primitives[node->primitivesOffset + i]->hit(r, t_min, t_max, hrec_temp, sampler);
                    if (prim_hrec && (!r.segment_absorption || !any_hit ||
                                      hrec_temp.OrderedDistance() <= rec.OrderedDistance())) {
                        any_hit = true;
                        rec = hrec_temp;
                        t_max = r.segment_absorption ? rec.DistanceUpperBound() : rec.t;
                    }
                }
                if (toVisitOffset == 0) {
                    break;
                }
                currentNodeIndex = nodesToVisit[--toVisitOffset];
            } else {
                // Put far BVH node on _nodesToVisit_ stack, advance to near node
                if (rayInvDirIsNeg(r, node->axis)) {
                    nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
                    currentNodeIndex = node->secondChildOffset;
                } else {
                    nodesToVisit[toVisitOffset++] = node->secondChildOffset;
                    currentNodeIndex = currentNodeIndex + 1;
                }
            }
        } else {
            if (toVisitOffset == 0) {
                break;
            }
            currentNodeIndex = nodesToVisit[--toVisitOffset];
        }
    }

    // bvhNodesVisited += nodesVisited;
    return any_hit;
}

bool BVHAggregate::HitP(const Ray& r, Float t_min, Float t_max, random_gen& rng) const {
    return traverseAnyPriorityBVH2(
        nodes.get(),
        r,
        t_min,
        t_max,
        [&](int primIndex, Float local_t_min, Float local_t_max) {
            return primitives[primIndex]->HitP(
                r, local_t_min, local_t_max, rng
            );
        }
    );
}

bool BVHAggregate::HitP(const Ray& r, Float t_min, Float t_max, Sampler* sampler) const {
    return traverseAnyPriorityBVH2(
        nodes.get(),
        r,
        t_min,
        t_max,
        [&](int primIndex, Float local_t_min, Float local_t_max) {
            return primitives[primIndex]->HitP(
                r, local_t_min, local_t_max, sampler
            );
        }
    );
}
#else

static_assert(SIMD_WIDTH == 4,
              "BVH4 traversal requires SIMD_WIDTH == 4.");

constexpr size_t kBVH4PriorityStaticCapacity = 4096;

namespace {

template <size_t StaticCapacity>
class BVH4PriorityFrontier {
public:
  void push(const BVHNodeEntry& entry) {
    if (using_heap_) {
      heap_.push(entry);
      return;
    }

    if (!static_.try_push(entry)) {
      spill_to_heap();
      heap_.push(entry);
    }
  }

  bool empty() const {
    return using_heap_ ? heap_.empty() : static_.empty();
  }

  BVHNodeEntry pop() {
    if (using_heap_) {
      ASSERT(!heap_.empty() && "Priority frontier underflow");
      BVHNodeEntry entry = heap_.top();
      heap_.pop();
      return entry;
    }

    return static_.pop();
  }

  float top_t() const {
    ASSERT(!empty() && "Priority frontier is empty");
    return using_heap_ ? heap_.top().tEnter : static_.top().tEnter;
  }

  bool using_heap() const {
    return using_heap_;
  }

  BVHNodeEntry next_child(const IVec4& offsets, const float* distances, int hitmask) {
    BVHNodeEntry children[4];
    size_t count = 0;
    for (int i = 0; i < 4; ++i) {
      if ((hitmask >> i) & 1) children[count++] = {offsets[i], distances[i]};
    }
    ASSERT(count > 0);

    // Follow an already-nearest single child without pushing and immediately
    // popping it. At capacity, retain the original spill and heap tie behavior.
    if (count == 1 && !using_heap_ && static_.size() < StaticCapacity &&
        (static_.empty() || children[0].tEnter <= static_.top().tEnter)) {
      return children[0];
    }
    if (!using_heap_ && static_.try_push_batch(children, count)) {
      return static_.pop();
    }
    // A batch that would overflow must spill at exactly the same insertion as
    // before. try_push_batch leaves the children untouched when it cannot fit.
    for (size_t i = 0; i < count; ++i) push(children[i]);
    return pop();
  }

private:
  void spill_to_heap() {
    ASSERT(!using_heap_);

    while (!static_.empty()) {
      heap_.push(static_.pop());
    }

    using_heap_ = true;
  }

  StaticPriorityQueue<StaticCapacity> static_;
  std::priority_queue<BVHNodeEntry> heap_;
  bool using_heap_ = false;
};

template <bool ContainmentOnly = false, typename IntersectPrimitive>
bool traverseClosestBVH4(
    const LinearBVHNode4* nodes4,
    const LinearBVHLeaf4* leaves4,
    int root4,
    const Ray& r,
    Float t_min,
    Float t_max,
    hit_record& rec,
    IntersectPrimitive&& intersectPrimitive
) {
  if (root4 == kBVH4Empty) {
    return false;
  }

  const RayBBox4 rbox(r);
  BVH4PriorityFrontier<kBVH4PriorityStaticCapacity> frontier;
  BVHNodeEntry entry{root4, -std::numeric_limits<float>::infinity()};

  bool any_hit = false;

  while (true) {
    if (isBVH4Leaf(entry.nodeIndex)) {
      const LinearBVHLeaf4& leaf = leaves4[bvh4LeafIndex(entry.nodeIndex)];
      for (int i = 0; i < leaf.nPrimitives; ++i) {
        hit_record tempRec;
        const int primIndex = leaf.primitivesOffset + i;

        if (intersectPrimitive(primIndex, t_min, t_max, tempRec) &&
            (!r.segment_absorption || !any_hit ||
             tempRec.OrderedDistance() <= rec.OrderedDistance())) {
          any_hit = true;
          rec = tempRec;
          t_max = r.segment_absorption ? rec.DistanceUpperBound() : rec.t;
        }
      }

      if (any_hit) {
        if (frontier.empty() || frontier.top_t() > t_max) {
          return true;
        }
      }

    } else {
      const LinearBVHNode4* node = &nodes4[entry.nodeIndex];
      IVec4 hits;
      FVec4 tEnters;
      rayBBoxIntersect4(rbox, node->bbox4, t_min, t_max, hits, tEnters);
      const IVec4 valid_hit =
          simd_and(hits, simd_not_equals_minus_one(node->childOffsets));
      int hitmask = simd_extract_hitmask(valid_hit);
      if constexpr (ContainmentOnly) {
        // Keep the ordinary closest-first ordering, but prune closed objects
        // that cannot contain the origin. Pad inclusively to retain rounded
        // contacts at faces/edges and transformed bounds. This restriction is
        // local to this aggregate: a candidate's mesh must still test ALL faces.
        for (int axis = 0; axis < 3; ++axis) {
          const Float padding = 8 * std::numeric_limits<Float>::epsilon() *
                                std::max(Float(1), std::abs(r.o[axis]));
          for (int child = 0; child < 4; ++child) {
            if (r.o[axis] < node->bbox4.corners[2 * axis][child] - padding ||
                r.o[axis] > node->bbox4.corners[2 * axis + 1][child] + padding)
              hitmask &= ~(1 << child);
          }
        }
      }
      if (hitmask != 0) {
        float tEntersArray[4];
        simd_extract_fvec4(tEnters, tEntersArray);
        entry = frontier.next_child(node->childOffsets, tEntersArray, hitmask);
        continue;
      }
    }
    if (frontier.empty()) break;
    entry = frontier.pop();
  }

  return any_hit;
}

template <typename HitPrimitive>
bool traverseAnyPriorityBVH4(
    const LinearBVHNode4* nodes4,
    const LinearBVHLeaf4* leaves4,
    int root4,
    const Ray& r,
    Float t_min,
    Float t_max,
    HitPrimitive&& hitPrimitive
) {
  // Keep any-hit traversal closest-first. Some primitive HitP paths can call
  // stochastic hit logic for alpha/transparency, so changing traversal order
  // would also change RNG/sampler consumption and rendered results.
  if (root4 == kBVH4Empty) {
    return false;
  }

  const RayBBox4 rbox(r);
  BVH4PriorityFrontier<kBVH4PriorityStaticCapacity> frontier;
  BVHNodeEntry entry{root4, -std::numeric_limits<float>::infinity()};

  while (true) {
    if (isBVH4Leaf(entry.nodeIndex)) {
      const LinearBVHLeaf4& leaf = leaves4[bvh4LeafIndex(entry.nodeIndex)];
      for (int i = 0; i < leaf.nPrimitives; ++i) {
        const int primIndex = leaf.primitivesOffset + i;

        if (hitPrimitive(primIndex, t_min, t_max)) {
          return true;
        }
      }

    } else {
      const LinearBVHNode4* node = &nodes4[entry.nodeIndex];
      IVec4 hits;
      FVec4 tEnters;
      rayBBoxIntersect4(rbox, node->bbox4, t_min, t_max, hits, tEnters);
      const IVec4 valid_hit =
          simd_and(hits, simd_not_equals_minus_one(node->childOffsets));
      const int hitmask = simd_extract_hitmask(valid_hit);
      if (hitmask != 0) {
        float tEntersArray[4];
        simd_extract_fvec4(tEnters, tEntersArray);
        entry = frontier.next_child(node->childOffsets, tEntersArray, hitmask);
        continue;
      }
    }
    if (frontier.empty()) break;
    entry = frontier.pop();
  }

  return false;
}

} // namespace

bool BVHAggregate::HitContainingObjects(const Ray& r, Float t_min, Float t_max,
                                       hit_record& rec, random_gen& rng,
                                       const std::vector<const hitable*>* excluded) const {
  return traverseClosestBVH4<true>(
      nodes4.get(), leaves4.get(), root4, r, t_min, t_max, rec,
      [&](int primIndex, Float lo, Float hi, hit_record& result) {
        // Leaves can contain several coincident-centroid objects, even when
        // maxPrimsInNode is one. Check each whole object's bounds as well.
        const auto& object = primitives[primIndex];
        if (excluded && std::find(excluded->begin(), excluded->end(), object.get()) != excluded->end())
          return false;
        aabb bounds;
        if (object->bounding_box(r.time(), r.time(), bounds)) {
          for (int axis = 0; axis < 3; ++axis) {
            const Float padding = 8 * std::numeric_limits<Float>::epsilon() *
                                  std::max(Float(1), std::abs(r.o[axis]));
            if (r.o[axis] < bounds.min()[axis] - padding ||
                r.o[axis] > bounds.max()[axis] + padding)
              return false;
          }
        }
        return object->hit(r, lo, hi, result, rng);
      });
}

const bool BVHAggregate::hit(
    const Ray& r,
    Float t_min,
    Float t_max,
    hit_record& rec,
    random_gen& rng
) const {
  return traverseClosestBVH4(
      nodes4.get(),
      leaves4.get(),
      root4,
      r,
      t_min,
      t_max,
      rec,
      [&](int primIndex, Float local_t_min, Float local_t_max,
          hit_record& tempRec) {
        return primitives[primIndex]->hit(
            r, local_t_min, local_t_max, tempRec, rng
        );
      }
  );
}

const bool BVHAggregate::hit(
    const Ray& r,
    Float t_min,
    Float t_max,
    hit_record& rec,
    Sampler* sampler
) const {
  return traverseClosestBVH4(
      nodes4.get(),
      leaves4.get(),
      root4,
      r,
      t_min,
      t_max,
      rec,
      [&](int primIndex, Float local_t_min, Float local_t_max,
          hit_record& tempRec) {
        return primitives[primIndex]->hit(
            r, local_t_min, local_t_max, tempRec, sampler
        );
      }
  );
}

bool BVHAggregate::HitP(
    const Ray& r,
    Float t_min,
    Float t_max,
    random_gen& rng
) const {
  return traverseAnyPriorityBVH4(
      nodes4.get(),
      leaves4.get(),
      root4,
      r,
      t_min,
      t_max,
      [&](int primIndex, Float local_t_min, Float local_t_max) {
        return primitives[primIndex]->HitP(
            r, local_t_min, local_t_max, rng
        );
      }
  );
}

bool BVHAggregate::HitP(
    const Ray& r,
    Float t_min,
    Float t_max,
    Sampler* sampler
) const {
  return traverseAnyPriorityBVH4(
      nodes4.get(),
      leaves4.get(),
      root4,
      r,
      t_min,
      t_max,
      [&](int primIndex, Float local_t_min, Float local_t_max) {
        return primitives[primIndex]->HitP(
            r, local_t_min, local_t_max, sampler
        );
      }
  );
}

#ifdef NOT_CRAN
context("BVH frontier batching") {
  test_that("child fast paths preserve sequential insertion order and heap spills") {
    bool same_order = true, saw_heap = false;
    int serial = 0;
    // Exercise every child mask, ties, infinities, and batches crossing capacity.
    for (int pending = 0; pending <= 18; ++pending) {
      for (int mask = 1; mask < 16; ++mask) {
        for (int pattern = 0; pattern < 5; ++pattern) {
          BVH4PriorityFrontier<16> reference, optimized;
          for (int i = 0; i < pending; ++i) {
            BVHNodeEntry entry{serial++, float(i % 3)};
            reference.push(entry);
            optimized.push(entry);
          }
          IVec4 offsets;
          float distances[4];
          for (int i = 0; i < 4; ++i) {
            offsets[i] = i % 2 ? bvh4LeafReference(serial++) : serial++;
            distances[i] = pattern == 0 ? 1.f : pattern == 1 ? float(i - 2) :
              pattern == 2 ? float(3 - i) : pattern == 3 ?
              -std::numeric_limits<float>::infinity() : std::numeric_limits<float>::infinity();
            if ((mask >> i) & 1) reference.push({offsets[i], distances[i]});
          }
          BVHNodeEntry actual = optimized.next_child(offsets, distances, mask);
          BVHNodeEntry expected = reference.pop();
          same_order &= actual.nodeIndex == expected.nodeIndex && actual.tEnter == expected.tEnter;
          same_order &= optimized.using_heap() == reference.using_heap();
          saw_heap |= optimized.using_heap();
          while (!reference.empty() && !optimized.empty()) {
            actual = optimized.pop();
            expected = reference.pop();
            same_order &= actual.nodeIndex == expected.nodeIndex && actual.tEnter == expected.tEnter;
          }
          same_order &= reference.empty() && optimized.empty();
        }
      }
    }
    expect_true(same_order);
    expect_true(saw_heap);
  }
}

namespace {
class BVHLeafProbe : public hitable {
public:
  BVHLeafProbe(int id, bool coincident, int& target, std::vector<int>& visits)
      : id(id), x(coincident ? 0 : 2 * id), target(target), visits(visits) {}
  const bool hit(const Ray& r, Float lo, Float hi, hit_record& rec,
                 random_gen& rng) const override {
    rng.unif_rand();
    return intersect(r, lo, hi, rec);
  }
  const bool hit(const Ray& r, Float lo, Float hi, hit_record& rec,
                 Sampler* sampler) const override {
    sampler->Get1D();
    return intersect(r, lo, hi, rec);
  }
  bool bounding_box(Float, Float, aabb& bounds) const override {
    bounds = aabb(point3f(x, -1, -1), point3f(x + 1, 1, 1));
    return true;
  }
  std::string GetName() const override { return "BVH leaf test probe"; }
  size_t GetSize() override { return sizeof(*this); }
  void hitable_info_bounds(Float, Float) const override {}
private:
  bool intersect(const Ray& r, Float lo, Float hi, hit_record& rec) const {
    visits.push_back(id);
    const Float distance = (x - r.o[0]) / r.d[0];
    if (id != target || distance < lo || distance > hi) return false;
    rec.t = distance;
    rec.shape = this;
    return true;
  }
  int id;
  Float x;
  int& target;
  std::vector<int>& visits;
};
}

context("BVH compact leaves") {
  test_that("containment probes prune distant whole objects, including root and coincident leaves") {
    random_gen rng(91);
    bool correct = true;
    for (bool coincident : {false, true}) {
      for (int count : {0, 1, 2, 17, 300}) {
        for (int maxLeaf : {1, 4}) {
          std::vector<int> visits;
          std::vector<std::shared_ptr<hitable>> probes;
          int target = -1;
          for (int i = 0; i < count; ++i)
            probes.push_back(std::make_shared<BVHLeafProbe>(i, coincident, target, visits));
          BVHAggregate bvh(probes, 0, 1, maxLeaf, true);
          for (int index : {0, count / 2, count}) {
            const Float x = coincident ? .5f : 2 * index + .5f;
            for (Float y : {0.f, 1.f, 2.f}) {
              visits.clear();
              hit_record rec;
              bool hit = bvh.HitContainingObjects(
                  Ray(point3f(x, y, 0), vec3f(1, 0, 0)), 0, FLT_MAX, rec, rng);
              std::vector<int> expected;
              if (y <= 1) {
                if (coincident)
                  for (int i = 0; i < count; ++i) expected.push_back(i);
                else if (index < count) expected.push_back(index);
              }
              std::sort(visits.begin(), visits.end());
              correct &= !hit && visits == expected;
              if (!probes.empty()) {
                std::vector<const hitable*> excluded{probes.front().get()};
                visits.clear();
                hit = bvh.HitContainingObjects(
                    Ray(point3f(x, y, 0), vec3f(1, 0, 0)), 0, FLT_MAX, rec, rng, &excluded);
                expected.erase(std::remove(expected.begin(), expected.end(), 0), expected.end());
                std::sort(visits.begin(), visits.end());
                correct &= !hit && visits == expected;
              }
            }
          }
        }
      }
    }
    expect_true(correct);
  }

  test_that("leaf references preserve the unused sentinel and full index range") {
    expect_false(isBVH4Leaf(kBVH4Empty));
    expect_false(isBVH4Leaf(0));
    for (int index : {0, 1, 255, 65535, std::numeric_limits<int>::max() - 1}) {
      const int reference = bvh4LeafReference(index);
      expect_true(isBVH4Leaf(reference));
      expect_true(bvh4LeafIndex(reference) == index);
    }
  }

  test_that("empty, root-leaf, large coincident, and mixed trees retain their hits") {
    Transform identity;
    Ray ray(point3f(-2, 0, 0), vec3f(1, 0, 0));
    random_gen rng(91);
    RandomSampler sampler(rng);
    bool correct = true;
    for (bool coincident : {false, true}) {
      for (int count : {0, 1, 2, 3, 7, 257, 300}) {
        for (int maxLeaf : {1, 4}) {
          std::vector<int> visits, expected;
          std::vector<std::shared_ptr<hitable>> probes;
          int target = -1;
          for (int i = 0; i < count; ++i) {
            probes.push_back(std::make_shared<BVHLeafProbe>(i, coincident, target, visits));
            expected.push_back(i);
          }
          for (bool placed : {false, true}) {
            std::unique_ptr<BVHAggregate> bvh;
            if (placed)
              bvh = std::make_unique<BVHAggregate>(probes, 0, 1, maxLeaf, true,
                                                  &identity, &identity, false);
            else
              bvh = std::make_unique<BVHAggregate>(probes, 0, 1, maxLeaf, true);
            for (bool sampled : {false, true}) {
              hit_record rec;
              target = -1;
              visits.clear();
              bool hit = sampled ? bvh->hit(ray, 0, FLT_MAX, rec, &sampler) :
                                   bvh->hit(ray, 0, FLT_MAX, rec, rng);
              correct &= !hit && visits == expected;
              visits.clear();
              hit = sampled ? bvh->HitP(ray, 0, FLT_MAX, &sampler) :
                              bvh->HitP(ray, 0, FLT_MAX, rng);
              correct &= !hit && visits == expected;
              if (count == 0) continue;

              target = count - 1;
              visits.clear();
              hit = sampled ? bvh->hit(ray, 0, FLT_MAX, rec, &sampler) :
                              bvh->hit(ray, 0, FLT_MAX, rec, rng);
              correct &= hit && rec.shape == probes.back().get() && visits == expected;
              visits.clear();
              hit = sampled ? bvh->HitP(ray, 0, FLT_MAX, &sampler) :
                              bvh->HitP(ray, 0, FLT_MAX, rng);
              correct &= hit && visits == expected;
            }
          }
        }
      }
    }
    expect_true(correct);
  }
}
#endif
#endif

#ifndef RAYSIMD
bool BVHAggregate::HitContainingObjects(const Ray& r, Float t_min, Float t_max,
                                       hit_record& rec, random_gen& rng,
                                       const std::vector<const hitable*>* excluded) const {
  // Scalar builds use the same whole-object filter without a point BVH walk.
  bool found = false;
  for (const auto& object : primitives) {
    if (excluded && std::find(excluded->begin(), excluded->end(), object.get()) != excluded->end())
      continue;
    aabb bounds;
    bool outside = false;
    if (object->bounding_box(r.time(), r.time(), bounds)) {
      for (int axis = 0; axis < 3; ++axis) {
        const Float padding = 8 * std::numeric_limits<Float>::epsilon() *
                              std::max(Float(1), std::abs(r.o[axis]));
        outside |= r.o[axis] < bounds.min()[axis] - padding ||
                   r.o[axis] > bounds.max()[axis] + padding;
      }
    }
    hit_record candidate;
    if (!outside && object->hit(r, t_min, t_max, candidate, rng) &&
        (!found || candidate.OrderedDistance() <= rec.OrderedDistance())) {
      rec = candidate;
      t_max = rec.DistanceUpperBound();
      found = true;
    }
  }
  return found;
}
#endif

Float BVHAggregate::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) {
  if (primitives.empty()) {
    return 0;
  }
  Float weight = 1.0 / primitives.size();
  Float sum = 0;
  for (const auto& object : primitives) {
    sum += weight*object->pdf_value(o,v, rng, time);
  }
  return(sum);
}

Float BVHAggregate::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) {
  if (primitives.empty()) {
    return 0;
  }
  Float weight = 1.0 / primitives.size();
  Float sum = 0;
  for (const auto& object : primitives) {
    sum += weight*object->pdf_value(o,v, sampler, time);
  }
  return(sum);
}

vec3f BVHAggregate::random(const point3f& o, random_gen& rng, Float time) {
  if (primitives.empty()) {
    return vec3f(0, 0, 0);
  }
  int index = int(rng.unif_rand() * primitives.size() * 0.99999999);
  return(primitives[index]->random(o, rng, time));
}

vec3f BVHAggregate::random(const point3f& o, Sampler* sampler, Float time) {
  if (primitives.empty()) {
    return vec3f(0, 0, 0);
  }
  int index = int(sampler->Get1D() * primitives.size() * 0.99999999);
  return(primitives[index]->random(o, sampler, time));
}

BVHBuildNode4* BVHAggregate::ConvertBVH2ToBVH4(BVHBuildNode* node, int* totalNodes4,
                                               int* totalLeaves4) {
    if (node == nullptr) {
        return nullptr;
    }

    // Skip nodes with no primitives and no children
    if (node->nPrimitives == 0 && node->children[0] == nullptr && node->children[1] == nullptr) {
        return nullptr;
    }

    BVHBuildNode4* newNode = new BVHBuildNode4();
    newNode->bounds = node->bounds;
    newNode->splitAxis = node->splitAxis;

    if (node->nPrimitives > 0) {
        // Leaf node
        newNode->nPrimitives = node->nPrimitives;
        newNode->firstPrimOffset = node->firstPrimOffset;
        newNode->nChildren = 0;
        for (int i = 0; i < 4; ++i) {
            newNode->children[i] = nullptr;
        }
    } else {
        // Interior node
        BVHBuildNode* potentialChildren[4];
        newNode->firstPrimOffset = -1;
        newNode->nPrimitives = 0;
        newNode->nChildren = 0;
        int nChildren = 0;

        // Process the two children of the BVH2 node
        for (int i = 0; i < 2; ++i) {
            BVHBuildNode* child = node->children[i];
            if (child == nullptr) {
                continue;
            }
            // Skip nodes with no primitives and no children
            if (child->nPrimitives == 0 && child->children[0] == nullptr && child->children[1] == nullptr) {
                continue;
            }
            if (child->nPrimitives > 0 || (child->children[0] == nullptr && child->children[1] == nullptr)) {
                // If the child is a leaf node or has no further children
                potentialChildren[nChildren++] = child;
            } else {
                // The child is an interior node, add its two children
                if (child->children[0] != nullptr) {
                    potentialChildren[nChildren++] = child->children[0];
                }
                if (child->children[1] != nullptr) {
                    potentialChildren[nChildren++] = child->children[1];
                }
            }
        }

        if (nChildren == 0) {
            // No valid children, delete the node and return nullptr
            delete newNode;
            return nullptr;
        }
        // Initialize the BVH4 node
        for (int i = 0; i < nChildren; ++i) {
            // Recursively convert the child nodes
            BVHBuildNode4* childNode = ConvertBVH2ToBVH4(potentialChildren[i], totalNodes4, totalLeaves4);
            if (childNode != nullptr) {
                newNode->children[newNode->nChildren++] = childNode;
                // Update the bounds
                newNode->bounds = surrounding_box(newNode->bounds, childNode->bounds);
            }
        }
        // Set any remaining child pointers to nullptr
        for (int i = newNode->nChildren; i < 4; ++i) {
            newNode->children[i] = nullptr;
        }

        if (newNode->nChildren == 0) {
            // No valid children, delete the node and return nullptr
            delete newNode;
            return nullptr;
        }
    }

    if (newNode->nPrimitives > 0) ++(*totalLeaves4);
    else ++(*totalNodes4);
    return newNode;
}


void BVHAggregate::buildBVH4(BVHBuildNode* root) {
    std::unique_ptr<BVHBuildNode> binaryRoot(root);
    std::unique_ptr<BVHBuildNode4> wideRoot(
        ConvertBVH2ToBVH4(root, &totalNodes4, &totalLeaves4));
    binaryRoot.reset();
    ASSERT(wideRoot && totalLeaves4 > 0 && "BVH4 conversion produced an empty tree");
    if (!wideRoot) return;

    // Allocate interiors and leaves separately. A one-leaf tree has no interior
    // allocation; its tagged root reference enters the ordinary leaf path.
    nodes4 = allocateBVH4NodeStorage(totalNodes4);
    if (totalNodes4 > 0) std::uninitialized_default_construct_n(nodes4.get(), totalNodes4);
    leaves4.reset(new LinearBVHLeaf4[totalLeaves4]);
    int offset = 0, leafOffset = 0;
    root4 = flattenBVH4(wideRoot.get(), &offset, &leafOffset);
    ASSERT(offset == totalNodes4 && leafOffset == totalLeaves4 &&
           "BVH4 flatten count mismatch");
#ifndef NDEBUG
    validateBVH4();
#endif
}

int BVHAggregate::flattenBVH4(BVHBuildNode4* node, int* offset, int* leafOffset) {
    if (!node) throw std::runtime_error("flattenBVH4 called with nullptr node");

    if (node->nChildren == 0) {
        ASSERT(node->nPrimitives > 0);
        const int index = (*leafOffset)++;
        leaves4[index] = {node->firstPrimOffset, node->nPrimitives};
        return bvh4LeafReference(index);
    }

    const int nodeOffset = (*offset)++;
    LinearBVHNode4& linearNode = nodes4[nodeOffset];
    linearNode.nChildren = node->nChildren;
    aabb childBoxes[4];
    for (int i = 0; i < node->nChildren; ++i) {
        linearNode.childOffsets[i] = flattenBVH4(node->children[i], offset, leafOffset);
        childBoxes[i] = node->children[i]->bounds;
    }
    for (int i = node->nChildren; i < 4; ++i) {
        linearNode.childOffsets[i] = kBVH4Empty;
        childBoxes[i] = aabb();
    }
    // Preserve each child's box and lane order exactly, including leaf boxes.
    linearNode.bbox4 = BBox4(childBoxes[0], childBoxes[1], childBoxes[2], childBoxes[3]);
    return nodeOffset;
}

void BVHAggregate::validateBVH4() const {
    if (root4 == kBVH4Empty) throw std::runtime_error("BVH4 tree is empty.");
    std::vector<int> pending{root4};
    std::vector<bool> visitedNodes(totalNodes4, false), visitedLeaves(totalLeaves4, false);
    size_t primitiveCount = 0, nodeCount = 0, leafCount = 0;
    while (!pending.empty()) {
        const int reference = pending.back();
        pending.pop_back();
        if (isBVH4Leaf(reference)) {
            const int index = bvh4LeafIndex(reference);
            if (index >= totalLeaves4 || visitedLeaves[index])
                throw std::runtime_error("Invalid or repeated BVH4 leaf reference.");
            visitedLeaves[index] = true;
            ++leafCount;
            const LinearBVHLeaf4& leaf = leaves4[index];
            if (leaf.primitivesOffset < 0 || leaf.nPrimitives <= 0 ||
                size_t(leaf.primitivesOffset) + size_t(leaf.nPrimitives) > primitives.size())
                throw std::runtime_error("Invalid BVH4 leaf primitive range.");
            primitiveCount += leaf.nPrimitives;
        } else {
            if (reference < 0 || reference >= totalNodes4 || visitedNodes[reference])
                throw std::runtime_error("Invalid or repeated BVH4 interior reference.");
            visitedNodes[reference] = true;
            ++nodeCount;
            const LinearBVHNode4& node = nodes4[reference];
            if (node.nChildren <= 0 || node.nChildren > 4)
                throw std::runtime_error("Invalid BVH4 child count.");
            for (int i = 0; i < node.nChildren; ++i) {
                if (node.childOffsets[i] == kBVH4Empty)
                    throw std::runtime_error("Empty BVH4 child inside the valid lanes.");
                pending.push_back(node.childOffsets[i]);
            }
            for (int i = node.nChildren; i < 4; ++i)
                if (node.childOffsets[i] != kBVH4Empty)
                    throw std::runtime_error("Nonempty BVH4 child outside the valid lanes.");
        }
    }
    if (nodeCount != size_t(totalNodes4) || leafCount != size_t(totalLeaves4) ||
        primitiveCount != primitives.size())
        throw std::runtime_error("BVH4 storage or primitive count mismatch.");
}
