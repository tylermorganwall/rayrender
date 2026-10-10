#ifndef RAYRENDER_HLBVH_H
#define RAYRENDER_HLBVH_H

#include <cstdint>
#include <memory>
#include <span>
#include <vector>

enum class BVHBuildMethod { SAH = 1, Equal = 2, HLBVH = 3, Metal = 4 };
struct BVHBuildOptions {
  BVHBuildMethod method = BVHBuildMethod::SAH;
  unsigned threads = 1;
  // Internal opt-in: every primitive's bounding_box and ShadowType must be
  // read-only, thread-safe, and callable without R APIs or lazy initialization.
  // Audited triangle-mesh builders enable this after loading their geometry.
  bool parallel_primitive_queries = false;
};

// Plain, explicitly sized host/device records. Bounds are conservatively rounded
// to float when Float is double; primitive intersection retains its native type.
struct HLBVHBounds {
  float lo[4], hi[4];
};
struct MortonPrimitive {
  uint32_t code, index;
};
struct HLBVHTreelet {
  uint32_t first, count;
};
struct HLBVHNode {
  HLBVHBounds bounds;
  uint32_t left, right, first, count, axis, padding[3];
};
static_assert(sizeof(HLBVHBounds) == 32 && sizeof(HLBVHNode) == 64);

struct HLBVHTree {
  std::vector<MortonPrimitive> order;
  std::vector<HLBVHTreelet> treelets;
  // CPU arenas own new[] storage; Metal arenas retain their shared MTLBuffer
  // through a custom deleter until upper-SAH stitching and packing finish.
  std::shared_ptr<HLBVHNode[]> nodes;
  uint32_t root = 0;
};

// PBRT v4's HLBVH: 30-bit Morton codes, radix sort, 12-bit treelet grouping,
// independent lower trees, then a 12-bucket SAH hierarchy over their roots.
HLBVHTree BuildHLBVH(std::span<const HLBVHBounds> bounds, unsigned max_leaf, BVHBuildOptions options);
bool MetalBVHAvailable();
// The Metal backend produces sorted order and complete lower treelets (including
// bounds). CPU SAH stitching and traversal-layout packing are shared by both.
void BuildMetalTreelets(std::span<const HLBVHBounds> bounds, const HLBVHBounds &centroids, unsigned max_leaf,
                        HLBVHTree &result);

#endif
