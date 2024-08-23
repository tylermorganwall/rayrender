#ifndef FLAT_BVH
#define FLAT_BVH

#include <vector>
#include <memory>
#include "bvh_node.h"
#include "ray.h"
#include "simd.h"
#include "aabb.h"

struct FlatBVHNode {
  aabb bounds;  // Regular AABB, not SIMD
  union {
    int32_t primitive_offset;    // Leaf
    int32_t second_child_offset; // Interior
  };
  uint16_t primitive_count;
  uint8_t axis;
  uint8_t pad[1];  // Ensure 32-byte alignment
};

class FlatBVH : public hitable {
public:
  FlatBVH() = default;
  void build(const std::shared_ptr<bvh_node>& root);
  virtual const bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng) const;
  std::string GetName() const {
    return(std::string("FlatBVH"));
  }
  void convertToSIMD();
  // size_t GetSize();
private:
  std::vector<FlatBVHNode> nodes;
  std::shared_ptr<bvh_node> root_node;
  std::vector<SimdAABB> simd_nodes;
  
  void flattenBVHNode(const bvh_node* node, int& offset);
};

#endif
