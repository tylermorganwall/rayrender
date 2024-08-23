#include "flat_bvh.h"

bool FlatBVH::bounding_box(Float t0, Float t1, aabb& b) const {
  root_node->bounding_box(t0,t1,b);
  return(true);
}

void FlatBVH::build(const std::shared_ptr<bvh_node>& root) {
  root_node = root;  // Store the root node
  // Count the total number of nodes
  std::pair<size_t, size_t> node_counts = root->CountNodeLeaf();
  size_t total_nodes = node_counts.first + node_counts.second;
  
  // Resize the nodes vector
  nodes.resize(total_nodes);
  
  // Flatten the tree
  int offset = 0;
  flattenBVHNode(root.get(), offset);
}

FlatBVH::FlatBVH(const std::shared_ptr<bvh_node>& root) {
  root_node = root;  // Store the root node
  
  // Count the total number of nodes
  std::pair<size_t, size_t> node_counts = root->CountNodeLeaf();
  size_t total_nodes = node_counts.first + node_counts.second;
  
  // Resize the nodes vector
  nodes.resize(total_nodes);
  
  // Flatten the tree
  int offset = 0;
  flattenBVHNode(root.get(), offset);
}


void FlatBVH::flattenBVHNode(const bvh_node* node, int& offset) {
  FlatBVHNode& flat_node = nodes[offset++];
  
  // Simply copy the AABB, no SIMD conversion yet
  flat_node.bounds = node->box;
  flat_node.axis = MaxDimension(node->box.Diag());
  
  if (node->left == node->right) {
    flat_node.primitive_count = 1;
    flat_node.primitive_offset = offset;
  } else {
    flat_node.primitive_count = 0;
    flattenBVHNode(dynamic_cast<const bvh_node*>(node->left.get()), offset);
    flat_node.second_child_offset = offset;
    flattenBVHNode(dynamic_cast<const bvh_node*>(node->right.get()), offset);
  }
}


void FlatBVH::convertToSIMD() {
  const size_t simd_width = SIMD_WIDTH;  // 4 for SSE, 8 for AVX
  size_t simd_node_count = (nodes.size() + simd_width - 1) / simd_width;
  // simd_nodes.resize(simd_node_count);
  simd_nodes.clear();
  for (size_t i = 0; i < simd_node_count; ++i) {
    std::array<aabb, simd_width> aabbs;
    for (size_t j = 0; j < simd_width; ++j) {
      size_t node_index = i * simd_width + j;
      if (node_index < nodes.size()) {
        aabbs[j] = nodes[node_index].bounds;
      } else {
        // Pad with empty AABBs if necessary
        aabbs[j] = aabb();
      }
    }
    // simd_nodes[i] = SimdAABB(aabbs);
    simd_nodes.push_back(SimdAABB(aabbs));
  }
}


const bool FlatBVH::hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng) const {
  int stack[64];
  int stack_ptr = 0;
  int node_index = 0;
  bool hit_anything = false;
  
  while (true) {
    const FlatBVHNode& node = nodes[node_index];
    if (node.bounds.hit(r, tmin, tmax, rng)) {
      if (node.primitive_count > 0) {
        // Leaf node, test primitive
        const bvh_node* original_node = root_node.get();
        for (int i = 0; i < node_index; ++i) {
          if (nodes[i].primitive_count == 0) {
            original_node = (r.direction()[nodes[i].axis] < 0) 
            ? dynamic_cast<const bvh_node*>(original_node->right.get())
              : dynamic_cast<const bvh_node*>(original_node->left.get());
          }
        }
        if (original_node->hit(r, tmin, tmax, rec, rng)) {
          hit_anything = true;
          tmax = rec.t;
        }
      } else {
        // Interior node, push far child and continue with near child
        int far_child = node_index + 1;
        int near_child = node.second_child_offset;
        
        if (r.direction()[node.axis] < 0.0f) {
          std::swap(near_child, far_child);
        }
        
        stack[stack_ptr++] = far_child;
        node_index = near_child;
        continue;
      }
    }
    
    if (stack_ptr == 0) break;
    node_index = stack[--stack_ptr];
  }
  
  return hit_anything;
}

const bool FlatBVH::hit(const ray& r, Float tmin, Float tmax, hit_record& rec, Sampler* sampler) const {
  int stack[64];
  int stack_ptr = 0;
  int node_index = 0;
  bool hit_anything = false;
  
  while (true) {
    const FlatBVHNode& node = nodes[node_index];
    if (node.bounds.hit(r, tmin, tmax, sampler)) {
      if (node.primitive_count > 0) {
        // Leaf node, test primitive
        const bvh_node* original_node = root_node.get();
        for (int i = 0; i < node_index; ++i) {
          if (nodes[i].primitive_count == 0) {
            original_node = (r.direction()[nodes[i].axis] < 0) 
            ? dynamic_cast<const bvh_node*>(original_node->right.get())
              : dynamic_cast<const bvh_node*>(original_node->left.get());
          }
        }
        if (original_node->hit(r, tmin, tmax, rec, sampler)) {
          hit_anything = true;
          tmax = rec.t;
        }
      } else {
        // Interior node, push far child and continue with near child
        int far_child = node_index + 1;
        int near_child = node.second_child_offset;
        
        if (r.direction()[node.axis] < 0.0f) {
          std::swap(near_child, far_child);
        }
        
        stack[stack_ptr++] = far_child;
        node_index = near_child;
        continue;
      }
    }
    
    if (stack_ptr == 0) break;
    node_index = stack[--stack_ptr];
  }
  
  return hit_anything;
}

size_t FlatBVH::GetSize()  {
  return(root_node->GetSize());
}