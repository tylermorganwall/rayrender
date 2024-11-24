#include "bvh.h"
#include "assert.h"
#include "mathinline.h"
#include <cmath>
#include <limits>
#include "raylog.h"
#include "aabb.h"
#include "simd.h"

#include <queue>
#include <vector>
#include <functional>

// Define a structure to hold node index and tEnter time
struct alignas(16) BVHNodeEntry {
    int nodeIndex;
    float tEnter;

    // Comparator for the priority queue (min-heap)
    bool operator<(const BVHNodeEntry& other) const {
        return tEnter > other.tEnter; // Inverted comparison for min-heap
    }
};

inline int floatToIntBits(float f) {
    int i;
    memcpy(&i, &f, sizeof(float));
    return i;
}

//Insertion sort is faster than a heap here
template <size_t MaxSize>
class StaticPriorityQueue {
public:
    StaticPriorityQueue() : size_(0) {}

    // Insert a new entry into the priority queue using insertion sort
    void push(const BVHNodeEntry& entry) {
        ASSERT(size_ < MaxSize && "Priority queue overflow");

        // Find the correct position to insert the new entry
        int i = size_ - 1;
        while (i >= 0 && data_[i].tEnter < entry.tEnter) {
            // Shift elements to make room for the new entry
            data_[i + 1] = data_[i];
            i--;
        }
        // Insert the new entry at the found position
        data_[i + 1] = entry;
        size_++;
    }

    // Remove and return the element with the smallest tEnter
    BVHNodeEntry pop() {
        ASSERT(size_ > 0 && "Priority queue underflow");
        // The smallest element is at the end of the array
        return data_[--size_];
    }

    // Return a const reference to the element with the smallest tEnter
    const BVHNodeEntry& top() const {
        ASSERT(size_ > 0 && "Priority queue is empty");
        // The smallest element is at the end
        return data_[size_ - 1];
    }

    bool empty() const {
        return size_ == 0;
    }

    size_t size() const {
        return size_;
    }

    void clear() {
        size_ = 0;
    }

private:
    BVHNodeEntry data_[MaxSize];
    size_t size_;
};

BVHAggregate::BVHAggregate(std::vector<std::shared_ptr<hitable> > prims,
        float t_min, float t_max, 
        int maxPrimsInNode, bool sah, 
        std::shared_ptr<Transform> ObjectToWorld, 
        std::shared_ptr<Transform> WorldToObject, 
        bool reverseOrientation) : 
        hitable(ObjectToWorld, WorldToObject, reverseOrientation), 
        maxPrimsInNode(std::min(255, maxPrimsInNode)),
        primitives(prims)
        { 
    SCOPED_CONTEXT("Initialization");
    SCOPED_TIMER_COUNTER("BVH Build");
    std::vector<BVHPrimitive> bvhPrimitives(primitives.size());
    for (size_t i = 0; i < primitives.size(); ++i) {
        aabb temp_box;
        primitives[i]->bounding_box(t_min, t_max, temp_box);
        bvhPrimitives[i] = BVHPrimitive(i, temp_box);
        scene_bounds = surrounding_box(scene_bounds, temp_box);
    }
    std::vector<std::shared_ptr<hitable> > orderedPrims(primitives.size());
    BVHBuildNode* root;
    std::atomic<int> totalNodes{0};
    std::atomic<int> orderedPrimsOffset{0};
    root = buildRecursive(std::span<BVHPrimitive>(bvhPrimitives),
                          &totalNodes, 
                          &orderedPrimsOffset, 
                          orderedPrims);
    primitives.swap(orderedPrims);
    bvhPrimitives.resize(0);
    bvhPrimitives.shrink_to_fit();

#ifndef RAYSIMD
    nodes.reset(new LinearBVHNode[totalNodes]);
    n_nodes = totalNodes;
    int offset = 0;
    flattenBVH(root, &offset);
#else
    totalNodes4 = 0;
    BVHBuildNode4* rootBVH4 = ConvertBVH2ToBVH4(root, &totalNodes4);
    nodes4.reset(new LinearBVHNode4[totalNodes4]);
    int offset4 = 0;
    flattenBVH4(rootBVH4, &offset4);
    // validateBVH4();
#endif
}


BVHAggregate::BVHAggregate(std::vector<std::shared_ptr<hitable> > prims,
                           float t_min, float t_max, 
                           int maxPrimsInNode, bool sah) :
                                maxPrimsInNode(std::min(255, maxPrimsInNode)),
                                primitives(prims) { 
    SCOPED_CONTEXT("Initialization");
    SCOPED_TIMER_COUNTER("BVH Build");
    std::vector<BVHPrimitive> bvhPrimitives(primitives.size());
    for (size_t i = 0; i < primitives.size(); ++i) {
        aabb temp_box;
        primitives[i]->bounding_box(t_min, t_max, temp_box);
        bvhPrimitives[i] = BVHPrimitive(i, temp_box);
        scene_bounds = surrounding_box(scene_bounds, temp_box);
    }

    std::vector<std::shared_ptr<hitable> > orderedPrims(primitives.size());
    BVHBuildNode* root;
    std::atomic<int> totalNodes{0};
    std::atomic<int> orderedPrimsOffset{0};
    root = buildRecursive(std::span<BVHPrimitive>(bvhPrimitives),
                          &totalNodes, 
                          &orderedPrimsOffset, 
                          orderedPrims);
    primitives.swap(orderedPrims);
    bvhPrimitives.resize(0);
    bvhPrimitives.shrink_to_fit();
#ifndef RAYSIMD
    nodes.reset(new LinearBVHNode[totalNodes]);
    n_nodes = totalNodes;
    int offset = 0;
    flattenBVH(root, &offset);
#else
    totalNodes4 = 0;
    BVHBuildNode4* rootBVH4 = ConvertBVH2ToBVH4(root, &totalNodes4);
    nodes4.reset(new LinearBVHNode4[totalNodes4]);
    int offset4 = 0;
    flattenBVH4(rootBVH4, &offset4);
    // validateBVH4();
#endif
}

bool BVHAggregate::bounding_box(Float t0, Float t1, aabb& box) const {
    box = scene_bounds;
    return(true);
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
    int nodeIndex = (*totalNodes)++;
    // Compute bounds of all primitives in BVH node
    aabb bounds;
    for (const auto &prim : bvhPrimitives) {
        bounds = surrounding_box(bounds, prim.bounds);
    }

    if ((bounds.surface_area() == 0 || bvhPrimitives.size() <= maxPrimsInNode) && !isRoot) {
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
            if (bvhPrimitives.size() <= 2) {
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
                if (bvhPrimitives.size() > maxPrimsInNode || minCost < leafCost) {
                    auto midIter = std::partition(
                        bvhPrimitives.begin(), bvhPrimitives.end(),
                        [=](const BVHPrimitive &bp) {
                            int b = nBuckets * centroidBounds.offset(bp.Centroid())[dim];
                            if (b == nBuckets) {
                                b = nBuckets - 1;
                            }
                            return b <= minCostSplitBucket;
                        });
                    mid = midIter - bvhPrimitives.begin();
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
        ASSERT(!node->children[0] && !node->children[1]);
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
const bool BVHAggregate::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const {
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
                    if (prim_hrec) {
                        any_hit = true;
                        rec = hrec_temp;
                        t_max = rec.t;
                    }
                }
                if (toVisitOffset == 0) {
                    break;
                }
                currentNodeIndex = nodesToVisit[--toVisitOffset];
            } else {
                // Put far BVH node on _nodesToVisit_ stack, advance to near node
                if (r.sign[node->axis]) {
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

const bool BVHAggregate::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const {
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
                    if (prim_hrec) {
                        any_hit = true;
                        rec = hrec_temp;
                        t_max = rec.t;
                    }
                }
                if (toVisitOffset == 0) {
                    break;
                }
                currentNodeIndex = nodesToVisit[--toVisitOffset];
            } else {
                // Put far BVH node on _nodesToVisit_ stack, advance to near node
                if (r.sign[node->axis]) {
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
#else

#define TOTAL_NODES_STATIC 4096
constexpr const int triggerLargeBVH = TOTAL_NODES_STATIC - SIMD_WIDTH;
#define RAY_SUPPORT_LARGE_BVH

const bool BVHAggregate::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const {
    if (!nodes4) {
        return false;
    }
    StaticPriorityQueue<TOTAL_NODES_STATIC> nodesToVisit;
    bool large_bvh = false;
    std::priority_queue<BVHNodeEntry> nodesToVisitLarge;

    nodesToVisit.push({0, -INFINITY});

    bool any_hit = false;
    
    while (!nodesToVisit.empty()) {
        if(nodesToVisit.size() > triggerLargeBVH) {
            large_bvh = true;
        }
        // Pop the node with the smallest tEnter
        BVHNodeEntry entry = nodesToVisit.top();
        nodesToVisit.pop();
        int currentNodeIndex = entry.nodeIndex;
        const LinearBVHNode4* node = &nodes4[currentNodeIndex];

        if (node->nPrimitives > 0) {
            // __builtin_prefetch(&primitives[node->primitivesOffset]);

            // Leaf node: test ray against primitives
            for (int i = 0; i < node->nPrimitives; ++i) {
                hit_record tempRec;
                bool hitPrimitive = primitives[node->primitivesOffset + i]->hit(r, t_min, t_max, tempRec, rng);
                if (hitPrimitive) {
                    any_hit = true;
                    rec = tempRec;
                    t_max = tempRec.t; // Update tMax to the closest hit
                }
            }
            if (any_hit) {
              if(nodesToVisit.empty()) {
                return true;
              }
              if(nodesToVisit.top().tEnter > t_max) {
                return true;
              }
            }
        } else {
            // Internal node
            const BBox4& bbox4 = node->bbox4;
            // Perform SIMD ray-box intersection

            IVec4 hits;
            FVec4 tEnters;
            
            rayBBoxIntersect4(r, bbox4, t_min, t_max, hits, tEnters);//, tExits);

            float tEntersArray[4];
            simd_extract_fvec4(tEnters, tEntersArray);
            
            const IVec4 valid_hit = simd_and(hits, simd_not_equals_minus_one(node->childOffsets));
            // int hitmask = simd_extract_hitmask(valid_hit);
            for (int i = 0; i < 4; ++i) {
                // __builtin_prefetch(&tEntersArray[i+1]);

                // const bool valid = (hitmask >> i) & 1;
                // const bool valid = valid_hit[i];

                if (valid_hit[i]) {
                    nodesToVisit.push({node->childOffsets[i], tEntersArray[i]});
                }
            }
        }
    }
    //If we exceed the static priority size, copy all elements to a normal priority queue and continue
    if(large_bvh) {
         while (!nodesToVisit.empty()) {
            nodesToVisitLarge.push(nodesToVisit.top());
            nodesToVisit.pop();
         }
        while (!nodesToVisitLarge.empty()) {
            // Pop the node with the smallest tEnter
            BVHNodeEntry entry = nodesToVisitLarge.top();
            nodesToVisitLarge.pop();
            int currentNodeIndex = entry.nodeIndex;
            const LinearBVHNode4* node = &nodes4[currentNodeIndex];

            if (node->nPrimitives > 0) {
                // Leaf node: test ray against primitives
                for (int i = 0; i < node->nPrimitives; ++i) {
                    hit_record tempRec;
                    bool hitPrimitive = primitives[node->primitivesOffset + i]->hit(r, t_min, t_max, tempRec, rng);
                    if (hitPrimitive) {
                        any_hit = true;
                        rec = tempRec;
                        t_max = tempRec.t; // Update tMax to the closest hit
                    }
                }
                if (any_hit && nodesToVisitLarge.top().tEnter > t_max) {
                    return true;
                }
            } else {
                // Internal node
                const BBox4& bbox4 = node->bbox4;
                // Perform SIMD ray-box intersection

                IVec4 hits;
                FVec4 tEnters;
                
                rayBBoxIntersect4(r, bbox4, t_min, t_max, hits, tEnters);//, tExits);

                const IVec4 valid_hit = simd_and(hits, simd_not_equals_minus_one(node->childOffsets));
                int hitmask = simd_extract_hitmask(valid_hit);

                float tEntersArray[4];
                simd_extract_fvec4(tEnters, tEntersArray);
                
                for (int i = 0; i < 4; ++i) {
                    const bool valid = (hitmask >> i) & 1;
                    if (valid) {
                        nodesToVisitLarge.push({node->childOffsets[i], tEntersArray[i]});
                    }
                }
            }
        }
    }

    return any_hit;
}

const bool BVHAggregate::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const {
    if (!nodes4) {
        return false;
    }
    StaticPriorityQueue<64> nodesToVisit;
    nodesToVisit.push({0, -INFINITY});

    bool any_hit = false;
    
    while (!nodesToVisit.empty()) {
        // Pop the node with the smallest tEnter
        BVHNodeEntry entry = nodesToVisit.top();
        nodesToVisit.pop();
        int currentNodeIndex = entry.nodeIndex;
        const LinearBVHNode4* node = &nodes4[currentNodeIndex];

        if (node->nPrimitives > 0) {
            // Leaf node: test ray against primitives
            for (int i = 0; i < node->nPrimitives; ++i) {
                hit_record tempRec;
                bool hitPrimitive = primitives[node->primitivesOffset + i]->hit(r, t_min, t_max, tempRec, sampler);
                if (hitPrimitive) {
                    any_hit = true;
                    rec = tempRec;
                    t_max = tempRec.t; // Update tMax to the closest hit
                }
            }
            if (any_hit && nodesToVisit.top().tEnter > t_max) {
                return true;
            }
        } else {
            // Internal node
            const BBox4& bbox4 = node->bbox4;
            // Perform SIMD ray-box intersection

            IVec4 hits;
            FVec4 tEnters;
            
            rayBBoxIntersect4(r, bbox4, t_min, t_max, hits, tEnters);//, tExits);

            const IVec4 valid_hit = simd_and(hits, simd_not_equals_minus_one(node->childOffsets));
            int hitmask = simd_extract_hitmask(valid_hit);

            float tEntersArray[4];
            simd_extract_fvec4(tEnters, tEntersArray);

            for (int i = 0; i < 4; ++i) {
                const bool valid = (hitmask >> i) & 1;
                if (valid) {
                    nodesToVisit.push({node->childOffsets[i], tEntersArray[i]});
                }
            }
        }
    }

    return any_hit;
}
#endif

Float BVHAggregate::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) {
  Float weight = 1.0 / primitives.size();
  Float sum = 0;
  for (const auto& object : primitives) {
    sum += weight*object->pdf_value(o,v, rng, time);
  }
  return(sum);
}

Float BVHAggregate::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) {
  Float weight = 1.0 / primitives.size();
  Float sum = 0;
  for (const auto& object : primitives) {
    sum += weight*object->pdf_value(o,v, sampler, time);
  }
  return(sum);
}

vec3f BVHAggregate::random(const point3f& o, random_gen& rng, Float time) {
  int index = int(rng.unif_rand() * primitives.size() * 0.99999999);
  return(primitives[index]->random(o, rng, time));
}

vec3f BVHAggregate::random(const point3f& o, Sampler* sampler, Float time) {
  int index = int(sampler->Get1D() * primitives.size() * 0.99999999);
  return(primitives[index]->random(o, sampler, time));
}

BVHBuildNode4* BVHAggregate::ConvertBVH2ToBVH4(BVHBuildNode* node, int* totalNodes4) {
    if (node == nullptr) {
        return nullptr;
    }

    // Skip nodes with no primitives and no children
    if (node->nPrimitives == 0 && node->children[0] == nullptr && node->children[1] == nullptr) {
        return nullptr;
    }

    // Increment the total node count for BVH4
    (*totalNodes4)++;

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
            BVHBuildNode4* childNode = ConvertBVH2ToBVH4(potentialChildren[i], totalNodes4);
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

    return newNode;
}


int BVHAggregate::flattenBVH4(BVHBuildNode4* node, int* offset) {
    if (node == nullptr) {
        throw std::runtime_error("flattenBVH4 called with nullptr node");
    }

    LinearBVHNode4* linearNode = &nodes4[*offset];
    int nodeOffset = (*offset)++;

    // linearNode->bounds = node->bounds;
    // linearNode->axis = node->splitAxis;
    linearNode->nChildren = node->nChildren;

    if (node->nChildren == 0) {
        // Leaf node
        linearNode->nPrimitives = node->nPrimitives;
        linearNode->primitivesOffset = node->firstPrimOffset;
        for (int i = 0; i < 4; ++i) {
            linearNode->childOffsets[i] = -1;
        }
        // Leaf nodes don't need bbox4; initialize to zero or leave as is
    } else {
        // Interior node
        linearNode->nPrimitives = 0;
        linearNode->primitivesOffset = -1;

        // Prepare arrays for child data
        aabb childBoxes[4];
        for (int i = 0; i < node->nChildren; ++i) {
            linearNode->childOffsets[i] = flattenBVH4(node->children[i], offset);
            childBoxes[i] = node->children[i]->bounds;
        }
        for (int i = node->nChildren; i < 4; ++i) {
            linearNode->childOffsets[i] = -1;
            childBoxes[i] = aabb(); // Empty bounding box
        }

        // Pack the child bounding boxes into bbox4
        linearNode->bbox4 = BBox4(childBoxes[0], childBoxes[1], childBoxes[2], childBoxes[3]);
    }
    return nodeOffset;
}

void BVHAggregate::validateBVH4() const {
    if (!nodes4) {
        throw std::runtime_error("BVH4 tree is empty.");
    }

    // Use a stack for traversal
    struct NodeInfo {
        int nodeIndex;
        int depth;
    };

    std::vector<NodeInfo> nodesToVisit;
    nodesToVisit.push_back({0, 0}); // Start with root node at index 0

    // int totalNodesVisited = 0;

    while (!nodesToVisit.empty()) {
        NodeInfo current = nodesToVisit.back();
        nodesToVisit.pop_back();

        if (current.nodeIndex < 0 || current.nodeIndex >= totalNodes4) {
            throw std::runtime_error("Invalid node index during traversal: " + std::to_string(current.nodeIndex));
        }

        const LinearBVHNode4* node = &nodes4[current.nodeIndex];
        // totalNodesVisited++;

        bool isLeaf = node->nPrimitives > 0;
        
        if (isLeaf) {
            // Leaf node validation
            if (node->nChildren != 0) {
                throw std::runtime_error("Invalid leaf node at index " + std::to_string(current.nodeIndex) + ": nChildren != 0");
            }
            if (node->primitivesOffset < 0) {
                throw std::runtime_error("Invalid leaf node at index " + std::to_string(current.nodeIndex) + ": primitivesOffset < 0");
            }
            if (node->nPrimitives <= 0) {
                throw std::runtime_error("Invalid leaf node at index " + std::to_string(current.nodeIndex) + ": nPrimitives <= 0");
            }
            if (node->primitivesOffset + node->nPrimitives > primitives.size()) {
                throw std::runtime_error("Invalid leaf node at index " + std::to_string(current.nodeIndex) + ": primitivesOffset + nPrimitives exceeds primitives.size()");
            }
        } else {
            // Interior node validation
            if (node->nChildren <= 0 || node->nChildren > 4) {
                throw std::runtime_error("Invalid interior node at index " + std::to_string(current.nodeIndex) + ": nChildren = " + std::to_string(static_cast<int>(node->nChildren)));
            }
            if (node->primitivesOffset != -1) {
                throw std::runtime_error("Invalid interior node at index " + std::to_string(current.nodeIndex) + ": primitivesOffset != -1");
            }
            if (node->nPrimitives != 0) {
                throw std::runtime_error("Invalid interior node at index " + std::to_string(current.nodeIndex) + ": nPrimitives != 0");
            }
            // Validate child offsets and push them onto the stack
            for (int i = 0; i < node->nChildren; ++i) {
                int childOffset = node->childOffsets[i];
                if (childOffset < 0 || childOffset >= totalNodes4) {
                    throw std::runtime_error("Invalid child offset at node " + std::to_string(current.nodeIndex) + ", child " + std::to_string(i) + ": " + std::to_string(childOffset));
                } else {
                    nodesToVisit.push_back({childOffset, current.depth + 1});
                }
            }
            // Ensure unused childOffsets are set to -1
            for (int i = node->nChildren; i < 4; ++i) {
                if (node->childOffsets[i] != -1) {
                    throw std::runtime_error("Invalid interior node at index " + std::to_string(current.nodeIndex) + ": unused childOffset[" + std::to_string(i) + "] != -1");
                }
            }
        }
    }
}