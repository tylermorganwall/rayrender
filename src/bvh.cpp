#include "bvh.h"

BVHAggregate::BVHAggregate(std::vector<std::shared_ptr<hitable> > prims,
        float t_min, float t_max, 
        int maxPrimsInNode, bool sah, 
        std::shared_ptr<Transform> ObjectToWorld, 
        std::shared_ptr<Transform> WorldToObject, 
        bool reverseOrientation) : 
        hitable(ObjectToWorld, WorldToObject, reverseOrientation), 
        maxPrimsInNode(std::min(255, maxPrimsInNode)),
        primitives(prims)//, splitMethod(splitMethod) 
        { 
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
    nodes.reset(new LinearBVHNode[totalNodes]);
    n_nodes = totalNodes;
    int offset = 0;
    flattenBVH(root, &offset);
    // transformToSimdFormat();
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
    ++*totalNodes;
    // Compute bounds of all primitives in BVH node
    aabb bounds;
    for (const auto &prim : bvhPrimitives) {
        bounds = surrounding_box(bounds, prim.bounds);
    }

    if (bounds.surface_area() == 0 || bvhPrimitives.size() == 1) {
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
            // Partition primitives based on _splitMethod_
            // switch (splitMethod) {
            // case SplitMethod::Middle: {
            //     // Partition primitives through node's midpoint
            //     Float pmid = (centroidBounds.pMin[dim] + centroidBounds.pMax[dim]) / 2;
            //     auto midIter = std::partition(bvhPrimitives.begin(), bvhPrimitives.end(),
            //                                   [dim, pmid](const BVHPrimitive &pi) {
            //                                       return pi.Centroid()[dim] < pmid;
            //                                   });
            //     mid = midIter - bvhPrimitives.begin();
            //     // For lots of prims with large overlapping bounding boxes, this
            //     // may fail to partition; in that case do not break and fall through
            //     // to EqualCounts.
            //     if (midIter != bvhPrimitives.begin() && midIter != bvhPrimitives.end())
            //         break;
            // }
            // case SplitMethod::EqualCounts: {
            //     // Partition primitives into equally sized subsets
            //     mid = bvhPrimitives.size() / 2;
            //     std::nth_element(bvhPrimitives.begin(), bvhPrimitives.begin() + mid,
            //                      bvhPrimitives.end(),
            //                      [dim](const BVHPrimitive &a, const BVHPrimitive &b) {
            //                          return a.Centroid()[dim] < b.Centroid()[dim];
            //                      });

            //     break;
            // }
            // case SplitMethod::SAH:
            // default: {
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

                // break;
            // }
            // }

            BVHBuildNode *children[2];
            // Recursively build BVHs for _children_
            // if (bvhPrimitives.size() > 128 * 1024) {
            //     // Recursively build child BVHs in parallel
            //     // ParallelFor(0, 2, [&](int i) {
            //         if (i == 0)
            //             children[0] = buildRecursive(
            //                 bvhPrimitives.subspan(0, mid), totalNodes,
            //                 orderedPrimsOffset, orderedPrims);
            //         else
            //             children[1] =
            //                 buildRecursive(threadAllocators, bvhPrimitives.subspan(mid),
            //                                totalNodes, orderedPrimsOffset, orderedPrims);
            //     // });

            // } else {
            // Recursively build child BVHs sequentially
            children[0] =
                buildRecursive(bvhPrimitives.subspan(0, mid),
                               totalNodes, orderedPrimsOffset, orderedPrims);
            children[1] =
                buildRecursive(bvhPrimitives.subspan(mid),
                                totalNodes, orderedPrimsOffset, orderedPrims);
            // }
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


#ifndef RAYSSE
const bool BVHAggregate::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const {
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
#endif

#if defined(__aarch64__)
inline uint32_t neonCompareAndMask(const float32x4_t& a, const float32x4_t& b) {
    uint32x4_t compResUint = vcleq_f32(a, b);
    static const int32x4_t shift = { 0, 1, 2, 3 };
    uint32x4_t tmp = vshrq_n_u32(compResUint, 31);
    return vaddvq_u32(vshlq_u32(tmp, shift));
}
#endif

void rayBBoxIntersect4(const ray& r,
                       const BBox4& bbox4,
                       IVec4& hits,
                       FVec4& tMins,
                       FVec4& tMaxs) {
#if defined(__x86_64__)
    __m128 origin_x = _mm_set1_ps(r.A.x());
    __m128 origin_y = _mm_set1_ps(r.A.y());
    __m128 origin_z = _mm_set1_ps(r.A.z());
    __m128 inv_dir_x = _mm_set1_ps(r.inv_dir.x());
    __m128 inv_dir_y = _mm_set1_ps(r.inv_dir.y());
    __m128 inv_dir_z = _mm_set1_ps(r.inv_dir.z());

    __m128 tMin = _mm_set1_ps(r.tMax);
    __m128 tMax = _mm_set1_ps(Infinity);

    const __m128* minCorner = bbox4.minCornerSSE();
    const __m128* maxCorner = bbox4.maxCornerSSE();

    for (int i = 0; i < 3; ++i) {
        __m128 bmin = r.sign[i] ? maxCorner[i] : minCorner[i];
        __m128 bmax = r.sign[i] ? minCorner[i] : maxCorner[i];
        __m128 origin = (i == 0) ? origin_x : ((i == 1) ? origin_y : origin_z);
        __m128 inv_dir = (i == 0) ? inv_dir_x : ((i == 1) ? inv_dir_y : inv_dir_z);

        __m128 t1 = _mm_mul_ps(_mm_sub_ps(bmin, origin), inv_dir);
        __m128 t2 = _mm_mul_ps(_mm_sub_ps(bmax, origin), inv_dir);

        tMin = _mm_max_ps(tMin, t1);
        tMax = _mm_min_ps(tMax, t2);
    }

    __m128 mask = _mm_cmple_ps(tMin, tMax);
    int hit = _mm_movemask_ps(mask);

    tMins = FVec4(tMin);
    tMaxs = FVec4(tMax);
    hits = IVec4(hit & 1, (hit >> 1) & 1, (hit >> 2) & 1, (hit >> 3) & 1);

#elif defined(__aarch64__)
    float32x4_t origin_x = vdupq_n_f32(r.A.x());
    float32x4_t origin_y = vdupq_n_f32(r.A.y());
    float32x4_t origin_z = vdupq_n_f32(r.A.z());
    float32x4_t inv_dir_x = vdupq_n_f32(r.inv_dir.x());
    float32x4_t inv_dir_y = vdupq_n_f32(r.inv_dir.y());
    float32x4_t inv_dir_z = vdupq_n_f32(r.inv_dir.z());

    float32x4_t tMin = vdupq_n_f32(r.tMax);
    float32x4_t tMax = vdupq_n_f32(Infinity);

    const float32x4_t* minCorner = bbox4.minCornerNeon();
    const float32x4_t* maxCorner = bbox4.maxCornerNeon();

    for (int i = 0; i < 3; ++i) {
        float32x4_t bmin = r.sign[i] ? maxCorner[i] : minCorner[i];
        float32x4_t bmax = r.sign[i] ? minCorner[i] : maxCorner[i];
        float32x4_t origin = (i == 0) ? origin_x : ((i == 1) ? origin_y : origin_z);
        float32x4_t inv_dir = (i == 0) ? inv_dir_x : ((i == 1) ? inv_dir_y : inv_dir_z);

        float32x4_t t1 = vmulq_f32(vsubq_f32(bmin, origin), inv_dir);
        float32x4_t t2 = vmulq_f32(vsubq_f32(bmax, origin), inv_dir);

        tMin = vmaxq_f32(tMin, t1);
        tMax = vminq_f32(tMax, t2);
    }

    uint32_t hit = neonCompareAndMask(tMin, tMax);

    tMins = FVec4(tMin);
    tMaxs = FVec4(tMax);
    hits = IVec4(hit & 1, (hit >> 1) & 1, (hit >> 2) & 1, (hit >> 3) & 1);

#else
    // Fallback implementation for other architectures
    for (int i = 0; i < 4; ++i) {
        Float tmin = r.tMax;
        Float tmax = Infinity;

        for (int j = 0; j < 3; ++j) {
            Float bmin = r.sign[j] ? bbox4.cornersFloatAlt[j+3][i] : bbox4.cornersFloatAlt[j][i];
            Float bmax = r.sign[j] ? bbox4.cornersFloatAlt[j][i] : bbox4.cornersFloatAlt[j+3][i];

            Float t1 = (bmin - r.A[j]) * r.inv_dir[j];
            Float t2 = (bmax - r.A[j]) * r.inv_dir[j];

            tmin = std::max(tmin, std::min(t1, t2));
            tmax = std::min(tmax, std::max(t1, t2));
        }

        hits[i] = tmin <= tmax;
        tMins[i] = tmin;
        tMaxs[i] = tmax;
    }
#endif
}

// std::vector<BVH4Node> transformToBVH4(const LinearBVHNode* bvh2Nodes) {
//     std::vector<BVH4Node> bvh4Nodes;
//     std::function<void(int, int&)> convertNode = [&](int nodeIdx, int& bvh4Idx) {
//         const LinearBVHNode& node = bvh2Nodes[nodeIdx];
        
//         if (node.nPrimitives > 0) {
//             // Leaf node
//             BVH4Node leaf;
//             leaf.bounds.setBBox(0, node.bounds.minCorner(), node.bounds.maxCorner());
//             leaf.primitiveOffset = node.primitivesOffset;
//             leaf.nPrimitives = node.nPrimitives;
//             leaf.nChildren = 0;
//             bvh4Nodes.push_back(leaf);
//             bvh4Idx++;
//         } else {
//             // Interior node
//             BVH4Node interior;
//             interior.nPrimitives = 0;
//             interior.nChildren = 0;
            
//             std::function<void(int, int)> gatherChildren = [&](int childIdx, int depth) {
//                 if (depth == 2 || bvh2Nodes[childIdx].nPrimitives > 0) {
//                     interior.bounds.setBBox(interior.nChildren, 
//                                             bvh2Nodes[childIdx].bounds.minCorner(),
//                                             bvh2Nodes[childIdx].bounds.maxCorner());
//                     interior.children[interior.nChildren] = bvh4Idx + 1;
//                     interior.nChildren++;
//                     convertNode(childIdx, bvh4Idx);
//                 } else {
//                     gatherChildren(bvh2Nodes[childIdx].secondChildOffset, depth + 1);
//                     gatherChildren(childIdx + 1, depth + 1);
//                 }
//             };
            
//             gatherChildren(nodeIdx + 1, 1);
//             gatherChildren(node.secondChildOffset, 1);
            
//             interior.secondChildOffset = bvh4Idx + 1;
//             bvh4Nodes.push_back(interior);
//             bvh4Idx++;
//         }
//     };
    
//     int bvh4Idx = 0;
//     convertNode(0, bvh4Idx);
    
//     return bvh4Nodes;
// }

#ifdef RAYSSE
const bool BVHAggregate::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const {
    if (simdNodes.empty()) return false;

    bool hit_anything = false;
    FVec4 origin(r.origin().x(), r.origin().y(), r.origin().z());
    FVec4 direction(r.direction().x(), r.direction().y(), r.direction().z());
    FVec4 invDirection(r.inv_dir.x(), r.inv_dir.y(), r.inv_dir.z());

    int toVisitOffset = 0;
    int nodesToVisit[64];
    nodesToVisit[toVisitOffset++] = 0;  // Start with root node

    while (toVisitOffset > 0) {
        int currentNodeIndex = nodesToVisit[--toVisitOffset];
        const BVH4Node& node = simdNodes[currentNodeIndex];

        IVec4 hits;
        FVec4 tMins, tMaxs;
        rayBBoxIntersect4(r, node.bounds, hits, tMins, tMaxs);

        // Check if any of the four boxes were hit using SIMD
        #if defined(__x86_64__)
        int hitMask = _mm_movemask_ps(_mm_castsi128_ps(hits.m128i));
        #elif defined(__aarch64__)
        int hitMask = vaddvq_u32(vreinterpretq_u32_s32(hits.i32x4));
        #else
        int hitMask = hits[0] | (hits[1] << 1) | (hits[2] << 2) | (hits[3] << 3);
        #endif

        if (hitMask != 0) {  // If any box was hit
            if (node.nPrimitives > 0) {
                // Leaf node
                for (int j = 0; j < node.nPrimitives; ++j) {
                    hit_record temp_rec;
                    if (primitives[node.primitiveOffset + j]->hit(r, t_min, t_max, temp_rec, rng)) {
                        hit_anything = true;
                        t_max = temp_rec.t;
                        rec = temp_rec;
                    }
                }
            } else {
                // Interior node
                for (int i = 0; i < node.nChildren; ++i) {
                    if (hitMask & (1 << i)) {
                        nodesToVisit[toVisitOffset++] = node.childrenIndices[i];
                    }
                }
            }
        }
    }

    return hit_anything;
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

#ifdef RAYSSE

// Add this function to your BVHAggregate class
void BVHAggregate::transformToSimdFormat() {
    size_t numNodes = n_nodes;
    size_t numSimdNodes = (numNodes + 3) / 4; // Round up to nearest multiple of 4
    std::vector<BBox4> simdNodes(numSimdNodes);

    for (size_t i = 0; i < numSimdNodes; ++i) {
        aabb a = (i * 4 < numNodes) ? nodes[i * 4].bounds : aabb();
        aabb b = (i * 4 + 1 < numNodes) ? nodes[i * 4 + 1].bounds : aabb();
        aabb c = (i * 4 + 2 < numNodes) ? nodes[i * 4 + 2].bounds : aabb();
        aabb d = (i * 4 + 3 < numNodes) ? nodes[i * 4 + 3].bounds : aabb();

        simdNodes[i] = BBox4(a, b, c, d);
    }

    // Store the SIMD nodes
    this->simdNodes = std::move(simdNodes);
}
#endif



// BVHBuildNode *BVHAggregate::buildRecursive(std::span<BVHPrimitive> bvhPrimitives,
//                                            std::atomic<int> *totalNodes,
//                                            std::atomic<int> *orderedPrimsOffset,
//                                            std::vector<std::shared_ptr<hitable> > &orderedPrims) {
//     BVHBuildNode *node = new BVHBuildNode();

//     // Initialize _BVHBuildNode_ for primitive range
//     ++*totalNodes;
//     // Compute bounds of all primitives in BVH node
//     aabb bounds;
//     for (const auto &prim : bvhPrimitives) {
//         bounds = surrounding_box(bounds, prim.bounds);
//     }

//     if (bounds.surface_area() == 0 || bvhPrimitives.size() == 1) {
//         // Create leaf _BVHBuildNode_
//         int firstPrimOffset = orderedPrimsOffset->fetch_add(bvhPrimitives.size());
//         for (size_t i = 0; i < bvhPrimitives.size(); ++i) {
//             int index = bvhPrimitives[i].primitiveIndex;
//             orderedPrims[firstPrimOffset + i] = primitives[index];
//         }
//         node->InitLeaf(firstPrimOffset, bvhPrimitives.size(), bounds);
//         return node;
//     } else {
//         // Compute bound of primitive centroids and choose split dimension _dim_
//         aabb centroidBounds;
//         for (const auto &prim : bvhPrimitives)
//             centroidBounds = surrounding_box(centroidBounds, prim.Centroid());
//         int dim = centroidBounds.MaxDimension();

//         // Partition primitives into two sets and build children
//         if (centroidBounds.max()[dim] == centroidBounds.min()[dim]) {
//             // Create leaf _BVHBuildNode_
//             int firstPrimOffset = orderedPrimsOffset->fetch_add(bvhPrimitives.size());
//             for (size_t i = 0; i < bvhPrimitives.size(); ++i) {
//                 int index = bvhPrimitives[i].primitiveIndex;
//                 orderedPrims[firstPrimOffset + i] = primitives[index];
//             }
//             node->InitLeaf(firstPrimOffset, bvhPrimitives.size(), bounds);
//             return node;

//         } else {
//             int mid = bvhPrimitives.size() / 2;
//             // Partition primitives using approximate SAH
//             if (bvhPrimitives.size() <= 2) {
//                 // Partition primitives into equally sized subsets
//                 mid = bvhPrimitives.size() / 2;
//                 std::nth_element(bvhPrimitives.begin(), bvhPrimitives.begin() + mid,
//                                     bvhPrimitives.end(),
//                                     [dim](const BVHPrimitive &a, const BVHPrimitive &b) {
//                                         return a.Centroid()[dim] < b.Centroid()[dim];
//                                     });

//             } else {
//                 // Allocate _BVHSplitBucket_ for SAH partition buckets
//                 constexpr int nBuckets = 12;
//                 BVHSplitBucket buckets[nBuckets];

//                 // Initialize _BVHSplitBucket_ for SAH partition buckets
//                 for (const auto &prim : bvhPrimitives) {
//                     int b = nBuckets * centroidBounds.offset(prim.Centroid())[dim];
//                     if (b == nBuckets) {
//                         b = nBuckets - 1;
//                     }
//                     // DCHECK_GE(b, 0);
//                     // DCHECK_LT(b, nBuckets);
//                     buckets[b].count++;
//                     buckets[b].bounds = surrounding_box(buckets[b].bounds, prim.bounds);
//                 }

//                 // Compute costs for splitting after each bucket
//                 constexpr int nSplits = nBuckets - 1;
//                 Float costs[nSplits] = {};
//                 // Partially initialize _costs_ using a forward scan over splits
//                 int countBelow = 0;
//                 aabb boundBelow;
//                 for (int i = 0; i < nSplits; ++i) {
//                     boundBelow = surrounding_box(boundBelow, buckets[i].bounds);
//                     countBelow += buckets[i].count;
//                     costs[i] += countBelow * boundBelow.surface_area();
//                 }

//                 // Finish initializing _costs_ using a backward scan over splits
//                 int countAbove = 0;
//                 aabb boundAbove;
//                 for (int i = nSplits; i >= 1; --i) {
//                     boundAbove = surrounding_box(boundAbove, buckets[i].bounds);
//                     countAbove += buckets[i].count;
//                     costs[i - 1] += countAbove * boundAbove.surface_area();
//                 }

//                 // Find bucket to split at that minimizes SAH metric
//                 int minCostSplitBucket = -1;
//                 Float minCost = INFINITY;
//                 for (int i = 0; i < nSplits; ++i) {
//                     // Compute cost for candidate split and update minimum if
//                     // necessary
//                     if (costs[i] < minCost) {
//                         minCost = costs[i];
//                         minCostSplitBucket = i;
//                     }
//                 }
//                 // Compute leaf cost and SAH split cost for chosen split
//                 Float leafCost = bvhPrimitives.size();
//                 minCost = 1.f / 2.f + minCost / bounds.surface_area();

//                 // Either create leaf or split primitives at selected SAH bucket
//                 if (bvhPrimitives.size() > maxPrimsInNode || minCost < leafCost) {
//                     auto midIter = std::partition(
//                         bvhPrimitives.begin(), bvhPrimitives.end(),
//                         [=](const BVHPrimitive &bp) {
//                             int b = nBuckets * centroidBounds.offset(bp.Centroid())[dim];
//                             if (b == nBuckets) {
//                                 b = nBuckets - 1;
//                             }
//                             return b <= minCostSplitBucket;
//                         });
//                     mid = midIter - bvhPrimitives.begin();
//                 } else {
//                     // Create leaf _BVHBuildNode_
//                     int firstPrimOffset =
//                         orderedPrimsOffset->fetch_add(bvhPrimitives.size());
//                     for (size_t i = 0; i < bvhPrimitives.size(); ++i) {
//                         int index = bvhPrimitives[i].primitiveIndex;
//                         orderedPrims[firstPrimOffset + i] = primitives[index];
//                     }
//                     node->InitLeaf(firstPrimOffset, bvhPrimitives.size(), bounds);
//                     return node;
//                 }
//             }
//             BVHBuildNode *children[2];
//             children[0] =
//                 buildRecursive(bvhPrimitives.subspan(0, mid),
//                                totalNodes, orderedPrimsOffset, orderedPrims);
//             children[1] =
//                 buildRecursive(bvhPrimitives.subspan(mid),
//                                 totalNodes, orderedPrimsOffset, orderedPrims);
//             // }
//             node->InitInterior(dim, children[0], children[1]);
//         }
//     }

//     return node;
// }

BVHBuildNode4* ConvertBVH2ToBVH4(BVHBuildNode* node) {
    if (node == nullptr) {
        return nullptr;
    }

    BVHBuildNode4* newNode = new BVHBuildNode4();
    newNode->bounds = node->bounds;

    if (node->nPrimitives > 0) {
        // Leaf node
        newNode->nPrimitives = node->nPrimitives;
        newNode->firstPrimOffset = node->firstPrimOffset;
        newNode->nChildren = 0;
        newNode->splitAxis = node->splitAxis;
        for (int i = 0; i < 4; ++i) {
            newNode->children[i] = nullptr;
        }
    } else {
        // Interior node
        BVHBuildNode* potentialChildren[4];
        int nChildren = 0;

        // Process the two children of the BVH2 node
        for (int i = 0; i < 2; ++i) {
            BVHBuildNode* child = node->children[i];
            if (child->nPrimitives > 0 || (child->children[0] == nullptr && child->children[1] == nullptr)) {
                // If the child is a leaf node or has no further children
                potentialChildren[nChildren++] = child;
            } else {
                // The child is an interior node, add its two children
                potentialChildren[nChildren++] = child->children[0];
                potentialChildren[nChildren++] = child->children[1];
            }
        }

        // Initialize the BVH4 node
        newNode->nChildren = nChildren;
        newNode->nPrimitives = 0;
        newNode->firstPrimOffset = -1;
        newNode->splitAxis = node->splitAxis;
        for (int i = 0; i < nChildren; ++i) {
            // Recursively convert the child nodes
            newNode->children[i] = ConvertBVH2ToBVH4(potentialChildren[i]);
            // Update the bounds
            newNode->bounds = surrounding_box(newNode->bounds, newNode->children[i]->bounds);
        }
        // Set any remaining child pointers to nullptr
        for (int i = nChildren; i < 4; ++i) {
            newNode->children[i] = nullptr;
        }
    }

    return newNode;
}

// int BVHAggregate::flattenBVH4(BVHBuildNode4* node, int* offset) {
//     LinearBVHNode4* linearNode = &nodes4[*offset];
//     int nodeOffset = (*offset)++;

//     linearNode->bounds = node->bounds;
//     linearNode->nPrimitives = node->nPrimitives;
//     linearNode->primitivesOffset = node->firstPrimOffset;
//     linearNode->axis = node->splitAxis;
//     linearNode->nChildren = node->nChildren;

//     if (node->nChildren == 0) {
//         // Leaf node, nothing more to do
//         for (int i = 0; i < 4; ++i) {
//             linearNode->childOffsets[i] = -1; // Invalid offset
//         }
//     } else {
//         // Interior node
//         for (int i = 0; i < node->nChildren; ++i) {
//             linearNode->childOffsets[i] = flattenBVH4(node->children[i], offset);
//         }
//         // Set any unused childOffsets to -1
//         for (int i = node->nChildren; i < 4; ++i) {
//             linearNode->childOffsets[i] = -1;
//         }
//     }

//     return nodeOffset;
// }


// bool BVHAggregate::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const {
//     if (!nodes4) {
//         return false;
//     }

//     // Initialize traversal stack
//     int nodesToVisit[64];
//     int toVisitOffset = 0;
//     int currentNodeIndex = 0;

//     bool any_hit = false;

//     while (true) {
//         const LinearBVHNode4* node = &nodes4[currentNodeIndex];

//         // Check ray against BVH node
//         if (node->bounds.hit(r, t_min, t_max, sampler)) {
//             if (node->nChildren == 0) {
//                 // Leaf node
//                 for (int i = 0; i < node->nPrimitives; ++i) {
//                     hit_record tempRec;
//                     bool hitPrimitive = primitives[node->primitivesOffset + i]->hit(r, t_min, t_max, tempRec, sampler);
//                     if (hitPrimitive) {
//                         any_hit = true;
//                         rec = tempRec;
//                         t_max = rec.t;
//                     }
//                 }
//             } else {
//                 // Interior node
//                 // Push child nodes onto stack
//                 for (int i = 0; i < node->nChildren; ++i) {
//                     int childOffset = node->childOffsets[i];
//                     if (childOffset != -1) {
//                         nodesToVisit[toVisitOffset++] = childOffset;
//                     }
//                 }
//             }
//         }

//         if (toVisitOffset == 0) {
//             break;
//         }
//         currentNodeIndex = nodesToVisit[--toVisitOffset];
//     }

//     return any_hit;
// }
