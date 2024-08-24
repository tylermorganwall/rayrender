#ifndef BVHH
#define BVHH
#include <memory>

struct BVHPrimitive {
    BVHPrimitive(size_t primitiveIndex, const aabb &bounds)
        : primitiveIndex(primitiveIndex), bounds(bounds) {}
    size_t primitiveIndex;
    aabb bounds;
    
    point3f Centroid() const { return .5f * bounds.min() + .5f * bounds.max(); }
};

struct BVHSplitBucket {
    int count = 0;
    Bounds3f bounds;
};


struct BVHBuildNode {
    // BVHBuildNode Public Methods
    void InitLeaf(int first, int n, const aabb &b) {
        firstPrimOffset = first;
        nPrimitives = n;
        bounds = b;
        children[0] = children[1] = nullptr;
        ++leafNodes;
        ++totalLeafNodes;
        totalPrimitives += n;
    }

    void InitInterior(int axis, BVHBuildNode *c0, BVHBuildNode *c1) {
        children[0] = c0;
        children[1] = c1;
        bounds = surrounding_box(c0->bounds, c1->bounds);
        splitAxis = axis;
        nPrimitives = 0;
        ++interiorNodes;
    }

    aabb bounds;
    BVHBuildNode *children[2];
    int splitAxis, firstPrimOffset, nPrimitives;
};

struct alignas(32) LinearBVHNode {
    aabb bounds;
    union {
        int primitivesOffset;    // leaf
        int secondChildOffset;   // interior
    };
    uint16_t nPrimitives;  // 0 -> interior node
    uint8_t axis;          // interior node: xyz
};

class BVHAggregate : hitable {
  public:
        BVHAggregate(std::vector<std::shared_ptr<hitable> >& p, 
                    int maxPrimsInNode = 1);
        
        // static BVHAggregate *Create(std::vector<Primitive> prims,
        //                         const ParameterDictionary &parameters);
        
        // Bounds3f Bounds() const;
        // pstd::optional<ShapeIntersection> Intersect(const Ray &ray, Float tMax) const;
        // bool IntersectP(const Ray &ray, Float tMax) const;

        virtual const bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const;
        virtual const bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const;

        virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
        
        Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0);
        Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0);
        vec3f random(const point3f& o, random_gen& rng, Float time = 0);
        vec3f random(const point3f& o, Sampler* sampler, Float time = 0);
        // void validate_bvh();
        
        std::string GetName() const {
            return(std::string("BVH"));
        }
        size_t GetSize();
        // std::pair<size_t,size_t> CountNodeLeaf();

  private:
       BVHBuildNode *buildRecursive(//ThreadLocal<Allocator> &threadAllocators,
                                    //pstd::span<BVHPrimitive> bvhPrimitives,
                                    std::span<bvh_node> bvhPrimitives,
                                    std::atomic<int> *totalNodes,
                                    std::atomic<int> *orderedPrimsOffset,
                                    std::vector<std::shared_ptr<hitable> > &orderedPrims);
    //    BVHBuildNode *buildHLBVH(Allocator alloc,
    //                             const std::vector<BVHPrimitive> &primitiveInfo,
    //                             std::atomic<int> *totalNodes,
    //                             std::vector<Primitive> &orderedPrims);
    //    BVHBuildNode *emitLBVH(BVHBuildNode *&buildNodes,
    //                           const std::vector<BVHPrimitive> &primitiveInfo,
    //                           MortonPrimitive *mortonPrims, int nPrimitives, int *totalNodes,
    //                           std::vector<Primitive> &orderedPrims,
    //                           std::atomic<int> *orderedPrimsOffset, int bitIndex);
    //    BVHBuildNode *buildUpperSAH(Allocator alloc,
    //                                std::vector<BVHBuildNode *> &treeletRoots, int start,
    //                                int end, std::atomic<int> *totalNodes) const;
       int flattenBVH(BVHBuildNode *node, int *offset);

       int maxPrimsInNode;
       std::vector<std::shared_ptr<hitable> > primitives;
    //    SplitMethod splitMethod;
       LinearBVHNode *nodes = nullptr;
       int totalNodes;

};

BVHAggregate::BVHAggregate(std::vector<std::shared_ptr<hitable> > prims,
        int maxPrimsInNode, bool sah) : 
        maxPrimsInNode(std::min(255, maxPrimsInNode)),
        primitives(std::move(prims))//, splitMethod(splitMethod) 
        { 
    std::vector<BVHPrimitive> bvhPrimitives(primitives.size());
    for (size_t i = 0; i < primitives.size(); ++i)
        bvhPrimitives[i] = BVHPrimitive(i, primitives[i].Bounds());

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
    nodes = new LinearBVHNode[totalNodes];
    int offset = 0;
    flattenBVH(root, &offset);
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

    if (bounds.SurfaceArea() == 0 || bvhPrimitives.size() == 1) {
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
        if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim]) {
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
                    int b = nBuckets * centroidBounds.Offset(prim.Centroid())[dim];
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
                    costs[i] += countBelow * boundBelow.SurfaceArea();
                }

                // Finish initializing _costs_ using a backward scan over splits
                int countAbove = 0;
                aabb boundAbove;
                for (int i = nSplits; i >= 1; --i) {
                    boundAbove = surrounding_box(boundAbove, buckets[i].bounds);
                    countAbove += buckets[i].count;
                    costs[i - 1] += countAbove * boundAbove.SurfaceArea();
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
                minCost = 1.f / 2.f + minCost / bounds.SurfaceArea();

                // Either create leaf or split primitives at selected SAH bucket
                if (bvhPrimitives.size() > maxPrimsInNode || minCost < leafCost) {
                    auto midIter = std::partition(
                        bvhPrimitives.begin(), bvhPrimitives.end(),
                        [=](const BVHPrimitive &bp) {
                            int b =
                                nBuckets * centroidBounds.Offset(bp.Centroid())[dim];
                            if (b == nBuckets)
                                b = nBuckets - 1;
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
            }

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

const bool bvh_node::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const {
    if (!nodes) {
        return false;
    }
    // Follow ray through BVH nodes to find primitive intersections
    int toVisitOffset = 0, currentNodeIndex = 0;
    int nodesToVisit[64];
    int nodesVisited = 0;
    while (true) {
        ++nodesVisited;
        const LinearBVHNode *node = &nodes[currentNodeIndex];
        // Check ray against BVH node
        if (node->bounds.hit(r, t_min, t_max, rec, rng)) {
            if (node->nPrimitives > 0) {
                // Intersect ray with primitives in leaf BVH node
                for (int i = 0; i < node->nPrimitives; ++i) {
                    // Check for intersection with primitive in BVH node
                    hit_record hrec_temp;
                    bool primSi = primitives[node->primitivesOffset + i]->hit(r, t_min, t_max, hrec_temp, rng);
                    if (primSi) {
                        rec = hrec_temp;
                        t_max = rec->t;
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

    bvhNodesVisited += nodesVisited;
    return si;
}

const bool bvh_node::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const {
    if (!nodes) {
        return false;
    }
    // Follow ray through BVH nodes to find primitive intersections
    int toVisitOffset = 0, currentNodeIndex = 0;
    int nodesToVisit[64];
    int nodesVisited = 0;
    while (true) {
        ++nodesVisited;
        const LinearBVHNode *node = &nodes[currentNodeIndex];
        // Check ray against BVH node
        if (node->bounds.hit(r, t_min, t_max, rec, sampler)) {
            if (node->nPrimitives > 0) {
                // Intersect ray with primitives in leaf BVH node
                for (int i = 0; i < node->nPrimitives; ++i) {
                    // Check for intersection with primitive in BVH node
                    hit_record hrec_temp;
                    bool primSi = primitives[node->primitivesOffset + i]->hit(r, t_min, t_max, hrec_temp, sampler);
                    if (primSi) {
                        rec = hrec_temp;
                        t_max = rec->t;
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

    bvhNodesVisited += nodesVisited;
    return si;
}

Float BVHAggregate::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) {
  return(0.5*left->pdf_value(o,v, rng, time) + 0.5*right->pdf_value(o,v, rng, time));
}

Float BVHAggregate::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) {
  return(0.5*left->pdf_value(o,v, sampler, time) + 0.5*right->pdf_value(o,v, sampler, time));
  
}

vec3f BVHAggregate::random(const point3f& o, random_gen& rng, Float time) {
  return(rng.unif_rand() > 0.5 ? left->random(o, rng, time) : right->random(o, rng, time));
}

vec3f BVHAggregate::random(const point3f& o, Sampler* sampler, Float time) {
  return(sampler->Get1D() > 0.5 ? left->random(o, sampler, time) : right->random(o, sampler, time));
  
}

#endif