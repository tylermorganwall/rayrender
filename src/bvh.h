#ifndef BVHH
#define BVHH
#include <memory>
#include <span>
#include "assert.h"
#include "hitable.h"
#include "simd.h"

struct BVHPrimitive {
    BVHPrimitive() : primitiveIndex(0), bounds(aabb()) {};
    BVHPrimitive(size_t primitiveIndex, const aabb &bounds)
        : primitiveIndex(primitiveIndex), bounds(bounds) {}
    size_t primitiveIndex;
    aabb bounds;
    
    point3f Centroid() const { return .5f * bounds.min() + .5f * bounds.max(); }
};

struct BVHSplitBucket {
    int count = 0;
    aabb bounds;
};

struct BVH4Node {
    BBox4 bounds;
    union {
        int primitiveOffset;    // Leaf
        int secondChildOffset;  // Interior
    };
    uint8_t nPrimitives;  // 0 -> interior node
    uint8_t nChildren;    // Number of valid children (1-4)
    uint8_t children[4];  // Indices of child nodes
};


struct BVHBuildNode {
    // BVHBuildNode Public Methods
    void InitLeaf(int first, int n, const aabb& b) {
        firstPrimOffset = first;
        nPrimitives = n;
        bounds = b;
        children[0] = children[1] = nullptr;

        // Debug print
        Rcpp::Rcout << "InitLeaf called: firstPrimOffset = " << firstPrimOffset
                    << ", nPrimitives = " << nPrimitives << "\n";
    }

    void InitInterior(int axis, BVHBuildNode* c0, BVHBuildNode* c1) {
        children[0] = c0;
        children[1] = c1;
        bounds = surrounding_box(c0->bounds, c1->bounds);
        splitAxis = axis;
        firstPrimOffset = -1; // Invalid for interior nodes
        nPrimitives = 0;      // Interior nodes contain no primitives directly

        // Debug print
        Rcpp::Rcout << "InitInterior called: splitAxis = " << splitAxis
                    << ", firstPrimOffset = " << firstPrimOffset
                    << ", nPrimitives = " << nPrimitives << "\n";
    }

    aabb bounds;
    BVHBuildNode *children[2];
    int splitAxis, firstPrimOffset, nPrimitives;
};

struct BVHBuildNode4 {
    aabb bounds;
    BVHBuildNode4* children[4]; // Up to 4 children
    int nChildren;              // Number of valid children (1 to 4)
    int splitAxis;              // Axis along which the node was split
    int firstPrimOffset;        // Offset in primitives array (for leaf nodes)
    int nPrimitives;            // Number of primitives in leaf node
};

struct LinearBVHNode {
    aabb bounds;
    union {
        int primitivesOffset;    // leaf
        int secondChildOffset;   // interior
    };
    uint16_t nPrimitives;  // 0 -> interior node
    uint8_t axis;          // interior node: xyz
};

struct LinearBVHNode4 {
    aabb bounds;          // Bounding box of this node
    uint8_t nPrimitives;  // 0 for interior nodes
    uint8_t axis;         // Split axis
    uint8_t nChildren;    // Number of children (up to 4)
    int primitivesOffset; // For leaf nodes
    int childOffsets[4];  // For interior nodes
    BBox4 bbox4;          // Packed bounding boxes of child nodes
};



class BVHAggregate : public hitable {
  public:
        BVHAggregate(std::vector<std::shared_ptr<hitable> > prims,
                    float t_min, float t_max, 
                    int maxPrimsInNode, bool sah, 
                    std::shared_ptr<Transform> ObjectToWorld, 
                    std::shared_ptr<Transform> WorldToObject, 
                    bool reverseOrientation);
        
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
        size_t GetSize() {
            return(0);
        };
        void transformToSimdFormat();
        aabb scene_bounds;
        std::vector<BVH4Node> simdNodes;
        int n_nodes;
        // std::pair<size_t,size_t> CountNodeLeaf();
  private:
       BVHBuildNode *buildRecursive(std::span<BVHPrimitive> bvhPrimitives,
                                    std::atomic<int> *totalNodes,
                                    std::atomic<int> *orderedPrimsOffset,
                                    std::vector<std::shared_ptr<hitable> > &orderedPrims);
       BVHBuildNode4* ConvertBVH2ToBVH4(BVHBuildNode* node, int* totalNodes4);
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
       int flattenBVH4(BVHBuildNode4* node, int* offset);
       void validateBVH4() const;

       int maxPrimsInNode;
       int totalNodes4;
       std::vector<std::shared_ptr<hitable> > primitives;
    //    SplitMethod splitMethod;
       std::unique_ptr<LinearBVHNode[]> nodes;
       std::unique_ptr<LinearBVHNode4[]> nodes4;
    //    int totalNodes;

};



#endif