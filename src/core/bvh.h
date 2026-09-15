#ifndef BVHH
#define BVHH
#include "../hitables/hitable.h"
#include "../math/aabb.h"
#include "../math/simd.h"
#include "../utils/assert.h"
#include <atomic>
#include <memory>
#include <span>

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

// struct BVH4Node {
//     BBox4 bounds;
//     union {
//         int primitiveOffset;    // Leaf
//         int secondChildOffset;  // Interior
//     };
//     uint8_t nPrimitives;  // 0 -> interior node
//     uint8_t nChildren;    // Number of valid children (1-4)
//     uint8_t children[4];  // Indices of child nodes
// };


struct BVHBuildNode {
    ~BVHBuildNode() {
        if(children[0]) {
            delete children[0];
        }
        if(children[1]) {
            delete children[1];
        }
    }
    // BVHBuildNode Public Methods
    void InitLeaf(int first, int n, const aabb& b) {
        firstPrimOffset = first;
        nPrimitives = n;
        bounds = b;
        children[0] = children[1] = nullptr;
    }

    void InitInterior(int axis, BVHBuildNode* c0, BVHBuildNode* c1) {
        children[0] = c0;
        children[1] = c1;
        bounds = surrounding_box(c0->bounds, c1->bounds);
        splitAxis = axis;
        firstPrimOffset = -1; // Invalid for interior nodes
        nPrimitives = 0;      // Interior nodes contain no primitives directly
    }

    aabb bounds;
    BVHBuildNode *children[2];
    int splitAxis, firstPrimOffset, nPrimitives;
};

struct BVHBuildNode4 {
    BVHBuildNode4()
        : bounds(),
          children{nullptr, nullptr, nullptr, nullptr},
          nChildren(0),
          splitAxis(0),
          firstPrimOffset(-1),
          nPrimitives(0) {}

    ~BVHBuildNode4() {
        for(int i = 0; i < 4; i++) {
            if(children[i]) {
                delete children[i];
            }
        }
    }
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

// Only interior nodes carry four child bounds. Nonnegative child references
// index this array; -1 is unused, and -2 - leafIndex refers to a compact leaf.
struct alignas(64) LinearBVHNode4 {
    IVec4 childOffsets;
    BBox4 bbox4;
    int nChildren;
};

// The parent already tested the leaf's bounds, so a leaf needs only its range
// in the existing ordered primitive array. Keep full counts for coincident
// primitives, whose leaves can exceed the requested maximum leaf size.
struct LinearBVHLeaf4 {
    int primitivesOffset;
    int nPrimitives;
};

static_assert(sizeof(LinearBVHLeaf4) == 8, "BVH4 leaves should occupy eight bytes");

class BVHAggregate : public hitable {
public:
    BVHAggregate(std::vector<std::shared_ptr<hitable> > prims,
                float t_min, float t_max, 
                int maxPrimsInNode, bool sah, 
                Transform* ObjectToWorld, 
                Transform* WorldToObject, 
                bool reverseOrientation);
                
    BVHAggregate(std::vector<std::shared_ptr<hitable> > prims,
                float t_min, float t_max, 
                int maxPrimsInNode, bool sah);
    
    // static BVHAggregate *Create(std::vector<Primitive> prims,
    //                         const ParameterDictionary &parameters);
    
    // Bounds3f Bounds() const;
    // std::optional<ShapeIntersection> Intersect(const Ray &ray, Float tMax) const;
    // bool IntersectP(const Ray &ray, Float tMax) const;

    virtual const bool hit(const Ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const;
    virtual const bool hit(const Ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const;
    virtual bool HitP(const Ray &r, Float t_min, Float t_max, random_gen& rng) const;
    virtual bool HitP(const Ray &r, Float t_min, Float t_max, Sampler* sampler) const;
    OpaqueShadowType ShadowType() const { return shadow_type; }
    bool OpaqueHit(const Ray&, Float, Float, random_gen&) const;
    // Mesh construction releases its temporary list; light sampling borrows
    // these same owned primitives instead of retaining a duplicate mesh list.
    std::span<const std::shared_ptr<hitable>> Primitives() const { return primitives; }

    virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
    
    Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0);
    Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0);
    vec3f random(const point3f& o, random_gen& rng, Float time = 0);
    vec3f random(const point3f& o, Sampler* sampler, Float time = 0);
    virtual void hitable_info_bounds(Float t0, Float t1) const {
      aabb box_top;
      bounding_box(t0, t1, box_top);
      Rcpp::Rcout << GetName() << ": " <<  box_top.min() << "-" << box_top.max() << "\n";
      for(size_t i = 0; i < primitives.size(); i++) {
        aabb box;
        primitives[i]->bounding_box(t0, t1, box);
        Rcpp::Rcout << "   " << primitives[i]->GetName() << ": " <<  box.min() << "-" << box.max() << "\n";
      }
    }
    // void validate_bvh();
    
    std::string GetName() const {
        return(std::string("BVH"));
    }
    size_t GetSize() {
        return(0);
    };
    void transformToSimdFormat();
    aabb scene_bounds;
    // std::vector<BVH4Node> simdNodes;
    int n_nodes;
    // std::pair<size_t,size_t> CountNodeLeaf();
private:
    void classifyOpaqueShadow();
    OpaqueShadowType shadow_type = OpaqueShadowType::Opaque;
    BVHBuildNode *buildRecursive(std::span<BVHPrimitive> bvhPrimitives,
                                std::atomic<int> *totalNodes,
                                std::atomic<int> *orderedPrimsOffset,
                                std::vector<std::shared_ptr<hitable> > &orderedPrims);
    BVHBuildNode4* ConvertBVH2ToBVH4(BVHBuildNode* node, int* totalNodes4,
                                   int* totalLeaves4);
    void buildBVH4(BVHBuildNode* root);
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
    int flattenBVH4(BVHBuildNode4* node, int* offset, int* leafOffset);
    void validateBVH4() const;

    int maxPrimsInNode;
    int totalNodes4 = 0;
    int totalLeaves4 = 0;
    int root4 = -1;
    std::vector<std::shared_ptr<hitable> > primitives;
    //    SplitMethod splitMethod;
       std::unique_ptr<LinearBVHNode[]> nodes;
       std::unique_ptr<LinearBVHNode4[]> nodes4;
       std::unique_ptr<LinearBVHLeaf4[]> leaves4;
    //    int totalNodes;

};



#endif
