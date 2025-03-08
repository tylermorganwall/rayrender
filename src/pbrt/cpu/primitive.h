// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#ifndef PBRT_CPU_PRIMITIVE_H
#define PBRT_CPU_PRIMITIVE_H

#include "pbrt.h"

#include "base/light.h"
#include "base/material.h"
#include "base/medium.h"
#include "base/shape.h"
#include "base/texture.h"
#include "util/stats.h"
#include "util/taggedptr.h"
#include "util/transform.h"

#include <memory>

namespace pbrt {

STAT_MEMORY_COUNTER("Memory/Primitives", primitiveMemory);

class SimplePrimitive;
class GeometricPrimitive;
class TransformedPrimitive;
class AnimatedPrimitive;
class BVHAggregate;
class KdTreeAggregate;

// Primitive Definition
class Primitive
    : public TaggedPointer<SimplePrimitive, GeometricPrimitive, TransformedPrimitive,
                           AnimatedPrimitive, BVHAggregate, KdTreeAggregate> {
  public:
    // Primitive Interface
    using TaggedPointer::TaggedPointer;

    Bounds3f Bounds() const;

    pstd::optional<ShapeIntersection> Intersect(const Ray &r,
                                                Float tMax = Infinity) const;
    bool IntersectP(const Ray &r, Float tMax = Infinity) const;
};

// GeometricPrimitive Definition
class GeometricPrimitive {
  public:
    // GeometricPrimitive Public Methods
    GeometricPrimitive(Shape shape, Material material, Light areaLight,
                       const MediumInterface &mediumInterface,
                       FloatTexture alpha = nullptr);
    Bounds3f Bounds() const;
    pstd::optional<ShapeIntersection> Intersect(const Ray &r, Float tMax) const;
    bool IntersectP(const Ray &r, Float tMax) const;

  private:
    // GeometricPrimitive Private Members
    Shape shape;
    Material material;
    Light areaLight;
    MediumInterface mediumInterface;
    FloatTexture alpha;
};

// SimplePrimitive Definition
class SimplePrimitive {
  public:
    // SimplePrimitive Public Methods
    Bounds3f Bounds() const;
    pstd::optional<ShapeIntersection> Intersect(const Ray &r, Float tMax) const;
    bool IntersectP(const Ray &r, Float tMax) const;
    SimplePrimitive(Shape shape, Material material);

  private:
    // SimplePrimitive Private Members
    Shape shape;
    Material material;
};

// TransformedPrimitive Definition
class TransformedPrimitive {
  public:
    // TransformedPrimitive Public Methods
    TransformedPrimitive(Primitive primitive, const Transform *renderFromPrimitive)
        : primitive(primitive), renderFromPrimitive(renderFromPrimitive) {
        primitiveMemory += sizeof(*this);
    }

    pstd::optional<ShapeIntersection> Intersect(const Ray &r, Float tMax) const;
    bool IntersectP(const Ray &r, Float tMax) const;

    Bounds3f Bounds() const { return (*renderFromPrimitive)(primitive.Bounds()); }

  private:
    // TransformedPrimitive Private Members
    Primitive primitive;
    const Transform *renderFromPrimitive;
};

// AnimatedPrimitive Definition
class AnimatedPrimitive {
  public:
    // AnimatedPrimitive Public Methods
    Bounds3f Bounds() const {
        return renderFromPrimitive.MotionBounds(primitive.Bounds());
    }

    AnimatedPrimitive(Primitive primitive, const AnimatedTransform &renderFromPrimitive);
    pstd::optional<ShapeIntersection> Intersect(const Ray &r, Float tMax) const;
    bool IntersectP(const Ray &r, Float tMax) const;

  private:
    // AnimatedPrimitive Private Members
    Primitive primitive;
    AnimatedTransform renderFromPrimitive;
};

}  // namespace pbrt

#endif  // PBRT_CPU_PRIMITIVE_H
