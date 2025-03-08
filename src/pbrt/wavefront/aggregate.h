// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#ifndef PBRT_WAVEFRONT_AGGREGATE_H
#define PBRT_WAVEFRONT_AGGREGATE_H

#include <pbrt/pbrt.h>

#include <pbrt/cpu/primitive.h>
#include <pbrt/lights.h>
#include <pbrt/materials.h>
#include <pbrt/scene.h>
#include <pbrt/util/containers.h>
#include <pbrt/util/pstd.h>
#include <pbrt/util/soa.h>
#include <pbrt/wavefront/integrator.h>
#include <pbrt/wavefront/workitems.h>

#include <map>
#include <string>

namespace pbrt {

class CPUAggregate : public WavefrontAggregate {
  public:
    CPUAggregate(BasicScene &scene, NamedTextures &textures,
                 const std::map<int, pstd::vector<Light> *> &shapeIndexToAreaLights,
                 const std::map<std::string, Medium> &media,
                 const std::map<std::string, pbrt::Material> &namedMaterials,
                 const std::vector<pbrt::Material> &materials);

    Bounds3f Bounds() const { return aggregate ? aggregate.Bounds() : Bounds3f(); }

    void IntersectClosest(int maxRays, const RayQueue *rayQueue,
                          EscapedRayQueue *escapedRayQueue,
                          HitAreaLightQueue *hitAreaLightQueue,
                          MaterialEvalQueue *basicEvalMaterialQueue,
                          MaterialEvalQueue *universalEvalMaterialQueue,
                          MediumSampleQueue *mediumSampleQueue,
                          RayQueue *nextRayQueue) const;

    void IntersectShadow(int maxRays, ShadowRayQueue *shadowRayQueue,
                         SOA<PixelSampleState> *pixelSampleState) const;

    void IntersectShadowTr(int maxRays, ShadowRayQueue *shadowRayQueue,
                           SOA<PixelSampleState> *pixelSampleState) const;

    void IntersectOneRandom(int maxRays,
                            SubsurfaceScatterQueue *subsurfaceScatterQueue) const;

  private:
    Primitive aggregate;
};

}  // namespace pbrt

#endif  // PBRT_WAVEFRONT_AGGREGATE_H
