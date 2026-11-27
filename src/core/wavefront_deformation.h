#ifndef RAYRENDER_WAVEFRONT_DEFORMATION_H
#define RAYRENDER_WAVEFRONT_DEFORMATION_H
#include "wavefront.h"
#include "../math/transform.h"
class mesh3d;
class triangle;

// Native pointers remain valid for the containing WavefrontContext's scene
// lifetime. Bind once in compiler order, including temporarily degenerate faces.
class WavefrontDeformation {
public:
  struct Mesh {
    const mesh3d *source;
    uint64_t revision;
    Transform placement;
    WFMeshRange range;
    uint32_t prototype = UINT32_MAX;
    Float orientation = 1;
    std::vector<const triangle *> faces;
  };
  std::vector<Mesh> meshes;
  std::vector<WFMeshRange> Update(WavefrontScene &);
};

// Update only geometric attributes; preserve UVs, material/light IDs, normal
// cone offsets and face identity. Returns false for a zero-area face.
bool UpdateWavefrontTriangle(const triangle &, const Transform &, WFTriangle &);
#if defined(NOT_CRAN) && defined(RAY_HAS_METAL_BVH)
std::vector<WFGeometryResult> ProbeMetalSessionGeometry(WavefrontSession &,
    const std::vector<WFGeometryProbe> &);
#endif
#endif
