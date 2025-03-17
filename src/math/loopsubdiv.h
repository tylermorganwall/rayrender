#ifndef LOOPSUBDIVH
#define LOOPSUBDIVH

#include "../math/hitable.h"
#include "../math/trianglemesh.h"

void LoopSubdivide(TriangleMesh* base_mesh,
                   const int nLevels,
                   bool verbose);

#endif  