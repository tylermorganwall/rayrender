#ifndef LOOPSUBDIVH
#define LOOPSUBDIVH

#include "hitable.h"
#include "trianglemesh.h"

void LoopSubdivide(TriangleMesh* base_mesh,
                   const int nLevels,
                   bool verbose);

#endif  