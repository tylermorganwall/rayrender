#ifndef LOOPSUBDIVH
#define LOOPSUBDIVH

#include "../hitables/hitable.h"
#include "../hitables/trianglemesh.h"

void LoopSubdivide(TriangleMesh* base_mesh,
                   const int nLevels,
                   bool verbose);

#endif  