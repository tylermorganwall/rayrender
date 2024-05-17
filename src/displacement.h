#ifndef DISPLACEMENTH
#define DISPLACEMENTH

#include "float.h"
#include <string>
struct TriangleMesh;


void DisplaceMesh(TriangleMesh* base_mesh,
                  std::string displacement_texture,
                  Float displacement_scale,
                  bool displacement_vector);

#endif