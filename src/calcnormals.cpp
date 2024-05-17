#include "vec3.h"
#include "normal.h"
#include "trianglemesh.h"
#include "assert.h"

// Calculate normal for a given face using three vertices
normal3f CalculateFaceNormal(const point3f& p0, const point3f& p1, const point3f& p2) {
  vec3f v0 = p1 - p0;
  vec3f v1 = p2 - p0;
  vec3f cross_product = cross(v0, v1);
  return unit_vector(cross_product);
}

void CalculateNormals(TriangleMesh *trianglemesh) {
  int vertexCount = trianglemesh->nVertices;
  int triangleCount = trianglemesh->nTriangles;
  
  // Initialize vertex normals to zero
  std::unique_ptr<normal3f[]> vertexNormals = std::make_unique<normal3f[]>(vertexCount);
  for (int i = 0; i < vertexCount; i++) {
    vertexNormals[i] = normal3f(0.f, 0.f, 0.f);
  }
  
  std::vector<int>& triangleArray = trianglemesh->vertexIndices;
  point3f* vertexArray = trianglemesh->p.get();
  
  // Calculate normals for each face and accumulate at each vertex
  for (int k = 0; k < triangleCount * 3; k += 3) {
    int i0 = triangleArray[k+0];
    int i1 = triangleArray[k+1];
    int i2 = triangleArray[k+2];
    ASSERT(i0 < vertexCount);
    ASSERT(i1 < vertexCount);
    ASSERT(i2 < vertexCount);
    
    const point3f& p0 = vertexArray[i0];
    const point3f& p1 = vertexArray[i1];
    const point3f& p2 = vertexArray[i2];
    
    normal3f faceNormal = CalculateFaceNormal(p0, p1, p2);
    
    if(any_is_nan(faceNormal)) {
      continue;
    }
    vertexNormals[i0] += faceNormal;
    vertexNormals[i1] += faceNormal;
    vertexNormals[i2] += faceNormal;
  }
  
  // Normalize the vertex normals
  for (int i = 0; i < vertexCount; i++) {
    vertexNormals[i] = unit_vector(vertexNormals[i]);
  }
  
  // Assign the calculated normals to the triangle mesh
  trianglemesh->has_normals = true;
  trianglemesh->nNormals = vertexCount;
  trianglemesh->n = std::move(vertexNormals);
  trianglemesh->normalIndices = trianglemesh->vertexIndices;
}