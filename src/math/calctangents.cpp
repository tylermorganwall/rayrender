#include "../math/calctangents.h"
#include "../math/vec3.h"
#include "../math/point3.h"
#include "../math/point2.h"
#include "../hitables/trianglemesh.h"
#include "../math/mathinline.h"
#include "../utils/assert.h"

vec3f Reject(vec3f t, normal3f n) {
  return(unit_vector(t - dot(t,n) * convert_to_vec3(n)));
}

void CalculateTangents(TriangleMesh *trianglemesh) {
  int vertexCount = trianglemesh->nVertices;
  
  int triangleCount = trianglemesh->nTriangles;

  std::unique_ptr<vec3f[]> tangent = std::make_unique<vec3f[]>(vertexCount); 
  std::unique_ptr<vec3f[]> bitangent = std::make_unique<vec3f[]>(vertexCount); 
  
  for (int i = 0; i < vertexCount; i++) {
    tangent[i] = vec3f(0.f, 0.f, 0.f);
    bitangent[i] = vec3f(0.f, 0.f, 0.f); 
  }
  
  std::vector<int>& triangleArray = trianglemesh->vertexIndices;
  point3f* vertexArray = trianglemesh->p.get();
  point2f* texcoordArray = trianglemesh->uv.get();
  normal3f* normalArray = trianglemesh->n.get();
  std::vector<bool>& tangent_right_handed = trianglemesh->tangent_right_handed;
  tangent_right_handed.clear();
  tangent_right_handed.resize(vertexCount);
  std::unique_ptr<vec3f[]> tangentArray = std::make_unique<vec3f[]>(vertexCount);

  // Calculate tangent and bitangent for each triangle and add to all three vertices.
  for (int k = 0; k < triangleCount * 3; k += 3) {
    int i0 = triangleArray[k+0];
    int i1 = triangleArray[k+1];
    int i2 = triangleArray[k+2];
    ASSERT(i0 < (int)trianglemesh->vertexIndices.size());
    ASSERT(i1 < (int)trianglemesh->vertexIndices.size());
    ASSERT(i2 < (int)trianglemesh->vertexIndices.size());
    ASSERT(i0 < vertexCount);
    ASSERT(i1 < vertexCount);
    ASSERT(i2 < vertexCount);
    
    const point3f& p0 = vertexArray[i0];
    const point3f& p1 = vertexArray[i1];
    const point3f& p2 = vertexArray[i2];
    const point2f& w0 = texcoordArray[i0];
    const point2f& w1 = texcoordArray[i1];
    const point2f& w2 = texcoordArray[i2];
    vec3f e1 = p1 - p0;
    vec3f e2 = p2 - p0;
    Float x1 = w1.xy.x - w0.xy.x;
    Float x2 = w2.xy.x - w0.xy.x;
    Float y1 = w1.xy.y - w0.xy.y; 
    Float y2 = w2.xy.y - w0.xy.y;
    Float inv_dop = DifferenceOfProducts(x1, y2, x2, y1);
    if(inv_dop == 0) {
      continue;
    }
    Float r = 1.f / inv_dop;
    vec3f t = (e1 * y2 - e2 * y1) * r;
    vec3f b = (e2 * x1 - e1 * x2) * r;
    tangent[i0] += t;
    tangent[i1] += t;
    tangent[i2] += t;
    bitangent[i0] += b;
    bitangent[i1] += b;
    bitangent[i2] += b;
  }
  // Orthonormalize each tangent and calculate the handedness.
  for (int i = 0; i < vertexCount; i++) {
    const vec3f& t = tangent[i];
    const vec3f& b = bitangent[i];
    const normal3f& n = normalArray[i];
    tangentArray[i] = Reject(t, n); 
    tangent_right_handed[i] = dot(cross(t, b), n) > 0.f;
  }
  trianglemesh->has_tangents = true;
  trianglemesh->nTangents = vertexCount;
  trianglemesh->t = std::move(tangentArray);
}