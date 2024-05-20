#include "displacement.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"
#define TINYEXR_IMPLEMENTATION
#define TINYEXR_USE_MINIZ 0
#define TINYEXR_USE_STB_ZLIB 1
#include "tinyobj/tinyexr.h"
#ifndef STBIMAGEH
#define STBIMAGEH
#include "stb/stb_image.h"
#endif
#include "assert.h"
#include "trianglemesh.h"
#include <filesystem> // C++17
namespace fs = std::filesystem;
#include "texturecache.h"


void DisplaceMesh(TriangleMesh* base_mesh,
                  std::string displacement_texture,
                  Float displacement_scale,
                  bool displacement_vector) {
  int nx,ny,nn;
  if(!base_mesh->has_tex) {
    throw std::runtime_error("Texcoords required for displacement mapping: no texcoords on mesh.");
  }
  if(base_mesh->nNormals != base_mesh->nVertices  || 
     base_mesh->nNormals != base_mesh->nTex ) {
    throw std::runtime_error("Number of normals (" + std::to_string(base_mesh->nNormals) +
                             ") and UV coords (" +  std::to_string(base_mesh->nTex) + 
                             ") in mesh must be exactly equal to number of vertices (" +  std::to_string(base_mesh->nVertices) + 
                             ")for displacement mapping.");
  }

  TextureCache temp_texcache;
  Float* displacement_data = temp_texcache.LookupFloat(displacement_texture, nx, ny, nn);
  auto tex = std::make_unique<image_texture_float>(displacement_data,
                                                   nx, ny, nn);
  
  if(!displacement_vector) {
    for(size_t i = 0; i < base_mesh->nVertices; i++) {
      const point2f& uv = base_mesh->uv[i];
      const point3f& pp = base_mesh->p[i];
      point3f disp = tex->value(uv.x(),uv.y(), pp);
      
      normal3f displace_n = displacement_scale * unit_vector(base_mesh->n[i]) * disp.x(); 
      base_mesh->p[i] += displace_n.convert_to_vec3();
    }
  } else {
    ASSERT(base_mesh->nVertices == base_mesh->nTex);
    ASSERT(base_mesh->nVertices == base_mesh->tangent_right_handed.size());
    
    for(size_t i = 0; i < base_mesh->nVertices; i++) {
      const point2f& uv = base_mesh->uv[i];
      const point3f& pp = base_mesh->p[i];
      point3f disp = unit_vector(tex->value(uv.x(),uv.y(), pp));
      vec3f tangent = unit_vector(base_mesh->t[i]);
      if(any_is_nan(tangent)) {
        continue;
      }
      vec3f n = unit_vector(base_mesh->n[i].convert_to_vec3());
      vec3f bitangent =  base_mesh->tangent_right_handed[i] ? -cross(n, tangent) : cross(n, tangent);
      
      vec3f displace_n = displacement_scale * (tangent * disp.x() + bitangent * disp.y() + n * disp.z());
      if(any_is_nan(displace_n)) {
        continue;
      }
      base_mesh->p[i] += displace_n;
    }
  }
  base_mesh->has_normals = false;
}