#ifndef TRIANGLEMESHH
#define TRIANGLEMESHH

#ifndef STBIMAGEH
#define STBIMAGEH
#include "stb_image.h"
#endif

#include "texture.h"
#include "transform.h"
#include "hitablelist.h"
#include "material.h"

inline char separator() {
#if defined _WIN32 || defined __CYGWIN__
  return '\\';
#else
  return '/';
#endif
}

struct TriangleMesh {
  // TriangleMesh Public Methods
  TriangleMesh(std::string inputfile, std::string basedir,
               std::shared_ptr<material> default_material, 
               std::shared_ptr<Transform> ObjectToWorld, 
               std::shared_ptr<Transform> WorldToObject, 
               bool reverseOrientation);
  ~TriangleMesh() {
    for(auto tex : obj_texture_data) {
      if(tex) stbi_image_free(tex);
    }
    for(auto bump : bump_texture_data) {
      if(bump) stbi_image_free(bump);
    }
  }
  // TriangleMesh Data
  size_t nTriangles, nVertices;
  std::vector<int> vertexIndices;
  std::vector<int> normalIndices;
  std::vector<int> texIndices;

  std::unique_ptr<point3f[]>  p;
  std::unique_ptr<normal3f[]> n;
  std::unique_ptr<vec3f[]>    s;
  std::unique_ptr<point2f[]>  uv;
  hitable_list triangles;
  
  // std::shared_ptr<alpha_texture> alphaMask;
  // std::vector<int> faceIndices;
  std::vector<std::shared_ptr<material> > mtl_materials;
  std::vector<size_t > face_material_id;
  
  //Texture Data (from MTL)
  std::vector<Float* > obj_texture_data;
  std::vector<Float* > bump_texture_data;
  std::vector<std::shared_ptr<bump_texture> > bump_textures;
  std::vector<std::shared_ptr<alpha_texture> > alpha_textures;
  size_t texture_size;
};

#endif