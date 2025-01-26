#ifndef TRIANGLEMESHH
#define TRIANGLEMESHH

#ifndef STBIMAGEH
#define STBIMAGEH
#include "stb/stb_image.h"
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

class TextureCache;

struct TriangleMesh {
  // TriangleMesh Public Methods
  TriangleMesh(std::string inputfile, std::string basedir,
               std::shared_ptr<material> default_material, 
               std::shared_ptr<alpha_texture> alpha_default,
               std::shared_ptr<bump_texture> bump_default,
               bool load_materials, bool load_textures, bool load_vertex_colors,  bool load_normals,
               bool verbose, Float scale, 
               bool calculate_consistent_normals, TextureCache& texCache,
               Transform* ObjectToWorld, 
               Transform* WorldToObject, 
               bool reverseOrientation);
  TriangleMesh(Rcpp::NumericMatrix vertices, 
               Rcpp::IntegerMatrix indices, 
               Rcpp::NumericMatrix normals, 
               Rcpp::NumericMatrix texcoords,
               Rcpp::NumericMatrix vertexcolors,
               unsigned char * mesh_texture_data,
               unsigned char * bump_texture_data,
               std::shared_ptr<alpha_texture> alpha,
               std::shared_ptr<bump_texture> bump,
               std::shared_ptr<material> default_material, 
               bool load_materials, bool load_textures, TextureCache& texCache,
               Transform* ObjectToWorld, 
               Transform* WorldToObject, 
               bool reverseOrientation);
  TriangleMesh(float* vertices, 
               int* indices, 
               float* normals, 
               float* texcoords,
               int numVerts, int numIndices,
               std::shared_ptr<alpha_texture> alpha,
               std::shared_ptr<bump_texture> bump,
               std::shared_ptr<material> default_material, 
               Transform* ObjectToWorld, 
               Transform* WorldToObject, 
               bool reverseOrientation);
  TriangleMesh(Rcpp::List raymesh, bool verbose, bool calculate_consistent_normals,
               bool override_material, bool flip_transmittance,
               std::shared_ptr<alpha_texture> alpha,
               std::shared_ptr<bump_texture> bump,
               TextureCache& texCache,
               std::shared_ptr<material> default_material, 
               Transform* ObjectToWorld, 
               Transform* WorldToObject, 
               bool reverseOrientation);
  
  ~TriangleMesh();
  size_t GetSize();
  void ValidateMesh();
    
  // TriangleMesh Data
  size_t nTriangles, nVertices, nNormals, nTex, nTangents;
  bool has_normals, has_tex, has_vertex_colors, has_consistent_normals, has_tangents;
  std::vector<int> vertexIndices;
  std::vector<int> normalIndices;
  std::vector<int> texIndices;

  std::unique_ptr<point3f[]>  p;
  std::unique_ptr<normal3f[]> n;
  std::unique_ptr<normal3f[]> face_n; //For consistent normals
  std::vector<Float> alpha_v; //For consistent normals
  
  std::unique_ptr<vec3f[]>    t; //tangent vector
  std::unique_ptr<point2f[]>  uv;
  std::unique_ptr<point3f[]>  vc;
  
  std::vector<std::shared_ptr<material> > mesh_materials;
  std::vector<Rcpp::List > ray_materials;
  
  std::vector<int > face_material_id;
  
  //Texture Data (from MTL)
  std::vector<unsigned char * > obj_texture_data;
  std::vector<unsigned char * > bump_texture_data;
  std::vector<std::shared_ptr<bump_texture> > bump_textures;
  std::vector<std::shared_ptr<alpha_texture> > alpha_textures;
  size_t texture_size;
  std::vector<bool> material_is_light;
  std::vector<bool> tangent_right_handed;
};

#endif
