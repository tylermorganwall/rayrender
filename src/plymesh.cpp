#define TINYPLY_IMPLEMENTATION

#include "miniply.h"
#include "plymesh.h"

inline char separator_ply() {
#if defined _WIN32 || defined __CYGWIN__
  return '\\';
#else
  return '/';
#endif
}


enum class Topology {
  Soup,   // Every 3 indices specify a triangle.
  Strip,  // Triangle strip, triangle i uses indices i, i-1 and i-2
  Fan,    // Triangle fan, triangle i uses indices, i, i-1 and 0.
};


struct TriMesh {
  // Per-vertex data
  float* pos          = nullptr; // has 3*numVerts elements.
  float* normal       = nullptr; // if non-null, has 3 * numVerts elements.
  float* uv           = nullptr; // if non-null, has 2 * numVerts elements.
  uint32_t numVerts   = 0;
  
  // Per-index data
  int* indices        = nullptr; // has numIndices elements.
  uint32_t numIndices = 0; // number of indices = 3 times the number of faces.
  
  Topology topology  = Topology::Soup; // How to interpret the indices.
  bool hasTerminator = false;          // Only applies when topology != Soup.
  int terminator     = -1;             // Value indicating the end of a strip or fan. Only applies when topology != Soup.
  
  ~TriMesh() {
    delete[] pos;
    delete[] normal;
    delete[] uv;
    delete[] indices;
  }
  
  bool all_indices_valid() const {
    bool checkTerminator = topology != Topology::Soup && hasTerminator && (terminator < 0 || terminator >= int(numVerts));
    for (uint32_t i = 0; i < numIndices; i++) {
      if (checkTerminator && indices[i] == terminator) {
        continue;
      }
      if (indices[i] < 0 || uint32_t(indices[i]) >= numVerts) {
        return false;
      }
    }
    return true;
  }
};


static TriMesh* parse_file_with_miniply(const char* filename, bool assumeTriangles) {
  miniply::PLYReader reader(filename);
  if (!reader.valid()) {
    Rcpp::Rcout << "Not valid reader \n";
    return nullptr;
  }
  
  uint32_t indexes[3];
  bool gotVerts = false, gotFaces = false;
  
  TriMesh* trimesh = new TriMesh();
  while (reader.has_element() && (!gotVerts || !gotFaces)) {
    if (reader.element_is(miniply::kPLYVertexElement) && reader.load_element() && reader.find_pos(indexes)) {
      trimesh->numVerts = reader.num_rows();
      trimesh->pos = new float[trimesh->numVerts * 3];
      reader.extract_properties(indexes, 3, miniply::PLYPropertyType::Float, trimesh->pos);
      if (reader.find_texcoord(indexes)) {
        trimesh->uv = new float[trimesh->numVerts * 2];
        reader.extract_properties(indexes, 2, miniply::PLYPropertyType::Float, trimesh->uv);
      }
      gotVerts = true;
    } else if (reader.element_is(miniply::kPLYFaceElement) && reader.load_element() && reader.find_indices(indexes)) {
      uint32_t propIdx = 1; 
      bool polys = reader.requires_triangulation(propIdx);
      if (polys && !gotVerts) {
        Rcpp::Rcout << "Error: need vertex positions to triangulate faces.\n";
        break;
      }
      if (polys) {
        trimesh->numIndices = reader.num_triangles(indexes[0]) * 3;
        trimesh->indices = new int[trimesh->numIndices];
        reader.extract_triangles(indexes[0], trimesh->pos, trimesh->numVerts, miniply::PLYPropertyType::Int, trimesh->indices);
      } else {
        trimesh->numIndices = reader.num_rows() * 3;
        trimesh->indices = new int[trimesh->numIndices];
        reader.extract_list_property(indexes[0], miniply::PLYPropertyType::Int, trimesh->indices);
      }
      gotFaces = true;
    }
    if (gotVerts && gotFaces) {
      break;
    }
    reader.next_element();
  }
  if (!gotVerts || !gotFaces) {
    std::string vert1 = gotVerts ? "" : "vertices ";
    std::string face1 = gotFaces ? "" : "faces";
    Rcpp::Rcout << "Failed to load: " << vert1 << face1 << "\n";
    delete trimesh;
    return nullptr;
  }
  
  return trimesh;
}


plymesh::plymesh(std::string inputfile, std::string basedir, std::shared_ptr<material> mat, 
            Float scale, Float shutteropen, Float shutterclose, int bvh_type, random_gen rng,
            std::shared_ptr<Transform> ObjectToWorld, std::shared_ptr<Transform> WorldToObject, bool reverseOrientation) :
  hitable(ObjectToWorld, WorldToObject, reverseOrientation) {
  TriMesh* tri = parse_file_with_miniply(inputfile.c_str(), false);
  mat_ptr = mat;
  
  if(tri == nullptr) {
    std::string err = inputfile;
    throw std::runtime_error("No mesh loaded: " + err);
  }
  bool has_normals = false;
  if(tri->normal != nullptr) {
    has_normals = true;
  }
  
  bool has_uv = false;
  if(tri->uv != nullptr) {
    has_uv = true;
  }
  int number_faces = tri->numIndices / 3;
  
  vec3f tris[3];
  vec3f normals[3];
  for (int i = 0; i < number_faces; i++) {
    bool tempnormal = false;
    int idx = 3*i;
    tris[0] = vec3f(tri->pos[3*tri->indices[idx  ]+0],
                   tri->pos[3*tri->indices[idx  ]+1],
                           tri->pos[3*tri->indices[idx  ]+2])*scale;
    tris[1] = vec3f(tri->pos[3*tri->indices[idx+1]+0],
                   tri->pos[3*tri->indices[idx+1]+1],
                           tri->pos[3*tri->indices[idx+1]+2])*scale;
    tris[2] = vec3f(tri->pos[3*tri->indices[idx+2]+0],
                   tri->pos[3*tri->indices[idx+2]+1],
                           tri->pos[3*tri->indices[idx+2]+2])*scale;
    if(has_normals) {
      tempnormal = true;
      normals[0] = vec3f(tri->normal[3*tri->indices[idx  ]+0],
                        tri->normal[3*tri->indices[idx  ]+1],
                                   tri->normal[3*tri->indices[idx  ]+2]);
      normals[1] = vec3f(tri->normal[3*tri->indices[idx+1]+0],
                        tri->normal[3*tri->indices[idx+1]+1],
                                   tri->normal[3*tri->indices[idx+1]+2]);
      normals[2] = vec3f(tri->normal[3*tri->indices[idx+2]+0],
                        tri->normal[3*tri->indices[idx+2]+1],
                                   tri->normal[3*tri->indices[idx+2]+2]);
    }
    if((normals[0].x() == 0 && normals[0].y() == 0 && normals[0].z() == 0) ||
       (normals[1].x() == 0 && normals[1].y() == 0 && normals[1].z() == 0) ||
       (normals[2].x() == 0 && normals[2].y() == 0 && normals[2].z() == 0)) {
      tempnormal = false;
    }
    
    if(has_normals && tempnormal) {
      triangles.add(std::make_shared<triangle>(tris[0],tris[1],tris[2],
                                       normals[0],normals[1],normals[2],
                                                                    false,
                                                                    mat_ptr, nullptr,  nullptr, 
                                                                    ObjectToWorld, WorldToObject, reverseOrientation));
    } else {
      triangles.add(std::make_shared<triangle>(tris[0],tris[1],tris[2], false, mat_ptr, nullptr, nullptr, 
                                               ObjectToWorld, WorldToObject, reverseOrientation));
    }
  }
  ply_mesh_bvh = std::make_shared<bvh_node>(triangles, shutteropen, shutterclose, bvh_type, rng);
  delete tri;
};


bool plymesh::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  return(ply_mesh_bvh->hit(r, t_min, t_max, rec, rng));
};

bool plymesh::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
  return(ply_mesh_bvh->hit(r, t_min, t_max, rec, sampler));
};

bool plymesh::bounding_box(Float t0, Float t1, aabb& box) const {
  return(ply_mesh_bvh->bounding_box(t0,t1,box));
};

