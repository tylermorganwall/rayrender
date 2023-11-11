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
    Rprintf("Not valid PLY reader for file: '%s' \n",filename);
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
      bool polys = reader.requires_triangulation(indexes[0]);
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
                 std::shared_ptr<alpha_texture> alpha, std::shared_ptr<bump_texture> bump,
                 Float scale, Float shutteropen, Float shutterclose, int bvh_type, random_gen rng,
                 std::shared_ptr<Transform> ObjectToWorld, std::shared_ptr<Transform> WorldToObject, bool reverseOrientation) :
  hitable(ObjectToWorld, WorldToObject, reverseOrientation) {
  TriMesh* tri = parse_file_with_miniply(inputfile.c_str(), false);

  if(tri == nullptr) {
    std::string err = inputfile;
    throw std::runtime_error("No mesh loaded: " + err);
  }

  mesh = std::unique_ptr<TriangleMesh>(new TriangleMesh(tri->pos, tri->indices, tri->normal, tri->uv, 
                                                        tri->numVerts, tri->numIndices,
                                                        alpha, bump, 
                                                        mat,
                                                        ObjectToWorld, WorldToObject, reverseOrientation));
  // mesh->ValidateMesh();
  size_t n = mesh->nTriangles * 3;
  
  for(size_t i = 0; i < n; i += 3) {
    triangles.add(std::make_shared<triangle>(mesh.get(), 
                                             &mesh->vertexIndices[i], 
                                             &mesh->normalIndices[i],
                                             &mesh->texIndices[i], i / 3,
                                             ObjectToWorld, WorldToObject, reverseOrientation));
  }
  ply_mesh_bvh = std::make_shared<bvh_node>(triangles, shutteropen, shutterclose, bvh_type, rng);
  // ply_mesh_bvh->validate_bvh();
  triangles.objects.clear();
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

std::pair<size_t,size_t> plymesh::CountNodeLeaf() {
  return(ply_mesh_bvh->CountNodeLeaf());
}
