#define TINYPLY_IMPLEMENTATION

#include "../ext/miniply/miniply.h"
#include "../hitables/plymesh.h"
#include "../math/loopsubdiv.h"
#include "../math/calcnormals.h"
#include "../utils/raylog.h"

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
  std::vector<int> ptex_face_indices;
  
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
      if (reader.find_normal(indexes)) {
        trimesh->normal = new float[trimesh->numVerts * 3];
        reader.extract_properties(indexes, 3, miniply::PLYPropertyType::Float, trimesh->normal);
      }
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
      const uint32_t face_property = reader.find_property("face_indices");
      if (face_property != miniply::kInvalidIndex) {
        std::vector<int> source_faces(reader.num_rows());
        if (!reader.extract_properties(&face_property, 1, miniply::PLYPropertyType::Int, source_faces.data())) {
          delete trimesh;
          Rcpp::stop("Cannot read PLY face_indices.");
        }
        const uint32_t* counts = reader.get_list_counts(indexes[0]);
        trimesh->ptex_face_indices.reserve(trimesh->numIndices / 3);
        for (size_t i = 0; i < source_faces.size(); ++i) {
          if (source_faces[i] < 0 || counts[i] < 3) {
            delete trimesh;
            Rcpp::stop("Invalid PLY face_indices or polygon size.");
          }
          // miniply emits each source polygon's n-2 triangles contiguously.
          trimesh->ptex_face_indices.insert(trimesh->ptex_face_indices.end(), counts[i] - 2, source_faces[i]);
        }
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
                 Float scale, int subdivision_levels, bool recalculate_normals,
                 bool verbose,
                 Float shutteropen, Float shutterclose, int bvh_type, random_gen rng,
                 Transform* ObjectToWorld, Transform* WorldToObject, bool reverseOrientation,
                 bool calculate_consistent_normals) :
  hitable(ObjectToWorld, WorldToObject, mat, reverseOrientation) {
  std::unique_ptr<TriMesh> tri(parse_file_with_miniply(inputfile.c_str(), false));

  if(tri == nullptr) {
    std::string err = inputfile;
    throw std::runtime_error("No mesh loaded: " + err);
  }

  mesh = std::unique_ptr<TriangleMesh>(new TriangleMesh(tri->pos, tri->indices, tri->normal, tri->uv, 
                                                        tri->numVerts, tri->numIndices,
                                                        alpha, bump, 
                                                        mat,
                                                        ObjectToWorld, WorldToObject, reverseOrientation,
                                                        calculate_consistent_normals));

  mesh->ptex_face_indices = std::move(tri->ptex_face_indices);
  // TriangleMesh copied the parser arrays; release them before refinement/BVH.
  tri.reset();

  //Loop subdivision automatically calculates new normals
  if(subdivision_levels > 1) {
    LoopSubdivide(mesh.get(),
                  subdivision_levels,
                  verbose);
  } else if (recalculate_normals) {
    CalculateNormals(mesh.get());
  }
  
  size_t n = mesh->nTriangles * 3;
  triangles.objects.reserve(mesh->nTriangles);
  
  // mesh->ValidateMesh();
  
  for(size_t i = 0; i < n; i += 3) {
    triangles.add(std::make_shared<triangle>(mesh.get(), 
                                             &mesh->vertexIndices[i], 
                                             &mesh->normalIndices[i],
                                             &mesh->texIndices[i], i / 3,
                                             ObjectToWorld, WorldToObject, reverseOrientation));
  }
  ply_mesh_bvh = std::make_shared<BVHAggregate>(std::move(triangles.objects), shutteropen, shutterclose, bvh_type, true);
  // ply_mesh_bvh->validate_bvh();
  // Transfer ownership to the BVH; retain no empty construction buffer.
  std::vector<std::shared_ptr<hitable>>().swap(triangles.objects);
};


const bool plymesh::hit(const Ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const {
  SCOPED_CONTEXT("MultiHit");
  SCOPED_TIMER_COUNTER("PlyMesh");
  return(ply_mesh_bvh->hit(r, t_min, t_max, rec, rng));
};

const bool plymesh::hit(const Ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const {
  SCOPED_CONTEXT("MultiHit");
  SCOPED_TIMER_COUNTER("PlyMesh");
  
  return(ply_mesh_bvh->hit(r, t_min, t_max, rec, sampler));
};

bool plymesh::HitP(const Ray& r, Float t_min, Float t_max,random_gen& rng) const {
  SCOPED_CONTEXT("MultiHit");
  SCOPED_TIMER_COUNTER("PlyMesh");
  return(ply_mesh_bvh->HitP(r, t_min, t_max, rng));
};

bool plymesh::HitP(const Ray& r, Float t_min, Float t_max, Sampler* sampler) const {
  SCOPED_CONTEXT("MultiHit");
  SCOPED_TIMER_COUNTER("PlyMesh");
  return(ply_mesh_bvh->HitP(r, t_min, t_max, sampler));
};

bool plymesh::bounding_box(Float t0, Float t1, aabb& box) const {
  return(ply_mesh_bvh->bounding_box(t0,t1,box));
};

std::pair<size_t,size_t> plymesh::CountNodeLeaf() {
  return(ply_mesh_bvh->CountNodeLeaf());
}
