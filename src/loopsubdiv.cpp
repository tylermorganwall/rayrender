#include "loopsubdiv.h"
#include <set>
#include <map>
#include "point3.h"
#include "Rcpp.h"
#include <queue>

struct SDFace;
struct SDVertex;

// LoopSubdiv Macros
#define NEXT(i) (((i) + 1) % 3)
#define PREV(i) (((i) + 2) % 3)

// LoopSubdiv Local Structures
struct SDVertex {
  // SDVertex Constructor
  SDVertex(const point3f &p = point3f(0, 0, 0),
           const point2f &uv = point2f(0, 0),
           const point3f &vc = point3f(0, 0, 0)) : 
    p(p), uv(uv), vc(vc) {}
  
  // SDVertex Methods
  int valence();
  void oneRing(point3f *p);
  void oneRing(point2f *p);
  point3f p;
  point2f uv;
  point3f vc;
  SDFace *startFace = nullptr;
  SDVertex *child = nullptr;
  bool regular = false, boundary = false;
  bool initialized = false; //This marks when the vertex has a startFace
};

struct SDFace {
  // SDFace Constructor
  SDFace() {
    for (int i = 0; i < 3; ++i) {
      v[i] = nullptr;
      f[i] = nullptr;
    }
    MatID = 0;
    vertices_initialized = false;
    faces_initialized = false;
    for (int i = 0; i < 4; ++i) children[i] = nullptr;
  }
  
  // SDFace Methods
  int vnum(SDVertex *vert) const {
    for (int i = 0; i < 3; ++i) {
      if(!vert) {
        Rcpp::stop("Vert not found");
      }
      if(!v[i]) {
        Rcpp::stop("Vert not initialized");
      }
      if (v[i] == vert) {
        return i; 
      }
    }
    Rcpp::stop( "Basic logic error in SDFace::vnum()");
    return -1;
  }
  SDFace *nextFace(SDVertex *vert) { return f[vnum(vert)]; }
  SDFace *prevFace(SDVertex *vert) { return f[PREV(vnum(vert))]; }
  SDVertex *nextVert(SDVertex *vert) { return v[NEXT(vnum(vert))]; }
  SDVertex *prevVert(SDVertex *vert) { return v[PREV(vnum(vert))]; }
  SDVertex *otherVert(SDVertex *v0, SDVertex *v1) {
    for (int i = 0; i < 3; ++i) {
      if (v[i] != v0 && v[i] != v1) {
        return v[i];
      }
    }
    Rcpp::stop( "Basic logic error in SDVertex::otherVert()");
    return nullptr;
  }
  SDVertex *v[3];
  SDFace *f[3];
  SDFace *children[4];
  int MatID;
  bool vertices_initialized; //This marks when the vertices are all initialized
  bool faces_initialized; //This marks when the vertices are all initialized
};

struct SDEdge {
  // SDEdge Constructor
  SDEdge(SDVertex *v0 = nullptr, SDVertex *v1 = nullptr) {
    v[0] = std::min(v0, v1);
    v[1] = std::max(v0, v1);
    f[0] = f[1] = nullptr;
    f0edgeNum = -1;
  }
  
  // SDEdge Comparison Function
  bool operator<(const SDEdge &e2) const {
    if (v[0] == e2.v[0]) return v[1] < e2.v[1];
    return v[0] < e2.v[0];
  }
  SDVertex *v[2];
  SDFace *f[2];
  int f0edgeNum;
};

// LoopSubdiv Local Declarations
template <typename T>
static T weightOneRing(SDVertex *vert, T vert_data, Float beta);

template <typename T>
static T weightBoundary(SDVertex *vert, T vert_data, Float beta);

// LoopSubdiv Inline Functions
inline int SDVertex::valence() {
  SDFace *f = startFace;
  if (!boundary) {
    // Compute valence of interior vertex
    int nf = 1;
    while ((f = f->nextFace(this)) != startFace) {
      ++nf;
    }
    return nf;
  } else {
    // Compute valence of boundary vertex
    int nf = 1; 

    while ((f = f->nextFace(this)) != nullptr) {
      ++nf;
    }
    f = startFace;
    while ((f = f->prevFace(this)) != nullptr) {
      ++nf;
    }
    return nf + 1;
  }
}

inline Float beta(int valence) {
  if (valence == 3) {
    return 3.f / 16.f;
  } else {
    return 3.f / (8.f * valence);
  }
}

inline Float loopGamma(int valence) {
  return 1.f / (valence + 3.f / (8.f * beta(valence)));
}


// SDVertex *calculateSharpVertexPosition(SDVertex *v) {
//   point3f sum(0.f, 0.f, 0.f);
//   float sharpWeight = 0.5; // Define how much sharpness influences
//   int count = 0;
//   
//   SDFace* f = v->startFace;
//   while(f)
//   
//   for (auto edge : v->connectedEdges()) {
//     if (edge->sharpness > 0) {
//       sum += edge->otherVertex(v)->p * sharpWeight;
//       count++;
//     }
//   }
//   
//   return (v->p * (1 - sharpWeight * count)) + sum;
// }

// LoopSubdiv Function Definitions
void LoopSubdivide(TriangleMesh* base_mesh,
                   const int nLevels,
                   bool verbose) {
  int nIndices = base_mesh->vertexIndices.size();
  const int *vertexIndices = base_mesh->vertexIndices.data();
  const int *texIndices    = base_mesh->texIndices.data();
  std::vector<std::unique_ptr<SDVertex> > vertex_storage;
  std::vector<std::unique_ptr<SDFace> > face_storage;
  
  //TODO: Need to check if UV indices equal vertex indices
  
  // const int *normalIndices = base_mesh->normalIndices.data();

  int nVertices = base_mesh->nVertices;
  const point3f *p = base_mesh->p.get();
  const point2f *uv = base_mesh->uv.get();
  const point3f *vc = base_mesh->vc.get();
  bool has_uv = base_mesh->has_tex;
  // bool has_normals = base_mesh->has_normals;
  bool has_vc = base_mesh->has_vertex_colors;
  
  std::vector<SDVertex *> vertices;
  std::vector<SDFace *> faces;
  // Allocate _LoopSubdiv_ vertices and faces
  std::vector<point2f> TexcoordsExpanded(nVertices);
  std::vector<int> vertexToTex(nVertices);
  if(has_uv) {
    for(size_t i = 0; i < base_mesh->vertexIndices.size(); ++i) {
      if(texIndices[i] != -1) {
        TexcoordsExpanded[vertexIndices[i]] = uv[texIndices[i]];
      }
    }
  }

  std::unique_ptr<SDVertex[]> verts(new SDVertex[nVertices]);
  for (int i = 0; i < nVertices; ++i) {
    verts[i] = SDVertex(p[i]);
    if(has_uv) {
      verts[i].uv = TexcoordsExpanded[i];
    }
    if(has_vc) {
      verts[i].vc = vc[i];
    }
    vertices.push_back(&verts[i]);
  }

  int nFaces = nIndices / 3;
  std::unique_ptr<SDFace[]> fs(new SDFace[nFaces]);
  for (int i = 0; i < nFaces; ++i) {
    fs[i].MatID = base_mesh->face_material_id[i];
    faces.push_back(&fs[i]);
  }

  // Set face to vertex pointers
  // rayrender addition: This also marks the vertices/faces as initialized
  const int *vp = vertexIndices;
  for (int i = 0; i < nFaces; ++i, vp += 3) {
    SDFace *f = faces[i];
    for (int j = 0; j < 3; ++j) {
      SDVertex *v = vertices[vp[j]];
      f->v[j] = v;
      v->initialized = true;
      v->startFace = f;
    }
    f->vertices_initialized = true;
  }

  // Set neighbor pointers in _faces_
  std::set<SDEdge> edges;
  for (int i = 0; i < nFaces; ++i) {
    SDFace *f = faces[i];
    for (int edgeNum = 0; edgeNum < 3; ++edgeNum) {
      // Update neighbor pointer for _edgeNum_
      int v0 = edgeNum, v1 = NEXT(edgeNum);
      SDEdge e(f->v[v0], f->v[v1]);
      if (edges.find(e) == edges.end()) {
        // Handle new edge
        e.f[0] = f;
        e.f0edgeNum = edgeNum;
        edges.insert(e);
      } else {
        // Handle previously seen edge
        e = *edges.find(e);
        e.f[0]->f[e.f0edgeNum] = f;
        f->f[edgeNum] = e.f[0];
        edges.erase(e);
      }
    }
    f->faces_initialized = true;
  }

  // Finish vertex initialization
  for (int i = 0; i < nVertices; ++i) {
    SDVertex *v = vertices[i];
    SDFace *f = v->startFace;
    // SDFace *prevf = nullptr;
    if(!f || !f->faces_initialized) {
      continue;
    }
    do {
      f = f->nextFace(v);
      Rcpp::checkUserInterrupt(); //Add logic here to auto detect loops due to flipped edges
    } while (f && f != v->startFace);
    v->boundary = (f == nullptr);
    if (!v->boundary && v->valence() == 6) {
      v->regular = true;
    } else if (v->boundary && v->valence() == 4) {
      v->regular = true;
    } else {
      v->regular = false;
    }
    v->initialized = true;
  }

  // Refine _LoopSubdiv_ into triangles
  std::vector<SDFace *> f = faces;
  std::vector<SDVertex *> v = vertices;
  for (int i = 0; i < nLevels; ++i) {
    if(verbose) {
      Rcpp::message(Rcpp::CharacterVector(std::string("* Subdividing mesh level " + 
        std::to_string(i+1) + "/" + std::to_string(nLevels))));
    }
    Rcpp::checkUserInterrupt();
    // Update _f_ and _v_ for next level of subdivision
    std::vector<SDFace *> newFaces;
    std::vector<SDVertex *> newVertices;
    
    // Allocate next level of children in mesh tree
    for (SDVertex *vertex : v) {
      if(vertex->initialized) {
        vertex_storage.push_back(std::make_unique<SDVertex>());
        vertex->child = vertex_storage.back().get();
        vertex->child->regular = vertex->regular;
        vertex->child->boundary = vertex->boundary;
        newVertices.push_back(vertex->child);
      }
    }
    for (SDFace *face : f) {
      for (int k = 0; k < 4; ++k) {
        face_storage.push_back(std::make_unique<SDFace>());
        face->children[k] = face_storage.back().get();
        face->children[k]->MatID = face->MatID;
        newFaces.push_back(face->children[k]);
      }
    }
    
    // Update vertex positions and create new edge vertices
    // Update vertex positions for even vertices
    for (SDVertex *vertex : v) {
      if(!vertex->initialized) {
        continue;
      }
      if (!vertex->boundary) {
        // Apply one-ring rule for even vertex
        if (vertex->regular) {
          vertex->child->p = weightOneRing<point3f>(vertex, vertex->p, 1.f / 16.f);
        } else {
          vertex->child->p = weightOneRing<point3f>(vertex, vertex->p, beta(vertex->valence()));
        }
        if(has_uv) {
          if (vertex->regular) {
            vertex->child->uv = weightOneRing<point2f>(vertex, vertex->uv, 1.f / 16.f);
          } else {
            vertex->child->uv = weightOneRing<point2f>(vertex, vertex->uv, beta(vertex->valence()));
          }
        }
        if(has_vc) {
          if (vertex->regular) {
            vertex->child->vc = weightOneRing<point3f>(vertex, vertex->vc, 1.f / 16.f);
          } else {
            vertex->child->vc = weightOneRing<point3f>(vertex, vertex->vc, beta(vertex->valence()));
          }
        }
      } else {
        
        // Apply boundary rule for even vertex
        vertex->child->p = weightBoundary<point3f>(vertex, vertex->p, 1.f / 8.f);
        if(has_uv) {
          vertex->child->uv = weightBoundary<point2f>(vertex, vertex->uv, 1.f / 8.f);
        }
        if(has_vc) {
          vertex->child->vc = weightBoundary<point3f>(vertex, vertex->vc, 1.f / 8.f);
        }
      }
    }
    
    // Compute new odd edge vertices
    std::map<SDEdge, SDVertex *> edgeVerts;
    for (SDFace *face : f) {
      for (int k = 0; k < 3; ++k) {
        // Compute odd vertex on _k_th edge
        SDEdge edge(face->v[k], face->v[NEXT(k)]);
        SDVertex *vert = edgeVerts[edge];
        if (!vert) {
          // Create and initialize new odd vertex
          vertex_storage.push_back(std::make_unique<SDVertex>());
          vert = vertex_storage.back().get();
          newVertices.push_back(vert);
          vert->regular = true;
          vert->boundary = (face->f[k] == nullptr);
          vert->startFace = face->children[3];
          vert->initialized = true;
          
          // Apply edge rules to compute new vertex position
          if (vert->boundary) {
            vert->p = 0.5f * edge.v[0]->p;
            vert->p += 0.5f * edge.v[1]->p;
            if(has_uv) {
              vert->uv = 0.5f * edge.v[0]->uv;
              vert->uv += 0.5f * edge.v[1]->uv;
            }
            if(has_vc) {
              vert->vc = 0.5f * edge.v[0]->vc;
              vert->vc += 0.5f * edge.v[1]->vc;
            }
          } else {
            vert->p = 3.f / 8.f * edge.v[0]->p;
            vert->p += 3.f / 8.f * edge.v[1]->p;
            vert->p += 1.f / 8.f * face->otherVert(edge.v[0], edge.v[1])->p;
            vert->p += 1.f / 8.f * face->f[k]->otherVert(edge.v[0], edge.v[1])->p;
            if(has_uv) {
              vert->uv = 3.f / 8.f * edge.v[0]->uv;
              vert->uv += 3.f / 8.f * edge.v[1]->uv;
              vert->uv += 1.f / 8.f * face->otherVert(edge.v[0], edge.v[1])->uv;
              vert->uv += 1.f / 8.f * face->f[k]->otherVert(edge.v[0], edge.v[1])->uv;
            }
            if(has_vc) {
              vert->vc = 3.f / 8.f * edge.v[0]->vc;
              vert->vc += 3.f / 8.f * edge.v[1]->vc;
              vert->vc += 1.f / 8.f * face->otherVert(edge.v[0], edge.v[1])->vc;
              vert->vc += 1.f / 8.f * face->f[k]->otherVert(edge.v[0], edge.v[1])->vc;
            }
          }
          edgeVerts[edge] = vert;
        }
      }
    }
    
    // Update new mesh topology
    
    // Update even vertex face pointers
    for (SDVertex *vertex : v) {
      if(vertex->initialized) {
        int vertNum = vertex->startFace->vnum(vertex);
        vertex->child->startFace = vertex->startFace->children[vertNum];
        vertex->child->initialized = true;
      }
    }
    
    // Update face neighbor pointers
    for (SDFace *face : f) {
      for (int j = 0; j < 3; ++j) {
        // Update children _f_ pointers for siblings
        face->children[3]->f[j] = face->children[NEXT(j)];
        face->children[j]->f[NEXT(j)] = face->children[3];
        
        // Update children _f_ pointers for neighbor children
        SDFace *f2 = face->f[j];
        face->children[j]->f[j] =
          f2 ? f2->children[f2->vnum(face->v[j])] : nullptr;
        f2 = face->f[PREV(j)];
        face->children[j]->f[PREV(j)] =
          f2 ? f2->children[f2->vnum(face->v[j])] : nullptr;
        //This helps avoid accessing floating triangles (e.g. which have no adjacent faces)
        if(!face->children[j]->f[j]) {
          continue;
        }
          
        face->children[j]->f[j]->faces_initialized = true;
      }
    }
    
    // Update face vertex pointers
    for (SDFace *face : f) {
      for (int j = 0; j < 3; ++j) {
        // Update child vertex pointer to new even vertex
        face->children[j]->v[j] = face->v[j]->child;

        // Update child vertex pointer to new odd vertex
        SDVertex *vert = edgeVerts[SDEdge(face->v[j], face->v[NEXT(j)])];
        face->children[j]->v[NEXT(j)] = vert;
        face->children[NEXT(j)]->v[j] = vert;
        face->children[3]->v[j] = vert;
        face->vertices_initialized = true;
      }
    }
    
    // Prepare for next level of subdivision
    f = newFaces;
    v = newVertices;
  }
  // Push vertices to limit surface
  std::unique_ptr<point3f[]> final_vertices(new point3f[v.size()]);
  for (size_t i = 0; i < v.size(); ++i) {
    if(v[i]->initialized) {
      if (v[i]->boundary) {
        final_vertices[i] = weightBoundary<point3f>(v[i], v[i]->p, 1.f / 5.f);
      } else {
        final_vertices[i] = weightOneRing<point3f>(v[i], v[i]->p, loopGamma(v[i]->valence()));
      }
    }
  }

  for (size_t i = 0; i < v.size(); ++i) {
    if(v[i]->initialized) {
      v[i]->p = final_vertices[i];
    }
  }

  // Compute UVs on limit surface
  std::unique_ptr<point2f[]> final_texcoords(new point2f[v.size()]);

  if(has_uv) {
    for (size_t i = 0; i < v.size(); ++i) {
      if(v[i]->initialized) {
        if (v[i]->boundary) {
          final_texcoords[i] = weightBoundary<point2f>(v[i], v[i]->uv, 1.f / 5.f);
        } else {
          final_texcoords[i] = weightOneRing<point2f>(v[i], v[i]->uv, loopGamma(v[i]->valence()));
        }
      }
    }
    
    for (size_t i = 0; i < v.size(); ++i) {
      if(v[i]->initialized) {
        v[i]->uv = final_texcoords[i];
      }
    }
  }

  if(has_vc) {
    // Compute UVs on limit surface
    std::unique_ptr<point3f[]> vcLimit(new point3f[v.size()]);
    for (size_t i = 0; i < v.size(); ++i) {
      if(v[i]->initialized) {
        if (v[i]->boundary) {
          vcLimit[i] = weightBoundary<point3f>(v[i], v[i]->vc, 1.f / 5.f);
        } else {
          vcLimit[i] = weightOneRing<point3f>(v[i], v[i]->vc, loopGamma(v[i]->valence()));
        }
      }
    }
    
    for (size_t i = 0; i < v.size(); ++i) {
      if(v[i]->initialized) {
        v[i]->vc = vcLimit[i];
      }
    }
  }

  
  //Visualize the faces as colors related to topology
  // Compute vertex tangents on limit surface
  std::unique_ptr<normal3f[]> final_normals = std::make_unique<normal3f[]>(v.size());
  // Ns.reserve(v.size());
  std::vector<point3f> pRing(16, point3f());
  int cntr = 0;
  for (SDVertex *vertex : v) {
    if(vertex->initialized) {
      vec3f S(0, 0, 0), T(0, 0, 0);
      int valence = vertex->valence();
      if (valence > (int)pRing.size()) {
        pRing.resize(valence);
      }
      vertex->oneRing(&pRing[0]);
      if (!vertex->boundary) {
        // Compute tangents of interior face
        for (int j = 0; j < valence; ++j) {
          S += std::cos(2 * M_PI * j / valence) * vec3f(pRing[j]);
          T += std::sin(2 * M_PI * j / valence) * vec3f(pRing[j]);
        }
      } else {
        // Compute tangents of boundary face
        S = pRing[valence - 1] - pRing[0];
        if (valence == 2) {
          T = vec3f(pRing[0] + pRing[1] - 2 * vertex->p);
        } else if (valence == 3) {
          T = pRing[1] - vertex->p;
        } else if (valence == 4) {   // regular
          T = vec3f(-1 * pRing[0] + 2 * pRing[1] + 2 * pRing[2] +
            -1 * pRing[3] + -2 * vertex->p);
        } else {
          Float theta = M_PI / float(valence - 1);
          T = vec3f(std::sin(theta) * (pRing[0] + pRing[valence - 1]));
          for (int k = 1; k < valence - 1; ++k) {
            Float wt = (2 * std::cos(theta) - 2) * std::sin((k)*theta);
            T += vec3f(wt * pRing[k]);
          }
          T = -T;
        }
      }
      final_normals[cntr] = (normal3f(-cross(S, T)));
    }
    cntr++;
  }
  
  // Create triangle mesh from subdivision mesh
  {
    size_t ntris = f.size();
    std::unique_ptr<int[]> indices(new int[3 * ntris]);
    std::vector<int> face_material_id(ntris);
    int *vp = indices.get();
    size_t totVerts = v.size();
    std::map<SDVertex *, int> usedVerts;
    for (size_t i = 0; i < totVerts; ++i) {
      if(v[i]->initialized) {
        usedVerts[v[i]] = i;
      } 
    }
    
    for (size_t i = 0; i < ntris; ++i) {
      face_material_id[i] = f[i]->MatID;
      for (int j = 0; j < 3; ++j) {
        *vp = usedVerts[f[i]->v[j]];
        ++vp;
      }
    }
    base_mesh->p = std::move(final_vertices);
    base_mesh->n = std::move(final_normals);
    
    if(base_mesh->nTex > 0) {
      base_mesh->has_tex = true;
      base_mesh->uv = std::move(final_texcoords);
    } else {
      uv = nullptr;
    }
    base_mesh->nVertices = v.size(); 
    base_mesh->nNormals = v.size();
    
    base_mesh->has_normals = true;
    base_mesh->nTex = base_mesh->has_tex ? v.size()  : 0;
    
    size_t numIndices = 3 * ntris; 

    base_mesh->vertexIndices.clear();
    base_mesh->normalIndices.clear();
    base_mesh->texIndices.clear();
    
    base_mesh->nTriangles = ntris;
    base_mesh->face_material_id = face_material_id;
    for (size_t s = 0; s < numIndices; s += 3) {
      base_mesh->vertexIndices.push_back(indices[s]);
      base_mesh->vertexIndices.push_back(indices[s+1]);
      base_mesh->vertexIndices.push_back(indices[s+2]);
      // if(base_mesh->has_normals) {
      base_mesh->normalIndices.push_back(indices[s]);
      base_mesh->normalIndices.push_back(indices[s+1]);
      base_mesh->normalIndices.push_back(indices[s+2]);
      // } else {
      //   base_mesh->normalIndices.push_back(-1);
      //   base_mesh->normalIndices.push_back(-1);
      //   base_mesh->normalIndices.push_back(-1);
      // }
      if(base_mesh->has_tex) {
        base_mesh->texIndices.push_back(indices[s+0]);
        base_mesh->texIndices.push_back(indices[s+1]);
        base_mesh->texIndices.push_back(indices[s+2]);
      } else {
        base_mesh->texIndices.push_back(-1);
        base_mesh->texIndices.push_back(-1);
        base_mesh->texIndices.push_back(-1);
      }
    }
    if(base_mesh->has_consistent_normals) {
      base_mesh->alpha_v.clear();
      base_mesh->face_n.reset(new normal3f[base_mesh->normalIndices.size() / 3]);
      std::map<int, std::priority_queue<Float> > alpha_values;
      for (size_t i = 0; i < base_mesh->normalIndices.size(); i += 3) {
        int idx_n1 = base_mesh->normalIndices[i];
        int idx_n2 = base_mesh->normalIndices[i+1];
        int idx_n3 = base_mesh->normalIndices[i+2];

        normal3f n1 = unit_vector(base_mesh->n[idx_n1]);
        normal3f n2 = unit_vector(base_mesh->n[idx_n2]);
        normal3f n3 = unit_vector(base_mesh->n[idx_n3]);
        
        normal3f face_normal = unit_vector(n1 + n2 + n3);
        base_mesh->face_n[i / 3] = face_normal;
        Float av1 = dot(n1,face_normal);
        Float av2 = dot(n2,face_normal);
        Float av3 = dot(n3,face_normal);
        alpha_values[idx_n1].push(-av1);
        alpha_values[idx_n2].push(-av2);
        alpha_values[idx_n3].push(-av3);
      }
      for (auto const& x : alpha_values) {
        base_mesh->alpha_v.push_back(-x.second.top());
      }
      for(size_t i = 0; i < base_mesh->alpha_v.size(); i++) {
        Float temp_av = clamp(base_mesh->alpha_v[i],-1,1);
        base_mesh->alpha_v[i] = std::acos(temp_av) * (1 + 0.03632 * (1 - temp_av) * (1 - temp_av));
      }
    }
  }
}

#define ALLOCA(TYPE, COUNT) (TYPE *) alloca((COUNT) * sizeof(TYPE))

template <typename T>
static T weightOneRing(SDVertex *vert, T vert_data, Float beta) {
  // Put _vert_ one-ring in _pRing_
  int valence = vert->valence();
  T *pRing = ALLOCA(T, valence);
  vert->oneRing(pRing);
  T p = (1 - valence * beta) * vert_data;
  for (int i = 0; i < valence; ++i) {
    p += beta * pRing[i];
  }
  return p;
}

//Need methods for each data type
void SDVertex::oneRing(point3f *p) {
  if (!boundary) {
    // Get one-ring vertices for interior vertex
    SDFace *face = startFace;
    do {
      *p++ = face->nextVert(this)->p;
      face = face->nextFace(this);
    } while (face != startFace);
  } else {
    // Get one-ring vertices for boundary vertex
    SDFace *face = startFace, *f2;
    while ((f2 = face->nextFace(this)) != nullptr) face = f2;
    *p++ = face->nextVert(this)->p;
    do {
      *p++ = face->prevVert(this)->p;
      face = face->prevFace(this);
    } while (face != nullptr);
  }
}

void SDVertex::oneRing(point2f *p) {
  if (!boundary) {
    // Get one-ring vertices for interior vertex
    SDFace *face = startFace;
    do {
      *p++ = face->nextVert(this)->uv;
      face = face->nextFace(this);
    } while (face != startFace);
  } else {
    // Get one-ring vertices for boundary vertex
    SDFace *face = startFace, *f2;
    while ((f2 = face->nextFace(this)) != nullptr) face = f2;
    *p++ = face->nextVert(this)->uv;
    do {
      *p++ = face->prevVert(this)->uv;
      face = face->prevFace(this);
    } while (face != nullptr);
  }
}

template<typename T>
static T weightBoundary(SDVertex *vert, T vert_data, Float beta) {
  // Put _vert_ one-ring in _pRing_
  int valence = vert->valence();
  T *pRing = ALLOCA(T, valence); //allocate on the stack
  vert->oneRing(pRing);
  T newValue = (1 - 2 * beta) * vert_data;
  newValue += beta * pRing[0];
  newValue += beta * pRing[valence - 1];
  return newValue;
}
