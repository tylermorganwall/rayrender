#include "trimesh.h"
#include "RProgress.h"
#include "bvh.h"
#include "loopsubdiv.h"
#include "displacement.h"
#include "calctangents.h"
#include "assert.h"
#include "calcnormals.h"
#include "raylog.h"

trimesh::trimesh(std::string inputfile, std::string basedir, Float scale, Float sigma,
                 std::shared_ptr<material> default_material, 
                 std::shared_ptr<alpha_texture> alpha_default, std::shared_ptr<bump_texture> bump_default,
                 bool load_materials, 
                 bool load_textures, bool load_vertex_colors,
                 bool importance_sample_lights, bool load_normals, bool calculate_consistent_normals,
                 int subdivision_levels, std::string displacement_texture, Float displacement,
                 bool displacement_vector, TextureCache& texCache, bool recalculate_normals,
                 hitable_list& imp_sample_objects, 
                 Float shutteropen, Float shutterclose, int bvh_type, random_gen rng, bool verbose,
                 Transform* ObjectToWorld, Transform* WorldToObject, bool reverseOrientation) : 
  hitable(ObjectToWorld, WorldToObject, default_material, reverseOrientation) {
  mesh = std::unique_ptr<TriangleMesh>(new TriangleMesh(inputfile, basedir, default_material, 
                                                        alpha_default, bump_default,
                                                        load_materials, load_textures, load_vertex_colors, load_normals, verbose,
                                                        scale, calculate_consistent_normals, texCache,
                                                        ObjectToWorld, WorldToObject, reverseOrientation));
  //Loop subdivision automatically calculates new normals
  if(subdivision_levels > 1) {
    LoopSubdivide(mesh.get(),
                  subdivision_levels,
                  verbose);
  } else if (recalculate_normals) {
    CalculateNormals(mesh.get());
  }
  if(displacement_texture.length() > 0) {
    if(mesh->nVertices != mesh->nNormals) {
      if(verbose) {
        Rcpp::message(Rcpp::CharacterVector("* Calculating mesh normals for displacement"));
      }
      CalculateNormals(mesh.get());
    }
    if(displacement_vector) {
      if(verbose) {
        Rcpp::message(Rcpp::CharacterVector("* Calculating mesh tangents for vector displacement"));
      }
      CalculateTangents(mesh.get());
    }
    if(verbose) {
      Rcpp::message(Rcpp::CharacterVector("* Displacing mesh"));
    }
    DisplaceMesh(mesh.get(),
                 displacement_texture,
                 displacement,
                 displacement_vector);
    // Calculate new normals
    if(recalculate_normals) {
      CalculateNormals(mesh.get());
    }
  }

  // mesh->ValidateMesh();
  size_t n = mesh->nTriangles * 3;
  for(size_t i = 0; i < n; i += 3) {
    triangles.add(std::make_shared<triangle>(mesh.get(), 
                                             &mesh->vertexIndices[i], 
                                             &mesh->normalIndices[i],
                                             &mesh->texIndices[i], i / 3,
                                             ObjectToWorld, WorldToObject, reverseOrientation));
   if(mesh->material_is_light[mesh->face_material_id[i / 3]] == 1 && importance_sample_lights) {
      imp_sample_objects.add(triangles.back());
    }
  }
  if(n > 0) {
    tri_mesh_bvh = std::make_shared<BVHAggregate>(triangles.objects, shutteropen, shutterclose, bvh_type, true);
    triangles.objects.clear();
  } else {
    throw std::runtime_error(inputfile + ": No triangles loaded.");
  }
}

const bool trimesh::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const {
  SCOPED_CONTEXT("MultiHit");
  SCOPED_TIMER_COUNTER("ObjMesh");
  
  return(tri_mesh_bvh->hit(r, t_min, t_max, rec, rng));
}

const bool trimesh::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const {
  SCOPED_CONTEXT("MultiHit");
  SCOPED_TIMER_COUNTER("ObjMesh");
  
  return(tri_mesh_bvh->hit(r, t_min, t_max, rec, sampler));
}

bool trimesh::HitP(const ray& r, Float t_min, Float t_max, random_gen& rng) const {
  SCOPED_CONTEXT("MultiHit");
  SCOPED_TIMER_COUNTER("ObjMesh");
  
  return(tri_mesh_bvh->HitP(r, t_min, t_max, rng));
}

bool trimesh::HitP(const ray& r, Float t_min, Float t_max, Sampler* sampler) const {
  SCOPED_CONTEXT("MultiHit");
  SCOPED_TIMER_COUNTER("ObjMesh");
  
  return(tri_mesh_bvh->HitP(r, t_min, t_max, sampler));
}

bool trimesh::bounding_box(Float t0, Float t1, aabb& box) const {
  return(tri_mesh_bvh->bounding_box(t0,t1,box));
}


Float trimesh::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) {
  return(triangles.pdf_value(o,v, rng, time));
}

Float trimesh::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) {
  return(triangles.pdf_value(o,v, sampler, time));
  
}

vec3f trimesh::random(const point3f& o, random_gen& rng, Float time) {
  return(triangles.random(o, rng, time));
}

vec3f trimesh::random(const point3f& o, Sampler* sampler, Float time) {
  return(triangles.random(o, sampler, time));
  
}

size_t trimesh::GetSize() {
  size_t total_size = tri_mesh_bvh->GetSize() + sizeof(*this);
  total_size += mesh->GetSize();
  return(total_size);
}

std::pair<size_t,size_t> trimesh::CountNodeLeaf() {
  return(tri_mesh_bvh->CountNodeLeaf());
}
