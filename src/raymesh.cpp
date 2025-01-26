#include "raymesh.h"
#include "RProgress.h"
#include "loopsubdiv.h"
#include "displacement.h"
#include "calcnormals.h"
#include "calctangents.h"
#ifndef STBIMAGEH
#define STBIMAGEH
#include "stb/stb_image.h"
#endif
#include "texturecache.h"
#include "trianglemesh.h"

#include "bvh.h"

#include "rng.h"
#include "triangle.h"
#include "raylog.h"

raymesh::raymesh(Rcpp::List raymesh_list, 
                 std::shared_ptr<material> default_material, 
                 std::shared_ptr<alpha_texture> alpha_mask, 
                 std::shared_ptr<bump_texture> bump_tex,
                 bool importance_sample_lights, 
                 bool calculate_consistent_normals,
                 bool override_material,
                 bool flip_transmittance,
                 int subdivision_levels,
                 std::string displacement_texture, Float displacement, bool displacement_vector,
                 TextureCache &texCache, bool recalculate_normals,
                 hitable_list& imp_sample_objects, 
                 bool verbose, 
                 Float shutteropen, Float shutterclose, int bvh_type, random_gen rng, 
                 Transform* ObjectToWorld, Transform* WorldToObject, bool reverseOrientation) : 
  hitable(ObjectToWorld, WorldToObject, default_material, reverseOrientation) {
  mesh = std::unique_ptr<TriangleMesh>(new TriangleMesh(raymesh_list, verbose, 
                                                        calculate_consistent_normals, override_material,
                                                        flip_transmittance,
                                                        alpha_mask,
                                                        bump_tex,
                                                        texCache,
                                                        default_material, 
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
    DisplaceMesh(mesh.get(),
                 displacement_texture,
                 displacement,
                 displacement_vector);
  }
  size_t n = mesh->nTriangles * 3;
  
#ifdef FULL_DEBUG
  mesh->ValidateMesh();
#endif
  for(size_t i = 0; i < n; i += 3) {
    triangles.add(std::make_shared<triangle>(mesh.get(), 
                                             &mesh->vertexIndices[i], 
                                             &mesh->normalIndices[i],
                                             &mesh->texIndices[i], i / 3,
                                             ObjectToWorld, WorldToObject, reverseOrientation));
    if(mesh->face_material_id[i / 3] < 0 || mesh->face_material_id[i / 3] >= (int)mesh->mesh_materials.size()) {
      throw std::runtime_error("Material ID out of range");
    }
    if(mesh->material_is_light[mesh->face_material_id[i / 3]] && importance_sample_lights) {
      imp_sample_objects.add(triangles.back());
    }
  }
  if(n > 0) {
    tri_mesh_bvh = std::make_shared<BVHAggregate>(triangles.objects, shutteropen, shutterclose, bvh_type, true);
#ifdef FULL_DEBUG
    tri_mesh_bvh->validate_bvh();
#endif
    triangles.objects.clear();
  } else {
    throw std::runtime_error("raymesh object not loaded (no triangles)");
  }
}

const bool raymesh::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const {
  SCOPED_CONTEXT("MultiHit");
  SCOPED_TIMER_COUNTER("RayMesh");
  
  return(tri_mesh_bvh->hit(r, t_min, t_max, rec, rng));
}

const bool raymesh::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const {
  SCOPED_CONTEXT("MultiHit");
  SCOPED_TIMER_COUNTER("RayMesh");
  
  return(tri_mesh_bvh->hit(r, t_min, t_max, rec, sampler));
}

bool raymesh::HitP(const ray& r, Float t_min, Float t_max, random_gen& rng) const {
  SCOPED_CONTEXT("MultiHit");
  SCOPED_TIMER_COUNTER("RayMesh");
  
  return(tri_mesh_bvh->HitP(r, t_min, t_max, rng));
}

bool raymesh::HitP(const ray& r, Float t_min, Float t_max, Sampler* sampler) const {
  SCOPED_CONTEXT("MultiHit");
  SCOPED_TIMER_COUNTER("RayMesh");
  
  return(tri_mesh_bvh->HitP(r, t_min, t_max, sampler));
}

bool raymesh::bounding_box(Float t0, Float t1, aabb& box) const {
  return(tri_mesh_bvh->bounding_box(t0,t1,box));
}


Float raymesh::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) {
  return(triangles.pdf_value(o,v, rng, time));
}

Float raymesh::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) {
  return(triangles.pdf_value(o,v, sampler, time));
  
}

vec3f raymesh::random(const point3f& o, random_gen& rng, Float time) {
  return(triangles.random(o, rng, time));
}

vec3f raymesh::random(const point3f& o, Sampler* sampler, Float time) {
  return(triangles.random(o, sampler, time));
  
}

size_t raymesh::GetSize() {
  size_t total_size = tri_mesh_bvh->GetSize() + sizeof(*this);
  total_size += mesh->GetSize();
  return(total_size);
}

std::pair<size_t,size_t> raymesh::CountNodeLeaf() {
  return(tri_mesh_bvh->CountNodeLeaf());
}
