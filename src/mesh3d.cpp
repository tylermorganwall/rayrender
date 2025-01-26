#include "mesh3d.h"
#include "texture.h"
#include "loopsubdiv.h"
#include "displacement.h"
#include "calcnormals.h"
#include "calctangents.h"
#ifndef STBIMAGEH
#define STBIMAGEH
#include "stb/stb_image.h"
#endif
#include "texturecache.h"

mesh3d::mesh3d(Rcpp::List mesh_info, std::shared_ptr<material> mat, 
               std::string displacement_texture, Float displacement, bool displacement_vector, 
               TextureCache &texCache, bool recalculate_normals,
               bool verbose, 
               Float shutteropen, Float shutterclose, int bvh_type, random_gen rng,
               Transform* ObjectToWorld, Transform* WorldToObject, bool reverseOrientation) :
  hitable(ObjectToWorld, WorldToObject, mat, reverseOrientation) {
  Rcpp::NumericMatrix vertices = Rcpp::as<Rcpp::NumericMatrix>(mesh_info["vertices"]);
  Rcpp::IntegerMatrix indices = Rcpp::as<Rcpp::IntegerMatrix>(mesh_info["indices"]);
  Rcpp::NumericMatrix norms = Rcpp::as<Rcpp::NumericMatrix>(mesh_info["normals"]);
  Rcpp::NumericMatrix txcoord = Rcpp::as<Rcpp::NumericMatrix>(mesh_info["texcoords"]);
  int subdivision_levels = Rcpp::as<int>(mesh_info["subdivision_levels"]);
  
  // float scale_mesh = Rcpp::as<float>(mesh_info["scale_mesh"]);
  
  
  std::string texture_location = Rcpp::as<std::string>(mesh_info["texture"]);
  std::string bump_text_location = Rcpp::as<std::string>(mesh_info["bump_texture"]);
  Float bump_intensity = Rcpp::as<Float>(mesh_info["bump_intensity"]);
  // int material_type = Rcpp::as<Float>(mesh_info["material_type"]);
  
  Rcpp::NumericMatrix colors = Rcpp::as<Rcpp::NumericMatrix>(mesh_info["color_vals"]);
  // int colortype = Rcpp::as<int>(mesh_info["color_type"]);
  
  int nx = 0, ny = 0,nn = 0;
  // bool has_texture = false;
  bool has_bump = false;
  unsigned char* mesh_material_data;
  unsigned char* bump_texture_data;
  if(strlen(texture_location.c_str()) > 0) {
    mesh_material_data = texCache.LookupChar(texture_location, nx, ny, nn, 4);
    nn = 4;
    // has_texture = nx != 0 && ny != 0 && nn != 0;
  } else {
    mesh_material_data = nullptr;
  }
  bool has_alpha = false;
  if(nn == 4) {
    for(int j = 0; j < nx - 1; j++) {
      for(int k = 0; k < ny - 1; k++) {
        if(mesh_material_data[4*j + 4*nx*k + 3] != 255) {
          has_alpha = true;
          break;
        }
      }
      if(has_alpha) {
        break;
      }
    }
  } 

  int nxb = 0, nyb = 0, nnb = 0;
  
  if(strlen(bump_text_location.c_str()) > 0) {
    bump_texture_data = texCache.LookupChar(bump_text_location, nxb, nyb, nnb, 4);
    has_bump = nxb != 0 && nyb != 0 && nnb != 0;
  } else {
    bump_texture_data = nullptr;
  }
  
  std::shared_ptr<alpha_texture> alpha = nullptr;
  std::shared_ptr<bump_texture> bump = nullptr;
  
  if(has_alpha) {
    alpha = std::make_shared<alpha_texture>(mesh_material_data, 
                                            nx, ny, nn);
  } 
  if(has_bump) {
    bump = std::make_shared<bump_texture>(bump_texture_data,
                                          nxb, nyb, nnb,
                                          bump_intensity);
  } 
  
  
  mesh = std::unique_ptr<TriangleMesh>(new TriangleMesh(vertices, indices, norms, txcoord, colors,
                                                        mesh_material_data, bump_texture_data,
                                                        alpha, bump,
                                                        mat, true, true, 
                                                        texCache,
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
    CalculateNormals(mesh.get());
  }
  size_t n = mesh->nTriangles * 3;
  // mesh->ValidateMesh();
  
  mesh->texture_size += sizeof(unsigned char) * (nx * ny * nn + nxb * nyb * nnb) ;
  
  for(size_t i = 0; i < n; i += 3) {
    triangles.add(std::make_shared<triangle>(mesh.get(), 
                                             &mesh->vertexIndices[i], 
                                             &mesh->normalIndices[i],
                                             &mesh->texIndices[i], i / 3,
                                             ObjectToWorld, WorldToObject, reverseOrientation));
  }
    
  mesh_bvh = std::make_shared<BVHAggregate>(triangles.objects, shutteropen, shutterclose, bvh_type, true);
  triangles.objects.clear();
}

const bool mesh3d::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const {
  return(mesh_bvh->hit(r, t_min, t_max, rec, rng));
};


const bool mesh3d::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const {
  return(mesh_bvh->hit(r, t_min, t_max, rec, sampler));
};

bool mesh3d::HitP(const ray& r, Float t_min, Float t_max, random_gen& rng) const {
  return(mesh_bvh->HitP(r, t_min, t_max, rng));
};

bool mesh3d::HitP(const ray& r, Float t_min, Float t_max, Sampler* sampler) const {
  return(mesh_bvh->HitP(r, t_min, t_max, sampler));
};

bool mesh3d::bounding_box(Float t0, Float t1, aabb& box) const {
  return(mesh_bvh->bounding_box(t0,t1,box));
};

std::pair<size_t,size_t> mesh3d::CountNodeLeaf() {
  return(mesh_bvh->CountNodeLeaf());
}
