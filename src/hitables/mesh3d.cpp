#include "../hitables/mesh3d.h"
#include "../materials/texture.h"
#include "../math/loopsubdiv.h"
#include "../math/displacement.h"
#include "../math/calcnormals.h"
#include "../math/calctangents.h"
#ifndef STBIMAGEH
#define STBIMAGEH
#include "../ext/stb/stb_image.h"
#endif
#include "../materials/texturecache.h"

mesh3d::mesh3d(Rcpp::List mesh_info, std::shared_ptr<material> mat, 
               std::string displacement_texture, Float displacement, bool displacement_vector, 
               TextureCache &texCache, bool recalculate_normals,
               bool verbose, 
               Float shutteropen, Float shutterclose, const BVHBuildOptions& bvh_type, random_gen rng,
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
  std::shared_ptr<const HeightImage> bump_texture_data;
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
    bump_texture_data = texCache.LookupHeight(bump_text_location, nxb, nyb, nnb);
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
  if(subdivision_levels > 0) {
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
  triangles.objects.reserve(mesh->nTriangles);
  // mesh->ValidateMesh();
  
  mesh->texture_size += sizeof(unsigned char) * (nx * ny * nn + nxb * nyb * nnb) ;
  
  for(size_t i = 0; i < n; i += 3) {
    triangles.add(std::make_shared<triangle>(mesh.get(), 
                                             &mesh->vertexIndices[i], 
                                             &mesh->normalIndices[i],
                                             &mesh->texIndices[i], i / 3,
                                             ObjectToWorld, WorldToObject, reverseOrientation));
  }
    
  // All triangle geometry/material tables are loaded; these queries only read them.
  auto mesh_build = bvh_type;
  mesh_build.parallel_primitive_queries = true;
  mesh_bvh = std::make_shared<BVHAggregate>(std::move(triangles.objects), shutteropen, shutterclose, 1, true, mesh_build);
  // Transfer ownership to the BVH; retain no empty construction buffer.
  std::vector<std::shared_ptr<hitable>>().swap(triangles.objects);
}

const bool mesh3d::hit(const Ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const {
  return(mesh_bvh->hit(r, t_min, t_max, rec, rng));
};


const bool mesh3d::hit(const Ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const {
  return(mesh_bvh->hit(r, t_min, t_max, rec, sampler));
};

bool mesh3d::HitP(const Ray& r, Float t_min, Float t_max, random_gen& rng) const {
  return(mesh_bvh->HitP(r, t_min, t_max, rng));
};

bool mesh3d::HitP(const Ray& r, Float t_min, Float t_max, Sampler* sampler) const {
  return(mesh_bvh->HitP(r, t_min, t_max, sampler));
};

bool mesh3d::bounding_box(Float t0, Float t1, aabb& box) const {
  return(mesh_bvh->bounding_box(t0,t1,box));
};

std::pair<size_t,size_t> mesh3d::CountNodeLeaf() {
  return(mesh_bvh->CountNodeLeaf());
}

void mesh3d::BindDeformation(std::shared_ptr<const animation::ModelInstance> instance,
                            size_t node, size_t primitive) {
  const auto &source=*instance->asset->nodes.at(node).meshes.at(primitive);
  if (source.positions.size()!=mesh->nVertices || (source.recalculate_normals ? source.positions.size() : source.normals.size())!=mesh->nNormals ||
      source.indices!=mesh->vertexIndices)
    throw std::runtime_error("Deforming meshes cannot change source topology or normal indexing.");
  deformation=std::move(instance);
  deformation_node=node; deformation_primitive=primitive;
}

std::function<void()> mesh3d::PrepareDeformation(double time) {
  auto pose=deformation->At(time);
  if (last_pose && !animation::ComparePoses(*last_pose,pose).pose) return {};
  auto surface=deformation->Deform(pose,deformation_node,deformation_primitive);
  // Mesh3d stores its existing static object/group transform in the vertices;
  // model node globals and outer placement remain in their instance wrappers.
  for (auto &position : surface.positions) position=(*ObjectToWorld)(position);
  for (auto &normal : surface.normals) normal=unit_vector((*ObjectToWorld)(normal));
  for (size_t i=0;i<surface.tangents.size();++i) {
    auto tangent=(*ObjectToWorld)(surface.tangents[i]);
    const auto &n=surface.normals[i];
    vec3f normal(n[0],n[1],n[2]);
    surface.tangents[i]=unit_vector(tangent-dot(tangent,normal)*normal);
    if (ObjectToWorld->SwapsHandedness()) surface.handedness[i]*=-1;
  }
  return [this,pose=std::move(pose),surface=std::move(surface)]() mutable {
    for (size_t i=0;i<surface.positions.size();++i) mesh->p[i]=surface.positions[i];
    for (size_t i=0;i<surface.normals.size();++i) mesh->n[i]=surface.normals[i];
    if (!surface.tangents.empty()) {
      if (!mesh->t) mesh->t=std::make_unique<vec3f[]>(surface.tangents.size());
      mesh->nTangents=surface.tangents.size(); mesh->has_tangents=true;
      mesh->tangent_right_handed.resize(mesh->nTangents);
      for (size_t i=0;i<surface.tangents.size();++i) {
        mesh->t[i]=surface.tangents[i];
        mesh->tangent_right_handed[i]=surface.handedness[i]>0;
      }
    }
    if (mesh->has_consistent_normals) mesh->CalculateConsistentNormals();
    last_pose=std::move(pose);
    ++deformation_revision;
  };
}

namespace {
// BVH random() chooses triangles uniformly, then samples each triangle's area.
// Convert that same mixture to solid-angle density using current hit geometry.
// Keep standalone triangle sampling behavior outside this mesh path unchanged.
template<class SamplerType>
Float mesh_area_pdf(const BVHAggregate &bvh, const point3f &origin,
                    const vec3f &direction, SamplerType &&sampler, Float time) {
  const auto primitives = bvh.Primitives();
  const Float length = direction.length();
  if (primitives.empty() || !(length > 0)) return 0;
  const Ray ray(origin,direction,time);
  Float sum = 0;
  for (const auto &primitive : primitives) {
    const auto &face = static_cast<const triangle&>(*primitive);
    hit_record hit;
    if (!face.hit(ray,.001,FLT_MAX,hit,sampler)) continue;
    const auto &positions = face.mesh->p;
    const vec3f area_normal = cross(positions[face.v[1]]-positions[face.v[0]],
                                   positions[face.v[2]]-positions[face.v[0]]);
    const Float projected_twice_area = std::abs(dot(area_normal,direction));
    if (projected_twice_area > 0)
      sum += 2*hit.t*hit.t*length*length*length/projected_twice_area;
  }
  return sum/Float(primitives.size());
}
}

Float mesh3d::pdf_value(const point3f &o, const vec3f &v, random_gen &rng, Float time) {
  return mesh_area_pdf(*mesh_bvh,o,v,rng,time);
}
Float mesh3d::pdf_value(const point3f &o, const vec3f &v, Sampler *sampler, Float time) {
  return mesh_area_pdf(*mesh_bvh,o,v,sampler,time);
}
