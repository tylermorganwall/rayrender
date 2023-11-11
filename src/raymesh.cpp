#include "raymesh.h"
#include "RProgress.h"

raymesh::raymesh(Rcpp::List raymesh_list, 
                 std::shared_ptr<material> default_material, 
                 std::shared_ptr<alpha_texture> alpha_mask, std::shared_ptr<bump_texture> bump_tex,
                 bool importance_sample_lights, 
                 bool calculate_consistent_normals,
                 bool override_material,
                 bool flip_transmittance,
                 hitable_list& imp_sample_objects, 
                 bool verbose, 
                 Float shutteropen, Float shutterclose, int bvh_type, random_gen rng, 
                 std::shared_ptr<Transform> ObjectToWorld, std::shared_ptr<Transform> WorldToObject, bool reverseOrientation) : 
  hitable(ObjectToWorld, WorldToObject, reverseOrientation) {
  mesh = std::unique_ptr<TriangleMesh>(new TriangleMesh(raymesh_list, verbose, 
                                                        calculate_consistent_normals, override_material,
                                                        flip_transmittance,
                                                        alpha_mask,
                                                        bump_tex,
                                                        default_material, 
                                                        ObjectToWorld, WorldToObject, reverseOrientation));
#ifdef FULL_DEBUG
  mesh->ValidateMesh();
#endif
  size_t n = mesh->nTriangles;
  for(size_t i = 0; i < 3*n; i += 3) {
    triangles.add(std::make_shared<triangle>(mesh.get(), 
                                             &mesh->vertexIndices[i], 
                                             &mesh->normalIndices[i],
                                             &mesh->texIndices[i], i / 3,
                                             ObjectToWorld, WorldToObject, reverseOrientation));
    if(mesh->face_material_id[i / 3] < 0 || mesh->face_material_id[i / 3] >= mesh->mesh_materials.size()) {
      throw std::runtime_error("Material ID out of range");
    }
    if(mesh->material_is_light[mesh->face_material_id[i / 3]] && importance_sample_lights) {
      imp_sample_objects.add(triangles.back());
    }
  }
  if(n > 0) {
    tri_mesh_bvh = std::make_shared<bvh_node>(triangles, shutteropen, shutterclose, bvh_type, rng);
#ifdef FULL_DEBUG
    tri_mesh_bvh->validate_bvh();
#endif
    triangles.objects.clear();
  } else {
    throw std::runtime_error("raymesh object not loaded (no triangles)");
  }
}

bool raymesh::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  return(tri_mesh_bvh->hit(r, t_min, t_max, rec, rng));
}

bool raymesh::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
  return(tri_mesh_bvh->hit(r, t_min, t_max, rec, sampler));
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
