#include "mesh3d.h"
#include "texture.h"

mesh3d::mesh3d(Rcpp::List mesh_info, std::shared_ptr<material> mat, 
       Float shutteropen, Float shutterclose, int bvh_type, random_gen rng,
       std::shared_ptr<Transform> ObjectToWorld, std::shared_ptr<Transform> WorldToObject, bool reverseOrientation,
       int prop_len, Rcpp::NumericVector tempvector, Rcpp::NumericVector temp_glossy, double sigma,
       double lightintensity) :
  hitable(ObjectToWorld, WorldToObject, reverseOrientation) {
  Rcpp::NumericMatrix vertices = Rcpp::as<Rcpp::NumericMatrix>(mesh_info["vertices"]);
  Rcpp::IntegerMatrix indices = Rcpp::as<Rcpp::IntegerMatrix>(mesh_info["indices"]);
  Rcpp::NumericMatrix norms = Rcpp::as<Rcpp::NumericMatrix>(mesh_info["normals"]);
  Rcpp::NumericMatrix txcoord = Rcpp::as<Rcpp::NumericMatrix>(mesh_info["texcoords"]);
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
    mesh_material_data = stbi_load(texture_location.c_str(), &nx, &ny, &nn, 0);
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
    bump_texture_data = stbi_load(bump_text_location.c_str(), &nxb, &nyb, &nnb, 0);
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
                                                        ObjectToWorld, WorldToObject, reverseOrientation));
  
  mesh->texture_size += sizeof(unsigned char) * (nx * ny * nn + nxb * nyb * nnb) ;
  
  size_t n = mesh->nTriangles * 3;
  for(size_t i = 0; i < n; i += 3) {
    triangles.add(std::make_shared<triangle>(mesh.get(), 
                                             &mesh->vertexIndices[i], 
                                             &mesh->normalIndices[i],
                                             &mesh->texIndices[i], i / 3,
                                             ObjectToWorld, WorldToObject, reverseOrientation));
  }
    
  mesh_bvh = std::make_shared<bvh_node>(triangles, shutteropen, shutterclose, bvh_type, rng);
  triangles.objects.clear();
}

bool mesh3d::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  return(mesh_bvh->hit(r, t_min, t_max, rec, rng));
};


bool mesh3d::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
  return(mesh_bvh->hit(r, t_min, t_max, rec, sampler));
};

bool mesh3d::bounding_box(Float t0, Float t1, aabb& box) const {
  return(mesh_bvh->bounding_box(t0,t1,box));
};

std::pair<size_t,size_t> mesh3d::CountNodeLeaf() {
  return(mesh_bvh->CountNodeLeaf());
}
