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
  float scale_mesh = Rcpp::as<float>(mesh_info["scale_mesh"]);
  
  std::string texture_location = Rcpp::as<std::string>(mesh_info["texture"]);
  std::string bump_text_location = Rcpp::as<std::string>(mesh_info["bump_texture"]);
  Float bump_intensity = Rcpp::as<Float>(mesh_info["bump_intensity"]);
  int material_type = Rcpp::as<Float>(mesh_info["material_type"]);
  
  
  Rcpp::NumericMatrix colors = Rcpp::as<Rcpp::NumericMatrix>(mesh_info["color_vals"]);
  int colortype = Rcpp::as<int>(mesh_info["color_type"]);
  
  int nx = 0,ny = 0,nn = 0;
  bool has_texture = false;
  bool has_bump = false;
  
  if(strlen(texture_location.c_str()) > 0) {
    mesh_materials = stbi_loadf(texture_location.c_str(), &nx, &ny, &nn, 0);
    has_texture = nx != 0 && ny != 0 && nn != 0;
  } else {
    mesh_materials = nullptr;
  }
  int nxb,nyb,nnb;
  
  if(strlen(bump_text_location.c_str()) > 0) {
    bump = stbi_loadf(bump_text_location.c_str(), &nxb, &nyb, &nnb, 0);
    has_bump = nxb != 0 && nyb != 0 && nnb != 0;
  } else {
    bump = nullptr;
  }
  mat_ptr = mat;
  int number_faces = indices.nrow();
  bool has_normals = norms.nrow() > 1;
  bool has_texcoords = txcoord.nrow() > 1;
  // std::shared_ptr<image_texture> tex = nullptr;
  bool single_tex = colortype != 3;
  
  for (int i = 0; i < number_faces; i++) {
    vec3f tris[3];
    vec3f normals[3];

    point2f tx[3];
    
    int idx[3] = {indices(i,0),indices(i,1),indices(i,2)};
    tris[0] = vec3f(vertices(idx[0],0),vertices(idx[0],1),vertices(idx[0],2))*scale_mesh;
    tris[1] = vec3f(vertices(idx[1],0),vertices(idx[1],1),vertices(idx[1],2))*scale_mesh;
    tris[2] = vec3f(vertices(idx[2],0),vertices(idx[2],1),vertices(idx[2],2))*scale_mesh;
    
    if(has_normals) {
      normals[0] = vec3f(norms(idx[0],0),norms(idx[0],1),norms(idx[0],2));
      normals[1] = vec3f(norms(idx[1],0),norms(idx[1],1),norms(idx[1],2));
      normals[2] = vec3f(norms(idx[2],0),norms(idx[2],1),norms(idx[2],2));
    }
    
    if(has_texcoords) {
      tx[0] = point2f(txcoord(idx[0],0),txcoord(idx[0],1));
      tx[1] = point2f(txcoord(idx[1],0),txcoord(idx[1],1));
      tx[2] = point2f(txcoord(idx[2],0),txcoord(idx[2],1));
    }
    std::shared_ptr<texture> tex = nullptr;
    if(colortype == 3 && has_texcoords && has_texture) {
      
      tex = std::make_shared<image_texture>(mesh_materials,
                                                     nx,ny,nn);
      MicrofacetDistribution *dist;
      switch(material_type) {
        case 1: {
          mat_ptr = std::make_shared<lambertian>(tex);
          break;
        }
        case 2: {
          mat_ptr = std::make_shared<metal>(tex,
                                            tempvector(3),
                                            point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)),
                                            point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
          break;
        }
        case 3: {
          mat_ptr = mat;
        }
        case 4: {
          mat_ptr = std::make_shared<orennayar>(tex, sigma);
          break;
        }
        case 5: {
          bool is_invisible = tempvector(3) == 1;
          mat_ptr = std::make_shared<diffuse_light>(tex, lightintensity, is_invisible);
          break;
        }
        case 6: {
          if(temp_glossy(0) == 1) {
            dist = new TrowbridgeReitzDistribution(temp_glossy(1), temp_glossy(2), nullptr, false,  true);
          } else {
            dist = new BeckmannDistribution(temp_glossy(1), temp_glossy(2), nullptr, false, true);
          }
          mat_ptr = std::make_shared<MicrofacetReflection>(tex, dist,
                                                             point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)),
                                                             point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
          break;
        }
        case 7: {
          if(temp_glossy(0) == 1) {
          dist = new TrowbridgeReitzDistribution(temp_glossy(1), temp_glossy(2), nullptr, false,  true);
        } else {
          dist = new BeckmannDistribution(temp_glossy(1), temp_glossy(2), nullptr, false, true);
        }
          mat_ptr =  std::make_shared<glossy>(tex, dist,
                                            point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)),
                                            point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
          break;
        }
        case 8: {
          Rprintf("`mesh3d()` object does not support spot lights--using `light()`");
          mat_ptr = std::make_shared<diffuse_light>(tex, lightintensity, false); //light intensity, invisible
          break;
        }
        case 9: {
          Rprintf("`mesh3d()` object does not support hair material--using `diffuse()`");
          mat_ptr = std::make_shared<lambertian>(tex);
          break;
        }
        case 10: {
          if(temp_glossy(0) == 1) {
            dist = new TrowbridgeReitzDistribution(temp_glossy(1), temp_glossy(2), nullptr, false,  true);
          } else {
            dist = new BeckmannDistribution(temp_glossy(1), temp_glossy(2), nullptr, false, true);
          }
          mat_ptr = std::make_shared<MicrofacetTransmission>(tex, dist,
                                                           point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)),
                                                           point3f(temp_glossy(6), temp_glossy(7), temp_glossy(8)));
          break;
        }
      }
    } else if(colortype == 2) {
      tex = std::make_shared<triangle_texture>(
        vec3f(colors(i,0),colors(i,1),colors(i,2)),
        vec3f(colors(i,0),colors(i,1),colors(i,2)),
        vec3f(colors(i,0),colors(i,1),colors(i,2)));
    } else if(colortype == 4) {
      tex = std::make_shared<triangle_texture>(
        vec3f(colors(idx[0],0),colors(idx[0],1),colors(idx[0],2)),
        vec3f(colors(idx[1],0),colors(idx[1],1),colors(idx[1],2)),
        vec3f(colors(idx[2],0),colors(idx[2],1),colors(idx[2],2)));
    } 
    std::shared_ptr<bump_texture> bump_tex = nullptr;
    if(has_bump) {
      bump_tex = std::make_shared<bump_texture>(bump, nxb, nyb, nnb,
                                                bump_intensity);
    }
    
    if(!has_texcoords) {
      if(has_normals) {
        triangles.add(std::make_shared<triangle>(tris[0],tris[1],tris[2],
                                         normals[0],normals[1],normals[2],
                                         single_tex,
                                         mat_ptr, nullptr,  bump_tex, 
                                         ObjectToWorld, WorldToObject, reverseOrientation));
      } else {
        triangles.add(std::make_shared<triangle>(tris[0],tris[1],tris[2], 
                                         single_tex, 
                                         mat_ptr, nullptr, bump_tex, 
                                         ObjectToWorld, WorldToObject, reverseOrientation));
      }
    } else {
      if(has_normals) {
        triangles.add(std::make_shared<triangle>(tris[0],tris[1],tris[2],
                                                 normals[0],normals[1],normals[2],
                                                 tx[0],tx[1],tx[2],                             
                                                 single_tex,
                                                 mat_ptr, nullptr,  bump_tex, 
                                                 ObjectToWorld, WorldToObject, reverseOrientation));
      } else {
        triangles.add(std::make_shared<triangle>(tris[0],tris[1],tris[2], 
                                                 tx[0],tx[1],tx[2],
                                                 single_tex, 
                                                 mat_ptr, nullptr, bump_tex, 
                                                 ObjectToWorld, WorldToObject, reverseOrientation));
      }
    }
  }
  mesh_bvh = std::make_shared<bvh_node>(triangles, shutteropen, shutterclose, bvh_type, rng);
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
