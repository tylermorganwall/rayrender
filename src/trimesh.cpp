#include "trimesh.h"

#define TINYOBJLOADER_IMPLEMENTATION
#include "tinyobj/tiny_obj_loader.h"


trimesh::trimesh(std::string inputfile, std::string basedir, Float scale, 
        Float shutteropen, Float shutterclose, int bvh_type, random_gen rng,
        std::shared_ptr<Transform> ObjectToWorld, std::shared_ptr<Transform> WorldToObject, bool reverseOrientation) : 
  hitable(ObjectToWorld, WorldToObject, reverseOrientation) {
  tinyobj::attrib_t attrib;
  std::vector<tinyobj::shape_t > shapes;
  std::vector<tinyobj::material_t > materials;
  std::string warn, err;
  mat_ptr = nullptr;
  
  bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, inputfile.c_str(), basedir.c_str());
  bool has_sep = true;
  if(strlen(basedir.c_str()) == 0) {
    has_sep = false;
  }
  if(ret) {
    int n = 0;
    for (size_t s = 0; s < shapes.size(); s++) {
      n += shapes[s].mesh.num_face_vertices.size();
    }
    std::shared_ptr<alpha_texture> alpha = nullptr;
    std::shared_ptr<bump_texture> bump = nullptr;
    
    obj_materials.reserve(materials.size()+1);
    bump_materials.reserve(materials.size()+1);
    bump_textures.reserve(materials.size()+1);
    
    std::vector<vec3f > diffuse_materials(materials.size()+1);
    std::vector<vec3f > specular_materials(materials.size()+1);
    std::vector<Float > ior_materials(materials.size()+1);
    std::vector<bool > has_diffuse(materials.size()+1);
    std::vector<bool > has_transparency(materials.size()+1);
    std::vector<bool > has_single_diffuse(materials.size()+1);
    std::vector<bool > has_alpha(materials.size()+1);
    std::vector<bool > has_bump(materials.size()+1);
    std::vector<Float > bump_intensity(materials.size()+1);
    
    std::vector<int > nx_mat(materials.size()+1);
    std::vector<int > ny_mat(materials.size()+1);
    std::vector<int > nn_mat(materials.size()+1);
    
    std::vector<int > nx_mat_bump(materials.size()+1);
    std::vector<int > ny_mat_bump(materials.size()+1);
    std::vector<int > nn_mat_bump(materials.size()+1);
    int nx, ny, nn;
    
    for (size_t i = 0; i < materials.size(); i++) {
      if(strlen(materials[i].diffuse_texname.c_str()) > 0) {
        if(has_sep) {
          obj_materials.push_back(stbi_loadf((basedir + separator() + materials[i].diffuse_texname).c_str(), &nx, &ny, &nn, 0));
        } else {
          obj_materials.push_back(stbi_loadf((materials[i].diffuse_texname).c_str(), &nx, &ny, &nn, 0));
        }
        if(nx == 0 || ny == 0 || nn == 0) {
          if(has_sep) {
            throw std::runtime_error("Could not find " + (basedir + separator() + materials[i].diffuse_texname));
          } else {
            throw std::runtime_error("Could not find " + materials[i].diffuse_texname);
          }
        }
        has_diffuse[i] = true;
        has_single_diffuse[i] = false;
        nx_mat[i] = nx;
        ny_mat[i] = ny;
        nn_mat[i] = nn;
        has_alpha[i] = false;
        if(nn == 4) {
          for(int j = 0; j < nx - 1; j++) {
            for(int k = 0; k < ny - 1; k++) {
              if(obj_materials[i][4*j + 4*nx*k + 3] != 1.0) {
                has_alpha[i] = true;
                break;
              }
            }
            if(has_alpha[i]) {
              break;
            }
          }
        } 
      } else if (sizeof(materials[i].diffuse) == 12) {
        obj_materials.push_back(nullptr);
        alpha_materials.push_back(nullptr);
        
        diffuse_materials[i] = vec3f(materials[i].diffuse[0],materials[i].diffuse[1],materials[i].diffuse[2]);
        has_diffuse[i] = true;
        has_alpha[i] = false;
        has_single_diffuse[i] = true;
      } else {
        obj_materials.push_back(nullptr);
        alpha_materials.push_back(nullptr);
        
        has_diffuse[i] = false;
        has_alpha[i] = false;
        has_single_diffuse[i] = false;
      }
      if(materials[i].dissolve < 1) {
        obj_materials.push_back(nullptr);
        alpha_materials.push_back(nullptr);
        
        specular_materials[i] = vec3f(materials[i].diffuse[0],materials[i].diffuse[1],materials[i].diffuse[2]);
        ior_materials[i] = materials[i].ior;
        has_alpha[i] = false;
        has_transparency[i] = true; 
      }
      if(strlen(materials[i].bump_texname.c_str()) > 0) {
        if(has_sep) {
          bump_materials[i] = stbi_loadf((basedir + separator() + materials[i].bump_texname).c_str(), &nx, &ny, &nn, 0);
        } else {
          bump_materials[i] = stbi_loadf(materials[i].bump_texname.c_str(), &nx, &ny, &nn, 0);
        }
        if(nx == 0 || ny == 0 || nn == 0) {
          if(has_sep) {
            throw std::runtime_error("Could not find " + basedir + separator() + materials[i].bump_texname);
          } else {
            throw std::runtime_error("Could not find " + materials[i].bump_texname);
          }
        }
        nx_mat_bump[i] = nx;
        ny_mat_bump[i] = ny;
        nn_mat_bump[i] = nn;
        bump_intensity[i] = materials[i].bump_texopt.bump_multiplier;
        has_bump[i] = true;
      } else {
        bump_materials.push_back(nullptr);
        bump_intensity[i] = 1.0f;
        has_bump[i] = false;
      }
    }
    bool has_normals = attrib.normals.size() > 0 ? true : false;
    
    vec3f tris[3];
    vec3f normals[3];
    Float tx[3];
    Float ty[3];
    for (size_t s = 0; s < shapes.size(); s++) {
      // Loop over faces(polygon)
      size_t index_offset = 0;
      for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
        bool tempnormal = false;
        
        // Loop over vertices in the face.
        for (size_t v = 0; v < 3; v++) {
          // access to vertex
          tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
          tris[v] = vec3f(attrib.vertices[3*idx.vertex_index+0],
                         attrib.vertices[3*idx.vertex_index+1],
                         attrib.vertices[3*idx.vertex_index+2])*scale;
          
          if(has_normals && idx.normal_index != -1) {
            tempnormal = true;
            normals[v] = vec3f(attrib.normals[3*idx.normal_index+0],
                              attrib.normals[3*idx.normal_index+1],
                              attrib.normals[3*idx.normal_index+2]);
          }
          if(idx.texcoord_index != -1) {
            tx[v] = attrib.texcoords[2*idx.texcoord_index+0];
            ty[v] = attrib.texcoords[2*idx.texcoord_index+1];
          }
        }
        index_offset += 3;
        if(std::isnan(tris[0].x()) || std::isnan(tris[0].y()) || std::isnan(tris[0].z()) ||
           std::isnan(tris[1].x()) || std::isnan(tris[1].y()) || std::isnan(tris[1].z()) ||
           std::isnan(tris[2].x()) || std::isnan(tris[2].y()) || std::isnan(tris[2].z())) {
          // triangles.pop_back();
          n--;
          continue;
        }
        // if(std::isnan(normals[0].x()) || std::isnan(normals[0].y()) || std::isnan(normals[0].z()) ||
        //    std::isnan(normals[1].x()) || std::isnan(normals[1].y()) || std::isnan(normals[1].z()) ||
        //    std::isnan(normals[2].x()) || std::isnan(normals[2].y()) || std::isnan(normals[2].z())) {
        //   // triangles.pop_back();
        //   tempnormal = false;
        // }
        
        // per-face material
        int material_num = shapes[s].mesh.material_ids[f];
        if(material_num > -1) {
          if(has_alpha[material_num]) {
            alpha = std::make_shared<alpha_texture>(obj_materials[material_num], 
                                      nx_mat[material_num], ny_mat[material_num], nn_mat[material_num], 
                                      vec3f(tx[0], tx[1], tx[2]), 
                                      vec3f(ty[0], ty[1], ty[2]));
            alpha_materials.push_back(alpha);
          } else {
            alpha = nullptr;
          }
          if(has_bump[material_num]) {
            bump = std::make_shared<bump_texture>(bump_materials[material_num],
                                    nx_mat_bump[material_num], ny_mat_bump[material_num], nn_mat_bump[material_num],
                                    vec3f(tx[0], tx[1], tx[2]),
                                    vec3f(ty[0], ty[1], ty[2]), bump_intensity[material_num]);
            bump_textures.push_back(bump);
          } else {
            bump = nullptr;
            bump_textures.push_back(nullptr);
          }
        } else {
          alpha = nullptr;
          bump = nullptr;
        }
        
        std::shared_ptr<material> tex = nullptr;
        
        if(has_normals && tempnormal) {
          if(material_num == -1) {
            tex = std::make_shared<lambertian>(std::make_shared<constant_texture>(vec3f(1,1,1)));
          } else {
            if(has_transparency[material_num]) {
              tex = std::make_shared<dielectric>(specular_materials[material_num], ior_materials[material_num], vec3f(0,0,0), 
                                   0);;
            } else if(has_diffuse[material_num]) {
              if(has_single_diffuse[material_num]) {
                tex = std::make_shared<lambertian>(std::make_shared<constant_texture>(diffuse_materials[material_num]));
              } else {
                tex = std::make_shared<lambertian>(std::make_shared<triangle_image_texture>(obj_materials[material_num],
                                                                nx_mat[material_num], ny_mat[material_num],nn_mat[material_num],
                                                                tx[0],ty[0], tx[1],ty[1],tx[2],ty[2]));
              }
            } else {
              tex = std::make_shared<lambertian>(std::make_shared<constant_texture>(vec3f(1,1,1)));
            }
          }
          triangles.add(std::make_shared<triangle>(tris[0],tris[1],tris[2], 
                                           normals[0], normals[1], normals[2], true, 
                                           tex, alpha, bump, 
                                           ObjectToWorld, WorldToObject, reverseOrientation));
        } else {
          if(material_num == -1) {
            tex = std::make_shared<lambertian>(std::make_shared<constant_texture>(diffuse_materials[material_num]));
          } else {
            if(has_transparency[material_num]) {
              tex = std::make_shared<dielectric>(specular_materials[material_num], ior_materials[material_num], vec3f(0,0,0),
                                   0);;
            } else if(has_diffuse[material_num]) {
              if(has_single_diffuse[material_num]) {
                tex = std::make_shared<lambertian>(std::make_shared<constant_texture>(diffuse_materials[material_num]));
              } else {
                tex = std::make_shared<lambertian>(std::make_shared<triangle_image_texture>(obj_materials[material_num],
                                                                nx_mat[material_num],ny_mat[material_num],nn_mat[material_num],
                                                                tx[0],ty[0],
                                                                tx[1],ty[1],
                                                                tx[2],ty[2]));
              }
            } else {
              tex = std::make_shared<lambertian>(std::make_shared<constant_texture>(vec3f(1,1,1)));
            }
          }
          triangles.add(std::make_shared<triangle>(tris[0],tris[1],tris[2], true, tex, alpha, bump, 
                                                   ObjectToWorld, WorldToObject, reverseOrientation));
        }
      }
    }
    tri_mesh_bvh = std::make_shared<bvh_node>(triangles, shutteropen, shutterclose, bvh_type, rng);
  } else {
    std::string mes = "Error reading " + inputfile + ": ";
    throw std::runtime_error(mes + warn + err);
  }
}

trimesh::trimesh(std::string inputfile, std::string basedir, Float scale, Float sigma,
        Float shutteropen, Float shutterclose, int bvh_type, random_gen rng,
        std::shared_ptr<Transform> ObjectToWorld, std::shared_ptr<Transform> WorldToObject, bool reverseOrientation) :
      hitable(ObjectToWorld, WorldToObject, reverseOrientation) {
  tinyobj::attrib_t attrib;
  std::vector<tinyobj::shape_t > shapes;
  std::vector<tinyobj::material_t > materials;
  std::string warn, err;
  mat_ptr = nullptr;
  
  bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, inputfile.c_str(), basedir.c_str());
  bool has_sep = true;
  if(strlen(basedir.c_str()) == 0) {
    has_sep = false;
  }
  if(ret) {
    int n = 0;
    for (size_t s = 0; s < shapes.size(); s++) {
      n += shapes[s].mesh.num_face_vertices.size();
    }
    std::shared_ptr<alpha_texture> alpha = nullptr;
    std::shared_ptr<bump_texture> bump = nullptr;
    
    obj_materials.reserve(materials.size()+1);
    bump_materials.reserve(materials.size()+1);
    
    std::vector<vec3f > diffuse_materials(materials.size()+1);
    std::vector<vec3f > specular_materials(materials.size()+1);
    std::vector<Float > ior_materials(materials.size()+1);
    std::vector<bool > has_diffuse(materials.size()+1);
    std::vector<bool > has_transparency(materials.size()+1);
    std::vector<bool > has_single_diffuse(materials.size()+1);
    std::vector<bool > has_alpha(materials.size()+1);
    std::vector<bool > has_bump(materials.size()+1);
    std::vector<Float > bump_intensity(materials.size()+1);
    
    std::vector<int > nx_mat(materials.size()+1);
    std::vector<int > ny_mat(materials.size()+1);
    std::vector<int > nn_mat(materials.size()+1);
    
    std::vector<int > nx_mat_bump(materials.size()+1);
    std::vector<int > ny_mat_bump(materials.size()+1);
    std::vector<int > nn_mat_bump(materials.size()+1);
    
    int nx,ny,nn;
    
    for (size_t i = 0; i < materials.size(); i++) {
      if(strlen(materials[i].diffuse_texname.c_str()) > 0) {
        if(has_sep) {
          obj_materials.push_back(stbi_loadf((basedir + separator() + materials[i].diffuse_texname).c_str(), &nx, &ny, &nn, 0));
        } else {
          obj_materials.push_back(stbi_loadf((materials[i].diffuse_texname).c_str(), &nx, &ny, &nn, 0));
        }
        if(nx == 0 || ny == 0 || nn == 0) {
          if(has_sep) {
            throw std::runtime_error("Could not find " + (basedir + separator() + materials[i].diffuse_texname));
          } else {
            throw std::runtime_error("Could not find " + materials[i].diffuse_texname);
          }
        }
        has_diffuse[i] = true;
        has_single_diffuse[i] = false;
        nx_mat[i] = nx;
        ny_mat[i] = ny;
        nn_mat[i] = nn;
        has_alpha[i] = false;
        if(nn == 4) {
          for(int j = 0; j < nx - 1; j++) {
            for(int k = 0; k < ny - 1; k++) {
              if(obj_materials[i][4*j + 4*nx*k + 3] != 1.0) {
                has_alpha[i] = true;
                break;
              }
            }
            if(has_alpha[i]) {
              break;
            }
          }
        } 
      } else if (sizeof(materials[i].diffuse) == 12) {
        obj_materials.push_back(nullptr);
        alpha_materials.push_back(nullptr);
        
        diffuse_materials[i] = vec3f(materials[i].diffuse[0],materials[i].diffuse[1],materials[i].diffuse[2]);
        has_diffuse[i] = true;
        has_single_diffuse[i] = true;
      } else {
        obj_materials.push_back(nullptr);
        alpha_materials.push_back(nullptr);
        
        has_diffuse[i] = false;
        has_single_diffuse[i] = false;
      }
      if(materials[i].dissolve < 1) {
        obj_materials.push_back(nullptr);
        alpha_materials.push_back(nullptr);
        
        specular_materials[i] = vec3f(materials[i].diffuse[0],materials[i].diffuse[1],materials[i].diffuse[2]);
        ior_materials[i] = materials[i].ior;
        has_transparency[i] = true; 
      }
      if(strlen(materials[i].bump_texname.c_str()) > 0) {
        bump_materials[i] = stbi_loadf((basedir + separator() + materials[i].bump_texname).c_str(), &nx, &ny, &nn, 0);
        if(has_sep) {
          bump_materials[i] = stbi_loadf((basedir + separator() + materials[i].bump_texname).c_str(), &nx, &ny, &nn, 0);
        } else {
          bump_materials[i] = stbi_loadf(materials[i].bump_texname.c_str(), &nx, &ny, &nn, 0);
        }
        
        if(nx == 0 || ny == 0 || nn == 0) {
          if(has_sep) {
            throw std::runtime_error("Could not find " + basedir + separator() + materials[i].bump_texname);
          } else {
            throw std::runtime_error("Could not find " + materials[i].bump_texname);
          }
        }
        nx_mat_bump[i] = nx;
        ny_mat_bump[i] = ny;
        nn_mat_bump[i] = nn;
        bump_intensity[i] = materials[i].bump_texopt.bump_multiplier;
        has_bump[i] = true;
      } else {
        bump_materials.push_back(nullptr);
        bump_intensity[i] = 1.0f;
        has_bump[i] = false;
      }
    }
    bool has_normals = attrib.normals.size() > 0 ? true : false;
    vec3f tris[3];
    vec3f normals[3];
    Float tx[3];
    Float ty[3];
    for (size_t s = 0; s < shapes.size(); s++) {
      // Loop over faces(polygon)
      size_t index_offset = 0;
      for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
        bool tempnormal = false;
        
        // Loop over vertices in the face.
        for (size_t v = 0; v < 3; v++) {
          // access to vertex
          tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
          tris[v] = vec3f(attrib.vertices[3*idx.vertex_index+0],
                         attrib.vertices[3*idx.vertex_index+1],
                                        attrib.vertices[3*idx.vertex_index+2])*scale;
          if(has_normals && idx.normal_index != -1) {
            tempnormal = true;
            normals[v] = vec3f(attrib.normals[3*idx.normal_index+0],
                              attrib.normals[3*idx.normal_index+1],
                                            attrib.normals[3*idx.normal_index+2]);
          }
          if(idx.texcoord_index != -1) {
            tx[v] = attrib.texcoords[2*idx.texcoord_index+0];
            ty[v] = attrib.texcoords[2*idx.texcoord_index+1];
          }
        }
        
        index_offset += 3;
        if(std::isnan(tris[0].x()) || std::isnan(tris[0].y()) || std::isnan(tris[0].z()) ||
           std::isnan(tris[1].x()) || std::isnan(tris[1].y()) || std::isnan(tris[1].z()) ||
           std::isnan(tris[2].x()) || std::isnan(tris[2].y()) || std::isnan(tris[2].z())) {
          // triangles.pop_back();
          n--;
          continue;
        }
        // per-face material
        int material_num = shapes[s].mesh.material_ids[f];
        if(material_num > -1) {
          if(has_alpha[material_num]) {
            alpha = std::make_shared<alpha_texture>(obj_materials[material_num], nx, ny, nn, 
                                      vec3f(tx[0], tx[1], tx[2]), 
                                      vec3f(ty[0], ty[1], ty[2]));
            alpha_materials.push_back(alpha);
          } else {
            alpha = nullptr;
            alpha_materials.push_back(nullptr);
          }
          if(has_bump[material_num]) {
            bump = std::make_shared<bump_texture>(bump_materials[material_num],
                                    nx_mat_bump[material_num], ny_mat_bump[material_num], nn_mat_bump[material_num],
                                                                                                     vec3f(tx[0], tx[1], tx[2]),
                                                                                                     vec3f(ty[0], ty[1], ty[2]), bump_intensity[material_num]);
            bump_textures.push_back(bump);
          } else {
            bump = nullptr;
            bump_textures.push_back(nullptr);
          }
        }
        
        std::shared_ptr<material> tex = nullptr;
        
        if(has_normals && tempnormal) {
          if(material_num == -1) {
            tex = std::make_shared<orennayar>(std::make_shared<constant_texture>(vec3f(1,1,1)), sigma);
          } else {
            if(has_transparency[material_num]) {
              tex = std::make_shared<dielectric>(specular_materials[material_num], ior_materials[material_num], vec3f(0,0,0), 
                                   0);;
            } else if(has_diffuse[material_num]) {
              if(has_single_diffuse[material_num]) {
                tex = std::make_shared<orennayar>(std::make_shared<constant_texture>(diffuse_materials[material_num]), sigma);
              } else {
                tex = std::make_shared<orennayar>(std::make_shared<triangle_image_texture>(obj_materials[material_num],
                                                               nx_mat[material_num], ny_mat[material_num],nn_mat[material_num],
                                                               tx[0],ty[0], tx[1],ty[1],tx[2],ty[2]), sigma);
              }
            } else {
              tex = std::make_shared<orennayar>(std::make_shared<constant_texture>(vec3f(1,1,1)), sigma);
            }
          }
          triangles.add(std::make_shared<triangle>(tris[0],tris[1],tris[2], 
                                           normals[0], normals[1], normals[2], true, 
                                           tex, alpha,  bump, 
                                           ObjectToWorld, WorldToObject, reverseOrientation));
        } else {
          if(material_num == -1) {
            tex = std::make_shared<orennayar>(std::make_shared<constant_texture>(diffuse_materials[material_num]), sigma);
          } else {
            if(has_transparency[material_num]) {
              tex = std::make_shared<dielectric>(specular_materials[material_num], ior_materials[material_num], vec3f(0,0,0), 
                                   0);;
            } else if(has_diffuse[material_num]) {
              if(has_single_diffuse[material_num]) {
                tex = std::make_shared<orennayar>(std::make_shared<constant_texture>(diffuse_materials[material_num]), sigma);
              } else {
                tex = std::make_shared<orennayar>(std::make_shared<triangle_image_texture>(obj_materials[material_num],
                                                               nx_mat[material_num],ny_mat[material_num],nn_mat[material_num],
                                                                                                               tx[0],ty[0],tx[1],ty[1],tx[2],ty[2]), sigma);
              }
            } else {
              tex = std::make_shared<orennayar>(std::make_shared<constant_texture>(vec3f(1,1,1)), sigma);
            }
          }
          triangles.add(std::make_shared<triangle>(tris[0],tris[1],tris[2], true, tex, alpha, bump, 
                                                   ObjectToWorld, WorldToObject, reverseOrientation));
        }
      }
    }
    tri_mesh_bvh = std::make_shared<bvh_node>(triangles, shutteropen, shutterclose, bvh_type, rng);
  } else {
    std::string mes = "Error reading " + inputfile + ": ";
    throw std::runtime_error(mes + warn + err);
  }
}

trimesh::trimesh(std::string inputfile, std::string basedir, std::shared_ptr<material> mat, 
        Float scale, Float shutteropen, Float shutterclose, int bvh_type, random_gen rng,
        std::shared_ptr<Transform> ObjectToWorld, std::shared_ptr<Transform> WorldToObject, bool reverseOrientation) :
    hitable(ObjectToWorld, WorldToObject, reverseOrientation) {
  tinyobj::attrib_t attrib;
  std::vector<tinyobj::shape_t > shapes;
  std::vector<tinyobj::material_t > materials;
  std::string warn, err;
  mat_ptr = mat;
  
  bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, inputfile.c_str(), basedir.c_str());
  bool has_sep = true;
  if(strlen(basedir.c_str()) == 0) {
    has_sep = false;
  }
  std::shared_ptr<alpha_texture> alpha = nullptr;
  std::shared_ptr<bump_texture> bump = nullptr;
  if(ret) {
    int n = 0;
    for (size_t s = 0; s < shapes.size(); s++) {
      n += shapes[s].mesh.num_face_vertices.size();
    }
    bool has_normals = attrib.normals.size() > 0 ? true : false;
    vec3f tris[3];
    vec3f normals[3];
    for (size_t s = 0; s < shapes.size(); s++) {
      
      size_t index_offset = 0;
      for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
        bool tempnormal = false;
        
        for (size_t v = 0; v < 3; v++) {
          tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
          
          tris[v] = vec3f(attrib.vertices[3*idx.vertex_index+0],
                         attrib.vertices[3*idx.vertex_index+1],
                                        attrib.vertices[3*idx.vertex_index+2])*scale;
          
          if(has_normals  && idx.normal_index != -1) {
            tempnormal = true;
            
            normals[v] = vec3f(attrib.normals[3*idx.normal_index+0],
                              attrib.normals[3*idx.normal_index+1],
                                            attrib.normals[3*idx.normal_index+2]);
          }
        }
        
        index_offset += 3;
        if(std::isnan(tris[0].x()) || std::isnan(tris[0].y()) || std::isnan(tris[0].z()) ||
           std::isnan(tris[1].x()) || std::isnan(tris[1].y()) || std::isnan(tris[1].z()) ||
           std::isnan(tris[2].x()) || std::isnan(tris[2].y()) || std::isnan(tris[2].z())) {
          // triangles.pop_back();
          n--;
          continue;
        }
        if((normals[0].x() == 0 && normals[0].y() == 0 && normals[0].z() == 0) ||
           (normals[1].x() == 0 && normals[1].y() == 0 && normals[1].z() == 0) ||
           (normals[2].x() == 0 && normals[2].y() == 0 && normals[2].z() == 0)) {
          has_normals = false;
        }
        
        if(has_normals && tempnormal) {
          triangles.add(std::make_shared<triangle>(tris[0],tris[1],tris[2],
                                           normals[0],normals[1],normals[2], 
                                           false,
                                           mat_ptr, alpha,  bump, 
                                           ObjectToWorld, WorldToObject, reverseOrientation));
        } else {
          triangles.add(std::make_shared<triangle>(tris[0],tris[1],tris[2], false, 
                                                   mat_ptr, alpha, bump, 
                                                   ObjectToWorld, WorldToObject, reverseOrientation));
        }
      }
    }
    tri_mesh_bvh = std::make_shared<bvh_node>(triangles, shutteropen, shutterclose, bvh_type, rng);
  } else {
    std::string mes = "Error reading " + inputfile + ": ";
    throw std::runtime_error(mes + warn + err);
  }
}

trimesh::trimesh(std::string inputfile, std::string basedir, float vertex_color_sigma,
        Float scale, bool is_vertex_color, Float shutteropen, Float shutterclose, int bvh_type, 
        random_gen rng,
        std::shared_ptr<Transform> ObjectToWorld, std::shared_ptr<Transform> WorldToObject, bool reverseOrientation) : hitable(ObjectToWorld, WorldToObject, reverseOrientation) {
  tinyobj::attrib_t attrib;
  std::vector<tinyobj::shape_t > shapes;
  std::vector<tinyobj::material_t > materials;
  std::string warn, err;
  mat_ptr = nullptr;
  
  bool is_lamb = vertex_color_sigma == 0 ? true : false;
  std::shared_ptr<alpha_texture> alpha = nullptr;
  std::shared_ptr<bump_texture> bump = nullptr;
  
  bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, inputfile.c_str(), basedir.c_str());
  bool has_sep = true;
  if(strlen(basedir.c_str()) == 0) {
    has_sep = false;
  }
  if(ret) {
    int n = 0;
    for (size_t s = 0; s < shapes.size(); s++) {
      n += shapes[s].mesh.num_face_vertices.size();
    }
    bool has_normals = attrib.normals.size() > 0 ? true : false;
    bool has_vertex_colors = attrib.colors.size() > 0 ? true : false;
    vec3f tris[3];
    vec3f normals[3];
    vec3f colors[3];
    for (size_t s = 0; s < shapes.size(); s++) {
      size_t index_offset = 0;
      for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
        bool tempnormal = false;
        
        colors[0] = vec3f(1,1,1);
        colors[1] = vec3f(1,1,1);
        colors[2] = vec3f(1,1,1);
        for (size_t v = 0; v < 3; v++) {
          tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
          
          tris[v] = vec3f(attrib.vertices[3*idx.vertex_index+0],
                         attrib.vertices[3*idx.vertex_index+1],
                                        attrib.vertices[3*idx.vertex_index+2])*scale;
          if(has_vertex_colors) {
            colors[v] = vec3f(attrib.colors[3*idx.vertex_index+0],
                             attrib.colors[3*idx.vertex_index+1],
                                          attrib.colors[3*idx.vertex_index+2]);
          }
          if(has_normals && idx.normal_index != -1) {
            normals[v] = vec3f(attrib.normals[3*idx.normal_index+0],
                              attrib.normals[3*idx.normal_index+1],
                                            attrib.normals[3*idx.normal_index+2]);
          }
        }
        
        index_offset += 3;
        if(std::isnan(tris[0].x()) || std::isnan(tris[0].y()) || std::isnan(tris[0].z()) ||
           std::isnan(tris[1].x()) || std::isnan(tris[1].y()) || std::isnan(tris[1].z()) ||
           std::isnan(tris[2].x()) || std::isnan(tris[2].y()) || std::isnan(tris[2].z())) {
          // triangles.pop_back();
          n--;
          continue;
        }
        std::shared_ptr<material> tex;
        if(is_lamb) {
          tex = std::make_shared<lambertian>(std::make_shared<triangle_texture>(colors[0],colors[1],colors[2]));
        } else {
          tex = std::make_shared<orennayar>(std::make_shared<triangle_texture>(colors[0],colors[1],colors[2]),
                              vertex_color_sigma);
        }
        
        if(has_normals && tempnormal) {
          triangles.add(std::make_shared<triangle>(tris[0],tris[1],tris[2],
                                           normals[0],normals[1],normals[2], 
                                          true,
                                          tex, alpha, bump, 
                                          ObjectToWorld, WorldToObject, reverseOrientation));
        } else {
          triangles.add(std::make_shared<triangle>(tris[0],tris[1],tris[2], true, tex, alpha, bump, 
                                                   ObjectToWorld, WorldToObject, reverseOrientation));
        }
      }
    }
    tri_mesh_bvh = std::make_shared<bvh_node>(triangles, shutteropen, shutterclose, bvh_type, rng);
  } else {
    std::string mes = "Error reading " + inputfile + ": ";
    throw std::runtime_error(mes + warn + err);
  }
}

bool trimesh::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  return(tri_mesh_bvh->hit(r, t_min, t_max, rec, rng));
}

bool trimesh::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
  return(tri_mesh_bvh->hit(r, t_min, t_max, rec, sampler));
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
