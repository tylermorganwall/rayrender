#ifndef TRIMESHH
#define TRIMESHH
#define TINYOBJLOADER_IMPLEMENTATION

#include "triangle.h"
#include "bvh_node.h"
#include "tinyobj/tiny_obj_loader.h"
#include "stb_image.h"
#include <Rcpp.h>

inline char separator() {
  #if defined _WIN32 || defined __CYGWIN__
    return '\\';
  #else
    return '/';
  #endif
}

class trimesh : public hitable {
public:
  trimesh() {}
  ~trimesh() {
    delete tri_mesh_bvh;
    for(auto mat : obj_materials) {
      if(mat) stbi_image_free(mat);
    }
    delete mat_ptr;
  }
  trimesh(std::string inputfile, std::string basedir, Float scale, 
          Float shutteropen, Float shutterclose, random_gen rng) {
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t > shapes;
    std::vector<tinyobj::material_t > materials;
    std::string warn, err;
    mat_ptr = nullptr;
    
    bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, inputfile.c_str(), basedir.c_str());
    // Rcpp::Rcerr << inputfile.c_str() << ": " << warn << " " << err << "\n";

    if(ret) {
      int n = 0;
      for (size_t s = 0; s < shapes.size(); s++) {
        n += shapes[s].mesh.num_face_vertices.size();
      }
      alpha_texture* alpha = nullptr;
      obj_materials.reserve(materials.size()+1);
      // std::vector<Float* > obj_materials(materials.size()+1);
      std::vector<vec3 > diffuse_materials(materials.size()+1);
      std::vector<vec3 > specular_materials(materials.size()+1);
      std::vector<Float > ior_materials(materials.size()+1);
      std::vector<bool > has_diffuse(materials.size()+1);
      std::vector<bool > has_transparency(materials.size()+1);
      std::vector<bool > has_single_diffuse(materials.size()+1);
      std::vector<bool > has_alpha(materials.size()+1);
      std::vector<vec3 > image_dim(materials.size()+1);
      std::vector<int > nx_mat(materials.size()+1);
      std::vector<int > ny_mat(materials.size()+1);
      std::vector<int > nn_mat(materials.size()+1);
      
      int nx, ny, nn;

      for (size_t i = 0; i < materials.size(); i++) {
        if(strlen(materials[i].diffuse_texname.c_str()) > 0) {
          obj_materials.push_back(stbi_loadf((basedir + separator() + materials[i].diffuse_texname).c_str(), &nx, &ny, &nn, 0));
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
          
          diffuse_materials[i] = vec3(materials[i].diffuse[0],materials[i].diffuse[1],materials[i].diffuse[2]);
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
          
          specular_materials[i] = vec3(materials[i].diffuse[0],materials[i].diffuse[1],materials[i].diffuse[2]);
          ior_materials[i] = materials[i].ior;
          has_alpha[i] = false;
          has_transparency[i] = true; 
        }
      }
      bool has_normals = attrib.normals.size() > 0 ? true : false;
      vec3 tris[3];
      vec3 normals[3];
      Float tx[3];
      Float ty[3];
      triangles.reserve(n+1);
      for (size_t s = 0; s < shapes.size(); s++) {
        bool tempnormal = false;
        // Loop over faces(polygon)
        size_t index_offset = 0;
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
          // Loop over vertices in the face.
          for (size_t v = 0; v < 3; v++) {
            // access to vertex
            tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
            tris[v] = vec3(attrib.vertices[3*idx.vertex_index+0],
                           attrib.vertices[3*idx.vertex_index+1],
                           attrib.vertices[3*idx.vertex_index+2])*scale;

            if(has_normals && idx.normal_index != -1) {
              tempnormal = true;
              normals[v] = vec3(attrib.normals[3*idx.normal_index+0],
                                attrib.normals[3*idx.normal_index+1],
                                attrib.normals[3*idx.normal_index+2]);
            }
            if(idx.texcoord_index != -1) {
              tx[v] = attrib.texcoords[2*idx.texcoord_index+0];
              ty[v] = attrib.texcoords[2*idx.texcoord_index+1];
              while(tx[v] < 0) tx[v] += 1;
              while(ty[v] < 0) ty[v] += 1;
              while(tx[v] > 1) tx[v] -= 1;
              while(ty[v] > 1) ty[v] -= 1;
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
              alpha = new alpha_texture(obj_materials[material_num], nx, ny, nn, 
                                        vec3(tx[0], tx[1], tx[2]), 
                                        vec3(ty[0], ty[1], ty[2]));
              alpha_materials.push_back(alpha);
            } else {
              alpha = nullptr;
            }
          } else {
            alpha = nullptr;
          }
          
          material *tex;
          
          if(has_normals && tempnormal) {
            if(material_num == -1) {
              tex = new lambertian(new constant_texture(vec3(1,1,1)));
            } else {
              if(has_transparency[material_num]) {
                tex = new dielectric(specular_materials[material_num], ior_materials[material_num], vec3(0,0,0), 
                                     0,  rng);;
              } else if(has_diffuse[material_num]) {
                if(has_single_diffuse[material_num]) {
                  tex = new lambertian(new constant_texture(diffuse_materials[material_num]));
                } else {
                  tex = new lambertian(new triangle_image_texture(obj_materials[material_num],
                                        nx_mat[material_num], ny_mat[material_num],nn_mat[material_num],
                                        tx[0],ty[0], tx[1],ty[1],tx[2],ty[2]));
                }
              } else {
                tex = new lambertian(new constant_texture(vec3(1,1,1)));
              }
            }
            triangles.push_back(new triangle(tris[0],tris[1],tris[2], normals[0], normals[1], normals[2], true, 
                                             tex, alpha));
          } else {
            if(material_num == -1) {
              tex = new lambertian(new constant_texture(diffuse_materials[material_num]));
            } else {
              if(has_transparency[material_num]) {
                tex = new dielectric(specular_materials[material_num], ior_materials[material_num], vec3(0,0,0),
                                     0,  rng);;
              } else if(has_diffuse[material_num]) {
                if(has_single_diffuse[material_num]) {
                  tex = new lambertian(new constant_texture(diffuse_materials[material_num]));
                } else {
                  tex = new lambertian(new triangle_image_texture(obj_materials[material_num],
                                                                  nx_mat[material_num],ny_mat[material_num],nn_mat[material_num],
                                                                  tx[0],ty[0],
                                                                  tx[1],ty[1],
                                                                  tx[2],ty[2]));
                }
              } else {
                tex = new lambertian(new constant_texture(vec3(1,1,1)));
              }
            }
            triangles.push_back(new triangle(tris[0],tris[1],tris[2], true, tex, alpha));
          }
        }
      }
      tri_mesh_bvh = new bvh_node(&triangles[0], n, shutteropen, shutterclose, rng);
    }
  };
  trimesh(std::string inputfile, std::string basedir, Float scale, Float sigma,
          Float shutteropen, Float shutterclose, random_gen rng) {
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t > shapes;
    std::vector<tinyobj::material_t > materials;
    std::string warn, err;
    mat_ptr = nullptr;
    
    bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, inputfile.c_str(), basedir.c_str());

    if(ret) {
      int n = 0;
      for (size_t s = 0; s < shapes.size(); s++) {
        n += shapes[s].mesh.num_face_vertices.size();
      }
      alpha_texture *alpha = nullptr;
      obj_materials.reserve(materials.size()+1);
      // std::vector<Float* > obj_materials(materials.size()+1);
      std::vector<vec3 > diffuse_materials(materials.size()+1);
      std::vector<vec3 > specular_materials(materials.size()+1);
      std::vector<Float > ior_materials(materials.size()+1);
      std::vector<bool > has_diffuse(materials.size()+1);
      std::vector<bool > has_transparency(materials.size()+1);
      std::vector<bool > has_single_diffuse(materials.size()+1);
      std::vector<bool > has_alpha(materials.size()+1);
      std::vector<vec3 > image_dim(materials.size()+1);
      std::vector<int > nx_mat(materials.size()+1);
      std::vector<int > ny_mat(materials.size()+1);
      std::vector<int > nn_mat(materials.size()+1);
      
      int nx,ny,nn;
      
      for (size_t i = 0; i < materials.size(); i++) {
        if(strlen(materials[i].diffuse_texname.c_str()) > 0) {
          obj_materials.push_back(stbi_loadf((basedir + separator() + materials[i].diffuse_texname).c_str(), &nx, &ny, &nn, 0));
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
          
          diffuse_materials[i] = vec3(materials[i].diffuse[0],materials[i].diffuse[1],materials[i].diffuse[2]);
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
          
          specular_materials[i] = vec3(materials[i].diffuse[0],materials[i].diffuse[1],materials[i].diffuse[2]);
          ior_materials[i] = materials[i].ior;
          has_transparency[i] = true; 
        }
      }
      bool has_normals = attrib.normals.size() > 0 ? true : false;
      vec3 tris[3];
      vec3 normals[3];
      Float tx[3];
      Float ty[3];
      triangles.reserve(n+1);
      for (size_t s = 0; s < shapes.size(); s++) {
        bool tempnormal = false;
        // Loop over faces(polygon)
        size_t index_offset = 0;
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
          // Loop over vertices in the face.
          for (size_t v = 0; v < 3; v++) {
            // access to vertex
            tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
            tris[v] = vec3(attrib.vertices[3*idx.vertex_index+0],
                           attrib.vertices[3*idx.vertex_index+1],
                           attrib.vertices[3*idx.vertex_index+2])*scale;
            if(has_normals && idx.normal_index != -1) {
              tempnormal = true;
              normals[v] = vec3(attrib.normals[3*idx.normal_index+0],
                                attrib.normals[3*idx.normal_index+1],
                                attrib.normals[3*idx.normal_index+2]);
            }
            if(idx.texcoord_index != -1) {
              tx[v] = attrib.texcoords[2*idx.texcoord_index+0];
              ty[v] = attrib.texcoords[2*idx.texcoord_index+1];
              while(tx[v] < 0) tx[v] += 1;
              while(ty[v] < 0) ty[v] += 1;
              while(tx[v] > 1) tx[v] -= 1;
              while(ty[v] > 1) ty[v] -= 1;
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
              alpha = new alpha_texture(obj_materials[material_num], nx, ny, nn, 
                                        vec3(tx[0], tx[1], tx[2]), 
                                        vec3(ty[0], ty[1], ty[2]));
              alpha_materials.push_back(alpha);
            } else {
              alpha = nullptr;
              alpha_materials.push_back(nullptr);
            }
          }
          
          material *tex;
          
          if(has_normals && tempnormal) {
            if(material_num == -1) {
              tex = new orennayar(new constant_texture(vec3(1,1,1)), sigma);
            } else {
              if(has_transparency[material_num]) {
                tex = new dielectric(specular_materials[material_num], ior_materials[material_num], vec3(0,0,0), 
                                     0,  rng);;
              } else if(has_diffuse[material_num]) {
                if(has_single_diffuse[material_num]) {
                  tex = new orennayar(new constant_texture(diffuse_materials[material_num]), sigma);
                } else {
                  tex = new orennayar(new triangle_image_texture(obj_materials[material_num],
                                                                  nx_mat[material_num], ny_mat[material_num],nn_mat[material_num],
                                                                  tx[0],ty[0], tx[1],ty[1],tx[2],ty[2]), sigma);
                }
              } else {
                tex = new orennayar(new constant_texture(vec3(1,1,1)), sigma);
              }
            }
            triangles.push_back(new triangle(tris[0],tris[1],tris[2], normals[0], normals[1], normals[2], true, 
                                             tex,alpha));
          } else {
            if(material_num == -1) {
              tex = new orennayar(new constant_texture(diffuse_materials[material_num]), sigma);
            } else {
              if(has_transparency[material_num]) {
                tex = new dielectric(specular_materials[material_num], ior_materials[material_num], vec3(0,0,0), 
                                     0,  rng);;
              } else if(has_diffuse[material_num]) {
                if(has_single_diffuse[material_num]) {
                  tex = new orennayar(new constant_texture(diffuse_materials[material_num]), sigma);
                } else {
                  tex = new orennayar(new triangle_image_texture(obj_materials[material_num],
                                                                  nx_mat[material_num],ny_mat[material_num],nn_mat[material_num],
                                                                  tx[0],ty[0],tx[1],ty[1],tx[2],ty[2]), sigma);
                }
              } else {
                tex = new orennayar(new constant_texture(vec3(1,1,1)), sigma);
              }
            }
            triangles.push_back(new triangle(tris[0],tris[1],tris[2], true, tex, alpha));
          }
        }
      }
      tri_mesh_bvh = new bvh_node(&triangles[0], n, shutteropen, shutterclose, rng);
    }
  };
  trimesh(std::string inputfile, std::string basedir, material *mat, 
          Float scale, Float shutteropen, Float shutterclose, random_gen rng) {
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t > shapes;
    std::vector<tinyobj::material_t > materials;
    std::string warn, err;
    mat_ptr = mat;

    bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, inputfile.c_str(), basedir.c_str());
    //need to define else
    if(ret) {
      int n = 0;
      for (size_t s = 0; s < shapes.size(); s++) {
        n += shapes[s].mesh.num_face_vertices.size();
      }
      bool has_normals = attrib.normals.size() > 0 ? true : false;
      vec3 tris[3];
      vec3 normals[3];
      triangles.reserve(n+1);
      for (size_t s = 0; s < shapes.size(); s++) {
        size_t index_offset = 0;
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
          for (size_t v = 0; v < 3; v++) {
            tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
            
            tris[v] = vec3(attrib.vertices[3*idx.vertex_index+0],
                           attrib.vertices[3*idx.vertex_index+1],
                           attrib.vertices[3*idx.vertex_index+2])*scale;

            if(has_normals) {
              normals[v] = vec3(attrib.normals[3*idx.normal_index+0],
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
          
          if(has_normals) {
            triangles.push_back(new triangle(tris[0],tris[1],tris[2],
                                             normals[0],normals[1],normals[2], 
                                             false,
                                             mat_ptr, nullptr));
          } else {
            triangles.push_back(new triangle(tris[0],tris[1],tris[2], false, mat_ptr, nullptr));
          }
        }
      }
      tri_mesh_bvh = new bvh_node(&triangles[0], n, shutteropen, shutterclose, rng);
    }
  };
  trimesh(std::string inputfile, std::string basedir, float vertex_color_sigma,
          Float scale, bool is_vertex_color, Float shutteropen, Float shutterclose, random_gen rng) {
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t > shapes;
    std::vector<tinyobj::material_t > materials;
    std::string warn, err;
    mat_ptr = nullptr;
    
    bool is_lamb = vertex_color_sigma == 0 ? true : false;
    
    bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, inputfile.c_str(), basedir.c_str());
    //need to define else
    if(ret) {
      int n = 0;
      for (size_t s = 0; s < shapes.size(); s++) {
        n += shapes[s].mesh.num_face_vertices.size();
      }
      bool has_normals = attrib.normals.size() > 0 ? true : false;
      bool has_vertex_colors = attrib.colors.size() > 0 ? true : false;
      vec3 tris[3];
      vec3 normals[3];
      vec3 colors[3];
      triangles.reserve(n+1);
      for (size_t s = 0; s < shapes.size(); s++) {
        size_t index_offset = 0;
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
          colors[0] = vec3(1,1,1);
          colors[1] = vec3(1,1,1);
          colors[2] = vec3(1,1,1);
          for (size_t v = 0; v < 3; v++) {
            tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
            
            tris[v] = vec3(attrib.vertices[3*idx.vertex_index+0],
                           attrib.vertices[3*idx.vertex_index+1],
                           attrib.vertices[3*idx.vertex_index+2])*scale;
            if(has_vertex_colors) {
              colors[v] = vec3(attrib.colors[3*idx.vertex_index+0],
                               attrib.colors[3*idx.vertex_index+1],
                               attrib.colors[3*idx.vertex_index+2]);
            }
            if(has_normals) {
              normals[v] = vec3(attrib.normals[3*idx.normal_index+0],
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
          material* tex;
          if(is_lamb) {
            tex = new lambertian(new triangle_texture(colors[0],colors[1],colors[2]));
          } else {
            tex = new orennayar(new triangle_texture(colors[0],colors[1],colors[2]),
                                vertex_color_sigma);
          }
          
          if(has_normals) {
            triangles.push_back(new triangle(tris[0],tris[1],tris[2],
                                             normals[0],normals[1],normals[2], 
                                             true,
                                             tex, nullptr));
          } else {
            triangles.push_back(new triangle(tris[0],tris[1],tris[2], true, tex, nullptr));
          }
        }
      }
      tri_mesh_bvh = new bvh_node(&triangles[0], n, shutteropen, shutterclose, rng);
    }
  };
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
    return(tri_mesh_bvh->hit(r, t_min, t_max, rec, rng));
  };
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const {
    return(tri_mesh_bvh->bounding_box(t0,t1,box));
  };
  bvh_node* tri_mesh_bvh;
  material *mat_ptr;
  std::vector<Float* > obj_materials;
  std::vector<alpha_texture* > alpha_materials;
  std::vector<hitable* > triangles;
};

#endif
