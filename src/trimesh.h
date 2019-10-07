#ifndef TRIMESHH
#define TRIMESHH
#define TINYOBJLOADER_IMPLEMENTATION

#include "triangle.h"
#include "bvh_node.h"
#include "tinyobj/tiny_obj_loader.h"
#include "stb_image.h"
#include <Rcpp.h>
#include <chrono>
#include <thread>

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
  trimesh(std::string inputfile, std::string basedir, Float scale, Float shutteropen, Float shutterclose, random_gen rng) {
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t > shapes;
    std::vector<tinyobj::material_t > materials;
    std::string warn, err;

    bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, inputfile.c_str(), basedir.c_str());
    // Rcpp::Rcerr << inputfile.c_str() << ": " << warn << " " << err << "\n";

    if(ret) {
      int n = 0;
      for (size_t s = 0; s < shapes.size(); s++) {
        n += shapes[s].mesh.num_face_vertices.size();
      }
      std::vector<unsigned char* > obj_materials(materials.size()+1);
      std::vector<vec3 > diffuse_materials(materials.size()+1);
      std::vector<vec3 > specular_materials(materials.size()+1);
      std::vector<Float > ior_materials(materials.size()+1);
      std::vector<bool > has_diffuse(materials.size()+1);
      std::vector<bool > has_transparency(materials.size()+1);
      std::vector<bool > has_single_diffuse(materials.size()+1);
      std::vector<vec3 > image_dim(materials.size()+1);
      std::vector<int > nx_mat(materials.size()+1);
      std::vector<int > ny_mat(materials.size()+1);
      std::vector<int > nn_mat(materials.size()+1);
      
      int nx,ny,nn;

      for (size_t i = 0; i < materials.size(); i++) {
        if(strlen(materials[i].diffuse_texname.c_str()) > 0) {
          obj_materials[i] = stbi_load((basedir + separator() + materials[i].diffuse_texname).c_str(), &nx, &ny, &nn, 0);
          has_diffuse[i] = true;
          has_single_diffuse[i] = false;
          nx_mat[i] = nx;
          ny_mat[i] = ny;
          nn_mat[i] = nn;
        } else if (sizeof(materials[i].diffuse) == 12) {
          diffuse_materials[i] = vec3(materials[i].diffuse[0],materials[i].diffuse[1],materials[i].diffuse[2]);
          has_diffuse[i] = true;
          has_single_diffuse[i] = true;
        } else {
          has_diffuse[i] = false;
          has_single_diffuse[i] = false;
        }
        if(materials[i].dissolve < 1) {
          specular_materials[i] = vec3(materials[i].diffuse[0],materials[i].diffuse[1],materials[i].diffuse[2]);
          ior_materials[i] = materials[i].ior;
          has_transparency[i] = true; 
        }
      }
      bool has_normals = attrib.normals.size() > 0 ? true : false;
      vec3 tris[3];
      vec3 normals[3];
      // vec3 colors[3];
      Float tx[3];
      Float ty[3];
      int currenttri=0;
      std::vector<hitable* > triangles(n+1);
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
            }
          }
          
          index_offset += 3;
          if(std::isnan(tris[0].x()) || std::isnan(tris[0].y()) || std::isnan(tris[0].z()) ||
             std::isnan(tris[1].x()) || std::isnan(tris[1].y()) || std::isnan(tris[1].z()) ||
             std::isnan(tris[2].x()) || std::isnan(tris[2].y()) || std::isnan(tris[2].z())) {
            triangles.pop_back();
            n--;
            continue;
          }
          // per-face material
          int material_num = shapes[s].mesh.material_ids[f];
          
          material *tex;
          
          if(has_normals && tempnormal) {
            if(material_num == -1) {
              tex = new lambertian(new constant_texture(vec3(1,1,1)));
            } else {
              if(has_transparency[material_num]) {
                tex = new dielectric(specular_materials[material_num], ior_materials[material_num], rng);;
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
            triangles[currenttri] = new triangle(tris[0],tris[1],tris[2], normals[0], normals[1], normals[2], tex);
            currenttri++;
          } else {
            if(material_num == -1) {
              tex = new lambertian(new constant_texture(diffuse_materials[material_num]));
            } else {
              if(has_transparency[material_num]) {
                tex = new dielectric(specular_materials[material_num], ior_materials[material_num], rng);;
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
            triangles[currenttri] = new triangle(tris[0],tris[1],tris[2], tex);
            currenttri++;
          }
        }
      }
      tri_mesh_bvh = bvh_node(&triangles[0], n, shutteropen, shutterclose, rng);
    }
  };
  trimesh(std::string inputfile, std::string basedir, material *mat, Float scale, Float shutteropen, Float shutterclose, random_gen rng) {
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t > shapes;
    std::vector<tinyobj::material_t > materials;
    std::string warn, err;

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
      // vec3 colors[3];
      std::vector<hitable* > triangles(n+1);
      int currenttri=0;
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
            triangles.pop_back();
            n--;
            continue;
          }
          if(has_normals) {
            triangles[currenttri] = new triangle(tris[0],tris[1],tris[2],
                                                 normals[0],normals[1],normals[2], 
                                                 mat);
          } else {
            triangles[currenttri] = new triangle(tris[0],tris[1],tris[2],mat);
          }
          currenttri++;
        }
      }
      tri_mesh_bvh = bvh_node(&triangles[0], n, shutteropen, shutterclose, rng);
    }
  };
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
    return(tri_mesh_bvh.hit(r, t_min, t_max, rec, rng));
  };
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const {
    return(tri_mesh_bvh.bounding_box(t0,t1,box));
  };
  bvh_node tri_mesh_bvh;
};

#endif
