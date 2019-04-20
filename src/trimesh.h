#ifndef TRIMESHH
#define TRIMESHH
#define TINYOBJLOADER_IMPLEMENTATION

#include "triangle.h"
#include "bvh_node.h"
#include "tinyobj/tiny_obj_loader.h"
#include "stb_image.h"
#include <Rcpp.h>

class trimesh : public hitable {
public:
  trimesh() {}
  trimesh(std::string inputfile, float scale, float shutteropen, float shutterclose, random_gen rng) {
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t > shapes;
    std::vector<tinyobj::material_t > materials;
    std::string warn, err;
    bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, inputfile.c_str());
    if(ret) {
      int n = 0;
      for (size_t s = 0; s < shapes.size(); s++) {
        n += shapes[s].mesh.num_face_vertices.size();
      }
      std::vector<unsigned char* > obj_materials(materials.size()+1);
      std::vector<bool > has_diffuse(materials.size()+1);
      std::vector<vec3 > image_dim(materials.size()+1);
      std::vector<int > nx_mat(materials.size()+1);
      std::vector<int > ny_mat(materials.size()+1);
      int nx,ny,nn;
      for (size_t i = 0; i < materials.size(); i++) {
        if(strlen(materials[i].diffuse_texname.c_str()) > 0) {
          obj_materials[i] = stbi_load(materials[i].diffuse_texname.c_str(), &nx, &ny, &nn, 0);
          // Rcpp::Rcout << materials[i].diffuse_texname << ": " <<  nx << " " << nx << " " << "\n";
          has_diffuse[i] = true;
          nx_mat[i] = nx;
          ny_mat[i] = ny;
        } else {
          has_diffuse[i] = false;
        }
      }
      bool has_normals = attrib.normals.size() > 0 ? true : false;
      vec3 tris[3];
      vec3 normals[3];
      vec3 colors[3];
      float tx[3];
      float ty[3];
      int currenttri=0;
      std::vector<hitable* > triangles(n+1);
      for (size_t s = 0; s < shapes.size(); s++) {
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
            if(has_normals) {
              normals[v] = vec3(attrib.normals[3*idx.normal_index+0],
                                attrib.normals[3*idx.normal_index+1],
                                attrib.normals[3*idx.normal_index+2]);
            }
            tx[v] = attrib.texcoords[2*idx.texcoord_index+0];
            // tx[v] = tx[v] < 0 ? -tx[v] : tx[v];
            ty[v] = attrib.texcoords[2*idx.texcoord_index+1];
            // ty[v] = ty[v] < 0 ? -ty[v] : ty[v];
            // Optional: vertex colors
            // colors[v] = vec3(attrib.colors[3*idx.vertex_index+0],
            //                  attrib.colors[3*idx.vertex_index+1],
            //                  attrib.colors[3*idx.vertex_index+2]);
          }
          
          index_offset += 3;
          int material_num = shapes[s].mesh.material_ids[f];
          
          // per-face material
          // shapes[s].mesh.material_ids[f];
          if(has_normals) {
            triangles[currenttri] = new triangle(tris[0],tris[1],tris[2],normals[0],normals[1],normals[2], 
                                        new lambertian(new triangle_image_texture(obj_materials[material_num],
                                                                                  nx_mat[material_num], ny_mat[material_num],
                                                                                  tx[0],ty[0],
                                                                                  tx[1],ty[1],
                                                                                  tx[2],ty[2])));
            currenttri++;
          } else {
            if(has_diffuse[material_num]) {
              triangles[currenttri] = new triangle(tris[0],tris[1],tris[2],
                                          new lambertian(new triangle_image_texture(obj_materials[material_num],
                                                                                    nx_mat[material_num],ny_mat[material_num],
                                                                                    tx[0],ty[0],
                                                                                    tx[1],ty[1],
                                                                                    tx[2],ty[2])));
            } else {
              triangles[currenttri] = new triangle(tris[0],tris[1],tris[2],new lambertian(new constant_texture(vec3(1,1,1))));
            }
            currenttri++;
          }
        }
      }
      tri_mesh_bvh = bvh_node(&triangles[0], n, shutteropen, shutterclose, rng);
    }
  };
  trimesh(std::string inputfile, material *mat, float scale, float shutteropen, float shutterclose, random_gen rng) {
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t > shapes;
    std::vector<tinyobj::material_t > materials;
    std::string warn, err;
    
    bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, inputfile.c_str());
    //need to define else
    if(ret) {
      int n = 0;
      for (size_t s = 0; s < shapes.size(); s++) {
        n += shapes[s].mesh.num_face_vertices.size();
      }
      bool has_normals = attrib.normals.size() > 0 ? true : false;
      vec3 tris[3];
      vec3 normals[3];
      vec3 colors[3];
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
  virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec, random_gen& rng) {
    return(tri_mesh_bvh.hit(r, t_min, t_max, rec, rng));
  };
  virtual bool bounding_box(float t0, float t1, aabb& box) const {
    return(tri_mesh_bvh.bounding_box(t0,t1,box));
  };
  bvh_node tri_mesh_bvh;
};

#endif