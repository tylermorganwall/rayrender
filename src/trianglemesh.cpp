#include "trianglemesh.h"

#ifndef TINYOBJLOADER_IMPLEMENTATION
#define TINYOBJLOADER_IMPLEMENTATION
#define TINYOBJLOADER_USE_MAPBOX_EARCUT
#include "tinyobj/tiny_obj_loader.h"
#endif

#include <queue>



inline Float luminance(point3f& color) {
  return(dot(color,point3f(0.2125,0.7154,0.0721)));
}

void LoadMtlMaterials(std::vector<std::shared_ptr<material> > &mtl_materials,
                      std::vector<tinyobj::material_t > &materials,
                      std::vector<unsigned char * > &obj_texture_data,
                      std::vector<unsigned char * > &bump_texture_data,
                      std::vector<std::shared_ptr<bump_texture> > &bump_textures,
                      std::vector<std::shared_ptr<alpha_texture> > &alpha_textures,
                      size_t &texture_size,
                      const std::string inputfile, const std::string basedir, bool has_sep,
                      std::shared_ptr<material> default_material, bool load_materials,
                      bool load_textures, bool verbose, std::vector<bool>& material_is_light) {
  mtl_materials.reserve(materials.size()+1);
  obj_texture_data.reserve(materials.size()+1);
  bump_texture_data.reserve(materials.size()+1);
  bump_textures.reserve(materials.size()+1);
  alpha_textures.reserve(materials.size()+1);
  material_is_light.reserve(materials.size()+1);
  
  //For default texture
  alpha_textures.push_back(nullptr);
  bump_textures.push_back(nullptr);
  
  std::vector<vec3f > diffuse_materials(materials.size()+1);
  std::vector<vec3f > specular_materials(materials.size()+1);
  std::vector<Float > ior_materials(materials.size()+1);
  std::vector<bool > has_diffuse(materials.size()+1, false);
  std::vector<bool > has_transparency(materials.size()+1, false);
  std::vector<bool > ior(materials.size()+1, 1.0);
  
  std::vector<bool > has_single_diffuse(materials.size()+1, false);
  std::vector<bool > has_alpha(materials.size()+1, false);
  std::vector<bool > has_bump(materials.size()+1, false);
  std::vector<Float > bump_intensity(materials.size()+1);
  
  std::vector<int > nx_mat(materials.size()+1);
  std::vector<int > ny_mat(materials.size()+1);
  std::vector<int > nn_mat(materials.size()+1);
  
  std::vector<int > nx_mat_bump(materials.size()+1);
  std::vector<int > ny_mat_bump(materials.size()+1);
  std::vector<int > nn_mat_bump(materials.size()+1);
  int nx, ny, nn;
  
  if(load_materials) {
    for (size_t i = 0; i < materials.size(); i++) {
      nx = 0; ny = 0; nn = 0;
      if(strlen(materials[i].diffuse_texname.c_str()) > 0 && load_textures) {
        int ok;
        std::replace(materials[i].diffuse_texname.begin(), materials[i].diffuse_texname.end(), '\\', separator());
        if(has_sep) {
          ok = stbi_info((basedir + separator() + materials[i].diffuse_texname).c_str(), &nx, &ny, &nn);
          obj_texture_data.push_back(stbi_load((basedir + separator() + materials[i].diffuse_texname).c_str(), &nx, &ny, &nn, 0));
        } else {
          ok = stbi_info((materials[i].diffuse_texname).c_str(), &nx, &ny, &nn);
          obj_texture_data.push_back(stbi_load((materials[i].diffuse_texname).c_str(), &nx, &ny, &nn, 0));
        }

        if(!obj_texture_data[i] || !ok) {
          REprintf("Load failed: %s\n", stbi_failure_reason());
          if(has_sep) {
            throw std::runtime_error("Loading failed of: " + (basedir + separator() + materials[i].diffuse_texname) + 
                                     "-- nx/ny/channels :"  + std::to_string(nx)  +  "/"  +  std::to_string(ny)  +  "/"  +  std::to_string(nn));
          } else {
            throw std::runtime_error("Loading failed of: " + materials[i].diffuse_texname + 
                                     "-- nx/ny/channels :" + std::to_string(nx)  +  "/"  +  std::to_string(ny)  +  "/"  +  std::to_string(nn));
          }
        }
        if(nx == 0 || ny == 0 || nn == 0) {
          if(has_sep) {
            throw std::runtime_error("Could not find " + (basedir + separator() + materials[i].diffuse_texname));
          } else {
            throw std::runtime_error("Could not find " + materials[i].diffuse_texname);
          }
        }
        if(verbose) {
          Rprintf("(%i/%i) Loading Material Texture %s (%i/%i/%i) \n", i+1, materials.size(),materials[i].name.c_str(),nx,ny,nn);
        }

        texture_size += sizeof(unsigned char) * nx * ny * nn;
        has_diffuse[i] = true;
        has_single_diffuse[i] = false;
        nx_mat[i] = nx;
        ny_mat[i] = ny;
        nn_mat[i] = nn;
        has_alpha[i] = false;
        if(nn == 4) {
          for(int j = 0; j < nx - 1; j++) {
            for(int k = 0; k < ny - 1; k++) {
              if(obj_texture_data[i][4*j + 4*nx*k + 3] != 255) {
                has_alpha[i] = true;
                break;
              }
            }
            if(has_alpha[i]) {
              break;
            }
          }
        } 
      } else if (sizeof(materials[i].diffuse) == 12 && materials[i].dissolve == 1) {
        obj_texture_data.push_back(nullptr);
        diffuse_materials[i] = vec3f(materials[i].diffuse[0],materials[i].diffuse[1],materials[i].diffuse[2]);
        has_diffuse[i] = true;
        has_alpha[i] = false;
        has_single_diffuse[i] = true;
      } else if(materials[i].dissolve < 1) {
        obj_texture_data.push_back(nullptr);
        specular_materials[i] = vec3f(materials[i].diffuse[0],materials[i].diffuse[1],materials[i].diffuse[2]);
        ior_materials[i] = materials[i].ior;
        has_alpha[i] = false;
        has_transparency[i] = true; 
      } else {
        obj_texture_data.push_back(nullptr);
        has_diffuse[i] = false;
        has_alpha[i] = false;
        has_single_diffuse[i] = false;
      }
      if(strlen(materials[i].bump_texname.c_str()) > 0 && load_textures) {
        std::replace(materials[i].bump_texname.begin(), materials[i].bump_texname.end(), '\\', separator());
        
        if(has_sep) {
          bump_texture_data[i] = stbi_load((basedir + separator() + materials[i].bump_texname).c_str(), &nx, &ny, &nn, 0);
        } else {
          bump_texture_data[i] = stbi_load(materials[i].bump_texname.c_str(), &nx, &ny, &nn, 0);
        }
        texture_size += sizeof(unsigned char) * nx * ny * nn;
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
        bump_texture_data.push_back(nullptr);
        bump_intensity[i] = 1.0f;
        has_bump[i] = false;
      }
      std::shared_ptr<alpha_texture> alpha = nullptr;
      std::shared_ptr<bump_texture> bump = nullptr;
      
      if(has_alpha[i]) {
        alpha = std::make_shared<alpha_texture>(obj_texture_data[i], 
                                                nx_mat[i], ny_mat[i], nn_mat[i]);
      } 
      if(has_bump[i]) {
        bump = std::make_shared<bump_texture>(bump_texture_data[i],
                                              nx_mat_bump[i], ny_mat_bump[i], nn_mat_bump[i],
                                              bump_intensity[i]);
      } 
      alpha_textures.push_back(alpha);
      bump_textures.push_back(bump);
    }
  } else {
    alpha_textures.push_back(nullptr);
    bump_textures.push_back(nullptr);
  }
  
  //First texture is default (when shapes[s].mesh.material_ids[f] == -1)
  mtl_materials.push_back(default_material);
  material_is_light.push_back(false);
  
  if(load_materials) {
    for(int material_num = 0; material_num < materials.size(); material_num++) {
      bool imp_sample_obj = false;
      std::shared_ptr<material> tex = nullptr;
      
      if(materials[material_num].illum <= 4) {
        point3f ke(materials[material_num].emission[0],
                   materials[material_num].emission[1],
                   materials[material_num].emission[2]);
        if(ke.x() != 0 || ke.y() != 0 || ke.z() != 0) {
          tex = std::make_shared<diffuse_light>(std::make_shared<constant_texture>(ke), 1.0, false);
          imp_sample_obj = true;
        } else {
          if(has_diffuse[material_num]) {
            if(has_single_diffuse[material_num]) {
              tex = std::make_shared<lambertian>(std::make_shared<constant_texture>(diffuse_materials[material_num]));
            } else {
              tex = std::make_shared<lambertian>(std::make_shared<image_texture_char>(obj_texture_data[material_num],
                                                                                      nx_mat[material_num], 
                                                                                      ny_mat[material_num],
                                                                                      nn_mat[material_num]));
            }
          } else {
            tex = default_material;
          }
        }
      } else {
        point3f spec = point3f(materials[material_num].specular[0],
                               materials[material_num].specular[1],
                               materials[material_num].specular[2]);
        if(materials[material_num].shininess == 1000) {
          tex = std::make_shared<metal>(std::make_shared<constant_texture>(spec),
                                        0., 
                                        point3f(0), 
                                        point3f(0));
        } else {
          point3f atten = point3f(materials[material_num].transmittance[0],
                                  materials[material_num].transmittance[1],
                                  materials[material_num].transmittance[2]);
          tex = std::make_shared<dielectric>(spec, 
                                             materials[material_num].ior, atten, 
                                             0);
        }
      }
      if(verbose) {
        Rprintf("(%i/%i) Loading Material %s (Imp Sample: %s) \n", material_num+1, 
                materials.size(),materials[material_num].name.c_str(), imp_sample_obj ? "true" : "false");
      }
      mtl_materials.push_back(tex);
      material_is_light.push_back(imp_sample_obj);
    }
  }
}

TriangleMesh::TriangleMesh(std::string inputfile, std::string basedir,
                           std::shared_ptr<material> default_material, 
                           bool load_materials, bool load_textures, bool load_vertex_colors, 
                           bool load_normals, bool verbose, Float scale, 
                           bool calculate_consistent_normals,
                           std::shared_ptr<Transform> ObjectToWorld, 
                           std::shared_ptr<Transform> WorldToObject, 
                           bool reverseOrientation) : nTriangles(0) {
  std::string warn, err;
  texture_size = 0;
  vertexIndices.clear();
  normalIndices.clear();
  texIndices.clear();
  face_material_id.clear();
  
  bool has_sep = true;
  bool ret = true;
  has_consistent_normals = calculate_consistent_normals;
  
  tinyobj::ObjReaderConfig reader_config;
  reader_config.mtl_search_path = basedir.c_str(); // Path to material files
  reader_config.vertex_color = load_vertex_colors;
  reader_config.triangulate = true;
  
  tinyobj::ObjReader reader;
  
  if (!reader.ParseFromFile(inputfile, reader_config)) {
    ret = false;
  }
  
  if (!reader.Warning().empty()) {
    Rcpp::Rcout << "TinyObjReader Warning: " << reader.Warning();
  }
  
  auto& attrib = reader.GetAttrib();
  auto& shapes = reader.GetShapes();
  auto& materials = reader.GetMaterials();
  has_normals = false;
  has_tex = false;
  
  has_vertex_colors = attrib.colors.size() > 0 ? true : false;
  if(strlen(basedir.c_str()) == 0) {
    has_sep = false;
  }
  if(ret) {
    nVertices = attrib.vertices.size();
    for (size_t s = 0; s < shapes.size(); s++) {
      nTriangles += shapes[s].mesh.indices.size();
      for(size_t m = 0; m < shapes[s].mesh.indices.size(); m++) {
        vertexIndices.push_back(shapes[s].mesh.indices[m].vertex_index);
        normalIndices.push_back(shapes[s].mesh.indices[m].normal_index);
        texIndices.push_back(shapes[s].mesh.indices[m].texcoord_index);
      }
      if(!has_vertex_colors) {
        if(load_materials) {
          for(size_t i = 0; i < shapes[s].mesh.material_ids.size(); i++ ) {
            face_material_id.push_back(shapes[s].mesh.material_ids[i] + 1);
          } 
        } else {
          for(size_t i = 0; i < shapes[s].mesh.material_ids.size(); i++ ) {
            face_material_id.push_back(0);
          } 
        }
      } 
    }
    
    nNormals = load_normals ? attrib.normals.size() : 0;
    nTex = !has_vertex_colors ? attrib.texcoords.size() : 0;
    p.reset(new point3f[nVertices / 3]);
    for (size_t i = 0; i < nVertices; i += 3) {
      p[i / 3] = (*ObjectToWorld)(point3f(attrib.vertices[i+0],
                                          attrib.vertices[i+1],
                                          attrib.vertices[i+2]) * scale);
    }
    
    
    if(nNormals > 0) {
      has_normals = true;
      n.reset(new normal3f[nNormals / 3]);
      for (size_t i = 0; i < nNormals; i += 3) {
        n[i / 3] = (*ObjectToWorld)(normal3f(attrib.normals[i+0],
                                             attrib.normals[i+1],
                                             attrib.normals[i+2]));
      }
      if(has_consistent_normals) {
        face_n.reset(new normal3f[normalIndices.size() / 3]);
        std::map<int, std::priority_queue<Float> > alpha_values;
        for (size_t i = 0; i < normalIndices.size(); i += 3) {
          int idx_n1 = normalIndices[i];
          int idx_n2 = normalIndices[i+1];
          int idx_n3 = normalIndices[i+2];
          
          normal3f n1 = unit_vector(n[idx_n1]);
          normal3f n2 = unit_vector(n[idx_n2]);
          normal3f n3 = unit_vector(n[idx_n3]);
          
          normal3f face_normal = unit_vector(n1 + n2 + n3);
          face_n[i / 3] = face_normal;
          Float av1 = dot(n1,face_normal);
          Float av2 = dot(n2,face_normal);
          Float av3 = dot(n3,face_normal);
          alpha_values[idx_n1].push(-av1);
          alpha_values[idx_n2].push(-av2);
          alpha_values[idx_n3].push(-av3);
        }
        for (auto const& x : alpha_values) {
          alpha_v.push_back(-x.second.top());
        }
        for(int i = 0; i < alpha_v.size(); i++) {
          Float temp_av = clamp(alpha_v[i],-1,1);
          alpha_v[i] = std::acos(temp_av) * (1 + 0.03632 * (1 - temp_av) * (1 - temp_av));
        }
      }
    } else {
      n = nullptr;
    }

    
    if(nTex > 0) {
      has_tex = true;
      uv.reset(new point2f[nTex / 2]);
      for (size_t i = 0; i < nTex; i += 2) {
        uv[i / 2] = point2f(attrib.texcoords[i+0],
                            attrib.texcoords[i+1]);
      }
    } else {
      uv = nullptr;
    }
    
    if(has_vertex_colors) {
      vc.reset(new point3f[nVertices / 3]);
      for (size_t i = 0; i < nVertices; i += 3) {
        vc[i / 3] = point3f(attrib.colors[i+0],
                            attrib.colors[i+1],
                            attrib.colors[i+2]);
      }
    } else {
      vc = nullptr;
    }
    
    if(!has_vertex_colors) {
      LoadMtlMaterials(mtl_materials, materials, obj_texture_data,
                       bump_texture_data, bump_textures, alpha_textures,
                       texture_size, inputfile, basedir, has_sep, default_material,
                       load_materials, load_textures, verbose, material_is_light);
    } else {
      mtl_materials.push_back(default_material);
      material_is_light.push_back(false);
      alpha_textures.push_back(nullptr);
      bump_textures.push_back(nullptr);
      for (size_t s = 0; s < vertexIndices.size(); s += 3) {
        std::shared_ptr<texture> tex = std::shared_ptr<triangle_texture>(
          new triangle_texture(vc[vertexIndices[s]],
                               vc[vertexIndices[s+1]],
                               vc[vertexIndices[s+2]]));
        mtl_materials.push_back(std::shared_ptr<material>(new lambertian(tex)));
        material_is_light.push_back(false);
        face_material_id.push_back(s / 3 + 1);
        alpha_textures.push_back(nullptr);
        bump_textures.push_back(nullptr);
      }
    }
  } else {
    std::string mes = "Error reading " + inputfile + ": ";
    throw std::runtime_error(mes + reader.Error());
  }
}

TriangleMesh::TriangleMesh(Rcpp::NumericMatrix vertices, 
                           Rcpp::IntegerMatrix indices, 
                           Rcpp::NumericMatrix normals, 
                           Rcpp::NumericMatrix texcoords,
                           Rcpp::NumericMatrix vertexcolors,
                           unsigned char * mesh_texture_data,
                           unsigned char * bump_texture_data_,
                           std::shared_ptr<alpha_texture> alpha,
                           std::shared_ptr<bump_texture> bump,
                           std::shared_ptr<material> default_material, 
                           bool load_materials, bool load_textures,
                           std::shared_ptr<Transform> ObjectToWorld, 
                           std::shared_ptr<Transform> WorldToObject, 
                           bool reverseOrientation) : nTriangles(0) {
  texture_size = 0;
  vertexIndices.clear();
  normalIndices.clear();
  texIndices.clear();
  face_material_id.clear();
  has_normals = false;
  has_tex = false;
  has_vertex_colors = vertexcolors.nrow() > 0;
  has_consistent_normals = false;
  
  nVertices = vertices.nrow();
  nNormals = normals.nrow();
  nTex = texcoords.nrow();
  int nVertexColors = vertexcolors.nrow();
  p.reset(new point3f[nVertices]);
  for (size_t i = 0; i < nVertices; i += 1) {
    p[i] = (*ObjectToWorld)(point3f(vertices(i,0),
                                    vertices(i,1),
                                    vertices(i,2)));
  }
  
  
  if(nNormals > 0) {
    has_normals = true;
    n.reset(new normal3f[nNormals]);
    for (size_t i = 0; i < nNormals; i++) {
      n[i] = (*ObjectToWorld)(normal3f(normals(i,0),
                                       normals(i,1),
                                       normals(i,2)));
    }
  } else {
    n = nullptr;
  }
  
  if(nTex > 0) {
    has_tex = true;
    uv.reset(new point2f[nTex]);
    for (size_t i = 0; i < nTex; i++) {
      uv[i] = point2f(texcoords(i,0),
                      texcoords(i,1));
    }
  } else {
    uv = nullptr;
  }
  if(has_vertex_colors) {
    vc.reset(new point3f[nVertexColors]);
    for (size_t i = 0; i < nVertexColors; i++) {
      vc[i] = point3f(vertexcolors(i,0),
                      vertexcolors(i,1),
                      vertexcolors(i,2));
    }
  } else {
    vc = nullptr;
  }
  
  nTriangles = 0;
  for (size_t s = 0; s < indices.nrow(); s++) {
    nTriangles++;
    vertexIndices.push_back(indices(s,0));
    vertexIndices.push_back(indices(s,1));
    vertexIndices.push_back(indices(s,2));
    if(has_normals) {
      normalIndices.push_back(indices(s,0));
      normalIndices.push_back(indices(s,1));
      normalIndices.push_back(indices(s,2));
    }
    if(has_tex) {
      texIndices.push_back(indices(s,0));
      texIndices.push_back(indices(s,1));
      texIndices.push_back(indices(s,2));
    }
  }
  
  
  //Material stuff
  if(!has_vertex_colors) {
    for (size_t s = 0; s < indices.nrow(); s++) {
      face_material_id.push_back(0);
    }
    mtl_materials.push_back(default_material);
    if(mesh_texture_data) {
      obj_texture_data.push_back(mesh_texture_data);
    }
    if(bump_texture_data_) {
      bump_texture_data.push_back(bump_texture_data_);
    }
    alpha_textures.push_back(alpha);
    bump_textures.push_back(bump);
  } else {
    mtl_materials.push_back(default_material);
    alpha_textures.push_back(nullptr);
    bump_textures.push_back(nullptr);
    for (size_t s = 0; s < vertexIndices.size(); s += 3) {
      std::shared_ptr<texture> tex = std::shared_ptr<triangle_texture>(
        new triangle_texture(vc[vertexIndices[s]],
                             vc[vertexIndices[s+1]],
                             vc[vertexIndices[s+2]]));
      mtl_materials.push_back(std::shared_ptr<material>(new lambertian(tex)));
      face_material_id.push_back(s / 3 + 1);
      alpha_textures.push_back(nullptr);
      bump_textures.push_back(nullptr);
    }
  }
}


TriangleMesh::TriangleMesh(float* vertices, 
                           int* indices, 
                           float* normals, 
                           float* texcoords,
                           int numVerts, int numIndices,
                           std::shared_ptr<alpha_texture> alpha,
                           std::shared_ptr<bump_texture> bump,
                           std::shared_ptr<material> default_material, 
                           std::shared_ptr<Transform> ObjectToWorld, 
                           std::shared_ptr<Transform> WorldToObject, 
                           bool reverseOrientation) : nTriangles(0) {
  texture_size = 0;
  vertexIndices.clear();
  normalIndices.clear();
  texIndices.clear();
  face_material_id.clear();
  has_normals = false;
  has_tex = false;
  has_consistent_normals = false;
  
  
  nVertices = numVerts * 3;
  nNormals = normals ? numVerts * 3 : 0;
  nTex = texcoords ? numVerts * 2 : 0;
  p.reset(new point3f[nVertices]);
  for (size_t i = 0; i < nVertices; i += 3) {
    p[i / 3] = (*ObjectToWorld)(point3f((Float)vertices[i+0],
                                        (Float)vertices[i+1],
                                        (Float)vertices[i+2]));
  }
  
  if(nNormals > 0) {
    has_normals = true;
    n.reset(new normal3f[nNormals]);
    for (size_t i = 0; i < nNormals; i += 3) {
      n[i / 3] = (*ObjectToWorld)(normal3f((Float)normals[i+0],
                                           (Float)normals[i+1],
                                           (Float)normals[i+2]));
    }
  } else {
    n = nullptr;
  }
  
  if(nTex > 0) {
    has_tex = true;
    uv.reset(new point2f[nTex]);
    for (size_t i = 0; i < nTex; i += 2) {
      uv[i / 2] = point2f((Float)texcoords[i+0],
                      (Float)texcoords[i+1]);
    }
  } else {
    uv = nullptr;
  }
  
  nTriangles = 0;
  for (size_t s = 0; s < numIndices; s += 3) {
    vertexIndices.push_back(indices[s]);
    vertexIndices.push_back(indices[s+1]);
    vertexIndices.push_back(indices[s+2]);
    
    // if(has_normals) {
      normalIndices.push_back(indices[s+0]);
      normalIndices.push_back(indices[s+1]);
      normalIndices.push_back(indices[s+2]);
    // }
    // if(has_tex) {
      texIndices.push_back(indices[s+0]);
      texIndices.push_back(indices[s+1]);
      texIndices.push_back(indices[s+2]);
    // }
    nTriangles++;
    face_material_id.push_back(0);
  }
  
  //Material stuff
  mtl_materials.push_back(default_material);
  alpha_textures.push_back(alpha);
  bump_textures.push_back(bump);
}

size_t TriangleMesh::GetSize() {
  size_t size = sizeof(*this);
  size += nTex / 2 * sizeof(point2f) + 
          nNormals / 3 * sizeof(normal3f) + 
          nVertices / 3 * sizeof(point3f);
  for(size_t i = 0; i < mtl_materials.size(); i++) {
    size += mtl_materials[i]->GetSize();
  }
  size += face_material_id.size()*sizeof(int);
  size += sizeof(unsigned char *) * bump_texture_data.size();
  size += sizeof(unsigned char *) * obj_texture_data.size();
  size += sizeof(std::shared_ptr<alpha_texture>) * alpha_textures.size();
  size += sizeof(std::shared_ptr<bump_texture>)  * bump_textures.size();
  size += sizeof(int) * vertexIndices.size();
  size += sizeof(int) * normalIndices.size();
  size += sizeof(int) * texIndices.size();
  size += texture_size;
  return(size);
}
