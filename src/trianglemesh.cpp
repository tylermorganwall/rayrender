#include "trianglemesh.h"

#ifndef TINYOBJLOADER_IMPLEMENTATION
#define TINYOBJLOADER_IMPLEMENTATION
#define TINYOBJLOADER_USE_MAPBOX_EARCUT
#include "tinyobj/tiny_obj_loader.h"
#endif

#include <queue>
#include "texturecache.h"

TriangleMesh::~TriangleMesh() {}


inline Float luminance(point3f& color) {
  return(dot(color,point3f(0.2125,0.7154,0.0721)));
}

void LoadRayMaterials(std::vector<std::shared_ptr<material> > &mesh_materials,
                      std::vector<Rcpp::List > &shape_materials,
                      std::vector<unsigned char * > &obj_texture_data,
                      std::vector<unsigned char * > &bump_texture_data,
                      std::vector<std::shared_ptr<bump_texture> > &bump_textures,
                      std::vector<std::shared_ptr<alpha_texture> > &alpha_textures,
                      std::shared_ptr<alpha_texture> alpha_default,
                      std::shared_ptr<bump_texture> bump_default,
                      size_t &texture_size,
                      std::shared_ptr<material> default_material,
                      bool override_material, bool flip_transmittance,
                      TextureCache &texCache,
                      bool verbose, std::vector<bool>& material_is_light) {
  int total_materials = 0;
  for(size_t i = 0; i < shape_materials.size(); i++) { 
    Rcpp::List single_material_batch = Rcpp::as<Rcpp::List>(shape_materials[0]);
    total_materials += single_material_batch.size();
  }
  //Need to ensure textures have full paths
  mesh_materials.reserve(total_materials+1);
  obj_texture_data.reserve(total_materials+1);
  bump_texture_data.reserve(total_materials+1);
  bump_textures.reserve(total_materials+1);
  alpha_textures.reserve(total_materials+1);
  material_is_light.reserve(total_materials+1);
  
  //For default texture
  if(override_material) {
    mesh_materials.push_back(default_material);
    //TODO: Do I need to check if the default material is a light and is importance sampled?
    material_is_light.push_back(false);
    alpha_textures.push_back(alpha_default);
    bump_textures.push_back(bump_default);
    return;
  } 
  
  std::vector<point3f > diffuse_materials(total_materials+1);
  std::vector<point3f > specular_materials(total_materials+1);
  std::vector<Float > ior_materials(total_materials+1);
  std::vector<bool > has_diffuse(total_materials+1, false);
  std::vector<bool > has_transparency(total_materials+1, false);
  std::vector<bool > ior(total_materials+1, 1.0);
  
  std::vector<bool > has_single_diffuse(total_materials+1, false);
  std::vector<bool > has_alpha(total_materials+1, false);
  std::vector<bool > has_bump(total_materials+1, false);
  std::vector<Float > bump_intensity(total_materials+1);
  
  std::vector<int > nx_mat(total_materials+1);
  std::vector<int > ny_mat(total_materials+1);
  std::vector<int > nn_mat(total_materials+1);
  
  std::vector<int > nx_mat_bump(total_materials+1);
  std::vector<int > ny_mat_bump(total_materials+1);
  std::vector<int > nn_mat_bump(total_materials+1);
  // int nx, ny, nn;
  
  int mat_num = 0;
  for(size_t ii = 0; ii < shape_materials.size(); ii++) {
    Rcpp::List materials = Rcpp::as<Rcpp::List>(shape_materials[ii]);
    for (size_t i = 0; i < (size_t)materials.size(); i++) {
      Rcpp::List single_material = Rcpp::as<Rcpp::List>(materials(i));
      std::string diffuse_texname = Rcpp::as<std::string>(single_material["diffuse_texname"]);
      std::string bump_texname = Rcpp::as<std::string>(single_material["bump_texname"]);
      std::string emissive_texname = Rcpp::as<std::string>(single_material["emissive_texname"]);
      std::string ambient_texname = Rcpp::as<std::string>(single_material["ambient_texname"]);
      std::string specular_texname = Rcpp::as<std::string>(single_material["specular_texname"]);
      std::string normal_texname = Rcpp::as<std::string>(single_material["normal_texname"]);
      Rcpp::NumericVector diffuse = Rcpp::as<Rcpp::NumericVector>(single_material["diffuse"]);
      Rcpp::NumericVector ambient = Rcpp::as<Rcpp::NumericVector>(single_material["ambient"]);
      Float bump_intensity_single = Rcpp::as<Float>(single_material["bump_intensity"]);
      Float dissolve = Rcpp::as<Float>(single_material["dissolve"]);
      Float ior = Rcpp::as<Float>(single_material["ior"]);

      if(strlen(diffuse_texname.c_str()) > 0) {
        int nx = 0; 
        int ny = 0; 
        int nn = 0;
        std::replace(diffuse_texname.begin(),diffuse_texname.end(), '\\', separator());
        obj_texture_data.push_back(texCache.LookupChar(diffuse_texname, nx, ny, nn, 4));
        if(nx == 0 || ny == 0 || nn == 0) {
          throw std::runtime_error("Could not find " + diffuse_texname);
        }

        if(verbose) {
          Rprintf("(%i/%i) Loading Material Texture %s (%i/%i/%i) \n", 
                  (int)i+1, (int)materials.size(),diffuse_texname.c_str(),nx,ny,nn);
        }
        
        texture_size += sizeof(unsigned char) * nx * ny * nn;
        has_diffuse[mat_num] = true;
        has_single_diffuse[mat_num] = false;
        
        has_alpha[mat_num] = false;
        if(nn == 4) {
          for(int j = 0; j < nx - 1; j++) {
            for(int k = 0; k < ny - 1; k++) {
              if(obj_texture_data[mat_num][4*j + 4*nx*k + 3] != 255) {
                has_alpha[mat_num] = true;
                break;
              }
            }
            if(has_alpha[mat_num]) {
              break;
            }
          }
        }
        //Always use 4 channels for RGB
        nn = 4;
        nx_mat[mat_num] = nx;
        ny_mat[mat_num] = ny;
        nn_mat[mat_num] = nn;
        if(has_alpha[mat_num]) {
          if(verbose) {
            Rprintf("(%i/%i) Material has alpha channel \n", (int)i+1, (int)materials.size());
          }
        }
      } else if (diffuse.size() == 3 && dissolve == 1) {
        obj_texture_data.push_back(nullptr);
        diffuse_materials[mat_num] = point3f(diffuse[0],diffuse[1],diffuse[2]);
        has_diffuse[mat_num] = true;
        has_alpha[mat_num] = false;
        has_single_diffuse[mat_num] = true;
      } else if(dissolve < 1) {
        obj_texture_data.push_back(nullptr);
        specular_materials[mat_num] = point3f(diffuse[0],diffuse[1],diffuse[2]);
        ior_materials[mat_num] = ior;
        has_alpha[mat_num] = false;
        has_transparency[mat_num] = true; 
      } else {
        obj_texture_data.push_back(nullptr);
        has_diffuse[mat_num] = false;
        has_alpha[mat_num] = false;
        has_single_diffuse[mat_num] = false;
      }
      if(strlen(bump_texname.c_str()) > 0) {
        int nxb = 0; 
        int nyb = 0; 
        int nnb = 0;
        std::replace(bump_texname.begin(),bump_texname.end(), '\\', separator());
  
        bump_texture_data[mat_num] = texCache.LookupChar(bump_texname, nxb, nyb, nnb, 1);
        // nn = 4;
        texture_size += sizeof(unsigned char) * nxb * nyb * nnb;
        if(nxb == 0 || nyb == 0 || nnb == 0) {
          throw std::runtime_error("Could not find " + bump_texname);
        }
        nx_mat_bump[mat_num] = nxb;
        ny_mat_bump[mat_num] = nyb;
        nn_mat_bump[mat_num] = 1;
        bump_intensity[mat_num] = bump_intensity_single;
        has_bump[mat_num] = true;
      } else {
        bump_texture_data.push_back(nullptr);
        bump_intensity[mat_num] = 1.0f;
        has_bump[mat_num] = false;
      }
      std::shared_ptr<alpha_texture> alpha = nullptr;
      std::shared_ptr<bump_texture> bump = nullptr;
      if(has_alpha[mat_num]) {
        alpha = std::make_shared<alpha_texture>(obj_texture_data[mat_num], 
                                                nx_mat[mat_num], ny_mat[mat_num], nn_mat[mat_num]);
      } 
  
      if(has_bump[mat_num]) {
        bump = std::make_shared<bump_texture>(bump_texture_data[mat_num],
                                              nx_mat_bump[mat_num], ny_mat_bump[mat_num], nn_mat_bump[mat_num],
                                              bump_intensity[mat_num]);
      }
      alpha_textures.push_back(alpha);
      bump_textures.push_back(bump);
      bool imp_sample_obj = false;
      std::shared_ptr<material> tex = nullptr;
      Rcpp::NumericVector emission = Rcpp::as<Rcpp::NumericVector>(single_material["emission"]);
      Rcpp::NumericVector specular = Rcpp::as<Rcpp::NumericVector>(single_material["specular"]);
      Rcpp::NumericVector transmittance = Rcpp::as<Rcpp::NumericVector>(single_material["transmittance"]);
      bool any_trans = transmittance(0) != 0 || transmittance(1) != 0 || transmittance(2) != 0;
      Float shininess = Rcpp::as<Float>(single_material["shininess"]);
      
      // int illum = Rcpp::as<int>(single_material["illum"]);
      
      if(dissolve == 1 && shininess != 1000 && !any_trans) {
        point3f ke(emission(0),
                   emission(1),
                   emission(2));
        if(ke.x() != 0 || ke.y() != 0 || ke.z() != 0) {
          tex = std::make_shared<diffuse_light>(std::make_shared<constant_texture>(ke), 1.0, false);
          imp_sample_obj = true;
        } else {
          if(has_diffuse[mat_num]) {
            if(has_single_diffuse[mat_num]) {
              tex = std::make_shared<lambertian>(std::make_shared<constant_texture>(diffuse_materials[mat_num]));
            } else {
              tex = std::make_shared<lambertian>(std::make_shared<image_texture_char>(obj_texture_data[mat_num],
                                                                                      nx_mat[mat_num], 
                                                                                      ny_mat[mat_num],
                                                                                      nn_mat[mat_num]));
            }
          } else {
            tex = default_material;
          }
        }
      } else {
        point3f spec = point3f(specular(0),
                               specular(1),
                               specular(2));
        if(shininess == 1000) {
          tex = std::make_shared<metal>(std::make_shared<constant_texture>(spec),
                                        0., 
                                        point3f(0), 
                                        point3f(0));
        } else {
          point3f atten;
          if(flip_transmittance) {
            atten = point3f(1.f-transmittance(0),
                            1.f-transmittance(1),
                            1.f-transmittance(2));
          } else {
            atten = point3f(transmittance(0),
                            transmittance(1),
                            transmittance(2));
          }
          tex = std::make_shared<dielectric>(spec, 
                                             ior, atten, 
                                             0);
        }
      }
      // if(verbose) {
      //   Rprintf("(%i/%i) Loading Material %s (Imp Sample: %s) \n", material_num+1, 
      //           materials.size(),materials[material_num].name.c_str(), imp_sample_obj ? "true" : "false");
      // }
      mesh_materials.push_back(tex);
      material_is_light.push_back(imp_sample_obj);
      mat_num++;
    }
  }
}

void LoadMtlMaterials(std::vector<std::shared_ptr<material> > &mesh_materials,
                      std::vector<tinyobj::material_t > &materials,
                      std::vector<unsigned char * > &obj_texture_data,
                      std::vector<unsigned char * > &bump_texture_data,
                      std::vector<std::shared_ptr<alpha_texture> > &alpha_textures,
                      std::vector<std::shared_ptr<bump_texture> > &bump_textures,
                      std::shared_ptr<alpha_texture> alpha_default,
                      std::shared_ptr<bump_texture> bump_default,
                      size_t &texture_size,
                      const std::string inputfile, const std::string basedir, bool has_sep,
                      std::shared_ptr<material> default_material, bool load_materials,
                      bool load_textures, 
                      TextureCache &texCache, bool verbose, std::vector<bool>& material_is_light) {
  mesh_materials.reserve(materials.size()+1);
  obj_texture_data.reserve(materials.size()+1);
  bump_texture_data.reserve(materials.size()+1);
  bump_textures.reserve(materials.size()+1);
  alpha_textures.reserve(materials.size()+1);
  material_is_light.reserve(materials.size()+1);

  
  std::vector<point3f > diffuse_materials(materials.size()+1);
  std::vector<point3f > specular_materials(materials.size()+1);
  std::vector<Float > ior_materials(materials.size()+1);
  std::vector<bool > has_diffuse_texture(materials.size()+1, false);
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

  if(load_materials) {
    //For default texture
    alpha_textures.push_back(alpha_default);
    bump_textures.push_back(bump_default);
    for (size_t i = 0; i < materials.size(); i++) {
      int nx = 0; 
      int ny = 0; 
      int nn = 0;
      if(strlen(materials[i].diffuse_texname.c_str()) > 0 && load_textures) {
        std::replace(materials[i].diffuse_texname.begin(), materials[i].diffuse_texname.end(), '\\', separator());
        if(has_sep) {
          obj_texture_data.push_back(texCache.LookupChar((basedir + 
            separator() + 
            materials[i].diffuse_texname), nx, ny, nn, 4));
          nn = 4;
        } else {
          obj_texture_data.push_back(texCache.LookupChar(materials[i].diffuse_texname, 
                                                         nx, ny, nn, 4));
          nn = 4;
        }
        
        if(nx == 0 || ny == 0) {
          if(has_sep) {
            throw std::runtime_error("Could not find " + (basedir + separator() + materials[i].diffuse_texname));
          } else {
            throw std::runtime_error("Could not find " + materials[i].diffuse_texname);
          }
        }
        if(verbose) {
          Rprintf("(%i/%i) Loading Material Texture %s (%i/%i/%i) \n", 
                  (int)i+1, (int)materials.size(),materials[i].name.c_str(),nx,ny,nn);
        }
        
        texture_size += sizeof(unsigned char) * nx * ny * nn;
        has_diffuse_texture[i] = true;
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
        diffuse_materials[i] = point3f(materials[i].diffuse[0],materials[i].diffuse[1],materials[i].diffuse[2]);
        has_diffuse_texture[i] = false;
        has_alpha[i] = false;
        has_single_diffuse[i] = true;
      } else if(materials[i].dissolve < 1) {
        obj_texture_data.push_back(nullptr);
        specular_materials[i] = point3f(materials[i].diffuse[0],materials[i].diffuse[1],materials[i].diffuse[2]);
        ior_materials[i] = materials[i].ior;
        has_alpha[i] = false;
        has_transparency[i] = true; 
      } else {
        obj_texture_data.push_back(nullptr);
        has_diffuse_texture[i] = false;
        has_alpha[i] = false;
        has_single_diffuse[i] = false;
      }
      if(strlen(materials[i].bump_texname.c_str()) > 0 && load_textures) {
        std::replace(materials[i].bump_texname.begin(), materials[i].bump_texname.end(), '\\', separator());
        
        if(has_sep) {
          bump_texture_data[i] = texCache.LookupChar(basedir + separator() + materials[i].bump_texname, 
                                                     nx, ny, nn, 4);
          nn = 4;
        } else {
          bump_texture_data[i] = texCache.LookupChar(materials[i].bump_texname, nx, ny, nn, 4);
          nn = 4;
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
    alpha_textures.push_back(alpha_default);
    bump_textures.push_back(bump_default);
  }
  
  //First texture is default (when shapes[s].mesh.material_ids[f] == -1)
  mesh_materials.push_back(default_material);
  material_is_light.push_back(false);
  
  if(load_materials) {
    for(size_t material_num = 0; material_num < materials.size(); material_num++) {
      bool imp_sample_obj = false;
      std::shared_ptr<material> tex = nullptr;
      // int illum = materials[material_num].illum;
      float shininess = materials[material_num].shininess;
      point3f spec = point3f(materials[material_num].specular[0],
                             materials[material_num].specular[1],
                             materials[material_num].specular[2]);
      float dissolve = materials[material_num].dissolve;
      
      point3f ke(materials[material_num].emission[0],
                 materials[material_num].emission[1],
                 materials[material_num].emission[2]);
      point3f atten = point3f(materials[material_num].transmittance[0],
                              materials[material_num].transmittance[1],
                              materials[material_num].transmittance[2]);
      bool any_transparent = dissolve < 1 && !has_diffuse_texture[material_num];
      bool is_metallic = shininess == 1000;
      bool is_glossy = shininess > 128;
      
      if(ke.x() != 0 || ke.y() != 0 || ke.z() != 0) { 
        //Any emitting material will be a light
        tex = std::make_shared<diffuse_light>(std::make_shared<constant_texture>(ke), 1.0, false);
        imp_sample_obj = true;
      } else if (any_transparent) {
        tex = std::make_shared<dielectric>(spec, 
                                           materials[material_num].ior, atten, 
                                           0);
      } else if (is_metallic) {
        tex = std::make_shared<metal>(std::make_shared<constant_texture>(spec),
                                      0., 
                                      point3f(0), 
                                      point3f(0));
      } else if (is_glossy) {
        if(has_diffuse_texture[material_num]) {
          float inv_shininess = 0.001;
          
          MicrofacetDistribution *dist = new TrowbridgeReitzDistribution(inv_shininess, 
                                                                         inv_shininess,
                                                                         nullptr, 
                                                                         false,  
                                                                         true);
          point3f Rd(0.f,0.f,0.f);
          point3f Rs(1.f,1.f,1.f);
          
          tex = std::make_shared<glossy>(std::make_shared<image_texture_char>(obj_texture_data[material_num],
                                                                              nx_mat[material_num], 
                                                                              ny_mat[material_num],
                                                                              nn_mat[material_num]),
                                                                              dist, Rd, Rs);
        } else if (has_single_diffuse[material_num]) {
          float inv_shininess = 0.001;
          MicrofacetDistribution *dist = new TrowbridgeReitzDistribution(inv_shininess, 
                                                                         inv_shininess,
                                                                         nullptr, 
                                                                         false,  
                                                                         true);
          point3f Rd = diffuse_materials[material_num];
          point3f Rs(1.f,1.f,1.f);
          
          tex = std::make_shared<glossy>(std::make_shared<constant_texture>(point3f(0,0,0)), 
                                         dist, Rd, Rs);
        } else {
          tex = default_material;
        }
      } else {
        if(has_diffuse_texture[material_num]) {
          tex = std::make_shared<lambertian>(std::make_shared<image_texture_char>(obj_texture_data[material_num],
                                                                                  nx_mat[material_num], 
                                                                                  ny_mat[material_num],
                                                                                  nn_mat[material_num]));
        } else if (has_single_diffuse[material_num]) {
          tex = std::make_shared<lambertian>(std::make_shared<constant_texture>(diffuse_materials[material_num]));
        } else {
          tex = default_material;
        }
      }
      if(verbose) {
        Rprintf("(%i/%i) Loading Material %s as %s (Imp Sample: %s) \n", 
                (int)material_num+1, 
                (int)materials.size(),
                materials[material_num].name.c_str(), 
                tex->GetName().c_str(),
                imp_sample_obj ? "true" : "false");
      }
      mesh_materials.push_back(tex);
      material_is_light.push_back(imp_sample_obj);
    }
  }
}

TriangleMesh::TriangleMesh(std::string inputfile, std::string basedir,
                           std::shared_ptr<material> default_material, 
                           std::shared_ptr<alpha_texture> alpha_default,
                           std::shared_ptr<bump_texture> bump_default,
                           bool load_materials, bool load_textures, bool load_vertex_colors, 
                           bool load_normals, bool verbose, Float scale, 
                           bool calculate_consistent_normals,
                           TextureCache& texCache,
                           Transform* ObjectToWorld, 
                           Transform* WorldToObject, 
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

  has_vertex_colors = attrib.colors.size() > 0 ? true : false;
  has_tex = attrib.texcoords.size() > 0 ? true : false;
  
  if(strlen(basedir.c_str()) == 0) {
    has_sep = false;
  }
  if(ret) {
    nVertices = attrib.vertices.size() / 3;
    for (size_t s = 0; s < shapes.size(); s++) {
      nTriangles += shapes[s].mesh.indices.size() / 3;
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
    nNormals = load_normals ? attrib.normals.size() / 3 : 0;
    nTex = !has_vertex_colors ? attrib.texcoords.size() / 2 : 0;
    p.reset(new point3f[nVertices]);
    for (size_t i = 0; i < nVertices * 3; i += 3) {
      p[i / 3] = (*ObjectToWorld)(point3f(attrib.vertices[i+0],
                                          attrib.vertices[i+1],
                                          attrib.vertices[i+2]) * scale);
    }
    
    
    if(nNormals > 0) {
      has_normals = true;
      n.reset(new normal3f[nNormals]);
      for (size_t i = 0; i < nNormals * 3; i += 3) {
        n[i / 3] = (*ObjectToWorld)(normal3f(attrib.normals[i+0],
                                             attrib.normals[i+1],
                                             attrib.normals[i+2]));
      }
      if(has_consistent_normals) {
        face_n.reset(new normal3f[normalIndices.size()]);
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
        for(size_t i = 0; i < alpha_v.size(); i++) {
          Float temp_av = clamp(alpha_v[i],-1,1);
          alpha_v[i] = std::acos(temp_av) * (1 + 0.03632 * (1 - temp_av) * (1 - temp_av));
        }
      }
    } else {
      n = nullptr;
    }

    
    if(nTex > 0) {
      has_tex = true;
      uv.reset(new point2f[nTex]);
      for (size_t i = 0; i < nTex * 2; i += 2) {
        uv[i / 2] = point2f(attrib.texcoords[i+0],
                            attrib.texcoords[i+1]);      }
    } else {
      uv = nullptr;
    }
    
    if(has_vertex_colors) {
      vc.reset(new point3f[nVertices]);
      for (size_t i = 0; i < nVertices * 3; i += 3) {
        vc[i / 3] = point3f(attrib.colors[i+0],
                            attrib.colors[i+1],
                            attrib.colors[i+2]);
      }
    } else {
      vc = nullptr;
    }
    
    if(!has_vertex_colors) {
      LoadMtlMaterials(mesh_materials, materials, obj_texture_data,
                       bump_texture_data, alpha_textures, bump_textures, 
                       alpha_default, bump_default,
                       texture_size, inputfile, basedir, has_sep, default_material,
                       load_materials, load_textures, 
                       texCache, 
                       verbose, material_is_light);
    } else {
      mesh_materials.push_back(default_material);
      material_is_light.push_back(false);
      alpha_textures.push_back(nullptr);
      bump_textures.push_back(nullptr);
      for (size_t s = 0; s < vertexIndices.size(); s += 3) {
        std::shared_ptr<texture> tex = std::shared_ptr<triangle_texture>(
          new triangle_texture(vc[vertexIndices[s]],
                               vc[vertexIndices[s+1]],
                               vc[vertexIndices[s+2]]));
        mesh_materials.push_back(std::shared_ptr<material>(new lambertian(tex)));
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
                           TextureCache& texCache,
                           Transform* ObjectToWorld, 
                           Transform* WorldToObject, 
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
  size_t nVertexColors = vertexcolors.nrow();
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
      for(size_t i = 0; i < alpha_v.size(); i++) {
        Float temp_av = clamp(alpha_v[i],-1,1);
        alpha_v[i] = std::acos(temp_av) * (1 + 0.03632 * (1 - temp_av) * (1 - temp_av));
      }
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
  for (size_t s = 0; s < static_cast<size_t>(indices.nrow()); s++) {
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
    for (size_t s = 0; s < static_cast<size_t>(indices.nrow()); s++) {
      face_material_id.push_back(0);
    }
    mesh_materials.push_back(default_material);
    if(mesh_texture_data) {
      obj_texture_data.push_back(mesh_texture_data);
    }
    if(bump_texture_data_) {
      bump_texture_data.push_back(bump_texture_data_);
    }
    alpha_textures.push_back(alpha);
    bump_textures.push_back(bump);
  } else {
    mesh_materials.push_back(default_material);
    alpha_textures.push_back(nullptr);
    bump_textures.push_back(nullptr);
    for (size_t s = 0; s < vertexIndices.size(); s += 3) {
      std::shared_ptr<texture> tex = std::shared_ptr<triangle_texture>(
        new triangle_texture(vc[vertexIndices[s]],
                             vc[vertexIndices[s+1]],
                             vc[vertexIndices[s+2]]));
      mesh_materials.push_back(std::shared_ptr<material>(new lambertian(tex)));
      face_material_id.push_back(s / 3 + 1);
      alpha_textures.push_back(nullptr);
      bump_textures.push_back(nullptr);
    }
  }
}


TriangleMesh::TriangleMesh(Rcpp::List raymesh, bool verbose, bool calculate_consistent_normals,
                           bool override_material, bool flip_transmittance,
                           std::shared_ptr<alpha_texture> alpha,
                           std::shared_ptr<bump_texture> bump,
                           TextureCache& texCache,
                           std::shared_ptr<material> default_material, 
                           Transform* ObjectToWorld, 
                           Transform* WorldToObject, 
                           bool reverseOrientation) : nTriangles(0) {
  Rcpp::List shape_container = Rcpp::as<Rcpp::List>(raymesh["shapes"]);
  has_vertex_colors = false;
  Rcpp::List vertex_raw = raymesh["vertices"];
  Rcpp::List normals_raw = raymesh["normals"];
  Rcpp::List tex_raw = raymesh["texcoords"];
  size_t number_shapes = shape_container.size();
  
  Rcpp::List materials_raw = Rcpp::as<Rcpp::List>(raymesh["materials"]);
  std::vector<Rcpp::List> materials;
  
  for(size_t i = 0; i < (size_t)materials_raw.size(); i++) {
    materials.push_back(Rcpp::as<Rcpp::List>(materials_raw(i)));
  }
  
  texture_size = 0;
  vertexIndices.clear();
  normalIndices.clear();
  nTriangles = 0;
  
  texIndices.clear();
  face_material_id.clear();
  has_normals = false;
  has_tex = false;
  has_consistent_normals = calculate_consistent_normals;
  nVertices = 0;
  nNormals = 0;
  nTex = 0;
  
  std::vector<Rcpp::NumericMatrix> vertices_shapes;
  std::vector<Rcpp::NumericMatrix> normals_shapes;
  std::vector<Rcpp::NumericMatrix> texcoords_shapes;
  
  vc = nullptr;
  
  for(size_t i = 0; i < number_shapes; i++) {
    vertices_shapes.push_back(Rcpp::as<Rcpp::NumericMatrix>(vertex_raw(i))); 
    normals_shapes.push_back(Rcpp::as<Rcpp::NumericMatrix>(normals_raw(i))); 
    texcoords_shapes.push_back(Rcpp::as<Rcpp::NumericMatrix>(tex_raw(i))); 
    nVertices += vertices_shapes[i].nrow();
    nNormals += normals_shapes[i].nrow();
    nTex += texcoords_shapes[i].nrow();
  }
  p.reset(new point3f[nVertices]);
  if(nNormals > 0) {
    n.reset(new normal3f[nNormals]);
    has_normals = true;
  } else {
    n = nullptr;
  }
  if(nTex > 0) {
    uv.reset(new point2f[nTex]);
    has_tex = true;
  } else {
    uv = nullptr;
  }
  size_t loaded_verts = 0;
  size_t loaded_norms = 0;
  size_t loaded_tex = 0;
  size_t max_mat_id = 0;
  
  bool any_normal_missing = !has_normals;
  
  for(size_t j = 0; j < number_shapes; j++) {
    Rcpp::List shape = Rcpp::as<Rcpp::List>(shape_container(j));
    Rcpp::NumericMatrix vertices = vertices_shapes[j];
    Rcpp::NumericMatrix normals = normals_shapes[j];
    Rcpp::NumericMatrix texcoords = texcoords_shapes[j];
    
    Rcpp::IntegerVector mat_ids = Rcpp::as<Rcpp::IntegerVector>(shape["material_ids"]); 
    
    Rcpp::IntegerMatrix indices = Rcpp::as<Rcpp::IntegerMatrix>(shape["indices"]); 
    Rcpp::IntegerMatrix tex_indices = Rcpp::as<Rcpp::IntegerMatrix>(shape["tex_indices"]); 
    Rcpp::IntegerMatrix norm_indices = Rcpp::as<Rcpp::IntegerMatrix>(shape["norm_indices"]); 
    Rcpp::LogicalVector has_vertex_tex = Rcpp::as<Rcpp::LogicalVector>(shape["has_vertex_tex"]); 
    Rcpp::LogicalVector has_vertex_normals = Rcpp::as<Rcpp::LogicalVector>(shape["has_vertex_normals"]); 
    
    size_t single_shape_verts = vertices.nrow();
    size_t single_shape_normals = normals.nrow();
    size_t single_shape_tex = texcoords.nrow();

    for (size_t ii = 0; ii < single_shape_verts; ii += 1) {
      p[ii + loaded_verts] = (*ObjectToWorld)(point3f(vertices(ii,0),
                                                      vertices(ii,1),
                                                      vertices(ii,2)));
    }

    if(single_shape_normals > 0) {
      for (size_t ii = 0; ii < single_shape_normals; ii++) {
        n[ii + loaded_norms] = (*ObjectToWorld)(normal3f(normals(ii,0),
                                                        normals(ii,1),
                                                        normals(ii,2)));
      }
    }

    if(single_shape_tex > 0) {
      for (size_t ii = 0; ii < single_shape_tex; ii++) {
        uv[ii + loaded_tex] = point2f(texcoords(ii,0),
                                      texcoords(ii,1));
      }
    }

    for (size_t s = 0; s < static_cast<size_t>(indices.nrow()); s++) {
      nTriangles++;
      vertexIndices.push_back(indices(s,0) + loaded_verts);
      vertexIndices.push_back(indices(s,1) + loaded_verts);
      vertexIndices.push_back(indices(s,2) + loaded_verts);
      if(has_normals && has_vertex_normals(s)) {
        normalIndices.push_back(norm_indices(s,0) + loaded_norms);
        normalIndices.push_back(norm_indices(s,1) + loaded_norms);
        normalIndices.push_back(norm_indices(s,2) + loaded_norms);
      } else {
        normalIndices.push_back(-1);
        normalIndices.push_back(-1);
        normalIndices.push_back(-1);
        any_normal_missing = true;
      }
      if(has_tex && has_vertex_tex(s)) {
        texIndices.push_back(tex_indices(s,0)+ loaded_tex);
        texIndices.push_back(tex_indices(s,1)+ loaded_tex);
        texIndices.push_back(tex_indices(s,2)+ loaded_tex);
      } else {
        texIndices.push_back(-1);
        texIndices.push_back(-1);
        texIndices.push_back(-1);
      }
    }
    loaded_verts += single_shape_verts;
    loaded_norms += single_shape_normals;
    loaded_tex += single_shape_tex;
    
    for (size_t s = 0; s < static_cast<size_t>(mat_ids.size()); s++) {
      //This deals with incrementing material_ids  so we don't have to flatten the raymesh in R
      int mat_id_tmp = !override_material ? mat_ids(s) : 0;
      face_material_id.push_back(mat_id_tmp + max_mat_id);
    }
    if(!override_material) {
      max_mat_id = *max_element(std::begin(face_material_id), std::end(face_material_id)) + 1;
    } else {
      max_mat_id = 0;
    }
  }
  if(has_consistent_normals && !any_normal_missing) {
    face_n.reset(new normal3f[normalIndices.size() / 3]);
    std::map<int, std::priority_queue<Float> > alpha_values;
    for (size_t ii = 0; ii < normalIndices.size(); ii += 3) {
      int idx_n1 = normalIndices[ii];
      int idx_n2 = normalIndices[ii+1];
      int idx_n3 = normalIndices[ii+2];
      
      normal3f n1 = unit_vector(n[idx_n1]);
      normal3f n2 = unit_vector(n[idx_n2]);
      normal3f n3 = unit_vector(n[idx_n3]);
      
      normal3f face_normal = unit_vector(n1 + n2 + n3);
      face_n[ii / 3] = face_normal;
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
    for(size_t ii = 0; ii < alpha_v.size(); ii++) {
      Float temp_av = clamp(alpha_v[ii],-1,1);
      alpha_v[ii] = std::acos(temp_av) * (1 + 0.03632 * (1 - temp_av) * (1 - temp_av));
    }
  }

  LoadRayMaterials(mesh_materials,
                   materials,
                   obj_texture_data,
                   bump_texture_data, 
                   bump_textures, alpha_textures,
                   alpha, bump,
                   texture_size, default_material,
                   override_material, flip_transmittance,
                   texCache,
                   verbose, material_is_light);
  
}

TriangleMesh::TriangleMesh(float* vertices, 
                           int* indices, 
                           float* normals, 
                           float* texcoords,
                           int numVerts, int numIndices, 
                           std::shared_ptr<alpha_texture> alpha,
                           std::shared_ptr<bump_texture> bump,
                           std::shared_ptr<material> default_material, 
                           Transform* ObjectToWorld, 
                           Transform* WorldToObject, 
                           bool reverseOrientation) : nTriangles(0) {
  has_vertex_colors = false;
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
  for (size_t s = 0; s < static_cast<size_t>(numIndices); s += 3) {
    vertexIndices.push_back(indices[s]);
    vertexIndices.push_back(indices[s+1]);
    vertexIndices.push_back(indices[s+2]);
    if(has_normals) {
      normalIndices.push_back(indices[s+0]);
      normalIndices.push_back(indices[s+1]);
      normalIndices.push_back(indices[s+2]);
    } else {
      normalIndices.push_back(-1);
      normalIndices.push_back(-1);
      normalIndices.push_back(-1);
    }
    if(has_tex) {
      texIndices.push_back(indices[s+0]);
      texIndices.push_back(indices[s+1]);
      texIndices.push_back(indices[s+2]);
    } else {
      texIndices.push_back(-1);
      texIndices.push_back(-1);
      texIndices.push_back(-1);
    }
    nTriangles++;
    face_material_id.push_back(0);
  }
  
  //Material stuff
  mesh_materials.push_back(default_material);
  alpha_textures.push_back(alpha);
  bump_textures.push_back(bump);
}

void TriangleMesh::ValidateMesh() {
  // 1. Check vertexIndices
  for (const int idx : vertexIndices) {
    if (idx < 0 || idx >= static_cast<int>(nVertices)) {
      throw std::runtime_error("Vertex index out of bounds");
    }
  }
  
  // 2. Check normalIndices
  if (has_normals) {
    for (const int idx : normalIndices) {
      if (idx < -1 || idx >= static_cast<int>(nNormals)) {
        throw std::runtime_error("Normal index out of bounds");
      }
    }
  }
  
  // 3. Check texIndices
  if (has_tex) {
    for (const int idx : texIndices) {
      if (idx < -1 || idx >= static_cast<int>(nTex)) {
        throw std::runtime_error("Texture index out of bounds");
      }
    }
  }
  
  // 4. Check for NaN or Inf in vertex data
  for (size_t i = 0; i < nVertices; ++i) {
    if (std::isnan(p[i].x()) || std::isnan(p[i].y()) || std::isnan(p[i].z()) ||
        std::isinf(p[i].x()) || std::isinf(p[i].y()) || std::isinf(p[i].z())) {
      throw std::runtime_error("Vertex data contains NaN or Inf values");
    }
  }

  // 5. Check for NaN or Inf in normal data
  if (has_normals) {
    for (size_t i = 0; i < nNormals; ++i) {
      if (std::isnan(n[i].x()) || std::isnan(n[i].y()) || std::isnan(n[i].z()) ||
          std::isinf(n[i].x()) || std::isinf(n[i].y()) || std::isinf(n[i].z())) {
        throw std::runtime_error("Normal data contains NaN or Inf values");
      }
    }
  }

  // 6. Check for NaN or Inf in texture coordinate data
  if (has_tex) {
    for (size_t i = 0; i < nTex; ++i) {
      if (std::isnan(uv[i].x()) || std::isnan(uv[i].y()) ||
          std::isinf(uv[i].x()) || std::isinf(uv[i].y())) {
        throw std::runtime_error("Texture coordinate data contains NaN or Inf values");
      }
    }
  }
}
 
size_t TriangleMesh::GetSize() {
  size_t size = sizeof(*this);
  size += nTex / 2 * sizeof(point2f) + 
          nNormals / 3 * sizeof(normal3f) + 
          nVertices / 3 * sizeof(point3f);
  for(size_t i = 0; i < mesh_materials.size(); i++) {
    size += mesh_materials[i]->GetSize();
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
