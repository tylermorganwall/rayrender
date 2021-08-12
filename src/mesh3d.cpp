#include "mesh3d.h"

mesh3d::mesh3d(Rcpp::List mesh_info, std::shared_ptr<material> mat, 
       Float shutteropen, Float shutterclose, int bvh_type, random_gen rng,
       std::shared_ptr<Transform> ObjectToWorld, std::shared_ptr<Transform> WorldToObject, bool reverseOrientation) :
  hitable(ObjectToWorld, WorldToObject, reverseOrientation) {
  Rcpp::NumericMatrix vertices = Rcpp::as<Rcpp::NumericMatrix>(mesh_info["vertices"]);
  Rcpp::IntegerMatrix indices = Rcpp::as<Rcpp::IntegerMatrix>(mesh_info["indices"]);
  Rcpp::NumericMatrix norms = Rcpp::as<Rcpp::NumericMatrix>(mesh_info["normals"]);
  Rcpp::NumericMatrix txcoord = Rcpp::as<Rcpp::NumericMatrix>(mesh_info["texcoords"]);
  float scale_mesh = Rcpp::as<float>(mesh_info["scale_mesh"]);
  
  std::string texture = Rcpp::as<std::string>(mesh_info["texture"]);
  
  Rcpp::NumericMatrix colors = Rcpp::as<Rcpp::NumericMatrix>(mesh_info["color_vals"]);
  int colortype = Rcpp::as<int>(mesh_info["color_type"]);
  
  int nx,ny,nn;
  bool has_texture = false;
  if(strlen(texture.c_str()) > 0) {
    mesh_materials = stbi_loadf(texture.c_str(), &nx, &ny, &nn, 0);
    has_texture = true;
  } else {
    mesh_materials = nullptr;
  }
  mat_ptr = mat;
  int number_faces = indices.nrow();
  bool has_normals = norms.nrow() > 1;
  bool has_texcoords = txcoord.nrow() > 1;
  std::shared_ptr<material> tex = nullptr;
  bool single_tex = colortype != 3;
  
  for (int i = 0; i < number_faces; i++) {
    bool tempnormal = false;
    vec3f tris[3];
    vec3f normals[3];

    vec2f tx[3];
    
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
      tx[0] = vec2f(txcoord(idx[0],0),txcoord(idx[0],1));
      tx[1] = vec2f(txcoord(idx[1],0),txcoord(idx[1],1));
      tx[2] = vec2f(txcoord(idx[2],0),txcoord(idx[2],1));
    }
    if(colortype == 1 && has_texcoords && has_texture) {
      tex = std::make_shared<lambertian>(std::make_shared<triangle_image_texture>(mesh_materials,
                                                      nx,ny,nn,
                                                      tx[0].x(),tx[0].y(),
                                                      tx[1].x(),tx[1].y(),
                                                      tx[2].x(),tx[2].y()));
    } else if(colortype == 2) {
      tex = std::make_shared<lambertian>(std::make_shared<triangle_texture>(
        vec3f(colors(i,0),colors(i,1),colors(i,2)),
        vec3f(colors(i,0),colors(i,1),colors(i,2)),
        vec3f(colors(i,0),colors(i,1),colors(i,2))));
    } else if(colortype == 4) {
      tex = std::make_shared<lambertian>(std::make_shared<triangle_texture>(
        vec3f(colors(idx[0],0),colors(idx[0],1),colors(idx[0],2)),
        vec3f(colors(idx[1],0),colors(idx[1],1),colors(idx[1],2)),
        vec3f(colors(idx[2],0),colors(idx[2],1),colors(idx[2],2))));
    } else {
      tex = mat_ptr;
    }
    
    if(has_normals) {
      triangles.add(std::make_shared<triangle>(tris[0],tris[1],tris[2],
                                       normals[0],normals[1],normals[2],
                                       single_tex,
                                       tex, nullptr,  nullptr, 
                                       ObjectToWorld, WorldToObject, reverseOrientation));
    } else {
      triangles.add(std::make_shared<triangle>(tris[0],tris[1],tris[2], 
                                       single_tex, tex, nullptr, nullptr, 
                                       ObjectToWorld, WorldToObject, reverseOrientation));
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
