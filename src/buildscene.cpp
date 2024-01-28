#include "buildscene.h"

Transform rotation_order_matrix(NumericVector temprotvec, NumericVector order_rotation) {
  Transform M;
  for(int i = 0; i < 3; i++) {
    if(order_rotation(i) == 1) {
      if(temprotvec(0) != 0) {
        M = RotateX(temprotvec(0)) * M ;
      }
    }
    if(order_rotation(i) == 2) {
      if(temprotvec(1) != 0) {
        M = RotateY(temprotvec(1)) * M;
      }
    }
    if(order_rotation(i) == 3) {
      if(temprotvec(2) != 0) {
        M = RotateZ(temprotvec(2)) * M;
      }
    }
  }
  return(M);
}

enum ShapeEnum {   
  SPHERE = 1,
  XY_RECT = 2, 
  XZ_RECT = 3,
  YZ_RECT = 4,
  BOX = 5, 
  OBJ = 6, 
  DISK = 7,
  CYLINDER = 8, 
  ELLIPSOID = 9, 
  CURVE = 10, 
  CSG_OBJECT = 11, 
  PLY = 12,
  MESH3D = 13,
  RAYMESH = 14
};

enum MaterialEnum {   
  DIFFUSE = 1,
  METAL = 2,
  DIELECTRIC = 3,
  ORENNAYER = 4, 
  LIGHT = 5, 
  MICROFACET = 6,
  GLOSSY = 7, 
  SPOTLIGHT = 8, 
  HAIR = 9, 
  MICROFACET_TRANSMISSION = 10
};


std::shared_ptr<material> LoadSingleMaterial(List SingleMaterial,
                                            Float* textures,
                                            int* nvec,
                                            unsigned char * rough_texture,  
                                            int* nvecr,
                                            NumericVector tricolorinfo) {
  MaterialEnum type = static_cast<MaterialEnum>(as<int>(SingleMaterial["type"]));
  NumericVector properties = as<NumericVector>(as<List>(SingleMaterial["properties"])(0));
  NumericVector checkercolor = as<NumericVector>(as<List>(SingleMaterial["checkercolor"])(0));
  bool ischeckered = checkercolor.size() == 4;
  
  NumericVector gradient_color = as<NumericVector>(as<List>(SingleMaterial["gradient_color"])(0));
  bool isgradient = gradient_color.size() == 3;
  
  bool gradient_transpose = as<bool>(SingleMaterial["gradient_transpose"]);
  bool is_world_gradient = as<bool>(SingleMaterial["world_gradient"]);
  NumericVector gradient_point_info = as<NumericVector>(as<List>(SingleMaterial["gradient_point_info"])(0));

  bool is_hsv_gradient = as<CharacterVector>(SingleMaterial["gradient_type"])(0) == "hsv";
  
  Float noise = as<Float>(SingleMaterial["noise"]);
  bool isnoise = noise != 0;
  Float noisephase = as<Float>(SingleMaterial["noisephase"]);
  Float noiseintensity = as<Float>(SingleMaterial["noiseintensity"]);
  NumericVector noisecolor = as<NumericVector>(as<List>(SingleMaterial["noisecolor"])(0));
  
  bool isimage = as<bool>(SingleMaterial["image"]);
  
  NumericVector image_repeat = as<NumericVector>(as<List>(SingleMaterial["image_repeat"])(0));
  LogicalVector alphaimage = as<LogicalVector>(SingleMaterial["alphaimage"]);
  
  Float lightintensity = as<Float>(SingleMaterial["lightintensity"]);
  Float sigma = as<Float>(SingleMaterial["sigma"]);
  NumericVector glossyinfo = as<NumericVector>(as<List>(SingleMaterial["glossyinfo"])(0));
  
  bool has_rough = as<bool>(SingleMaterial["roughness_texture"]);
  
  std::shared_ptr<material> mat = nullptr;
  std::shared_ptr<texture> material_texture;
  
  bool is_tri_color = tricolorinfo.size() == 9;
  
  std::shared_ptr<roughness_texture> roughness;
  if(has_rough) {
    roughness = std::make_shared<roughness_texture>(rough_texture, 
                                                    nvecr[0], nvecr[1], nvecr[2]);
  }
  
  if(isimage) {
    material_texture = std::make_shared<image_texture_float>(textures, nvec[0], nvec[1], nvec[2], 
                                                             image_repeat[0], image_repeat[1], 1.0);
  } else if (isnoise) {
    material_texture = std::make_shared<noise_texture>(noise,point3f(properties(0),properties(1),properties(2)),
                                                       point3f(noisecolor(0),noisecolor(1),noisecolor(2)),
                                                       noisephase, noiseintensity);
  } else if (ischeckered) {
    material_texture = std::make_shared<checker_texture>(std::make_shared<constant_texture>(point3f(checkercolor(0),checkercolor(1),checkercolor(2))),
                                                         std::make_shared<constant_texture>(point3f(properties(0),properties(1),properties(2))),
                                                         checkercolor(3));
  } else if (isgradient && !is_world_gradient) {
    material_texture = std::make_shared<gradient_texture>(point3f(properties(0),properties(1),properties(2)),
                                                          point3f(gradient_color(0),gradient_color(1),gradient_color(2)),
                                                          gradient_transpose, is_hsv_gradient);
  } else if (is_tri_color) {
    material_texture = std::make_shared<triangle_texture>(point3f(tricolorinfo(0),tricolorinfo(1),tricolorinfo(2)),
                                                          point3f(tricolorinfo(3),tricolorinfo(4),tricolorinfo(5)),
                                                          point3f(tricolorinfo(6),tricolorinfo(7),tricolorinfo(8)));
  } else if (is_world_gradient) {
    material_texture = std::make_shared<world_gradient_texture>(point3f(gradient_point_info(0),gradient_point_info(1),gradient_point_info(2)),
                                                                point3f(gradient_point_info(3),gradient_point_info(4),gradient_point_info(5)),
                                                                point3f(properties(0),properties(1),properties(2)),
                                                                point3f(gradient_color(0),gradient_color(1),gradient_color(2)),
                                                                is_hsv_gradient);
  } else {
    material_texture = std::make_shared<constant_texture>(point3f(properties(0),properties(1),properties(2)));
  }
  
  bool is_invisible;
  switch (type) {
    case DIFFUSE: {
      mat = std::make_shared<lambertian>(material_texture);
      break;
    }
    case METAL: {
      Float fuzz = properties(3);
      point3f eta(glossyinfo(3), glossyinfo(4), glossyinfo(5));
      point3f kappa(glossyinfo(6),glossyinfo(7),glossyinfo(8));
      mat = std::make_shared<metal>(material_texture, fuzz, eta, kappa);
      break;
    }
    case DIELECTRIC: {
      point3f tint(properties(0),properties(1),properties(2));
      Float ior = properties(3);
      point3f atten(properties(4),properties(5),properties(6));
      int priority = static_cast<int>(properties(7));
      mat = std::make_shared<dielectric>(tint, ior, atten, priority);
      break;
    }
    case ORENNAYER: {
      mat = std::make_shared<orennayar>(material_texture, sigma);
      break;
    }
    case LIGHT: {
      is_invisible = properties(3) == 1;
      mat = std::make_shared<diffuse_light>(material_texture, lightintensity, is_invisible);
      break;
    }
    case MICROFACET: {
      MicrofacetDistribution *dist;
      if(glossyinfo(0) == 1) {
        dist = new TrowbridgeReitzDistribution(glossyinfo(1), glossyinfo(2),roughness, has_rough,  true);
      } else {
        dist = new BeckmannDistribution(glossyinfo(1), glossyinfo(2),roughness, has_rough, true);
      }
      point3f eta(glossyinfo(3), glossyinfo(4), glossyinfo(5));
      point3f kappa(glossyinfo(6), glossyinfo(7), glossyinfo(8));
      mat = std::make_shared<MicrofacetReflection>(material_texture, dist, eta, kappa);
      break;
    }
    case GLOSSY: {
      MicrofacetDistribution *dist;
      if(glossyinfo(0) == 1) {
        dist = new TrowbridgeReitzDistribution(glossyinfo(1), glossyinfo(2),roughness, has_rough,  true);
      } else {
        dist = new BeckmannDistribution(glossyinfo(1), glossyinfo(2),roughness, has_rough, true);
      }
      point3f eta(glossyinfo(3), glossyinfo(4), glossyinfo(5));
      point3f kappa(glossyinfo(6),glossyinfo(7),glossyinfo(8));
      mat = std::make_shared<glossy>(material_texture, dist, eta, kappa);
      break;
    }
    case SPOTLIGHT: {
      is_invisible = properties(8) == 1;
      vec3f spotlight_direction(properties(3), properties(4), properties(5));
      Float spotlight_width = properties(6);
      Float spotlight_falloff = properties(7);
      mat = std::make_shared<spot_light>(material_texture, spotlight_direction, 
                                         spotlight_width, spotlight_falloff, lightintensity, is_invisible);
      break;
    }
    case HAIR: {
      point3f sigma_a(properties(0),properties(1),properties(2));
      Float eta = properties(3);
      Float beta_m = properties(4);
      Float beta_n = properties(5);
      Float alpha = properties(6);
      mat = std::make_shared<hair>(sigma_a, eta, beta_m, beta_n, alpha);
      break;
    }
    case MICROFACET_TRANSMISSION: {
      MicrofacetDistribution *dist;
      if(glossyinfo(0) == 1) {
        dist = new TrowbridgeReitzDistribution(glossyinfo(1), glossyinfo(2),roughness, has_rough,  true);
      } else {
        dist = new BeckmannDistribution(glossyinfo(1), glossyinfo(2),roughness, has_rough, true);
      }
      point3f eta(glossyinfo(3), glossyinfo(4), glossyinfo(5));
      point3f kappa(glossyinfo(6), glossyinfo(7), glossyinfo(8));
      mat = std::make_shared<MicrofacetTransmission>(material_texture, dist, eta, kappa);
      break;
    }
    default: {
      mat = std::make_shared<lambertian>(material_texture);
      break;
    }
  }
  return(mat);
}


std::shared_ptr<hitable> build_scene(List& scene,
                                     IntegerVector& shape,
                                     List& position_list,
                                     Float shutteropen, 
                                     Float shutterclose,
                                     std::vector<Float* >& textures, 
                                     std::vector<int* >& nvec,
                                     std::vector<unsigned char * >& alpha_textures, 
                                     std::vector<int* >& nveca,
                                     std::vector<unsigned char * >& bump_textures, 
                                     std::vector<int* >& nvecb,
                                     std::vector<unsigned char * >& roughness_textures,  
                                     std::vector<int* >& nvecr,
                                     std::vector<std::shared_ptr<material> >* shared_materials, 
                                     int bvh_type,
                                     TransformCache& transformCache, 
                                     hitable_list& imp_sample_objects,
                                     bool verbose,
                                     random_gen& rng) {
  hitable_list list;
  NumericMatrix IdentityMat(4,4);
  IdentityMat.fill_diag(1);
  
  List RayMaterials = as<List>(scene["material"]);
  List ShapeInfo = as<List>(scene["shape_info"]);
  List Transforms = as<List>(scene["transforms"]);
  List AnimationInfo = as<List>(scene["animation_info"]);
  size_t n = ShapeInfo.length();

  NumericVector x = position_list["xvec"];
  NumericVector y = position_list["yvec"];
  NumericVector z = position_list["zvec"];
  
  std::vector<std::shared_ptr<alpha_texture> > alpha(n);
  std::vector<std::shared_ptr<bump_texture> > bump(n);
  std::vector<std::shared_ptr<roughness_texture> > roughness(n);
  
  
  for(size_t i = 0; i < n; i++) {
    List SingleShape = ShapeInfo(i);
    List SingleMaterial = RayMaterials(i);
    List SingleTransform = Transforms(i);
    List SingleAnimation = AnimationInfo(i);
    MaterialEnum material_type = static_cast<MaterialEnum>(as<int>(SingleMaterial["type"]));
    
    /// Load material for shape
    NumericVector tricolorinfo = as<NumericVector>(as<List>(SingleShape["tricolorinfo"])(0));
    bool isvolume = as<bool>(SingleMaterial["fog"]);
    Float fog_density = as<Float>(SingleMaterial["fogdensity"]);
    Float importance_sample = as<bool>(SingleMaterial["implicit_sample"]);
    
    IntegerVector material_id_vec = as<IntegerVector>(SingleShape["material_id"]);
    bool is_shared_mat = Rcpp::IntegerVector::is_na(material_id_vec(0));
    int material_id = material_id_vec(0);
    std::shared_ptr<material> shape_material;
    
    if(is_shared_mat && shared_materials->size() > static_cast<size_t>(material_id - 1)) {
      shape_material = shared_materials->at(material_id-1);
    } else {
      shape_material = LoadSingleMaterial(SingleMaterial,
                                          textures[i],
                                          nvec[i],
                                          roughness_textures[i],  
                                          nvecr[i],
                                          tricolorinfo);
    }
    if(is_shared_mat && shared_materials->size() < static_cast<size_t>(material_id)) {
      shared_materials->push_back(shape_material);
    }
    
    /// Done loading material 
  
    /// Load group/animation transforms
    // List GroupTransformContainer = as<List>(SingleTransform["group_transform"])(0);

    NumericMatrix GroupTransformMat = as<NumericMatrix>(as<List>(SingleTransform["group_transform"])(0));
    NumericVector scales = as<NumericVector>(as<List>(SingleTransform["scale"])(0));
    NumericVector order_rotation = as<NumericVector>(as<List>(SingleTransform["order_rotation"])(0));
    NumericVector angle = as<NumericVector>(as<List>(SingleTransform["angle"])(0));
    bool is_grouped = GroupTransformMat.nrow() == 4;
                                              
    NumericMatrix StartTransformAnimationMat = as<NumericMatrix>(as<List>(SingleAnimation["start_transform_animation"])(0));
    NumericMatrix EndTransformAnimationMat   = as<NumericMatrix>(as<List>(SingleAnimation["end_transform_animation"])(0));
    Float start_time = as<Float>(SingleAnimation["start_time"]);
    Float end_time = as<Float>(SingleAnimation["end_time"]);
    bool is_animated = StartTransformAnimationMat.nrow() == 4 && EndTransformAnimationMat.nrow() == 4;
    
    // Get default values if not transformed/animated
    if(!is_grouped) {
      GroupTransformMat = IdentityMat;
    }
    if(!is_animated) {
      StartTransformAnimationMat = IdentityMat;
      EndTransformAnimationMat = IdentityMat;
    }
    
    //Generate texture
    std::shared_ptr<material> tex = nullptr;
    Float bump_intensity = as<Float>(SingleMaterial["bump_intensity"]);
    NumericVector image_repeat = as<NumericVector>(as<List>(SingleMaterial["image_repeat"])(0));
    bool has_alpha = as<bool>(SingleMaterial["alphaimage"]);
    bool has_bump = as<bool>(SingleMaterial["bump_texture"]);
    
    if(has_alpha) {
      alpha[i] = std::make_shared<alpha_texture>(alpha_textures[i], 
                                                 nveca[i][0], nveca[i][1], nveca[i][2]);
    }
    if(has_bump) {
      bump[i] = std::make_shared<bump_texture>(bump_textures[i], 
                                               nvecb[i][0], nvecb[i][1], nvecb[i][2], 
                                               bump_intensity, image_repeat[0], image_repeat[1]);
    }
    //Generate center vector
    ShapeEnum shape_type = static_cast<ShapeEnum>(shape(i));
    vec3f center = vec3f(x(i), y(i), z(i));
    
    Transform GroupTransform(GroupTransformMat);
    Transform TempM = GroupTransform * Translate(center) *
      rotation_order_matrix(angle, order_rotation) * 
      Scale(scales[0], scales[1], scales[2]);
    std::shared_ptr<Transform> ObjToWorld = transformCache.Lookup(TempM);
    std::shared_ptr<Transform> WorldToObj = transformCache.Lookup(TempM.GetInverseMatrix());
    
    
    Transform AnimationStartTransform(StartTransformAnimationMat);
    Transform AnimationEndTransform(EndTransformAnimationMat);
    std::shared_ptr<Transform> StartAnim = transformCache.Lookup(AnimationStartTransform);
    std::shared_ptr<Transform> EndAnim = transformCache.Lookup(AnimationEndTransform);
    
    AnimatedTransform Animate(StartAnim, start_time, EndAnim, end_time);
    
    List shape_properties = as<List>(SingleShape["shape_properties"]);
    bool is_flipped = as<bool>(SingleShape["flipped"]);
    NumericVector fog_color = as<NumericVector>(as<List>(SingleMaterial["properties"])(0));
    NumericVector noise_color = as<NumericVector>(as<List>(SingleMaterial["noisecolor"])(0));
    
    
    // New code
    std::shared_ptr<hitable> entry;
    switch(shape_type) {
      case SPHERE: {
        Float radius = as<Float>(shape_properties["radius"]);
        entry = std::make_shared<sphere>(radius, shape_material, alpha[i], bump[i],
                                         ObjToWorld, WorldToObj, is_flipped);
        if(isvolume) {
          entry = std::make_shared<constant_medium>(entry, fog_density, 
                                                    std::make_shared<constant_texture>(point3f(fog_color(0),fog_color(1),fog_color(2))));
        }
        if(is_animated) {
          entry = std::make_shared<AnimatedHitable>(entry, Animate);
        }
        list.add(entry);
        break;
      }
      case XY_RECT: {
        NumericVector widths =  as<NumericVector>(shape_properties["rectinfo"]);
        entry = std::make_shared<xy_rect>(-widths(0)/2, widths(0)/2,
                                          -widths(1)/2, widths(1)/2,
                                          0, shape_material, alpha[i], bump[i], 
                                          ObjToWorld,WorldToObj, is_flipped);
        if(is_animated) {
          entry = std::make_shared<AnimatedHitable>(entry, Animate);
        }
        list.add(entry);
        break;
      }
      case XZ_RECT: {
        NumericVector widths =  as<NumericVector>(shape_properties["rectinfo"]);
        entry = std::make_shared<xz_rect>(-widths(0)/2, widths(0)/2,
                                          -widths(1)/2, widths(1)/2,
                                          0, shape_material, alpha[i], bump[i], 
                                          ObjToWorld,WorldToObj, is_flipped);
        if(is_animated) {
          entry = std::make_shared<AnimatedHitable>(entry, Animate);
        }
        list.add(entry);
        break;
      }
      case YZ_RECT: {
        NumericVector widths =  as<NumericVector>(shape_properties["rectinfo"]);
        entry = std::make_shared<yz_rect>(-widths(0)/2, widths(0)/2,
                                          -widths(1)/2, widths(1)/2,
                                          0, shape_material, alpha[i], bump[i], 
                                          ObjToWorld,WorldToObj, is_flipped);
        if(is_animated) {
          entry = std::make_shared<AnimatedHitable>(entry, Animate);
        }
        list.add(entry);
        break;
      }
      case BOX: {
        NumericVector widths =  as<NumericVector>(shape_properties["boxinfo"]);
        
        entry = std::make_shared<box>(-vec3f(widths(0),widths(1),widths(2))/2, 
                                       vec3f(widths(0),widths(1),widths(2))/2, 
                                       shape_material, alpha[i], bump[i],
                                       ObjToWorld,WorldToObj, is_flipped);
        
        if(isvolume) {
          Float noise = as<Float>(SingleMaterial["noise"]);
          bool isnoise = noise != 0;
          if(!isnoise) {
            entry = std::make_shared<constant_medium>(entry, fog_density, 
                                                      std::make_shared<constant_texture>(point3f(fog_color(0),fog_color(1),fog_color(2))));
          } else {
            Float noisephase = as<Float>(SingleMaterial["noisephase"]);
            Float noiseintensity = as<Float>(SingleMaterial["noiseintensity"]);
            
            entry = std::make_shared<constant_medium>(entry, fog_density, 
                                                      std::make_shared<noise_texture>(noise,
                                                                                      point3f(fog_color(0),fog_color(1),fog_color(2)),
                                                                                      point3f(noise_color(0),noise_color(1),noise_color(2)),
                                                                                      noisephase, noiseintensity));
          }
        } 
        if(is_animated) {
          entry = std::make_shared<AnimatedHitable>(entry, Animate);
        }
        list.add(entry);
        break;
      }
      case OBJ: {
        std::string objfilename = Rcpp::as<std::string>(SingleShape["fileinfo"]);
        std::string objbasename =  Rcpp::as<std::string>(shape_properties["basename"]);
        Float scale_obj = as<Float>(shape_properties["scale_obj"]);
        bool load_textures = as<bool>(shape_properties["load_textures"]);
        bool load_material = as<bool>(shape_properties["load_material"]);
        bool load_vertex_colors = as<bool>(shape_properties["vertex_colors"]);
        bool importance_sample_lights = as<bool>(shape_properties["importance_sample_lights"]);
        bool load_normals = as<bool>(shape_properties["load_normals"]);
        bool calculate_consistent_normals = as<bool>(shape_properties["calculate_consistent_normals"]);
        Float sigma_obj = as<Float>(SingleMaterial["sigma"]);
        entry = std::make_shared<trimesh>(objfilename, objbasename, 
                                          scale_obj, sigma_obj, shape_material, 
                                          load_material, load_textures, load_vertex_colors,
                                          importance_sample_lights, load_normals, calculate_consistent_normals,
                                          imp_sample_objects,
                                          shutteropen, shutterclose, bvh_type, rng, verbose,
                                          ObjToWorld,WorldToObj, is_flipped);
        if(isvolume) {
          entry = std::make_shared<constant_medium>(entry, fog_density, 
                                                    std::make_shared<constant_texture>(point3f(fog_color(0),fog_color(1),fog_color(2))));
        }
        if(is_animated) {
          entry = std::make_shared<AnimatedHitable>(entry, Animate);
        }
        list.add(entry);
        break;
      }
      case DISK: {
        Float radius = as<Float>(shape_properties["radius"]);
        Float inner_radius = as<Float>(shape_properties["inner_radius"]);
        
        entry = std::make_shared<disk>(vec3f(0,0,0), radius, inner_radius, shape_material, alpha[i], bump[i],
                                       ObjToWorld,WorldToObj, is_flipped);
        if(is_animated) {
          entry = std::make_shared<AnimatedHitable>(entry, Animate);
        }
        list.add(entry);
        break;
      }
      case CYLINDER: {
        Float radius =  Rcpp::as<Float>(shape_properties["radius"]);
        Float length =  Rcpp::as<Float>(shape_properties["length"]);
        Float phi_min =  Rcpp::as<Float>(shape_properties["phi_min"]);
        Float phi_max =  Rcpp::as<Float>(shape_properties["phi_max"]);
        bool has_cap =  Rcpp::as<bool>(shape_properties["has_cap"]);
        
        bool has_caps_not_light = material_type != LIGHT && material_type != SPOTLIGHT;
        bool has_cap_option = has_cap && has_caps_not_light;
        entry = std::make_shared<cylinder>(radius, length, 
                                           phi_min, phi_max, has_cap_option,
                                           shape_material, alpha[i], bump[i],
                                           ObjToWorld,WorldToObj, is_flipped);
        if(is_animated) {
          entry = std::make_shared<AnimatedHitable>(entry, Animate);
        }
        list.add(entry);
        break;
      }
      case ELLIPSOID: {
        Float a =  Rcpp::as<Float>(shape_properties["a"]);
        Float b =  Rcpp::as<Float>(shape_properties["b"]);
        Float c =  Rcpp::as<Float>(shape_properties["c"]);
        entry = std::make_shared<ellipsoid>(vec3f(0,0,0), 1, 
                                            vec3f(a,b,c),
                                            shape_material, alpha[i], bump[i],
                                            ObjToWorld,WorldToObj, is_flipped);
        if(isvolume) {
          entry = std::make_shared<constant_medium>(entry, fog_density, 
                                                    std::make_shared<constant_texture>(point3f(fog_color(0),fog_color(1),fog_color(2))));
        }
        if(is_animated) {
          entry = std::make_shared<AnimatedHitable>(entry, Animate);
        }
        list.add(entry);
        break;
      }
      case CURVE: {
        NumericVector p1 = as<NumericVector>(shape_properties["p1"]);
        NumericVector p2 = as<NumericVector>(shape_properties["p2"]);
        NumericVector p3 = as<NumericVector>(shape_properties["p3"]);
        NumericVector p4 = as<NumericVector>(shape_properties["p4"]);
        Float width      = as<Float>(shape_properties["width"]);
        Float width_end  = as<Float>(shape_properties["width_end"]);
        Float u_min      = as<Float>(shape_properties["u_min"]);
        Float u_max      = as<Float>(shape_properties["u_max"]);
        CurveType curvetype    = static_cast<CurveType>(Rcpp::as<int>(shape_properties["curvetype"]));
        NumericVector normal     = as<NumericVector>(shape_properties["normal"]);
        NumericVector normal_end = as<NumericVector>(shape_properties["normal_end"]);
        
        vec3f p[4];
        vec3f n[2];
        p[0] = vec3f(p1(0),p1(1),p1(2));
        p[1] = vec3f(p2(0),p2(1),p2(2));
        p[2] = vec3f(p3(0),p3(1),p3(2));
        p[3] = vec3f(p4(0),p4(1),p4(2));
        
        n[0] = vec3f(normal(0),normal(1),normal(2));
        n[1] = vec3f(normal_end(0),normal_end(1),normal_end(2));
        
        std::shared_ptr<CurveCommon> curve_data = std::make_shared<CurveCommon>(p, width, width_end, curvetype, n);
        entry = std::make_shared<curve>(u_min, u_max, curve_data, shape_material,
                                        ObjToWorld, WorldToObj, is_flipped);
        if(is_animated) {
          entry = std::make_shared<AnimatedHitable>(entry, Animate);
        }
        
        list.add(entry);
        break;
      }
      case CSG_OBJECT: {
        List csg_object = as<List>(SingleShape["csg_object"])(0);
        std::shared_ptr<ImplicitShape> shapes = parse_csg(csg_object);
        entry = std::make_shared<csg>(shape_material, shapes,
                                      ObjToWorld,WorldToObj, is_flipped);
        if(is_animated) {
          entry = std::make_shared<AnimatedHitable>(entry, Animate);
        }
        list.add(entry);
        break;
      }
      case PLY: {
        std::string plyfilename = Rcpp::as<std::string>(SingleShape["fileinfo"]);
        std::string plybasename =  Rcpp::as<std::string>(shape_properties["basename"]);
        Float scale_ply = as<Float>(shape_properties["scale_ply"]);
        entry = std::make_shared<plymesh>(plyfilename, plybasename, 
                                          shape_material, alpha[i], bump[i], 
                                          scale_ply,
                                          shutteropen, shutterclose, bvh_type, rng,
                                          ObjToWorld,WorldToObj, is_flipped);
        if(entry == nullptr) {
          continue;
        }
        if(isvolume) {
          entry = std::make_shared<constant_medium>(entry, fog_density, 
                                                    std::make_shared<constant_texture>(point3f(fog_color(0),
                                                                                               fog_color(1),
                                                                                               fog_color(2))));
        }
        if(is_animated) {
          entry = std::make_shared<AnimatedHitable>(entry, Animate);
        }
        list.add(entry);
        break;
      }
      case MESH3D: {
        List mesh_entry = as<List>(SingleShape["mesh_info"])(0);
        entry = std::make_shared<mesh3d>(mesh_entry, shape_material,
                                         shutteropen, shutterclose, bvh_type, rng,
                                         ObjToWorld,WorldToObj, is_flipped);
        if(is_animated) {
          entry = std::make_shared<AnimatedHitable>(entry, Animate);
        }
        list.add(entry);
        break;
      }
      case RAYMESH: {
        List raymesh_object = as<List>(SingleShape["mesh_info"])(0);
        bool importance_sample_lights = as<bool>(shape_properties["importance_sample_lights"]);
        bool calculate_consistent_normals = as<bool>(shape_properties["calculate_consistent_normals"]);
        bool override_material = as<bool>(shape_properties["override_material"]);
        bool flip_transmittance = as<bool>(shape_properties["flip_transmittance"]);
        //importance sample lights--need to change
        ////calculate consistent normals--need to change
        entry = std::make_shared<raymesh>(raymesh_object,
                                          shape_material,
                                          alpha[i], bump[i], 
                                          importance_sample_lights, 
                                          calculate_consistent_normals, 
                                          override_material,
                                          flip_transmittance,
                                          imp_sample_objects, 
                                          verbose, 
                                          shutteropen, shutterclose, bvh_type, rng, 
                                          ObjToWorld, WorldToObj, is_flipped);
        if(is_animated) {
          entry = std::make_shared<AnimatedHitable>(entry, Animate);
        }
        list.add(entry);
        break;
      }
    }
    if(importance_sample) {
      imp_sample_objects.add(entry);
    }
  }
  auto world_bvh = std::make_shared<bvh_node>(list, shutteropen, shutterclose, bvh_type, rng);
#ifdef FULL_DEBUG
  world_bvh->validate_bvh();
#endif
  // auto nodeleaf = world_bvh->CountNodeLeaf();
  // Rcpp::Rcout << "Node/Leaf: " << nodeleaf.first << " " << nodeleaf.second << " " << world_bvh->GetSize() << "\n";
  return(world_bvh);
}

