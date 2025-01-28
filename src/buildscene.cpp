#include "buildscene.h"
#include "hitable.h"
#include "sphere.h"
#include "hitablelist.h"
#include "perlin.h"
#include "texture.h"
#include "rectangle.h"
#include "box.h"
#include "constant.h"
#include "triangle.h"
#include "pdf.h"
#include "trimesh.h"
#include "disk.h"
#include "cylinder.h"
#include "ellipsoid.h"
#include "curve.h"
#include "csg.h"
#include "plymesh.h"
#include "mesh3d.h"
#include "raymesh.h"
#include "instance.h"
#include "transform.h"
#include "transformcache.h"
#include "texturecache.h"
#include "raylog.h"
#include "bvh.h"

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
  RAYMESH = 14,
  INSTANCE = 15
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

void LoadTexture(std::string image_file,
                 std::string alpha_file,
                 std::string bump_file,
                 std::string roughness_file,
                 std::vector<Float* >& textures,
                 std::vector<unsigned char * >& alpha_textures,
                 std::vector<unsigned char * >& bump_textures,
                 std::vector<unsigned char * >& roughness_textures,
                 int* nvec,
                 int* nveca,
                 int* nvecb,
                 int* nvecr,
                 NumericVector glossy_info,
                 bool has_image, bool has_alpha, bool has_bump, bool has_roughness,
                 TextureCache& texCache) {
  
  if(has_image) {
    int nx, ny, nn;
    Float* texture_data = texCache.LookupFloat(image_file, nx, ny, nn, 4);
    nn = 4;
    // texture_bytes += nx * ny * nn;
    textures.push_back(texture_data);
    nvec[0] = nx;
    nvec[1] = ny;
    nvec[2] = nn;
  } else {
    textures.push_back(nullptr);
  }
  if(has_alpha) {
    int nxa, nya, nna;
    unsigned char* alpha_data = texCache.LookupChar(alpha_file, nxa, nya, nna, 4);
    // texture_bytes += nxa * nya * nna;
    nna = 4;
    alpha_textures.push_back(alpha_data);
    nveca[0] = nxa;
    nveca[1] = nya;
    nveca[2] = nna;
  } else {
    alpha_textures.push_back(nullptr);
  }
  if(has_bump) {
    int nxb, nyb, nnb;
    unsigned char* bump_data = texCache.LookupChar(bump_file, nxb, nyb, nnb, 1);

    nnb = 1;
    bump_textures.push_back(bump_data);
    nvecb[0] = nxb;
    nvecb[1] = nyb;
    nvecb[2] = nnb;
  } else {
    bump_textures.push_back(nullptr);
  }
  if(has_roughness) {
    int nxr, nyr, nnr;
    unsigned char* roughness_data = texCache.LookupChar(roughness_file, nxr, nyr, nnr, 3);
    nnr = 3;
    // texture_bytes += nxr * nyr * nnr;
    
    Float min = glossy_info(9), max = glossy_info(10);
    Float rough_range = max-min;
    Float maxr = 0, minr = 1;
    for(int ii = 0; ii < nxr; ii++) {
      for(int jj = 0; jj < nyr; jj++) {
        Float temp_rough = roughness_data[nnr*ii + nnr*nxr*jj];
        maxr = maxr < temp_rough ? temp_rough : maxr;
        minr = minr > temp_rough ? temp_rough : minr;
        if(nnr > 1) {
          temp_rough = roughness_data[nnr*ii + nnr*nxr*jj+1];
          maxr = maxr < temp_rough ? temp_rough : maxr;
          minr = minr > temp_rough ? temp_rough : minr;
        }
      }
    }
    Float data_range = maxr-minr;
    for(int ii = 0; ii < nxr; ii++) {
      for(int jj = 0; jj < nyr; jj++) {
        if(!glossy_info(11)) {
          roughness_data[nnr*ii + nnr*nxr*jj] =
            (roughness_data[nnr*ii + nnr*nxr*jj]-minr)/data_range * rough_range + min;
          if(nnr > 1) {
            roughness_data[nnr*ii + nnr*nxr*jj+1] =
              (roughness_data[nnr*ii + nnr*nxr*jj+1]-minr)/data_range * rough_range + min;
          }
        } else {
          roughness_data[nnr*ii + nnr*nxr*jj] =
            (1.0-(roughness_data[nnr*ii + nnr*nxr*jj]-minr)/data_range) * rough_range + min;
          if(nnr > 1) {
            roughness_data[nnr*ii + nnr*nxr*jj+1] =
              (1.0-(roughness_data[nnr*ii + nnr*nxr*jj+1]-minr)/data_range) * rough_range + min;
          }
        }
      }
    }
    roughness_textures.push_back(roughness_data);
    nvecr[0] = nxr;
    nvecr[1] = nyr;
    nvecr[2] = nnr;
  } else {
    roughness_textures.push_back(nullptr);
  }
} 


std::shared_ptr<material> LoadSingleMaterial(List SingleMaterial,
                                             TextureCache& texCache,
                                             std::vector<Float* >& textures,
                                             std::vector<unsigned char * >& alpha_textures,
                                             std::vector<unsigned char * >& bump_textures,
                                             std::vector<unsigned char * >& roughness_textures,
                                             int* nvec,
                                             int* nveca,
                                             int* nvecb,
                                             int* nvecr,
                                             bool& has_image,
                                             bool& has_alpha,
                                             bool& has_bump,
                                             bool& has_roughness,
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
  

  NumericVector image_repeat = as<NumericVector>(as<List>(SingleMaterial["image_repeat"])(0));

  Float lightintensity = as<Float>(SingleMaterial["lightintensity"]);
  Float sigma = as<Float>(SingleMaterial["sigma"]);
  NumericVector glossyinfo = as<NumericVector>(as<List>(SingleMaterial["glossyinfo"])(0));
  
  std::string image_file = as<std::string>(SingleMaterial["image"]);
  std::string alpha_file = as<std::string>(SingleMaterial["alphaimage"]);
  std::string bump_file = as<std::string>(SingleMaterial["bump_texture"]);
  std::string roughness_file = as<std::string>(SingleMaterial["roughness_texture"]);
  has_image = !image_file.empty();
  has_alpha = !alpha_file.empty();
  has_bump = !bump_file.empty();
  has_roughness = !roughness_file.empty();
  
  LoadTexture(image_file, 
              alpha_file, 
              bump_file, 
              roughness_file,
              textures,
              alpha_textures,
              bump_textures,
              roughness_textures,
              nvec,
              nveca,
              nvecb,
              nvecr,
              glossyinfo,
              has_image,
              has_alpha,
              has_bump,
              has_roughness,
              texCache);
  
  std::shared_ptr<material> mat = nullptr;
  std::shared_ptr<texture> material_texture;
  
  bool is_tri_color = tricolorinfo.size() == 9;
  
  std::shared_ptr<roughness_texture> roughness;
  if(has_roughness) {
    roughness = std::make_shared<roughness_texture>(roughness_textures.back(), 
                                                    nvecr[0], nvecr[1], nvecr[2]);
  }
  
  if(has_image) {
    material_texture = std::make_shared<image_texture_float>(textures.back(), nvec[0], nvec[1], nvec[2], 
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
        dist = new TrowbridgeReitzDistribution(glossyinfo(1), glossyinfo(2),roughness, has_roughness,  true);
      } else {
        dist = new BeckmannDistribution(glossyinfo(1), glossyinfo(2),roughness, has_roughness, true);
      }
      point3f eta(glossyinfo(3), glossyinfo(4), glossyinfo(5));
      point3f kappa(glossyinfo(6), glossyinfo(7), glossyinfo(8));
      mat = std::make_shared<MicrofacetReflection>(material_texture, dist, eta, kappa);
      break;
    }
    case GLOSSY: {
      MicrofacetDistribution *dist;
      if(glossyinfo(0) == 1) {
        dist = new TrowbridgeReitzDistribution(glossyinfo(1), glossyinfo(2),roughness, has_roughness,  true);
      } else {
        dist = new BeckmannDistribution(glossyinfo(1), glossyinfo(2),roughness, has_roughness, true);
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
        dist = new TrowbridgeReitzDistribution(glossyinfo(1), glossyinfo(2),roughness, has_roughness,  true);
      } else {
        dist = new BeckmannDistribution(glossyinfo(1), glossyinfo(2),roughness, has_roughness, true);
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
                                     Float shutteropen, 
                                     Float shutterclose,
                                     std::vector<Float* >& textures, 
                                     std::vector<unsigned char * >& alpha_textures, 
                                     std::vector<unsigned char * >& bump_textures, 
                                     std::vector<unsigned char * >& roughness_textures,  
                                     std::vector<std::shared_ptr<material> >* shared_materials, 
                                     std::vector<std::shared_ptr<alpha_texture> >& alpha,
                                     std::vector<std::shared_ptr<bump_texture> >& bump,
                                     std::vector<std::shared_ptr<roughness_texture> >& roughness,
                                     int bvh_type,
                                     TransformCache& transformCache, 
                                     TextureCache& texCache,
                                     hitable_list& imp_sample_objects,
                                     std::vector<std::shared_ptr<hitable> >& instanced_objects,
                                     std::vector<std::shared_ptr<hitable_list> >& instance_importance_sampled,
                                     std::vector<int>& texture_idx,
                                     bool verbose, 
                                     random_gen& rng) {
  auto nvec  = std::make_unique<int[]>(3);
  auto nveca = std::make_unique<int[]>(3);
  auto nvecb = std::make_unique<int[]>(3);
  auto nvecr = std::make_unique<int[]>(3);
  int init_texture_size = texture_idx.size();
  
  hitable_list list;
  NumericMatrix IdentityMat(4,4);
  IdentityMat.fill_diag(1);
  Transform IdentityTransform(IdentityMat);
  Transform* Iden = transformCache.Lookup(IdentityTransform);
  
  List RayMaterials = as<List>(scene["material"]);
  List ShapeInfo = as<List>(scene["shape_info"]);
  List Transforms = as<List>(scene["transforms"]);
  List AnimationInfo = as<List>(scene["animation_info"]);
  size_t n = ShapeInfo.length();

  NumericVector x = scene["x"];
  NumericVector y = scene["y"];
  NumericVector z = scene["z"];
  
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
    bool is_shared_mat = !Rcpp::IntegerVector::is_na(material_id_vec(0));
    int material_id = material_id_vec(0);
    std::shared_ptr<material> shape_material;
    bool has_image = false;
    bool has_alpha = false;
    bool has_bump = false;
    bool has_roughness = false;
    
    if(is_shared_mat && shared_materials->size() > static_cast<size_t>(material_id - 1)) {
      shape_material = shared_materials->at(material_id-1);
      texture_idx.push_back(material_id-1);
    } else {
      shape_material = LoadSingleMaterial(SingleMaterial,
                                          texCache,
                                          textures,
                                          alpha_textures,
                                          bump_textures,
                                          roughness_textures,
                                          nvec.get(),
                                          nveca.get(),
                                          nvecb.get(),
                                          nvecr.get(),
                                          has_image,
                                          has_alpha,
                                          has_bump,
                                          has_roughness,
                                          tricolorinfo);
      texture_idx.push_back(init_texture_size + i);
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
    
    
    //`mat_idx` selects the index of the texture, in case there's a shared material.
    int mat_idx = texture_idx.back();
    if(has_alpha) {
      alpha.push_back(std::make_shared<alpha_texture>(alpha_textures[mat_idx], nveca[0], nveca[1], nveca[2]));
    } else {
      alpha.push_back(nullptr);
    }
    if(has_bump) {
      bump.push_back(std::make_shared<bump_texture>(bump_textures[mat_idx], nvecb[0], nvecb[1], nvecb[2], 
                                                    bump_intensity, image_repeat[0], image_repeat[1]));
    } else {
      bump.push_back(nullptr);
    }
    //Generate center vector
    ShapeEnum shape_type = static_cast<ShapeEnum>(shape(i));
    vec3f center = vec3f(x(i), y(i), z(i));
    
    Transform GroupTransform(GroupTransformMat);
    Transform TempM = GroupTransform * Translate(center) *
      rotation_order_matrix(angle, order_rotation) * 
      Scale(scales[0], scales[1], scales[2]);
    Transform* ObjToWorld = transformCache.Lookup(TempM);
    Transform* WorldToObj = transformCache.Lookup(TempM.GetInverseMatrix());
    
    
    Transform AnimationStartTransform(StartTransformAnimationMat);
    Transform AnimationEndTransform(EndTransformAnimationMat);
    Transform* StartAnim = transformCache.Lookup(AnimationStartTransform);
    Transform* EndAnim = transformCache.Lookup(AnimationEndTransform);
    
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
        entry = std::make_shared<sphere>(radius, shape_material, alpha[mat_idx], bump[mat_idx],
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
                                          0, shape_material, alpha[mat_idx], bump[mat_idx], 
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
                                          0, shape_material, alpha[mat_idx], bump[mat_idx], 
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
                                          0, shape_material, alpha[mat_idx], bump[mat_idx], 
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
                                       shape_material, alpha[mat_idx], bump[mat_idx],
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
        int subdivision_levels = as<int>(shape_properties["subdivision_levels"]);
        std::string displacement_texture = as<std::string>(shape_properties["displacement_texture"]);
        Float displacement_intensity = as<Float>(shape_properties["displacement_intensity"]);
        bool is_vector_displacement = as<bool>(shape_properties["displacement_vector"]);
        bool recalculate_normals = as<bool>(shape_properties["recalculate_normals"]);
        
        
        Float sigma_obj = as<Float>(SingleMaterial["sigma"]);

        
        entry = std::make_shared<trimesh>(objfilename, objbasename, 
                                          scale_obj, sigma_obj, shape_material, 
                                          alpha[mat_idx], bump[mat_idx], 
                                          load_material, load_textures, load_vertex_colors,
                                          importance_sample_lights, load_normals, calculate_consistent_normals,
                                          subdivision_levels, 
                                          displacement_texture,
                                          displacement_intensity,
                                          is_vector_displacement,
                                          texCache,
                                          recalculate_normals,
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
        
        entry = std::make_shared<disk>(vec3f(0,0,0), radius, inner_radius, shape_material, alpha[mat_idx], bump[mat_idx],
                                       ObjToWorld, WorldToObj, is_flipped);
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
                                           shape_material, alpha[mat_idx], bump[mat_idx],
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
        entry = std::make_shared<ellipsoid>(point3f(0), 1, 
                                            vec3f(a,b,c),
                                            shape_material, alpha[mat_idx], bump[mat_idx],
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
        
        point3f p[4];
        vec3f n[2];
        p[0] = point3f(p1(0),p1(1),p1(2));
        p[1] = point3f(p2(0),p2(1),p2(2));
        p[2] = point3f(p3(0),p3(1),p3(2));
        p[3] = point3f(p4(0),p4(1),p4(2));
        
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
        int subdivision_levels = as<int>(shape_properties["subdivision_levels"]);
        bool recalculate_normals = as<bool>(shape_properties["recalculate_normals"]);
        
        entry = std::make_shared<plymesh>(plyfilename, plybasename, 
                                          shape_material, alpha[mat_idx], bump[mat_idx], 
                                          scale_ply, subdivision_levels, recalculate_normals,
                                          verbose,
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
        std::string displacement_texture = as<std::string>(mesh_entry["displacement_texture"]);
        Float displacement_intensity = as<Float>(mesh_entry["displacement_intensity"]);
        bool is_vector_displacement = as<bool>(mesh_entry["displacement_vector"]);
        bool recalculate_normals = as<bool>(mesh_entry["recalculate_normals"]);
        
        entry = std::make_shared<mesh3d>(mesh_entry, shape_material,
                                         displacement_texture,
                                         displacement_intensity,
                                         is_vector_displacement,
                                         texCache, recalculate_normals,
                                         verbose,
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
        int subdivision_levels = as<int>(shape_properties["subdivision_levels"]);
        std::string displacement_texture = as<std::string>(shape_properties["displacement_texture"]);
        Float displacement_intensity = as<Float>(shape_properties["displacement_intensity"]);
        bool is_vector_displacement = as<bool>(shape_properties["displacement_vector"]);
        bool recalculate_normals = as<bool>(shape_properties["recalculate_normals"]);
        
        //importance sample lights--need to change
        ////calculate consistent normals--need to change
        entry = std::make_shared<raymesh>(raymesh_object,
                                          shape_material,
                                          alpha[mat_idx], bump[mat_idx], 
                                          importance_sample_lights, 
                                          calculate_consistent_normals, 
                                          override_material,
                                          flip_transmittance,
                                          subdivision_levels,
                                          displacement_texture,
                                          displacement_intensity,
                                          is_vector_displacement,
                                          texCache, recalculate_normals,
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
      case INSTANCE: {
        List original_scene = as<List>(shape_properties["original_scene"])(0);
        NumericVector x_values = as<NumericVector>(shape_properties["x_values"]);
        NumericVector y_values = as<NumericVector>(shape_properties["y_values"]);
        NumericVector z_values = as<NumericVector>(shape_properties["z_values"]);
        NumericVector angle_x = as<NumericVector>(shape_properties["angle_x"]);
        NumericVector angle_y = as<NumericVector>(shape_properties["angle_y"]);
        NumericVector angle_z = as<NumericVector>(shape_properties["angle_z"]);
        NumericVector scale_x = as<NumericVector>(shape_properties["scale_x"]);
        NumericVector scale_y = as<NumericVector>(shape_properties["scale_y"]);
        NumericVector scale_z = as<NumericVector>(shape_properties["scale_z"]);
        IntegerVector shape_vec = as<IntegerVector>(original_scene["shape"]);
        auto instance_importance_sample_list = std::make_shared<hitable_list>();
        std::shared_ptr<hitable> instance_scene = build_scene(original_scene,
                                                              shape_vec,
                                                              shutteropen,
                                                              shutterclose,
                                                              textures, 
                                                              alpha_textures, 
                                                              bump_textures, 
                                                              roughness_textures,  
                                                              shared_materials,
                                                              alpha, bump, roughness,
                                                              bvh_type,
                                                              transformCache,
                                                              texCache,
                                                              (*instance_importance_sample_list),
                                                              instanced_objects,
                                                              instance_importance_sampled,
                                                              texture_idx,
                                                              false,
                                                              rng);
        instanced_objects.push_back(instance_scene);
        bool any_importance_sampled = instance_importance_sample_list->size() > 0;
        if(any_importance_sampled) {
          instance_importance_sampled.push_back(instance_importance_sample_list);
        }
        
        for(size_t ii = 0; ii < (size_t)x_values.size(); ii++) {
          vec3f center_instance = vec3f(x_values(ii), y_values(ii), z_values(ii));
          NumericVector angle_instance = {angle_x(ii), angle_y(ii), angle_z(ii)};

          Transform InstanceTransform = GroupTransform * 
            Translate(center_instance) * 
            rotation_order_matrix(angle_instance, order_rotation) * 
            Scale(scale_x(ii), scale_y(ii), scale_z(ii));
          Transform* ObjToWorldInst = transformCache.Lookup(InstanceTransform);
          Transform* WorldToObjInst = transformCache.Lookup(InstanceTransform.GetInverseMatrix());
          list.add(std::make_shared<instance>(instance_scene.get(),
                                              ObjToWorldInst, 
                                              WorldToObjInst,
                                              instance_importance_sample_list.get()));
          if(any_importance_sampled) {
            imp_sample_objects.add(list.back());
          }
        }
        break;
      }
    }
    if(importance_sample && shape_type != INSTANCE) {
      imp_sample_objects.add(entry);
    }
  }
  std::shared_ptr<BVHAggregate> world_bvh = std::make_shared<BVHAggregate>(list.objects, shutteropen, shutterclose, 1, true,
                                                   Iden,Iden,false );
  
  #ifdef FULL_DEBUG
  world_bvh->validate_bvh();
#endif
  // auto nodeleaf = world_bvh->CountNodeLeaf();
  // Rcpp::Rcout << "Node/Leaf: " << nodeleaf.first << " " << nodeleaf.second << " " << world_bvh->GetSize() << "\n";
  return(world_bvh);
}

