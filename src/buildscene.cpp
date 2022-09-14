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


std::shared_ptr<hitable> build_scene(IntegerVector& type, 
                     NumericVector& radius, IntegerVector& shape,
                     List& position_list,
                     List& properties, 
                     int n, Float shutteropen, Float shutterclose,
                     LogicalVector& ischeckered, List& checkercolors, 
                     List gradient_info,
                     NumericVector& noise, LogicalVector& isnoise,
                     NumericVector& noisephase, NumericVector& noiseintensity, List noisecolorlist,
                     List& angle, 
                     LogicalVector& isimage, LogicalVector has_alpha,
                     std::vector<unsigned char * >& alpha_textures, std::vector<int* >& nveca,
                     std::vector<Float* >& textures, std::vector<int* >& nvec,
                     LogicalVector has_bump,
                     std::vector<unsigned char * >& bump_textures, std::vector<int* >& nvecb,
                     NumericVector& bump_intensity,
                     std::vector<unsigned char * >& roughness_textures,  std::vector<int* >& nvecr,
                     LogicalVector has_roughness,
                     NumericVector& lightintensity,
                     LogicalVector& isflipped,
                     LogicalVector& isvolume, NumericVector& voldensity,
                     List& order_rotation_list, 
                     LogicalVector& isgrouped, List& group_transform,
                     LogicalVector& tri_normal_bools, LogicalVector& is_tri_color, List& tri_color_vert, 
                     CharacterVector& fileinfo, CharacterVector& filebasedir,
                     List& scale_list, NumericVector& sigma,  List &glossyinfo,
                     IntegerVector& shared_id_mat, LogicalVector& is_shared_mat,
                     std::vector<std::shared_ptr<material> >* shared_materials, List& image_repeat_list,
                     List& csg_info, List& mesh_list, int bvh_type,
                     TransformCache& transformCache, List& animation_info, 
                     LogicalVector& implicit_sample, hitable_list& imp_sample_objects,
                     bool verbose,
                     random_gen& rng) {
  hitable_list list;
  NumericMatrix IdentityMat(4,4);
  IdentityMat.fill_diag(1);
  LogicalVector isgradient = gradient_info["isgradient"];
  List gradient_colors = gradient_info["gradient_colors"];
  LogicalVector gradient_trans = gradient_info["gradient_trans"];
  List gradient_control_points = gradient_info["gradient_control_points"];
  LogicalVector is_world_gradient = gradient_info["is_world_gradient"];
  LogicalVector gradient_is_hsv = gradient_info["type"];
  
  NumericVector x = position_list["xvec"];
  NumericVector y = position_list["yvec"];
  NumericVector z = position_list["zvec"];
  
  List csg_list = csg_info["csg"];
  
  //Animation
  LogicalVector has_animation        = animation_info["animation_bool"];
  List start_transform_animation     = animation_info["start_transform_animation"];
  List end_transform_animation       = animation_info["end_transform_animation"];
  NumericVector animation_start_time = animation_info["animation_start_time"];
  NumericVector animation_end_time   = animation_info["animation_end_time"];     
  
  NumericVector tempvector;
  NumericVector tempchecker;
  NumericVector tempgradient;
  NumericVector tempnoisecolor;
  NumericVector temprotvec;
  NumericVector order_rotation;
  NumericMatrix temp_group_transform;
  NumericVector temp_tri_color;
  NumericVector temp_scales;
  NumericVector temp_glossy;
  NumericVector temp_repeat;
  NumericVector temp_gradient_control;
  NumericMatrix temp_animation_transform_start;
  NumericMatrix temp_animation_transform_end;
  
  
  
  int prop_len;
  
  List templist;
  vec3f center(x(0), y(0), z(0));
  std::vector<std::shared_ptr<alpha_texture> > alpha(n);
  std::vector<std::shared_ptr<bump_texture> > bump(n);
  std::vector<std::shared_ptr<roughness_texture> > roughness(n);
  
  for(int i = 0; i < n; i++) {
    tempvector = as<NumericVector>(properties(i));
    tempgradient = as<NumericVector>(gradient_colors(i));
    tempchecker = as<NumericVector>(checkercolors(i));
    tempnoisecolor = as<NumericVector>(noisecolorlist(i));
    temprotvec = as<NumericVector>(angle(i));
    temp_tri_color = as<NumericVector>(tri_color_vert(i));
    order_rotation = as<NumericVector>(order_rotation_list(i));
    temp_scales = as<NumericVector>(scale_list(i));
    temp_glossy = as<NumericVector>(glossyinfo(i));
    temp_repeat = as<NumericVector>(image_repeat_list(i));
    temp_gradient_control = as<NumericVector>(gradient_control_points(i));

    if(isgrouped(i)) {
      temp_group_transform = as<NumericMatrix>(group_transform(i));
    } else {
      temp_group_transform = IdentityMat;
    }
    if(has_animation(i)) {
      temp_animation_transform_start = as<NumericMatrix>(start_transform_animation(i));
      temp_animation_transform_end = as<NumericMatrix>(end_transform_animation(i));
    } else {
      temp_animation_transform_start = IdentityMat;
      temp_animation_transform_end = IdentityMat;
    }
    prop_len=2;
    //Generate texture
    std::shared_ptr<material> tex = nullptr;
    if(has_alpha(i)) {
      alpha[i] = std::make_shared<alpha_texture>(alpha_textures[i], 
                                                 nveca[i][0], nveca[i][1], nveca[i][2]);
    }
    if(has_bump(i)) {
      bump[i] = std::make_shared<bump_texture>(bump_textures[i], 
                                               nvecb[i][0], nvecb[i][1], nvecb[i][2], 
                                               bump_intensity(i),temp_repeat[0], temp_repeat[1]);
    }
    if(has_roughness(i)) {
      roughness[i] = std::make_shared<roughness_texture>(roughness_textures[i], 
                                                         nvecr[i][0], nvecr[i][1], nvecr[i][2]);
    }
    if(type(i) == 2) {
      prop_len = 3;
    } else if (type(i) == 3) {
      prop_len = 7;
    } else if (type(i) == 5) {
      prop_len = 3;
    } else if (type(i) == 8) {
      prop_len = 8;
    } else if (type(i) == 9) {
      prop_len = 6;
    }

    if(is_shared_mat(i) && shared_materials->size() > static_cast<int>(shared_id_mat(i)-1)) {
      tex = shared_materials->at(shared_id_mat(i)-1);
    } else {
      if(type(i) == 1) {
        if(isimage(i)) {
          tex = std::make_shared<lambertian>(std::make_shared<image_texture_float>(textures[i],nvec[i][0],nvec[i][1],nvec[i][2], 
                                                 temp_repeat[0], temp_repeat[1], 1.0));
        } else if (isnoise(i)) {
          tex = std::make_shared<lambertian>(std::make_shared<noise_texture>(noise(i),point3f(tempvector(0),tempvector(1),tempvector(2)),
                                                 point3f(tempnoisecolor(0),tempnoisecolor(1),tempnoisecolor(2)),
                                                 noisephase(i), noiseintensity(i)));
        } else if (ischeckered(i)) {
          tex = std::make_shared<lambertian>(std::make_shared<checker_texture>(
            std::make_shared<constant_texture>(point3f(tempchecker(0),tempchecker(1),tempchecker(2))),
            std::make_shared<constant_texture>(point3f(tempvector(0),tempvector(1),tempvector(2))),
                                             tempchecker(3)));
        } else if (isgradient(i) && !is_world_gradient(i)) {
          tex = std::make_shared<lambertian>(std::make_shared<gradient_texture>(point3f(tempvector(0),tempvector(1),tempvector(2)),
                                                                                point3f(tempgradient(0),tempgradient(1),tempgradient(2)),
                                                    gradient_trans(i), gradient_is_hsv(i)));
        } else if (is_tri_color(i)) {
          tex = std::make_shared<lambertian>(std::make_shared<triangle_texture>(point3f(temp_tri_color(0),temp_tri_color(1),temp_tri_color(2)),
                                                                                point3f(temp_tri_color(3),temp_tri_color(4),temp_tri_color(5)),
                                                                                point3f(temp_tri_color(6),temp_tri_color(7),temp_tri_color(8))));
        } else if (is_world_gradient(i)) {
          tex = std::make_shared<lambertian>(std::make_shared<world_gradient_texture>(point3f(temp_gradient_control(0),temp_gradient_control(1),temp_gradient_control(2)),
                                                                                      point3f(temp_gradient_control(3),temp_gradient_control(4),temp_gradient_control(5)),
                                                                                      point3f(tempvector(0),tempvector(1),tempvector(2)),
                                                                                      point3f(tempgradient(0),tempgradient(1),tempgradient(2)),
                                                          gradient_is_hsv(i)));
        } else {
          tex = std::make_shared<lambertian>(std::make_shared<constant_texture>(point3f(tempvector(0),tempvector(1),tempvector(2))));
        }
      } else if (type(i) == 2) {
        if(isimage(i)) {
          tex = std::make_shared<metal>(std::make_shared<image_texture_float>(textures[i],nvec[i][0],nvec[i][1],nvec[i][2], 
                                            temp_repeat[0], temp_repeat[1], 1.0),
                          tempvector(3), 
                          point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                          point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        } else if (isnoise(i)) {
          tex = std::make_shared<metal>(std::make_shared<noise_texture>(noise(i),point3f(tempvector(0),tempvector(1),tempvector(2)),
                                                                        point3f(tempnoisecolor(0),tempnoisecolor(1),tempnoisecolor(2)),
                                                 noisephase(i), noiseintensity(i)),
                          tempvector(3), 
                          point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                          point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        } else if (ischeckered(i)) {
          tex = std::make_shared<metal>(std::make_shared<checker_texture>(std::make_shared<constant_texture>(point3f(tempchecker(0),tempchecker(1),tempchecker(2))),
                                              std::make_shared<constant_texture>(point3f(tempvector(0),tempvector(1),tempvector(2))),
                                              tempchecker(3)),
                          tempvector(3), 
                          point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                          point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        } else if (isgradient(i) && !is_world_gradient(i)) {
          tex = std::make_shared<metal>(std::make_shared<gradient_texture>(point3f(tempvector(0),tempvector(1),tempvector(2)),
                                                                           point3f(tempgradient(0),tempgradient(1),tempgradient(2)),
                                                    gradient_trans(i), gradient_is_hsv(i)),
                          tempvector(3), 
                          point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                          point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        } else if (is_tri_color(i)) {
          tex = std::make_shared<metal>(std::make_shared<triangle_texture>(point3f(temp_tri_color(0),temp_tri_color(1),temp_tri_color(2)),
                                                                           point3f(temp_tri_color(3),temp_tri_color(4),temp_tri_color(5)),
                                                                           point3f(temp_tri_color(6),temp_tri_color(7),temp_tri_color(8))),
                          tempvector(3), 
                          point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                          point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        } else if (is_world_gradient(i)) {
          tex = std::make_shared<metal>(std::make_shared<world_gradient_texture>(point3f(temp_gradient_control(0),temp_gradient_control(1),temp_gradient_control(2)),
                                                                                 point3f(temp_gradient_control(3),temp_gradient_control(4),temp_gradient_control(5)),
                                                     point3f(tempvector(0),tempvector(1),tempvector(2)),
                                                     point3f(tempgradient(0),tempgradient(1),tempgradient(2)),
                                                     gradient_is_hsv(i)),
                         tempvector(3), 
                         point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                         point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        } else {
          tex = std::make_shared<metal>(std::make_shared<constant_texture>(point3f(tempvector(0),tempvector(1),tempvector(2))),
                          tempvector(3), 
                          point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                          point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        }
      } else if (type(i) == 3) {
        tex = std::make_shared<dielectric>(point3f(tempvector(0),tempvector(1),tempvector(2)), tempvector(3), 
                                           point3f(tempvector(4),tempvector(5),tempvector(6)), 
                             tempvector(7));
      } else if (type(i) == 4) {
        if(isimage(i)) {
          tex = std::make_shared<orennayar>(std::make_shared<image_texture_float>(textures[i],nvec[i][0],nvec[i][1],nvec[i][2], 
                                                temp_repeat[0], temp_repeat[1], 1.0), sigma(i));
        } else if (isnoise(i)) {
          tex = std::make_shared<orennayar>(std::make_shared<noise_texture>(noise(i),point3f(tempvector(0),tempvector(1),tempvector(2)),
                                                 point3f(tempnoisecolor(0),tempnoisecolor(1),tempnoisecolor(2)),
                                                 noisephase(i), noiseintensity(i)), sigma(i));
        } else if (ischeckered(i)) {
          tex = std::make_shared<orennayar>(std::make_shared<checker_texture>(std::make_shared<constant_texture>(point3f(tempchecker(0),tempchecker(1),tempchecker(2))),
                                                   std::make_shared<constant_texture>(point3f(tempvector(0),tempvector(1),tempvector(2))),tempchecker(3)), 
                                                   sigma(i));
        } else if (isgradient(i) && !is_world_gradient(i)) {
          tex = std::make_shared<orennayar>(std::make_shared<gradient_texture>(point3f(tempvector(0),tempvector(1),tempvector(2)),
                                                    point3f(tempgradient(0),tempgradient(1),tempgradient(2)),
                                                    gradient_trans(i), gradient_is_hsv(i)), sigma(i));
        } else if (is_tri_color(i)) {
          tex = std::make_shared<orennayar>(std::make_shared<triangle_texture>(point3f(temp_tri_color(0),temp_tri_color(1),temp_tri_color(2)),
                                                    point3f(temp_tri_color(3),temp_tri_color(4),temp_tri_color(5)),
                                                    point3f(temp_tri_color(6),temp_tri_color(7),temp_tri_color(8))), sigma(i) );
        } else if (is_world_gradient(i)) {
          tex = std::make_shared<orennayar>(std::make_shared<world_gradient_texture>(point3f(temp_gradient_control(0),temp_gradient_control(1),temp_gradient_control(2)),
                                                         point3f(temp_gradient_control(3),temp_gradient_control(4),temp_gradient_control(5)),
                                                         point3f(tempvector(0),tempvector(1),tempvector(2)),
                                                         point3f(tempgradient(0),tempgradient(1),tempgradient(2)),
                                                         gradient_is_hsv(i)), sigma(i));
        } else {
          tex = std::make_shared<orennayar>(std::make_shared<constant_texture>(point3f(tempvector(0),tempvector(1),tempvector(2))), sigma(i)); //marked as small definite loss in valgrind memcheck
        }
      } else if (type(i) == 5) {
        std::shared_ptr<texture> light_tex = nullptr;
        if(isimage(i)) {
          light_tex = std::make_shared<image_texture_float>(textures[i],nvec[i][0],nvec[i][1],nvec[i][2], 
                                        temp_repeat[0], temp_repeat[1], 1.0);
        } else if (isnoise(i)) {
          light_tex = std::make_shared<noise_texture>(noise(i),point3f(tempvector(0),tempvector(1),tempvector(2)),
                                                point3f(tempnoisecolor(0),tempnoisecolor(1),tempnoisecolor(2)),
                                                noisephase(i), noiseintensity(i));
        } else if (ischeckered(i)) {
          light_tex = std::make_shared<checker_texture>(std::make_shared<constant_texture>(point3f(tempchecker(0),tempchecker(1),tempchecker(2))),
                                                  std::make_shared<constant_texture>(point3f(tempvector(0),tempvector(1),tempvector(2))),tempchecker(3));
        } else if (isgradient(i) && !is_world_gradient(i)) {
          light_tex = std::make_shared<gradient_texture>(point3f(tempvector(0),tempvector(1),tempvector(2)),
                                                   point3f(tempgradient(0),tempgradient(1),tempgradient(2)),
                                                   gradient_trans(i), gradient_is_hsv(i));
        } else if (is_tri_color(i)) {
          light_tex = std::make_shared<triangle_texture>(point3f(temp_tri_color(0),temp_tri_color(1),temp_tri_color(2)),
                                                   point3f(temp_tri_color(3),temp_tri_color(4),temp_tri_color(5)),
                                                   point3f(temp_tri_color(6),temp_tri_color(7),temp_tri_color(8)));
        } else if (is_world_gradient(i)) {
          light_tex = std::make_shared<world_gradient_texture>(point3f(temp_gradient_control(0),temp_gradient_control(1),temp_gradient_control(2)),
                                                         point3f(temp_gradient_control(3),temp_gradient_control(4),temp_gradient_control(5)),
                                                         point3f(tempvector(0),tempvector(1),tempvector(2)),
                                                         point3f(tempgradient(0),tempgradient(1),tempgradient(2)),
                                                         gradient_is_hsv(i));
        } else {
          light_tex = std::make_shared<constant_texture>(point3f(tempvector(0),tempvector(1),tempvector(2)));
        }
        bool is_invisible = tempvector(3) == 1;
        tex = std::make_shared<diffuse_light>(light_tex, lightintensity(i), is_invisible);
      } else if (type(i) == 6) {
        MicrofacetDistribution *dist;
        if(temp_glossy(0) == 1) {
          dist = new TrowbridgeReitzDistribution(temp_glossy(1), temp_glossy(2),roughness[i], has_roughness(i),  true);
        } else {
          dist = new BeckmannDistribution(temp_glossy(1), temp_glossy(2),roughness[i], has_roughness(i), true);
        }
        if(isimage(i)) {
          tex = std::make_shared<MicrofacetReflection>(std::make_shared<image_texture_float>(textures[i],nvec[i][0],nvec[i][1],nvec[i][2], 
                                                           temp_repeat[0], temp_repeat[1], 1.0), dist, 
                                         point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                                         point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        } else if (isnoise(i)) {
          tex = std::make_shared<MicrofacetReflection>(std::make_shared<noise_texture>(noise(i),point3f(tempvector(0),tempvector(1),tempvector(2)),
                                                point3f(tempnoisecolor(0),tempnoisecolor(1),tempnoisecolor(2)),
                                                noisephase(i), noiseintensity(i)), dist, 
                                                point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                                                point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        } else if (ischeckered(i)) {
          tex = std::make_shared<MicrofacetReflection>(std::make_shared<checker_texture>(std::make_shared<constant_texture>(point3f(tempchecker(0),tempchecker(1),tempchecker(2))),
                                                  std::make_shared<constant_texture>(point3f(tempvector(0),tempvector(1),tempvector(2))),tempchecker(3)), 
                                                  dist, 
                                                  point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                                                  point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        }  else if (isgradient(i) && !is_world_gradient(i)) {
          tex = std::make_shared<MicrofacetReflection>(std::make_shared<gradient_texture>(point3f(tempvector(0),tempvector(1),tempvector(2)),
                                                   point3f(tempgradient(0),tempgradient(1),tempgradient(2)),
                                                   gradient_trans(i), gradient_is_hsv(i)), dist, 
                                                   point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                                                   point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        } else if (is_tri_color(i)) {
          tex = std::make_shared<MicrofacetReflection>(std::make_shared<triangle_texture>(point3f(temp_tri_color(0),temp_tri_color(1),temp_tri_color(2)),
                                                   point3f(temp_tri_color(3),temp_tri_color(4),temp_tri_color(5)),
                                                   point3f(temp_tri_color(6),temp_tri_color(7),temp_tri_color(8))), dist, 
                                                   point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                                                   point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        } else if (is_world_gradient(i)) {
          tex = std::make_shared<MicrofacetReflection>(std::make_shared<world_gradient_texture>(point3f(temp_gradient_control(0),temp_gradient_control(1),temp_gradient_control(2)),
                                                                    point3f(temp_gradient_control(3),temp_gradient_control(4),temp_gradient_control(5)),
                                                                    point3f(tempvector(0),tempvector(1),tempvector(2)),
                                                                    point3f(tempgradient(0),tempgradient(1),tempgradient(2)), gradient_is_hsv(i)), dist, 
                                                                    point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                                                                    point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        } else {
          tex = std::make_shared<MicrofacetReflection>(std::make_shared<constant_texture>(point3f(tempvector(0),tempvector(1),tempvector(2))), dist, 
                                         point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                                         point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        }
      } if (type(i) == 7) {
        MicrofacetDistribution *dist;
        if(temp_glossy(0) == 1) {
          dist = new TrowbridgeReitzDistribution(temp_glossy(1), temp_glossy(2),roughness[i], has_roughness(i), true);
        } else {
          dist = new BeckmannDistribution(temp_glossy(1), temp_glossy(2),roughness[i], has_roughness(i), true);
        }
        if(isimage(i)) {
          tex = std::make_shared<glossy>(std::make_shared<image_texture_float>(textures[i],nvec[i][0],nvec[i][1],nvec[i][2], 
                                             temp_repeat[0], temp_repeat[1], 1.0), dist, 
                           point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                           point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        } else if (isnoise(i)) {
          tex = std::make_shared<glossy>(std::make_shared<noise_texture>(noise(i),point3f(tempvector(0),tempvector(1),tempvector(2)),
                                             point3f(tempnoisecolor(0),tempnoisecolor(1),tempnoisecolor(2)),
                                             noisephase(i), noiseintensity(i)), dist, 
                           point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                           point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        } else if (ischeckered(i)) {
          tex = std::make_shared<glossy>(std::make_shared<checker_texture>(std::make_shared<constant_texture>(point3f(tempchecker(0),tempchecker(1),tempchecker(2))),
                                               std::make_shared<constant_texture>(point3f(tempvector(0),tempvector(1),tempvector(2))),tempchecker(3)), 
                           dist, 
                           point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                           point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        }  else if (isgradient(i) && !is_world_gradient(i)) {
          tex = std::make_shared<glossy>(std::make_shared<gradient_texture>(point3f(tempvector(0),tempvector(1),tempvector(2)),
                                                     point3f(tempgradient(0),tempgradient(1),tempgradient(2)),
                                                     gradient_trans(i), gradient_is_hsv(i)), dist, 
                           point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                           point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        } else if (is_tri_color(i)) {
          tex = std::make_shared<glossy>(std::make_shared<triangle_texture>(point3f(temp_tri_color(0),temp_tri_color(1),temp_tri_color(2)),
                                                point3f(temp_tri_color(3),temp_tri_color(4),temp_tri_color(5)),
                                                point3f(temp_tri_color(6),temp_tri_color(7),temp_tri_color(8))), dist, 
                           point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                           point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        } else if (is_world_gradient(i)) {
          tex = std::make_shared<glossy>(std::make_shared<world_gradient_texture>(point3f(temp_gradient_control(0),temp_gradient_control(1),temp_gradient_control(2)),
                                                      point3f(temp_gradient_control(3),temp_gradient_control(4),temp_gradient_control(5)),
                                                      point3f(tempvector(0),tempvector(1),tempvector(2)),
                                                      point3f(tempgradient(0),tempgradient(1),tempgradient(2)), gradient_is_hsv(i)), dist, 
                          point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                          point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        } else {
          tex = std::make_shared<glossy>(std::make_shared<constant_texture>(point3f(tempvector(0),tempvector(1),tempvector(2))), dist, 
                           point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                           point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        }
      } else if (type(i) == 8) {
        bool is_invisible = tempvector(8) == 1;
        tex = std::make_shared<spot_light>(std::make_shared<constant_texture>(point3f(tempvector(0),tempvector(1),tempvector(2))*lightintensity(i)),
                             vec3f(tempvector(3),tempvector(4),tempvector(5)), tempvector(6),tempvector(7), is_invisible);
      } else if (type(i) == 9) {
        tex = std::make_shared<hair>(point3f(tempvector(0),tempvector(1),tempvector(2)), 
                       tempvector(3),tempvector(4),tempvector(5),tempvector(6));
      } else if (type(i) == 10) {
        MicrofacetDistribution *dist;
        if(temp_glossy(0) == 1) {
          dist = new TrowbridgeReitzDistribution(temp_glossy(1), temp_glossy(2),roughness[i], has_roughness(i), true);
        } else {
          dist = new BeckmannDistribution(temp_glossy(1), temp_glossy(2),roughness[i], has_roughness(i), true);
        }
        if(isimage(i)) {
          tex = std::make_shared<MicrofacetTransmission>(std::make_shared<image_texture_float>(textures[i],nvec[i][0],nvec[i][1],nvec[i][2], 
                                                                                       temp_repeat[0], temp_repeat[1], 1.0), dist, 
                                                                                       point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                                                                                       point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        } else if (isnoise(i)) {
          tex = std::make_shared<MicrofacetTransmission>(std::make_shared<noise_texture>(noise(i),point3f(tempvector(0),tempvector(1),tempvector(2)),
                                                                                       point3f(tempnoisecolor(0),tempnoisecolor(1),tempnoisecolor(2)),
                                                                                       noisephase(i), noiseintensity(i)), dist, 
                                                                                       point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                                                                                       point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        } else if (ischeckered(i)) {
          tex = std::make_shared<MicrofacetTransmission>(std::make_shared<checker_texture>(std::make_shared<constant_texture>(point3f(tempchecker(0),tempchecker(1),tempchecker(2))),
                                                                                         std::make_shared<constant_texture>(point3f(tempvector(0),tempvector(1),tempvector(2))),tempchecker(3)), 
                                                                                         dist, 
                                                                                         point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                                                                                         point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        }  else if (isgradient(i) && !is_world_gradient(i)) {
          tex = std::make_shared<MicrofacetTransmission>(std::make_shared<gradient_texture>(point3f(tempvector(0),tempvector(1),tempvector(2)),
                                                                                          point3f(tempgradient(0),tempgradient(1),tempgradient(2)),
                                                                                          gradient_trans(i), gradient_is_hsv(i)), dist, 
                                                                                          point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                                                                                          point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        } else if (is_tri_color(i)) {
          tex = std::make_shared<MicrofacetTransmission>(std::make_shared<triangle_texture>(point3f(temp_tri_color(0),temp_tri_color(1),temp_tri_color(2)),
                                                                                          point3f(temp_tri_color(3),temp_tri_color(4),temp_tri_color(5)),
                                                                                          point3f(temp_tri_color(6),temp_tri_color(7),temp_tri_color(8))), dist, 
                                                                                          point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                                                                                          point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        } else if (is_world_gradient(i)) {
          tex = std::make_shared<MicrofacetTransmission>(std::make_shared<world_gradient_texture>(point3f(temp_gradient_control(0),temp_gradient_control(1),temp_gradient_control(2)),
                                                                                                point3f(temp_gradient_control(3),temp_gradient_control(4),temp_gradient_control(5)),
                                                                                                point3f(tempvector(0),tempvector(1),tempvector(2)),
                                                                                                point3f(tempgradient(0),tempgradient(1),tempgradient(2)), gradient_is_hsv(i)), dist, 
                                                                                                point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                                                                                                point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        } else {
          tex = std::make_shared<MicrofacetTransmission>(std::make_shared<constant_texture>(point3f(tempvector(0),tempvector(1),tempvector(2))), dist, 
                                                       point3f(temp_glossy(3), temp_glossy(4), temp_glossy(5)), 
                                                       point3f(temp_glossy(6),temp_glossy(7),temp_glossy(8)));
        }
      }
    }
    if(is_shared_mat(i) && shared_materials->size() < static_cast<size_t>(shared_id_mat(i)) ) {
      shared_materials->push_back(tex);
    }
    //Generate center vector
    if(shape(i) == 1) {
      center =  vec3f(x(i), y(i), z(i));
    } else if(shape(i) == 2) {
      center =  vec3f(tempvector(prop_len+1),tempvector(prop_len+3),tempvector(prop_len+5));
    } else if(shape(i) == 3) {
      center =  vec3f(tempvector(prop_len+1),tempvector(prop_len+5), tempvector(prop_len+3));
    } else if(shape(i) == 4) {
      center =  vec3f(tempvector(prop_len+5),tempvector(prop_len+1),tempvector(prop_len+3));
    } else if(shape(i) == 5) {
      center =  vec3f(x(i), y(i), z(i));
    } else if(shape(i) == 6) {
      center =  vec3f(0, 0, 0);
    } else if(shape(i) == 7) {
      center = vec3f(x(i), y(i), z(i));
    } else if(shape(i) == 8) {
      center = vec3f(x(i), y(i), z(i));
    } else if(shape(i) == 9) {
      center = vec3f(x(i), y(i), z(i));
    } else if(shape(i) == 10) {
      center = vec3f(x(i), y(i), z(i));
    } else if(shape(i) == 11) {
      center = vec3f(x(i), y(i), z(i));
    } else if(shape(i) == 12) {
      center = vec3f(x(i), y(i), z(i));
    } else if(shape(i) == 13) {
      center = vec3f(x(i), y(i), z(i));
    } else if(shape(i) == 14) {
      center = vec3f(x(i), y(i), z(i));
    } else if(shape(i) == 15) {
      center = vec3f(x(i), y(i), z(i));
    } else if(shape(i) == 16) {
      center = vec3f(x(i), y(i), z(i));
    } else if(shape(i) == 17) {
      center = vec3f(x(i), y(i), z(i));
    }
    
    Transform GroupTransform(temp_group_transform);
    Transform TempM = GroupTransform * Translate(center) *
      rotation_order_matrix(temprotvec, order_rotation) * 
      Scale(temp_scales[0], temp_scales[1], temp_scales[2]);
    std::shared_ptr<Transform> ObjToWorld = transformCache.Lookup(TempM);
    std::shared_ptr<Transform> WorldToObj = transformCache.Lookup(TempM.GetInverseMatrix());
    
    
    Transform AnimationStartTransform(temp_animation_transform_start);
    Transform AnimationEndTransform(temp_animation_transform_end);
    std::shared_ptr<Transform> StartAnim = transformCache.Lookup(AnimationStartTransform);
    std::shared_ptr<Transform> EndAnim = transformCache.Lookup(AnimationEndTransform);
    
    AnimatedTransform Animate(StartAnim, animation_start_time(i), 
                              EndAnim, animation_end_time(i));
    
    //Generate objects
    std::shared_ptr<hitable> entry;
    if (shape(i) == 1) {
      entry = std::make_shared<sphere>(radius(i), tex, alpha[i], bump[i],
                                       ObjToWorld, WorldToObj, isflipped(i));
      if(isvolume(i)) {
        entry = std::make_shared<constant_medium>(entry, voldensity(i), 
                                                  std::make_shared<constant_texture>(point3f(tempvector(0),tempvector(1),tempvector(2))));
      }
      if(has_animation(i)) {
        entry = std::make_shared<AnimatedHitable>(entry, Animate);
      }
      list.add(entry);
    } else if (shape(i)  == 2) {
      entry = std::make_shared<xy_rect>(-tempvector(prop_len+2)/2,tempvector(prop_len+2)/2,
                                   -tempvector(prop_len+4)/2,tempvector(prop_len+4)/2,
                                   0, tex, alpha[i], bump[i], 
                                   ObjToWorld,WorldToObj, isflipped(i));
      if(has_animation(i)) {
        entry = std::make_shared<AnimatedHitable>(entry, Animate);
      }
      list.add(entry);
    } else if (shape(i)  == 3) {
      entry = std::make_shared<xz_rect>(-tempvector(prop_len+2)/2,tempvector(prop_len+2)/2,
                                   -tempvector(prop_len+4)/2,tempvector(prop_len+4)/2,
                                   0, tex, alpha[i], bump[i],
                                   ObjToWorld,WorldToObj, isflipped(i));
      if(has_animation(i)) {
        entry = std::make_shared<AnimatedHitable>(entry, Animate);
      }
      list.add(entry);
    } else if (shape(i)  == 4) {
      entry = std::make_shared<yz_rect>(-tempvector(prop_len+2)/2,tempvector(prop_len+2)/2,
                                   -tempvector(prop_len+4)/2,tempvector(prop_len+4)/2,
                                   0, tex, alpha[i], bump[i], 
                                   ObjToWorld,WorldToObj, isflipped(i));
      if(has_animation(i)) {
        entry = std::make_shared<AnimatedHitable>(entry, Animate);
      }
      list.add(entry);
    } else if (shape(i)  == 5) {
      entry = std::make_shared<box>(-vec3f(tempvector(prop_len+1),tempvector(prop_len+2),tempvector(prop_len+3))/2, 
                               vec3f(tempvector(prop_len+1),tempvector(prop_len+2),tempvector(prop_len+3))/2, 
                               tex, alpha[i], bump[i],
                               ObjToWorld,WorldToObj, isflipped(i));
      
      if(isvolume(i)) {
        if(!isnoise(i)) {
          entry = std::make_shared<constant_medium>(entry, voldensity(i), 
                                                    std::make_shared<constant_texture>(point3f(tempvector(0),tempvector(1),tempvector(2))));
        } else {
          entry = std::make_shared<constant_medium>(entry, voldensity(i), 
                                                   std::make_shared<noise_texture>(noise(i),
                                                        point3f(tempvector(0),tempvector(1),tempvector(2)),
                                                        point3f(tempnoisecolor(0),tempnoisecolor(1),tempnoisecolor(2)),
                                                        noisephase(i), noiseintensity(i)));
        }
      } 
      if(has_animation(i)) {
        entry = std::make_shared<AnimatedHitable>(entry, Animate);
      }
      list.add(entry);
    } else if (shape(i) == 8) {
      std::string objfilename = Rcpp::as<std::string>(fileinfo(i));
      std::string objbasedirname = Rcpp::as<std::string>(filebasedir(i));
      entry = std::make_shared<trimesh>(objfilename, objbasedirname, 
                                        tempvector(prop_len+1), sigma(i), tex, 
                                        (bool)tempvector(prop_len+3), (bool)tempvector(prop_len+2), (bool)tempvector(prop_len+4),
                                        (bool)tempvector(prop_len+5), (bool)tempvector(prop_len+6), (bool)tempvector(prop_len+7),
                                        imp_sample_objects,
                                        shutteropen, shutterclose, bvh_type, rng, verbose,
                                        ObjToWorld,WorldToObj, isflipped(i));
      if(isvolume(i)) {
        entry = std::make_shared<constant_medium>(entry, voldensity(i), 
                                                  std::make_shared<constant_texture>(point3f(tempvector(0),tempvector(1),tempvector(2))));
      }
      if(has_animation(i)) {
        entry = std::make_shared<AnimatedHitable>(entry, Animate);
      }
      list.add(entry);
    } else if (shape(i) == 9) {
      entry = std::make_shared<disk>(vec3f(0,0,0), radius(i), tempvector(prop_len+1), tex, alpha[i], bump[i],
                                     ObjToWorld,WorldToObj, isflipped(i));
      if(has_animation(i)) {
        entry = std::make_shared<AnimatedHitable>(entry, Animate);
      }
      list.add(entry);
    } else if (shape(i) == 10) {
      bool has_caps = type(i) != 5 && type(i) != 8;
      bool has_cap_option = tempvector(prop_len+4) == 1;
      entry = std::make_shared<cylinder>(radius(i), tempvector(prop_len+1), 
                                    tempvector(prop_len+2), tempvector(prop_len+3),has_caps && has_cap_option,
                                    tex, alpha[i], bump[i],
                                    ObjToWorld,WorldToObj, isflipped(i));
      if(has_animation(i)) {
        entry = std::make_shared<AnimatedHitable>(entry, Animate);
      }
      list.add(entry);
    } else if (shape(i) == 11) {
      entry = std::make_shared<ellipsoid>(vec3f(0,0,0), radius(i), 
                                      vec3f(tempvector(prop_len + 1), tempvector(prop_len + 2), tempvector(prop_len + 3)),
                                      tex, alpha[i], bump[i],
                                      ObjToWorld,WorldToObj, isflipped(i));
      if(isvolume(i)) {
        entry = std::make_shared<constant_medium>(entry, voldensity(i), 
                                                  std::make_shared<constant_texture>(point3f(tempvector(0),tempvector(1),tempvector(2))));
      }
      if(has_animation(i)) {
        entry = std::make_shared<AnimatedHitable>(entry, Animate);
      }
      list.add(entry);
    } else if (shape(i) == 13) {
      entry = std::make_shared<cone>(radius(i), tempvector(prop_len+1), tex, alpha[i], bump[i],
                                                              ObjToWorld,WorldToObj, isflipped(i));
      if(has_animation(i)) {
        entry = std::make_shared<AnimatedHitable>(entry, Animate);
      }
      list.add(entry);
    } else if (shape(i) == 14) {
      int pl = prop_len;
      vec3f p[4];
      vec3f n[2];
      CurveType type_curve;
      if(tempvector(pl+17) == 1) {
        type_curve = CurveType::Flat;
      } else if (tempvector(pl+17) == 2) {
        type_curve = CurveType::Cylinder;
      } else {
        type_curve = CurveType::Ribbon;
      }
      p[0] = vec3f(tempvector(pl+1),tempvector(pl+2),tempvector(pl+3));
      p[1] = vec3f(tempvector(pl+4),tempvector(pl+5),tempvector(pl+6));
      p[2] = vec3f(tempvector(pl+7),tempvector(pl+8),tempvector(pl+9));
      p[3] = vec3f(tempvector(pl+10),tempvector(pl+11),tempvector(pl+12));
      
      n[0] = vec3f(tempvector(pl+18), tempvector(pl+19), tempvector(pl+20));
      n[1] = vec3f(tempvector(pl+21), tempvector(pl+22), tempvector(pl+23));
      
      std::shared_ptr<CurveCommon> curve_data = std::make_shared<CurveCommon>(p, tempvector(pl+13), tempvector(pl+14), type_curve, n);
      entry = std::make_shared<curve>(tempvector(pl+15),tempvector(pl+16), curve_data, tex,
                                                               ObjToWorld, WorldToObj, isflipped(i));
      if(has_animation(i)) {
        entry = std::make_shared<AnimatedHitable>(entry, Animate);
      }
      list.add(entry);
    } else if (shape(i) == 15) {
      List temp_csg = csg_list(i);
      std::shared_ptr<ImplicitShape> shapes = parse_csg(temp_csg);
      entry = std::make_shared<csg>(tex, shapes,
                                                             ObjToWorld,WorldToObj, isflipped(i));
      if(has_animation(i)) {
        entry = std::make_shared<AnimatedHitable>(entry, Animate);
      }
      list.add(entry);
    } else if (shape(i) == 16) {
      std::string objfilename = Rcpp::as<std::string>(fileinfo(i));
      std::string objbasedirname = Rcpp::as<std::string>(filebasedir(i));
      entry = std::make_shared<plymesh>(objfilename, objbasedirname, 
                                        tex, alpha[i], bump[i], 
                                        tempvector(prop_len+1),
                                        shutteropen, shutterclose, bvh_type, rng,
                                        ObjToWorld,WorldToObj, isflipped(i));
      if(entry == nullptr) {
        continue;
      }
      if(isvolume(i)) {
        entry = std::make_shared<constant_medium>(entry, voldensity(i), 
                                                  std::make_shared<constant_texture>(point3f(tempvector(0),tempvector(1),tempvector(2))));
      }
      if(has_animation(i)) {
        entry = std::make_shared<AnimatedHitable>(entry, Animate);
      }
      list.add(entry);
    } else if (shape(i) == 17) {
      List mesh_entry = mesh_list(i);
      entry = std::make_shared<mesh3d>(mesh_entry, tex,
                                       shutteropen, shutterclose, bvh_type, rng,
                                       ObjToWorld,WorldToObj, isflipped(i),
                                       prop_len, tempvector, temp_glossy, sigma(i),
                                       lightintensity(i));
      if(has_animation(i)) {
        entry = std::make_shared<AnimatedHitable>(entry, Animate);
      }
      list.add(entry);
    }
    if(implicit_sample(i)) {
      imp_sample_objects.add(entry);
    }
  }
  auto world_bvh = std::make_shared<bvh_node>(list, shutteropen, shutterclose, bvh_type, rng);
  // auto nodeleaf = world_bvh->CountNodeLeaf();
  // Rcpp::Rcout << "Node/Leaf: " << nodeleaf.first << " " << nodeleaf.second << " " << world_bvh->GetSize() << "\n";
  return(world_bvh);
}

