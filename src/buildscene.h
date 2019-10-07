#ifndef BUILDSCENEH
#define BUILDSCENEH

#include "hitable.h"
#include "sphere.h"
#include "hitablelist.h"
#include "bvh_node.h"
#include "perlin.h"
#include "texture.h"
#include "xyrect.h"
#include "box.h"
#include "constant.h"
#include "triangle.h"
#include "pdf.h"
#include "trimesh.h"
#include "disk.h"
#include "cylinder.h"
#include "ellipsoid.h"
#include <Rcpp.h>
using namespace Rcpp;


hitable* rotation_order(hitable* entry, NumericVector temprotvec, NumericVector order_rotation) {
  for(int i = 0; i < 3; i++) {
    if(order_rotation(i) == 1) {
      if(temprotvec(0) != 0) {
        entry = new rotate_x(entry,temprotvec(0));
      }
    }
    if(order_rotation(i) == 2) {
      if(temprotvec(1) != 0) {
        entry = new rotate_y(entry,temprotvec(1));
      }
    }
    if(order_rotation(i) == 3) {
      if(temprotvec(2) != 0) {
        entry = new rotate_z(entry,temprotvec(2));
      }
    }
  }
  return(entry);
}


hitable *build_scene(IntegerVector& type, 
                     NumericVector& radius, IntegerVector& shape,
                     NumericVector& x, NumericVector& y, NumericVector& z,
                     List& properties, List& velocity, LogicalVector& moving,
                     int n, Float shutteropen, Float shutterclose,
                     LogicalVector& ischeckered, List& checkercolors, 
                     NumericVector& noise, LogicalVector& isnoise,
                     NumericVector& noisephase, NumericVector& noiseintensity, List noisecolorlist,
                     List& angle, 
                     LogicalVector& isimage, CharacterVector& filelocation,
                     LogicalVector& islight, NumericVector& lightintensity,
                     LogicalVector& isflipped,
                     LogicalVector& isvolume, NumericVector& voldensity,
                     List& order_rotation_list, 
                     LogicalVector& isgrouped, List& group_pivot, List& group_translate,
                     List& group_angle, List& group_order_rotation, List& group_scale,
                     LogicalVector& tri_normal_bools, LogicalVector& is_tri_color, List& tri_color_vert, 
                     CharacterVector& fileinfo, CharacterVector& filebasedir,
                     List& scale_list, random_gen& rng) {
  hitable **list = new hitable*[n+1];
  NumericVector tempvector;
  NumericVector tempchecker;
  NumericVector tempvel;
  NumericVector tempnoisecolor;
  NumericVector temprotvec;
  NumericVector order_rotation;
  NumericVector temp_gpivot;
  NumericVector temp_gtrans;
  NumericVector temp_gorder;
  NumericVector temp_gangle;
  NumericVector temp_gscale;
  NumericVector temp_tri_color;
  NumericVector temp_scales;
  vec3 gpivot;
  vec3 gtrans;
  vec3 gorder;
  vec3 gangle;
  int prop_len;
  
  List templist;
  vec3 center(x(0), y(0), z(0));
  vec3 vel(x(0), y(0), z(0));
  for(int i = 0; i < n; i++) {
    tempvector = as<NumericVector>(properties(i));
    tempchecker = as<NumericVector>(checkercolors(i));
    tempvel = as<NumericVector>(velocity(i));
    tempnoisecolor = as<NumericVector>(noisecolorlist(i));
    temprotvec = as<NumericVector>(angle(i));
    temp_tri_color = as<NumericVector>(tri_color_vert(i));
    order_rotation = as<NumericVector>(order_rotation_list(i));
    temp_scales = as<NumericVector>(scale_list(i));
    bool is_scaled = false;
    bool is_group_scaled = false;
    if(temp_scales[0] != 1 || temp_scales[1] != 1 || temp_scales[2] != 1) {
      is_scaled = true;
    }
    if(isgrouped(i)) {
      temp_gpivot = as<NumericVector>(group_pivot(i));
      temp_gangle = as<NumericVector>(group_angle(i));
      temp_gorder = as<NumericVector>(group_order_rotation(i));
      temp_gtrans = as<NumericVector>(group_translate(i));
      temp_gscale = as<NumericVector>(group_scale(i));
      if(temp_gscale[0] != 1 || temp_gscale[1] != 1 || temp_gscale[2] != 1) {
        is_group_scaled = true;
      }
      gpivot = vec3(temp_gpivot(0),temp_gpivot(1),temp_gpivot(2));
      gtrans = vec3(temp_gtrans(0),temp_gtrans(1),temp_gtrans(2));
    } else {
      gpivot = vec3(0,0,0); 
      gtrans = vec3(0,0,0); 
    }
    prop_len=2;
    vel = vec3(tempvel(0),tempvel(1),tempvel(2));
    //Generate texture
    material *tex;
    if(type(i) == 1) {
      if(isimage(i)) {
        int nx, ny, nn;
        unsigned char *tex_data = stbi_load(filelocation(i), &nx, &ny, &nn, 4);
        tex = new lambertian(new image_texture(tex_data,nx,ny,nn));
      } else if (islight(i)) {
        tex = new diffuse_light(new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2))*lightintensity(i)) );
      } else if (isnoise(i)) {
        tex = new lambertian(new noise_texture(noise(i),vec3(tempvector(0),tempvector(1),tempvector(2)),
                                               vec3(tempnoisecolor(0),tempnoisecolor(1),tempnoisecolor(2)),
                                               noisephase(i), noiseintensity(i)));
      } else if (ischeckered(i)) {
        tex = new lambertian(new checker_texture(new constant_texture(vec3(tempchecker(0),tempchecker(1),tempchecker(2))),
                                                 new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2))),tempchecker(3)));
      } else if (is_tri_color(i)) {
        tex = new lambertian(new triangle_texture(vec3(temp_tri_color(0),temp_tri_color(1),temp_tri_color(2)),
                                                  vec3(temp_tri_color(3),temp_tri_color(4),temp_tri_color(5)),
                                                  vec3(temp_tri_color(6),temp_tri_color(7),temp_tri_color(8))) );
      } else {
        tex = new lambertian(new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2))) );
      }
    } 
    else if (type(i) == 2) {
      tex = new metal(vec3(tempvector(0),tempvector(1),tempvector(2)),tempvector(3));
      prop_len = 3;
    } else {
      tex = new dielectric(vec3(tempvector(0),tempvector(1),tempvector(2)),tempvector(3), rng);
      prop_len = 3;
    }
    //Generate center vector
    if(shape(i) == 1) {
      center =  vec3(x(i), y(i), z(i));
    } else if(shape(i) == 2) {
      center =  vec3(tempvector(prop_len+1),tempvector(prop_len+3),tempvector(prop_len+5));
    } else if(shape(i) == 3) {
      center =  vec3(tempvector(prop_len+1),tempvector(prop_len+5), tempvector(prop_len+3));
    } else if(shape(i) == 4) {
      center =  vec3(tempvector(prop_len+5),tempvector(prop_len+1),tempvector(prop_len+3));
    } else if(shape(i) == 5) {
      center =  vec3(x(i), y(i), z(i));
    } else if(shape(i) == 6) {
      center =  vec3(0, 0, 0);
    } else if(shape(i) == 7) {
      center = vec3(x(i), y(i), z(i));
    } else if(shape(i) == 8) {
      center = vec3(x(i), y(i), z(i));
    } else if(shape(i) == 9) {
      center = vec3(x(i), y(i), z(i));
    } else if(shape(i) == 10) {
      center = vec3(x(i), y(i), z(i));
    } else if(shape(i) == 11) {
      center = vec3(x(i), y(i), z(i));
    }
    //Generate objects
    if (shape(i) == 1) {
      hitable *entry = new sphere(vec3(0,0,0), radius(i), tex);
      if(is_scaled) {
        entry = new scale(entry, vec3(temp_scales[0], temp_scales[1], temp_scales[2]));
      }
      entry = rotation_order(entry, temprotvec, order_rotation);
      if(isgrouped(i)) {
        entry = new translate(entry, center - gpivot);
        if(is_group_scaled) {
          entry = new scale(entry, vec3(temp_gscale[0], temp_gscale[1], temp_gscale[2]));
        }
        entry = rotation_order(entry, temp_gangle, temp_gorder);
        entry = new translate(entry, -center + gpivot );
      }
      if(!moving(i)) {
        entry = new translate(entry, center + gtrans + vel * shutteropen);
      } else {
        entry = new moving_sphere(center + vel * shutteropen, 
                                  center + vel * shutterclose, 
                                  shutteropen, shutterclose, radius(i), tex);
      }
      if(isflipped(i)) {
        entry = new flip_normals(entry);
      }
      if(isvolume(i)) {
        list[i] = new constant_medium(entry, voldensity(i), new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2))));
      } else {
        list[i] = entry;
      }
    } else if (shape(i)  == 2) {
      hitable *entry = new xy_rect(-tempvector(prop_len+2)/2,tempvector(prop_len+2)/2,
                                   -tempvector(prop_len+4)/2,tempvector(prop_len+4)/2,
                                   0, tex);
      if(is_scaled) {
        entry = new scale(entry, vec3(temp_scales[0], temp_scales[1], temp_scales[2]));
      }
      entry = rotation_order(entry, temprotvec, order_rotation);
      if(isgrouped(i)) {
        entry = new translate(entry, center - gpivot);
        if(is_group_scaled) {
          entry = new scale(entry, vec3(temp_gscale[0], temp_gscale[1], temp_gscale[2]));
        }
        entry = rotation_order(entry, temp_gangle, temp_gorder);
        entry = new translate(entry, -center + gpivot );
      }
      entry = new translate(entry,center + gtrans + vel * shutteropen);
      if(isflipped(i)) {
        list[i] = new flip_normals(entry);
      } else {
        list[i] = entry;
      }
    } else if (shape(i)  == 3) {
      hitable *entry = new xz_rect(-tempvector(prop_len+2)/2,tempvector(prop_len+2)/2,
                                   -tempvector(prop_len+4)/2,tempvector(prop_len+4)/2,
                                   0, tex);
      if(is_scaled) {
        entry = new scale(entry, vec3(temp_scales[0], temp_scales[1], temp_scales[2]));
      }
      entry = rotation_order(entry, temprotvec, order_rotation);
      if(isgrouped(i)) {
        entry = new translate(entry, center - gpivot);
        if(is_group_scaled) {
          entry = new scale(entry, vec3(temp_gscale[0], temp_gscale[1], temp_gscale[2]));
        }
        entry = rotation_order(entry, temp_gangle, temp_gorder);
        entry = new translate(entry, -center + gpivot );
      }
      entry = new translate(entry,center + gtrans + vel * shutteropen);
      if(isflipped(i)) {
        list[i] = new flip_normals(entry);
      } else {
        list[i] = entry;
      }
    } else if (shape(i)  == 4) {
      hitable *entry = new yz_rect(-tempvector(prop_len+2)/2,tempvector(prop_len+2)/2,
                                   -tempvector(prop_len+4)/2,tempvector(prop_len+4)/2,
                                   0, tex);
      if(is_scaled) {
        entry = new scale(entry, vec3(temp_scales[0], temp_scales[1], temp_scales[2]));
      }
      entry = rotation_order(entry, temprotvec, order_rotation);
      if(isgrouped(i)) {
        entry = new translate(entry, center - gpivot);
        if(is_group_scaled) {
          entry = new scale(entry, vec3(temp_gscale[0], temp_gscale[1], temp_gscale[2]));
        }
        entry = rotation_order(entry, temp_gangle, temp_gorder);
        entry = new translate(entry, -center + gpivot );
      }
      entry = new translate(entry,center + gtrans + vel * shutteropen);
      if(isflipped(i)) {
        list[i] = new flip_normals(entry);
      } else {
        list[i] = entry;
      }
    } else if (shape(i)  == 5) {
      hitable *entry = new box(-vec3(tempvector(prop_len+1),tempvector(prop_len+2),tempvector(prop_len+3))/2, 
                               vec3(tempvector(prop_len+1),tempvector(prop_len+2),tempvector(prop_len+3))/2, 
                               tex);
      if(is_scaled) {
        entry = new scale(entry, vec3(temp_scales[0], temp_scales[1], temp_scales[2]));
      }
      entry = rotation_order(entry, temprotvec, order_rotation);
      if(isgrouped(i)) {
        entry = new translate(entry, center - gpivot);
        if(is_group_scaled) {
          entry = new scale(entry, vec3(temp_gscale[0], temp_gscale[1], temp_gscale[2]));
        }
        entry = rotation_order(entry, temp_gangle, temp_gorder);
        entry = new translate(entry, -center + gpivot );
      }
      entry = new translate(entry, center + gtrans + vel * shutteropen);
      if(isvolume(i)) {
        if(!isnoise(i)) {
          list[i] = new constant_medium(entry,voldensity(i), new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2))));
        } else {
          list[i] = new constant_medium(entry,voldensity(i), 
                                        new noise_texture(noise(i),
                                                          vec3(tempvector(0),tempvector(1),tempvector(2)),
                                                          vec3(tempnoisecolor(0),tempnoisecolor(1),tempnoisecolor(2)),
                                                          noisephase(i), noiseintensity(i)));
        }
      } else {
        list[i] = entry;
      }
    } else if (shape(i)  == 6) {
      hitable *entry;
      if(tri_normal_bools(i)) {
        entry= new triangle(vec3(tempvector(prop_len+1),tempvector(prop_len+2),tempvector(prop_len+3)),
                            vec3(tempvector(prop_len+4),tempvector(prop_len+5),tempvector(prop_len+6)),
                            vec3(tempvector(prop_len+7),tempvector(prop_len+8),tempvector(prop_len+9)),
                            vec3(tempvector(prop_len+10),tempvector(prop_len+11),tempvector(prop_len+12)),
                            vec3(tempvector(prop_len+13),tempvector(prop_len+14),tempvector(prop_len+15)),
                            vec3(tempvector(prop_len+16),tempvector(prop_len+17),tempvector(prop_len+18)),
                            tex);
      } else {
        entry= new triangle(vec3(tempvector(prop_len+1),tempvector(prop_len+2),tempvector(prop_len+3)),
                            vec3(tempvector(prop_len+4),tempvector(prop_len+5),tempvector(prop_len+6)),
                            vec3(tempvector(prop_len+7),tempvector(prop_len+8),tempvector(prop_len+9)),
                            tex);
      }
      if(is_scaled) {
        entry = new scale(entry, vec3(temp_scales[0], temp_scales[1], temp_scales[2]));
      }
      entry = rotation_order(entry, temprotvec, order_rotation);
      if(isgrouped(i)) {
        entry = new translate(entry, center - gpivot);
        if(is_group_scaled) {
          entry = new scale(entry, vec3(temp_gscale[0], temp_gscale[1], temp_gscale[2]));
        }
        entry = rotation_order(entry, temp_gangle, temp_gorder);
        entry = new translate(entry, -center + gpivot );
      }
      entry = new translate(entry, center + gtrans + vel * shutteropen);
      if(isflipped(i)) {
        list[i] = new flip_normals(entry);
      } else {
        list[i] = entry;
      }
    } else if (shape(i) == 7) {
      hitable *entry;
      std::string objfilename = Rcpp::as<std::string>(fileinfo(i));
      std::string objbasedirname = Rcpp::as<std::string>(filebasedir(i));
      entry = new trimesh(objfilename, objbasedirname, 
                         tex,
                         tempvector(prop_len+1),
                         shutteropen, shutterclose, rng);
      if(is_scaled) {
        entry = new scale(entry, vec3(temp_scales[0], temp_scales[1], temp_scales[2]));
      }
      entry = rotation_order(entry, temprotvec, order_rotation);
      if(isgrouped(i)) {
        entry = new translate(entry, center - gpivot);
        if(is_group_scaled) {
          entry = new scale(entry, vec3(temp_gscale[0], temp_gscale[1], temp_gscale[2]));
        }
        entry = rotation_order(entry, temp_gangle, temp_gorder);
        entry = new translate(entry, -center + gpivot );
      }
      entry = new translate(entry, center + gtrans + vel * shutteropen);
      if(isflipped(i)) {
        entry = new flip_normals(entry);
      } 
      if(isvolume(i)) {
        list[i] = new constant_medium(entry, voldensity(i), 
                                      new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2))));
      } else {
        list[i] = entry;
      }
    } else if (shape(i) == 8) {
      hitable *entry;
      std::string objfilename = Rcpp::as<std::string>(fileinfo(i));
      std::string objbasedirname = Rcpp::as<std::string>(filebasedir(i));
      entry = new trimesh(objfilename, objbasedirname, 
                          tempvector(prop_len+1),
                          shutteropen, shutterclose, rng);
      if(is_scaled) {
        entry = new scale(entry, vec3(temp_scales[0], temp_scales[1], temp_scales[2]));
      }
      entry = rotation_order(entry, temprotvec, order_rotation);
      if(isgrouped(i)) {
        entry = new translate(entry, center - gpivot);
        if(is_group_scaled) {
          entry = new scale(entry, vec3(temp_gscale[0], temp_gscale[1], temp_gscale[2]));
        }
        entry = rotation_order(entry, temp_gangle, temp_gorder);
        entry = new translate(entry, -center + gpivot );
      }
      entry = new translate(entry, center + gtrans + vel * shutteropen);
      if(isflipped(i)) {
        entry = new flip_normals(entry);
      } 
      if(isvolume(i)) {
        list[i] = new constant_medium(entry, voldensity(i), 
                                      new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2))));
      } else {
        list[i] = entry;
      }
    } else if (shape(i) == 9) {
      hitable *entry;
      entry = new disk(vec3(0,0,0), radius(i), tempvector(prop_len+1), tex);
      if(is_scaled) {
        entry = new scale(entry, vec3(temp_scales[0], temp_scales[1], temp_scales[2]));
      }
      entry = rotation_order(entry, temprotvec, order_rotation);
      if(isgrouped(i)) {
        entry = new translate(entry, center - gpivot);
        if(is_group_scaled) {
          entry = new scale(entry, vec3(temp_gscale[0], temp_gscale[1], temp_gscale[2]));
        }
        entry = rotation_order(entry, temp_gangle, temp_gorder);
        entry = new translate(entry, -center + gpivot );
      }
      entry = new translate(entry, center + gtrans + vel * shutteropen);
      if(isflipped(i)) {
        list[i] = new flip_normals(entry);
      } else {
        list[i] = entry;
      }
    } else if (shape(i) == 10) {
      hitable *entry = new cylinder(radius(i), tempvector(prop_len+1), 
                                    tempvector(prop_len+2), tempvector(prop_len+3), tex);
      if(is_scaled) {
        entry = new scale(entry, vec3(temp_scales[0], temp_scales[1], temp_scales[2]));
      }
      entry = rotation_order(entry, temprotvec, order_rotation);
      if(isgrouped(i)) {
        entry = new translate(entry, center - gpivot);
        if(is_group_scaled) {
          entry = new scale(entry, vec3(temp_gscale[0], temp_gscale[1], temp_gscale[2]));
        }
        entry = rotation_order(entry, temp_gangle, temp_gorder);
        entry = new translate(entry, -center + gpivot );
      }
      entry = new translate(entry, center + gtrans + vel * shutteropen);
      if(isflipped(i)) {
        entry = new flip_normals(entry);
      }
      list[i] = entry;
    } else if (shape(i) == 11) {
      hitable *entry = new ellipsoid(vec3(0,0,0), radius(i), 
                                     vec3(tempvector(prop_len + 1), tempvector(prop_len + 2), tempvector(prop_len + 3)),
                                     tex);
      if(is_scaled) {
        entry = new scale(entry, vec3(temp_scales[0], temp_scales[1], temp_scales[2]));
      }
      entry = rotation_order(entry, temprotvec, order_rotation);
      if(isgrouped(i)) {
        entry = new translate(entry, center - gpivot);
        if(is_group_scaled) {
          entry = new scale(entry, vec3(temp_gscale[0], temp_gscale[1], temp_gscale[2]));
        }
        entry = rotation_order(entry, temp_gangle, temp_gorder);
        entry = new translate(entry, -center + gpivot );
      }
      entry = new translate(entry, center + gtrans + vel * shutteropen);
      if(isflipped(i)) {
        entry = new flip_normals(entry);
      }
      if(isvolume(i)) {
        list[i] = new constant_medium(entry, voldensity(i), 
                                      new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2))));
      } else {
        list[i] = entry;
      }
    }
  }
  return(new bvh_node(list, n, shutteropen, shutterclose, rng));
}

hitable* build_imp_sample(IntegerVector& type, 
                          NumericVector& radius, IntegerVector& shape,
                          NumericVector& x, NumericVector& y, NumericVector& z,
                          List& properties, List& velocity, 
                          int n, Float shutteropen, Float shutterclose, 
                          List& angle, int i, List& order_rotation_list,
                          LogicalVector& isgrouped, 
                          List& group_pivot, List& group_translate,
                          List& group_angle, List& group_order_rotation, List& group_scale,
                          CharacterVector& fileinfo, CharacterVector& filebasedir,
                          List& scale_list, random_gen& rng) {
  NumericVector tempvector;
  NumericVector temprotvec;
  NumericVector tempvel;
  NumericVector order_rotation;
  NumericVector temp_gpivot;
  NumericVector temp_gtrans;
  NumericVector temp_gorder;
  NumericVector temp_gangle;
  NumericVector temp_gscale;
  NumericVector temp_scales;
  vec3 gpivot;
  vec3 gtrans;
  vec3 gorder;
  vec3 gangle;
  
  List templist;
  vec3 center(x(0), y(0), z(0));
  vec3 vel(x(0), y(0), z(0));
  int prop_len;
  tempvector = as<NumericVector>(properties(i));
  tempvel = as<NumericVector>(velocity(i));
  vel = vec3(tempvel(0),tempvel(1),tempvel(2));
  temprotvec =  as<NumericVector>(angle(i));
  order_rotation = as<NumericVector>(order_rotation_list(i));
  temp_scales = as<NumericVector>(scale_list(i));
  bool is_scaled = false;
  bool is_group_scaled = false;
  if(temp_scales[0] != 1 || temp_scales[1] != 1 || temp_scales[2] != 1) {
    is_scaled = true;
  }
  prop_len=2;
  center =  vec3(x(i), y(i), z(i));
  if(isgrouped(i)) {
    temp_gpivot = as<NumericVector>(group_pivot(i));
    temp_gangle = as<NumericVector>(group_angle(i));
    temp_gorder = as<NumericVector>(group_order_rotation(i));
    temp_gtrans = as<NumericVector>(group_translate(i));
    temp_gscale = as<NumericVector>(group_scale(i));
    if(temp_gscale[0] != 1 || temp_gscale[1] != 1 || temp_gscale[2] != 1) {
      is_group_scaled = true;
    }
    gpivot = vec3(temp_gpivot(0),temp_gpivot(1),temp_gpivot(2));
    gtrans = vec3(temp_gtrans(0),temp_gtrans(1),temp_gtrans(2));
  } else{
    gpivot = vec3(0,0,0); 
    gtrans = vec3(0,0,0); 
  }
  if(type(i) != 1) {
    prop_len = 3;
  }
  
  if(shape(i) == 1) {
    center =  vec3(x(i), y(i), z(i));
  } else if(shape(i) == 2) {
    center =  vec3(tempvector(prop_len+1),tempvector(prop_len+3),tempvector(prop_len+5));
  } else if(shape(i) == 3) {
    center =  vec3(tempvector(prop_len+1),tempvector(prop_len+5), tempvector(prop_len+3));
  } else if(shape(i) == 4) {
    center =  vec3(tempvector(prop_len+5),tempvector(prop_len+1),tempvector(prop_len+3));
  } else if(shape(i) == 5) {
    center =  vec3(x(i), y(i), z(i));
  } else if(shape(i) == 6) {
    center =  vec3(x(i), y(i), z(i));
  } else if(shape(i) == 7) {
    center = vec3(x(i), y(i), z(i));
  } else if(shape(i) == 8) {
    center = vec3(x(i), y(i), z(i));
  } else if(shape(i) == 9) {
    center = vec3(x(i), y(i), z(i));
  } else if(shape(i) == 10) {
    center = vec3(x(i), y(i), z(i));
  } else if(shape(i) == 11) {
    center = vec3(x(i), y(i), z(i));
  }
  
  if(shape(i) == 1) {
    hitable *entry = new sphere(vec3(0,0,0), radius(i), 0);
    if(is_scaled) {
      entry = new scale(entry, vec3(temp_scales[0], temp_scales[1], temp_scales[2]));
    }
    if(isgrouped(i)) {
      entry = new translate(entry, center - gpivot);
      if(is_group_scaled) {
        entry = new scale(entry, vec3(temp_gscale[0], temp_gscale[1], temp_gscale[2]));
      }
      entry = rotation_order(entry, temp_gangle, temp_gorder);
      entry = new translate(entry, -center + gpivot );
    }
    entry = new translate(entry, center + gtrans + vel * shutteropen);
    return(entry);
  } else if (shape(i) == 2) {
    hitable *entry = new xy_rect(-tempvector(prop_len+2)/2,tempvector(prop_len+2)/2,
                                 -tempvector(prop_len+4)/2,tempvector(prop_len+4)/2,
                                 0, 0);
    if(is_scaled) {
      entry = new scale(entry, vec3(temp_scales[0], temp_scales[1], temp_scales[2]));
    }
    entry = rotation_order(entry, temprotvec, order_rotation);
    if(isgrouped(i)) {
      entry = new translate(entry, center - gpivot);
      if(is_group_scaled) {
        entry = new scale(entry, vec3(temp_gscale[0], temp_gscale[1], temp_gscale[2]));
      }
      entry = rotation_order(entry, temp_gangle, temp_gorder);
      entry = new translate(entry, -center + gpivot );
    }
    entry = new translate(entry, center + gtrans + vel * shutteropen);
    return(entry);
  } else if (shape(i) == 3) {
    hitable *entry = new xz_rect(-tempvector(prop_len+2)/2,tempvector(prop_len+2)/2,
                                 -tempvector(prop_len+4)/2,tempvector(prop_len+4)/2,
                                 0, 0);
    if(is_scaled) {
      entry = new scale(entry, vec3(temp_scales[0], temp_scales[1], temp_scales[2]));
    }
    entry = rotation_order(entry, temprotvec, order_rotation);
    if(isgrouped(i)) {
      entry = new translate(entry, center - gpivot);
      if(is_group_scaled) {
        entry = new scale(entry, vec3(temp_gscale[0], temp_gscale[1], temp_gscale[2]));
      }
      entry = rotation_order(entry, temp_gangle, temp_gorder);
      entry = new translate(entry, -center + gpivot );
    }
    entry = new translate(entry, center + gtrans + vel * shutteropen);
    return(entry);
  } else if (shape(i) == 4) {
    hitable *entry = new yz_rect(-tempvector(prop_len+2)/2,tempvector(prop_len+2)/2,
                                 -tempvector(prop_len+4)/2,tempvector(prop_len+4)/2,
                                 0, 0);
    if(is_scaled) {
      entry = new scale(entry, vec3(temp_scales[0], temp_scales[1], temp_scales[2]));
    }
    entry = rotation_order(entry, temprotvec, order_rotation);
    if(isgrouped(i)) {
      entry = new translate(entry, center - gpivot);
      if(is_group_scaled) {
        entry = new scale(entry, vec3(temp_gscale[0], temp_gscale[1], temp_gscale[2]));
      }
      entry = rotation_order(entry, temp_gangle, temp_gorder);
      entry = new translate(entry, -center + gpivot );
    }
    entry = new translate(entry, center + gtrans + vel * shutteropen);
    return(entry);
  } else if (shape(i) == 5) {
    hitable *entry = new box(-vec3(tempvector(prop_len+1),tempvector(prop_len+2),tempvector(prop_len+3))/2, 
                             vec3(tempvector(prop_len+1),tempvector(prop_len+2),tempvector(prop_len+3))/2, 
                             0);
    if(is_scaled) {
      entry = new scale(entry, vec3(temp_scales[0], temp_scales[1], temp_scales[2]));
    }
    entry = rotation_order(entry, temprotvec, order_rotation);
    if(isgrouped(i)) {
      entry = new translate(entry, center - gpivot);
      if(is_group_scaled) {
        entry = new scale(entry, vec3(temp_gscale[0], temp_gscale[1], temp_gscale[2]));
      }
      entry = rotation_order(entry, temp_gangle, temp_gorder);
      entry = new translate(entry, -center + gpivot );
    }
    entry = new translate(entry, center + gtrans + vel * shutteropen);
    return(entry);
  } else if (shape(i) == 6) {
    hitable *entry = new triangle(vec3(tempvector(prop_len+1),tempvector(prop_len+2),tempvector(prop_len+3)),
                                  vec3(tempvector(prop_len+4),tempvector(prop_len+5),tempvector(prop_len+6)),
                                  vec3(tempvector(prop_len+7),tempvector(prop_len+8),tempvector(prop_len+9)),
                                  0);
    if(is_scaled) {
      entry = new scale(entry, vec3(temp_scales[0], temp_scales[1], temp_scales[2]));
    }
    entry = rotation_order(entry, temprotvec, order_rotation);
    if(isgrouped(i)) {
      entry = new translate(entry, center - gpivot);
      if(is_group_scaled) {
        entry = new scale(entry, vec3(temp_gscale[0], temp_gscale[1], temp_gscale[2]));
      }
      entry = rotation_order(entry, temp_gangle, temp_gorder);
      entry = new translate(entry, -center + gpivot );
    }
    entry = new translate(entry, center + gtrans + vel * shutteropen);
    return(entry);
  } else if (shape(i) == 7 || shape(i) == 8) {
    hitable *entry;
    std::string objfilename = Rcpp::as<std::string>(fileinfo(i));
    std::string objbasedirname = Rcpp::as<std::string>(filebasedir(i));
    entry = new trimesh(objfilename, objbasedirname,
                        0,
                        tempvector(prop_len+1),
                        shutteropen, shutterclose, rng);
    if(is_scaled) {
      entry = new scale(entry, vec3(temp_scales[0], temp_scales[1], temp_scales[2]));
    }
    entry = rotation_order(entry, temprotvec, order_rotation);
    if(isgrouped(i)) {
      entry = new translate(entry, center - gpivot);
      if(is_group_scaled) {
        entry = new scale(entry, vec3(temp_gscale[0], temp_gscale[1], temp_gscale[2]));
      }
      entry = rotation_order(entry, temp_gangle, temp_gorder);
      entry = new translate(entry, -center + gpivot );
    }
    entry = new translate(entry, center + gtrans + vel * shutteropen);
    return(entry);
  } else if (shape(i) == 9) {
    hitable *entry;
    entry = new disk(vec3(0,0,0), radius(i), tempvector(prop_len+1), 0);
    if(is_scaled) {
      entry = new scale(entry, vec3(temp_scales[0], temp_scales[1], temp_scales[2]));
    }
    entry = rotation_order(entry, temprotvec, order_rotation);
    if(isgrouped(i)) {
      entry = new translate(entry, center - gpivot);
      if(is_group_scaled) {
        entry = new scale(entry, vec3(temp_gscale[0], temp_gscale[1], temp_gscale[2]));
      }
      entry = rotation_order(entry, temp_gangle, temp_gorder);
      entry = new translate(entry, -center + gpivot );
    }
    entry = new translate(entry, center + gtrans + vel * shutteropen);
    return(entry);
  } else if (shape(i) == 10) {
    hitable *entry = new cylinder(radius(i), tempvector(prop_len+1), 
                                  tempvector(prop_len+2), tempvector(prop_len+3), 0);
    if(is_scaled) {
      entry = new scale(entry, vec3(temp_scales[0], temp_scales[1], temp_scales[2]));
    }
    entry = rotation_order(entry, temprotvec, order_rotation);
    if(isgrouped(i)) {
      entry = new translate(entry, center - gpivot);
      if(is_group_scaled) {
        entry = new scale(entry, vec3(temp_gscale[0], temp_gscale[1], temp_gscale[2]));
      }
      entry = rotation_order(entry, temp_gangle, temp_gorder);
      entry = new translate(entry, -center + gpivot );
    }
    entry = new translate(entry, center + gtrans + vel * shutteropen);
    return(entry);
  } else {
    hitable *entry = new ellipsoid(vec3(0,0,0), radius(i), 
                                   vec3(tempvector(prop_len + 1), tempvector(prop_len + 2), tempvector(prop_len + 3)),
                                   0);
    if(is_scaled) {
      entry = new scale(entry, vec3(temp_scales[0], temp_scales[1], temp_scales[2]));
    }
    entry = rotation_order(entry, temprotvec, order_rotation);
    if(isgrouped(i)) {
      entry = new translate(entry, center - gpivot);
      if(is_group_scaled) {
        entry = new scale(entry, vec3(temp_gscale[0], temp_gscale[1], temp_gscale[2]));
      }
      entry = rotation_order(entry, temp_gangle, temp_gorder);
      entry = new translate(entry, -center + gpivot );
    }
    entry = new translate(entry, center + gtrans + vel * shutteropen);
    return(entry);
  }
}

#endif
