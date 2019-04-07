#define STB_IMAGE_IMPLEMENTATION 
#include "vec3.h"
#include "stb_image.h"
#include <Rcpp.h>
#include "sphere.h"
#include "hitablelist.h"
#include "camera.h"
#include "float.h"
#include "bvh_node.h"
#include "perlin.h"
#include "texture.h"
#include "xyrect.h"
#include "box.h"
#include "constant.h"
#include "triangle.h"
#include "pdf.h"
#include <RcppParallel.h>
#include "rng.h"
using namespace Rcpp;


inline vec3 de_nan(const vec3& c) {
  vec3 temp = c;
  if(std::isnan(c[0])) temp.e[0] = 0;
  if(std::isnan(c[1])) temp.e[1] = 0;
  if(std::isnan(c[2])) temp.e[2] = 0;
  return(temp);
}

inline vec3 clamp(const vec3& c, float clampval) {
  vec3 temp = c;
  if(c[0] > clampval) temp.e[0] = clampval;
  if(c[1] > clampval) temp.e[1] = clampval;
  if(c[2] > clampval) temp.e[2] = clampval;
  return(temp);
}

vec3 color(const ray& r, hitable *world, hitable *hlist, int depth, random_gen& rng) {
  hit_record hrec;
  if(world->hit(r, 0.001, FLT_MAX, hrec, rng)) {
    scatter_record srec;
    vec3 emitted = hrec.mat_ptr->emitted(r, hrec, hrec.u, hrec.v,hrec.p);
    float pdf_val;
    if(depth < 50 && hrec.mat_ptr->scatter(r, hrec, srec, rng)) {
      if(srec.is_specular) {
        return(srec.attenuation * 
               color(srec.specular_ray, world, 
                     hlist, depth + 1, rng));
      }
      hitable_pdf p_imp(hlist, hrec.p);
      mixture_pdf p(&p_imp, srec.pdf_ptr);
      ray scattered = ray(hrec.p, p.generate(rng), r.time());
      pdf_val = p.value(scattered.direction(), rng);
      return(emitted + srec.attenuation * 
             hrec.mat_ptr->scattering_pdf(r, hrec, scattered) *  
             color(scattered, world, 
                  hlist, depth + 1, rng) / pdf_val);
    } else {
      return(emitted);
    }
  } else {
    return(vec3(0,0,0));
  }
}

vec3 color_amb(const ray& r, hitable *world, hitable *hlist, int depth,
           const vec3& backgroundhigh, const vec3& backgroundlow, random_gen& rng) {
  hit_record hrec;
  if(world->hit(r, 0.001, FLT_MAX, hrec, rng)) {
    scatter_record srec;
    vec3 emitted = hrec.mat_ptr->emitted(r, hrec, hrec.u, hrec.v,hrec.p);
    float pdf_val;
    if(depth < 50 && hrec.mat_ptr->scatter(r, hrec, srec, rng)) {
      if(srec.is_specular) {
        return(srec.attenuation * 
               color_amb(srec.specular_ray, world, hlist, depth + 1, 
                         backgroundhigh,backgroundlow, rng));
      }
      hitable_pdf p_imp(hlist, hrec.p);
      mixture_pdf p(&p_imp, srec.pdf_ptr);
      ray scattered = ray(hrec.p, p.generate(rng), r.time());
      pdf_val = p.value(scattered.direction(), rng);
      return(emitted + srec.attenuation * 
             hrec.mat_ptr->scattering_pdf(r, hrec, scattered) *  
             color_amb(scattered, world, hlist, depth + 1, 
                       backgroundhigh,backgroundlow, rng) / pdf_val);
    } else {
      vec3 unit_direction = unit_vector(r.direction());
      float t = 0.5 * (unit_direction.y() + 1.0);
      return(emitted+(1.0 - t) * backgroundlow + t * backgroundhigh);
    }
  } else {
    vec3 unit_direction = unit_vector(r.direction());
    float t = 0.5 * (unit_direction.y() + 1.0);
    return (1.0 - t) * backgroundlow + t * backgroundhigh;
  }
}

vec3 color_uniform(const ray& r, hitable *world, int depth, random_gen& rng) {
  hit_record hrec;
  if(world->hit(r, 0.001, FLT_MAX, hrec, rng)) {
    scatter_record srec;
    vec3 emitted = hrec.mat_ptr->emitted(r, hrec, hrec.u, hrec.v,hrec.p);
    float pdf_val;
    if(depth < 50 && hrec.mat_ptr->scatter(r, hrec, srec, rng)) {
      if(srec.is_specular) {
        return(srec.attenuation * color_uniform(srec.specular_ray, world, depth + 1, rng));
      }
      cosine_pdf p(hrec.normal);
      ray scattered = ray(hrec.p, p.generate(rng), r.time());
      pdf_val = p.value(scattered.direction(), rng);
      return(emitted + srec.attenuation * hrec.mat_ptr->scattering_pdf(r, hrec, scattered) *  
             color_uniform(scattered, world, depth + 1, rng) / pdf_val);
    } else {
      return(emitted);
    }
  } else {
    return(vec3(0,0,0));
  }
}

vec3 color_amb_uniform(const ray& r, hitable *world, int depth,
               const vec3& backgroundhigh, const vec3& backgroundlow, random_gen& rng) {
  hit_record hrec;
  if(world->hit(r, 0.001, FLT_MAX, hrec, rng)) {
    scatter_record srec;

    vec3 emitted = hrec.mat_ptr->emitted(r, hrec, hrec.u, hrec.v,hrec.p);
    float pdf_val;
    if(depth < 50 && hrec.mat_ptr->scatter(r, hrec, srec, rng)) {
      if(srec.is_specular) {
        return(srec.attenuation * 
               color_amb_uniform(srec.specular_ray, world, depth + 1, 
                                 backgroundhigh,backgroundlow, rng));
      }
      cosine_pdf p(hrec.normal);
      ray scattered = ray(hrec.p, p.generate(rng), r.time());

      pdf_val = p.value(scattered.direction(), rng);
      return(emitted + srec.attenuation * 
             hrec.mat_ptr->scattering_pdf(r, hrec, scattered) *  
             color_amb_uniform(scattered, world, depth + 1, 
                               backgroundhigh,backgroundlow, rng) / pdf_val);
    } else {
      vec3 unit_direction = unit_vector(r.direction());
      float t = 0.5 * (unit_direction.y() + 1.0);
      return(emitted+(1.0 - t) * backgroundlow + t * backgroundhigh);
    }
  } else {
    vec3 unit_direction = unit_vector(r.direction());
    float t = 0.5 * (unit_direction.y() + 1.0);
    return (1.0 - t) * backgroundlow + t * backgroundhigh;
  }
}

struct Colorworker : public RcppParallel::Worker {
  Colorworker(NumericMatrix outputr, NumericMatrix outputg, NumericMatrix outputb,
              bool ambient_light, int nx, int ny, int ns, camera cam, vec3 backgroundhigh, vec3 backgroundlow, 
              hitable *world, hitable *hlist, int numbertosample, float clampval)
  : outputr(outputr), outputg(outputg), outputb(outputb), ambient_light(ambient_light),
    nx(nx), ny(ny), ns(ns), cam(cam),
    backgroundhigh(backgroundhigh), backgroundlow(backgroundlow), world(world), hlist(hlist), 
    numbertosample(numbertosample), clampval(clampval) {}
  void operator()(std::size_t begin, std::size_t end) {
    random_gen rng;
    for(std::size_t j = begin; j < end; j++) {
      for(int i = 0; i < nx; i++) {
        vec3 col(0,0,0);
        for(int s = 0; s < ns; s++) {
          float u = float(i + rng.unif_rand()) / float(nx);
          float v = float(j + rng.unif_rand()) / float(ny);
          ray r = cam.get_ray(u,v);
          if(numbertosample) {
            if(ambient_light) {
              col += clamp(de_nan(color_amb(r, world, hlist, 0, 
                                            backgroundhigh, backgroundlow, rng)),clampval);
            } else {
              col += clamp(de_nan(color(r, world, hlist, 0, rng)),clampval);
            }
          } else {
            if(ambient_light) {
              col += clamp(de_nan(color_amb_uniform(r, world, 0, backgroundhigh, backgroundlow, rng)),clampval);
            } else {
              col += clamp(de_nan(color_uniform(r, world, 0, rng)),clampval);
            }
          }
        }
        col /= float(ns);
        outputr(i,j) = pow(col[0],1/2.2);
        outputg(i,j) = pow(col[1],1/2.2);
        outputb(i,j) = pow(col[2],1/2.2);
      }
    }
  }

  RcppParallel::RMatrix<double> outputr, outputg, outputb;
  bool ambient_light;
  int nx, ny, ns;
  camera cam;
  vec3 backgroundhigh, backgroundlow;
  hitable *world;
  hitable *hlist;
  int numbertosample;
  float clampval;
};

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
                        int n, float shutteropen, float shutterclose,
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
                        List& group_angle, List& group_order_rotation,
                        random_gen& rng, 
                        LogicalVector tri_normal_bools) {
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
    order_rotation = as<NumericVector>(order_rotation_list(i));
    if(isgrouped(i)) {
      temp_gpivot = as<NumericVector>(group_pivot(i));
      temp_gangle = as<NumericVector>(group_angle(i));
      temp_gorder = as<NumericVector>(group_order_rotation(i));
      temp_gtrans = as<NumericVector>(group_translate(i));
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
        unsigned char *tex_data = stbi_load(filelocation(i), &nx, &ny, &nn, 0);
        tex = new lambertian(new image_texture(tex_data,nx,ny));
      } else if (islight(i)) {
        tex = new diffuse_light(new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2))*lightintensity(i)) );
      } else if (isnoise(i)) {
        tex = new lambertian(new noise_texture(noise(i),vec3(tempvector(0),tempvector(1),tempvector(2)),
                                                    vec3(tempnoisecolor(0),tempnoisecolor(1),tempnoisecolor(2)),
                                                    noisephase(i), noiseintensity(i), rng));
      } else if (ischeckered(i)) {
        tex = new lambertian(new checker_texture(new constant_texture(vec3(tempchecker(0),tempchecker(1),tempchecker(2))),
                                                      new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2))),tempchecker(3)));
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
    }
    //Generate objects
    if (shape(i) == 1) {
      hitable *entry = new sphere(vec3(0,0,0), radius(i), tex);
      entry = rotation_order(entry, temprotvec, order_rotation);
      if(isgrouped(i)) {
        entry = new translate(entry, center - gpivot);
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
      entry = rotation_order(entry, temprotvec, order_rotation);
      if(isgrouped(i)) {
        entry = new translate(entry, center - gpivot);
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
      entry = rotation_order(entry, temprotvec, order_rotation);
      if(isgrouped(i)) {
        entry = new translate(entry, center - gpivot);
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
      entry = rotation_order(entry, temprotvec, order_rotation);
      if(isgrouped(i)) {
        entry = new translate(entry, center - gpivot);
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
      entry = rotation_order(entry, temprotvec, order_rotation);
      if(isgrouped(i)) {
        entry = new translate(entry, center - gpivot);
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
                                                          noisephase(i), noiseintensity(i), rng));
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
      entry = rotation_order(entry, temprotvec, order_rotation);
      if(isgrouped(i)) {
        entry = new translate(entry, center - gpivot);
        entry = rotation_order(entry, temp_gangle, temp_gorder);
        entry = new translate(entry, -center + gpivot );
      }
      entry = new translate(entry, center + gtrans + vel * shutteropen);
      if(isflipped(i)) {
        list[i] = new flip_normals(entry);
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
                         int n, float shutteropen, float shutterclose, 
                         List& angle, int i, List& order_rotation_list,
                         LogicalVector& isgrouped, 
                         List& group_pivot, List& group_translate,
                         List& group_angle, List& group_order_rotation,
                         random_gen& rng) {
  NumericVector tempvector;
  NumericVector temprotvec;
  NumericVector tempvel;
  NumericVector order_rotation;
  NumericVector temp_gpivot;
  NumericVector temp_gtrans;
  NumericVector temp_gorder;
  NumericVector temp_gangle;
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
  prop_len=2;
  center =  vec3(x(i), y(i), z(i));
  if(isgrouped(i)) {
    temp_gpivot = as<NumericVector>(group_pivot(i));
    temp_gangle = as<NumericVector>(group_angle(i));
    temp_gorder = as<NumericVector>(group_order_rotation(i));
    temp_gtrans = as<NumericVector>(group_translate(i));
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
  }
  
  if(shape(i) == 1) {
    hitable *entry = new sphere(vec3(0,0,0), radius(i), 0);
    if(isgrouped(i)) {
      entry = new translate(entry, center - gpivot);
      entry = rotation_order(entry, temp_gangle, temp_gorder);
      entry = new translate(entry, -center + gpivot );
    }
    return(new translate(entry, center + gtrans + vel * shutteropen));
  } else if (shape(i) == 2) {
    hitable *entry = new xy_rect(-tempvector(prop_len+2)/2,tempvector(prop_len+2)/2,
                                 -tempvector(prop_len+4)/2,tempvector(prop_len+4)/2,
                                 0, 0);
    entry = rotation_order(entry, temprotvec, order_rotation);
    if(isgrouped(i)) {
      entry = new translate(entry, center - gpivot);
      entry = rotation_order(entry, temp_gangle, temp_gorder);
      entry = new translate(entry, -center + gpivot );
    }
    return(new translate(entry,center + gtrans + vel * shutteropen));
  } else if (shape(i) == 3) {
    hitable *entry = new xz_rect(-tempvector(prop_len+2)/2,tempvector(prop_len+2)/2,
                                 -tempvector(prop_len+4)/2,tempvector(prop_len+4)/2,
                                 0, 0);
    entry = rotation_order(entry, temprotvec, order_rotation);
    if(isgrouped(i)) {
      entry = new translate(entry, center - gpivot);
      entry = rotation_order(entry, temp_gangle, temp_gorder);
      entry = new translate(entry, -center + gpivot );
    }
    return(new translate(entry,center + gtrans + vel * shutteropen));
  } else if (shape(i) == 4) {
    hitable *entry = new yz_rect(-tempvector(prop_len+2)/2,tempvector(prop_len+2)/2,
                                 -tempvector(prop_len+4)/2,tempvector(prop_len+4)/2,
                                 0, 0);
    entry = rotation_order(entry, temprotvec, order_rotation);
    if(isgrouped(i)) {
      entry = new translate(entry, center - gpivot);
      entry = rotation_order(entry, temp_gangle, temp_gorder);
      entry = new translate(entry, -center + gpivot );
    }
    return(new translate(entry,center + gtrans + vel * shutteropen));
  } else if (shape(i) == 5) {
    hitable *entry = new box(-vec3(tempvector(prop_len+1),tempvector(prop_len+2),tempvector(prop_len+3))/2, 
                             vec3(tempvector(prop_len+1),tempvector(prop_len+2),tempvector(prop_len+3))/2, 
                             0);
    entry = rotation_order(entry, temprotvec, order_rotation);
    if(isgrouped(i)) {
      entry = new translate(entry, center - gpivot);
      entry = rotation_order(entry, temp_gangle, temp_gorder);
      entry = new translate(entry, -center + gpivot );
    }
    return(new translate(entry,center + gtrans + vel * shutteropen));
  } else {
    hitable *entry = new triangle(vec3(tempvector(prop_len+1),tempvector(prop_len+2),tempvector(prop_len+3)),
                                  vec3(tempvector(prop_len+4),tempvector(prop_len+5),tempvector(prop_len+6)),
                                  vec3(tempvector(prop_len+7),tempvector(prop_len+8),tempvector(prop_len+9)),
                                  0);
    entry = rotation_order(entry, temprotvec, order_rotation);
    if(isgrouped(i)) {
      entry = new translate(entry, center - gpivot);
      entry = rotation_order(entry, temp_gangle, temp_gorder);
      entry = new translate(entry, -center + gpivot );
    }
    return(new translate(entry, center + gtrans + vel * shutteropen));
  }
}

// [[Rcpp::export]]
List render_scene_rcpp(int nx, int ny, int ns, float fov, bool ambient_light,
                      NumericVector lookfromvec, NumericVector lookatvec, 
                      float aperture, NumericVector camera_up,
                      IntegerVector type, 
                      NumericVector radius, IntegerVector shape,
                      NumericVector x, NumericVector y, NumericVector z,
                      List properties, List velocity, LogicalVector moving,
                      int n,
                      NumericVector& bghigh, NumericVector& bglow,
                      float shutteropen, float shutterclose,
                      LogicalVector ischeckered, List checkercolors, 
                      NumericVector noise, LogicalVector isnoise,
                      NumericVector& noisephase, NumericVector& noiseintensity, List noisecolorlist,
                      List& angle,
                      LogicalVector& isimage, CharacterVector& filelocation,
                      LogicalVector& islight, NumericVector& lightintensity,
                      LogicalVector& isflipped, float focus_distance,
                      LogicalVector& isvolume, NumericVector& voldensity,
                      bool parallel, LogicalVector& implicit_sample, List& order_rotation_list,
                      float clampval,
                      LogicalVector& isgrouped, List& group_pivot, List& group_translate,
                      List& group_angle, List& group_order_rotation, 
                      LogicalVector tri_normal_bools) {
  NumericMatrix routput(nx,ny);
  NumericMatrix goutput(nx,ny);
  NumericMatrix boutput(nx,ny);
  vec3 lookfrom(lookfromvec[0],lookfromvec[1],lookfromvec[2]);
  vec3 lookat(lookatvec[0],lookatvec[1],lookatvec[2]);
  vec3 backgroundhigh(bghigh[0],bghigh[1],bghigh[2]);
  vec3 backgroundlow(bglow[0],bglow[1],bglow[2]);
  float dist_to_focus = focus_distance;
  random_gen rng;
  camera cam(lookfrom, lookat, vec3(camera_up(0),camera_up(1),camera_up(2)), fov, float(nx)/float(ny), 
             aperture, dist_to_focus,
             shutteropen, shutterclose, rng);
  hitable *world = build_scene(type, radius, shape, x, y, z, 
                                  properties, velocity, moving,
                                  n,shutteropen,shutterclose,
                                  ischeckered, checkercolors, 
                                  noise, isnoise,noisephase,noiseintensity, noisecolorlist,
                                  angle, 
                                  isimage, filelocation,
                                  islight, lightintensity,
                                  isflipped,
                                  isvolume, voldensity, order_rotation_list, 
                                  isgrouped, group_pivot, group_translate,
                                  group_angle, group_order_rotation, rng, tri_normal_bools);
  int numbertosample = 0;
  for(int i = 0; i < implicit_sample.size(); i++) {
    if(implicit_sample(i)) {
      numbertosample++;
    }
  }

  std::vector<hitable* > implicit_sample_vector(numbertosample);
  int counter = 0;
  for(int i = 0; i < n; i++)  {
    if(implicit_sample(i)) {
      implicit_sample_vector[counter] = build_imp_sample(type, radius, shape, x, y, z,
                               properties, velocity,
                               n, shutteropen, shutterclose,
                               angle, i, order_rotation_list,
                               isgrouped, group_pivot, group_translate,
                               group_angle, group_order_rotation,
                               rng);
      counter++;
    }
  }
  hitable_list hlist(&implicit_sample_vector[0],numbertosample);
  
  if(!parallel) {
    for(int j = ny - 1; j >= 0; j--) {
      for(int i = 0; i < nx; i++) {
        vec3 col(0,0,0);
        for(int s = 0; s < ns; s++) {
          float u = float(i + rng.unif_rand()) / float(nx);
          float v = float(j + rng.unif_rand()) / float(ny);
          ray r = cam.get_ray(u,v);
          if(numbertosample) {
            if(ambient_light) {
              col += clamp(de_nan(color_amb(r, world, &hlist, 0, 
                                            backgroundhigh, backgroundlow, rng)),clampval);
            } else {
              col += clamp(de_nan(color(r, world, &hlist, 0, rng)),clampval);
            }
          } else {
            if(ambient_light) {
              col += clamp(de_nan(color_amb_uniform(r, world, 0, 
                                                    backgroundhigh, backgroundlow, rng)),clampval);
            } else {
              col += clamp(de_nan(color_uniform(r, world, 0, rng)),clampval);
            }
          }
        }
        col /= float(ns);
        routput(i,j) = pow(col[0],1/2.2);
        goutput(i,j) = pow(col[1],1/2.2);
        boutput(i,j) = pow(col[2],1/2.2);
      }
    }
  } else {
    Colorworker color_worker(routput, goutput, boutput,
                      ambient_light, nx, ny, ns,
                      cam, backgroundhigh, backgroundlow, world, &hlist, 
                      numbertosample, clampval);
    RcppParallel::parallelFor(0, ny, color_worker);
  }
  return(List::create(_["r"] = routput, _["g"] = goutput, _["b"] = boutput));
}

