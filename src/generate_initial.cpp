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
#include <RcppParallel.h>
using namespace Rcpp;

vec3 color(const ray& r, hitable *world, int depth) {
  hit_record rec;
  if(world->hit(r, 0.001, FLT_MAX, rec)) {
    ray scattered;
    vec3 albedo;
    vec3 emitted = rec.mat_ptr->emitted(rec.u,rec.v,rec.p);
    float pdf;
    if(depth < 50 && rec.mat_ptr->scatter(r, rec, albedo, scattered, pdf)) {
      return(emitted + albedo * color(scattered, world, depth + 1));
    } else {
      return(emitted);
    }
  } else {
    return(vec3(0,0,0));
  }
}

vec3 color_amb(const ray& r, hitable *world, int depth, 
           const vec3& backgroundhigh, const vec3& backgroundlow) {
  hit_record rec;
  if(world->hit(r, 0.001, FLT_MAX, rec)) {
    ray scattered;
    vec3 albedo;
    vec3 emitted = rec.mat_ptr->emitted(rec.u,rec.v,rec.p);
    float pdf;
    if(depth < 50 && rec.mat_ptr->scatter(r, rec, albedo, scattered, pdf)) {
      return(emitted + albedo * color_amb(scattered, world, depth + 1, backgroundhigh, backgroundlow));
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
              bool ambient_light, int nx, int ny, int ns, camera cam, vec3 backgroundhigh, vec3 backgroundlow, hitable *world)
  : outputr(outputr), outputg(outputg), outputb(outputb), ambient_light(ambient_light),
    nx(nx), ny(ny), ns(ns), cam(cam),
    backgroundhigh(backgroundhigh), backgroundlow(backgroundlow), world(world) {}
  void operator()(std::size_t begin, std::size_t end) {
    srand(end);
    for(int j = begin; j < end; j++) {
      for(int i = 0; i < nx; i++) {
        vec3 col(0,0,0);
        for(int s = 0; s < ns; s++) {
          float u = float(i + drand48()) / float(nx);
          float v = float(j + drand48()) / float(ny);
          ray r = cam.get_ray(u,v);
          if(ambient_light) {
            col += color_amb(r, world, 0, backgroundhigh, backgroundlow);
          } else {
            col += color(r, world, 0);
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
};


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
                        LogicalVector& isvolume, List& fogcolor, NumericVector& voldensity) {
  hitable **list = new hitable*[n+1];
  NumericVector tempvector;
  NumericVector tempchecker;
  NumericVector tempvel;
  NumericVector tempfog;
  NumericVector tempnoisecolor;
  NumericVector temprotvec;
  int prop_len;

  List templist;
  vec3 center(x(0), y(0), z(0));
  vec3 vel(x(0), y(0), z(0));
  for(int i = 0; i < n; i++) {
    tempvector = as<NumericVector>(properties(i));
    tempchecker = as<NumericVector>(checkercolors(i));
    tempvel = as<NumericVector>(velocity(i));
    tempfog = as<NumericVector>(fogcolor(i));
    tempnoisecolor = as<NumericVector>(noisecolorlist(i));
    temprotvec = as<NumericVector>(angle(i));
    prop_len=2;
    
    center =  vec3(x(i), y(i), z(i));
    vel = vec3(tempvel(0),tempvel(1),tempvel(2));
    if (shape(i) == 1) {
      material *sphere_tex;
      if(type(i) == 1) {
        if(isimage(i)) {
          // int nx, ny, nn;
          // unsigned char *tex_data = stbi_load(filelocation(i), &nx, &ny, &nn, 0);
          // sphere_tex = new lambertian(new image_texture(tex_data,nx,ny));
        } else if (islight(i)) {
          sphere_tex = new diffuse_light(new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2))*lightintensity(i)) );
        } else if (isnoise(i)) {
          // sphere_tex = new lambertian(new noise_texture(noise(i),vec3(tempvector(0),tempvector(1),tempvector(2)), 
          //                                               vec3(tempnoisecolor(0),tempnoisecolor(1),tempnoisecolor(2)),
          //                                               noisephase(i), noiseintensity(i)));
        } else if (ischeckered(i)) {
          // sphere_tex = new lambertian(new checker_texture(new constant_texture(vec3(tempchecker(0),tempchecker(1),tempchecker(2))),
          //                                               new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2))),tempchecker(3))); 
        } else {
          sphere_tex = new lambertian(new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2))) );
        }
      } 
      // else if (type(i) == 2) {
      //   sphere_tex = new metal(vec3(tempvector(0),tempvector(1),tempvector(2)),tempvector(3));
      // } else {
      //   sphere_tex = new dielectric(vec3(tempvector(0),tempvector(1),tempvector(2)),tempvector(3));
      // }
      hitable *entry = new sphere(vec3(0,0,0), radius(i), sphere_tex);
      if(temprotvec(0) != 0) {
        entry = new rotate_x(entry,temprotvec(0));
      }
      if(temprotvec(1) != 0) {
        entry = new rotate_y(entry,temprotvec(1));
      }
      if(temprotvec(2) != 0) {
        entry = new rotate_z(entry,temprotvec(2));
      }
      if(!moving(i)) {
        entry = new translate(entry, center + vel * shutteropen);
      } else {
        entry = new moving_sphere(center + vel * shutteropen, 
                                  center + vel * shutterclose, 
                                  shutteropen, shutterclose, radius(i), sphere_tex);
      }
      if(isvolume(i)) {
        list[i] = new constant_medium(entry,voldensity(i), new constant_texture(vec3(tempfog(0),tempfog(1),tempfog(2))));
      } else {
        list[i] = entry;
      }
    } else if (shape(i)  == 2) {
      material *rect_tex;
      if(type(i) == 1) {
        if(isimage(i)) {
          // int nx, ny, nn;
          // unsigned char *tex_data = stbi_load(filelocation(i), &nx, &ny, &nn, 0);
          // rect_tex = new lambertian(new image_texture(tex_data,nx,ny));
        } else if (islight(i)) {
          rect_tex = new diffuse_light(new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2))*lightintensity(i)) );
        } else if (isnoise(i)) {
          // rect_tex = new lambertian(new noise_texture(noise(i),vec3(tempvector(0),tempvector(1),tempvector(2)), 
          //                                               vec3(tempnoisecolor(0),tempnoisecolor(1),tempnoisecolor(2)),
          //                                               noisephase(i), noiseintensity(i)));
        } else if (ischeckered(i)) {
          // rect_tex = new lambertian(new checker_texture(new constant_texture(vec3(tempchecker(0),tempchecker(1),tempchecker(2))),
          //                                                 new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2))),tempchecker(3))); 
        } else {
          rect_tex = new lambertian(new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2))) );
        }
      } 
      // else if (type(i) == 2) {
      //   rect_tex = new metal(vec3(tempvector(0),tempvector(1),tempvector(2)),tempvector(3));
      //   prop_len = 3;
      // } else {
      //   rect_tex = new dielectric(vec3(tempvector(0),tempvector(1),tempvector(2)),tempvector(3));
      //   prop_len = 3;
      // }
      hitable *entry = new xy_rect(-tempvector(prop_len+2)/2,tempvector(prop_len+2)/2,
                                                              -tempvector(prop_len+4)/2,tempvector(prop_len+4)/2,
                                                              0, rect_tex);
      if(temprotvec(0) != 0) {
        entry = new rotate_x(entry,temprotvec(0));
      }
      if(temprotvec(1) != 0) {
        entry = new rotate_y(entry,temprotvec(1));
      }
      if(temprotvec(2) != 0) {
        entry = new rotate_z(entry,temprotvec(2));
      }
      entry = new translate(entry,vec3(tempvector(prop_len+1),tempvector(prop_len+3),tempvector(prop_len+5)) + vel * shutteropen);
      if(isflipped(i)) {
        list[i] = new flip_normals(entry);
      } else {
        list[i] = entry;
      }
    } else if (shape(i)  == 3) {
      material *rect_tex;
      if(type(i) == 1) {
        if(isimage(i)) {
          // int nx, ny, nn;
          // unsigned char *tex_data = stbi_load(filelocation(i), &nx, &ny, &nn, 0);
          // rect_tex = new lambertian(new image_texture(tex_data,nx,ny));
        } else if (islight(i)) {
          rect_tex = new diffuse_light(new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2))*lightintensity(i)) );
        } else if (isnoise(i)) {
          // rect_tex = new lambertian(new noise_texture(noise(i),vec3(tempvector(0),tempvector(1),tempvector(2)), 
          //                                             vec3(tempnoisecolor(0),tempnoisecolor(1),tempnoisecolor(2)),
          //                                             noisephase(i), noiseintensity(i)));
        } else if (ischeckered(i)) {
          // rect_tex = new lambertian(new checker_texture(new constant_texture(vec3(tempchecker(0),tempchecker(1),tempchecker(2))),
          //                                               new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2))),tempchecker(3))); 
        } else {
          rect_tex = new lambertian(new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2))) );
        }
      } 
      // else if (type(i) == 2) {
      //   rect_tex = new metal(vec3(tempvector(0),tempvector(1),tempvector(2)),tempvector(3));
      //   prop_len = 3;
      // } else {
      //   rect_tex = new dielectric(vec3(tempvector(0),tempvector(1),tempvector(2)),tempvector(3));
      //   prop_len = 3;
      // }
      hitable *entry = new xz_rect(-tempvector(prop_len+2)/2,tempvector(prop_len+2)/2,
                                   -tempvector(prop_len+4)/2,tempvector(prop_len+4)/2,
                                   0, rect_tex);
      if(temprotvec(0) != 0) {
        entry = new rotate_x(entry,temprotvec(0));
      }
      if(temprotvec(1) != 0) {
        entry = new rotate_y(entry,temprotvec(1));
      }
      if(temprotvec(2) != 0) {
        entry = new rotate_z(entry,temprotvec(2));
      }
      entry = new translate(entry,vec3(tempvector(prop_len+1),tempvector(prop_len+5), tempvector(prop_len+3)) + vel * shutteropen);
      if(isflipped(i)) {
        list[i] = new flip_normals(entry);
      } else {
        list[i] = entry;
      }
    } else if (shape(i)  == 4) {
      material *rect_tex;
      if(type(i) == 1) {
        if(isimage(i)) {
          // int nx, ny, nn;
          // unsigned char *tex_data = stbi_load(filelocation(i), &nx, &ny, &nn, 0);
          // rect_tex = new lambertian(new image_texture(tex_data,nx,ny));
        } else if (islight(i)) {
          rect_tex = new diffuse_light(new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2))*lightintensity(i)) );
        } else if (isnoise(i)) {
          // rect_tex = new lambertian(new noise_texture(noise(i),vec3(tempvector(0),tempvector(1),tempvector(2)), 
          //                                             vec3(tempnoisecolor(0),tempnoisecolor(1),tempnoisecolor(2)),
          //                                             noisephase(i), noiseintensity(i)));
        } else if (ischeckered(i)) {
          // rect_tex = new lambertian(new checker_texture(new constant_texture(vec3(tempchecker(0),tempchecker(1),tempchecker(2))),
          //                                               new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2))),tempchecker(3))); 
        } else {
          rect_tex = new lambertian(new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2))) );
        }
      } 
      // else if (type(i) == 2) {
      //   rect_tex = new metal(vec3(tempvector(0),tempvector(1),tempvector(2)),tempvector(3));
      //   prop_len = 3;
      // } else {
      //   rect_tex = new dielectric(vec3(tempvector(0),tempvector(1),tempvector(2)),tempvector(3));
      //   prop_len = 3;
      // }
      hitable *entry = new yz_rect(-tempvector(prop_len+2)/2,tempvector(prop_len+2)/2,
                                   -tempvector(prop_len+4)/2,tempvector(prop_len+4)/2,
                                   0, rect_tex);
      if(temprotvec(0) != 0) {
        entry = new rotate_x(entry,temprotvec(0));
      }
      if(temprotvec(1) != 0) {
        entry = new rotate_y(entry,temprotvec(1));
      }
      if(temprotvec(2) != 0) {
        entry = new rotate_z(entry,temprotvec(2));
      }
      entry = new translate(entry,vec3(tempvector(prop_len+5),tempvector(prop_len+1),tempvector(prop_len+3)) + vel * shutteropen);
      if(isflipped(i)) {
        list[i] = new flip_normals(entry);
      } else {
        list[i] = entry;
      }
    } else if (shape(i)  == 5) {
      material *rect_tex;
      if(type(i) == 1) {
        if(isimage(i)) {
          // int nx, ny, nn;
          // unsigned char *tex_data = stbi_load(filelocation(i), &nx, &ny, &nn, 0);
          // rect_tex = new lambertian(new image_texture(tex_data,nx,ny));
        } else if (islight(i)) {
          rect_tex = new diffuse_light(new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2))*lightintensity(i)) );
        } else if (isnoise(i)) {
          // rect_tex = new lambertian(new noise_texture(noise(i),
          //                                             vec3(tempvector(0),tempvector(1),tempvector(2)), 
          //                                             vec3(tempnoisecolor(0),tempnoisecolor(1),tempnoisecolor(2)),
          //                                             noisephase(i), noiseintensity(i)));
        } else if (ischeckered(i)) {
          // rect_tex = new lambertian(new checker_texture(new constant_texture(vec3(tempchecker(0),tempchecker(1),tempchecker(2))),
          //                                 new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2))),tempchecker(3))); 
        } else {
          rect_tex = new lambertian(new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2))) );
        }
      }
      // else if (type(i) == 2) {
      //   rect_tex = new metal(vec3(tempvector(0),tempvector(1),tempvector(2)),tempvector(3));
      //   prop_len = 3;
      // } else {
      //   rect_tex = new dielectric(vec3(tempvector(0),tempvector(1),tempvector(2)),tempvector(3));
      //   prop_len = 3;
      // }
      hitable *entry = new box(-vec3(tempvector(prop_len+1),tempvector(prop_len+2),tempvector(prop_len+3))/2, 
                                vec3(tempvector(prop_len+1),tempvector(prop_len+2),tempvector(prop_len+3))/2, 
                                rect_tex);
      if(temprotvec(0) != 0) {
        entry = new rotate_x(entry,temprotvec(0));
      }
      if(temprotvec(1) != 0) {
        entry = new rotate_y(entry,temprotvec(1));
      }
      if(temprotvec(2) != 0) {
        entry = new rotate_z(entry,temprotvec(2));
      }
      entry = new translate(entry,center + vel * shutteropen);
      if(isvolume(i)) {
        if(!isnoise(i)) {
          list[i] = new constant_medium(entry,voldensity(i), new constant_texture(vec3(tempfog(0),tempfog(1),tempfog(2))));
        } else {
          list[i] = new constant_medium(entry,voldensity(i), 
                                        new noise_texture(noise(i),
                                                          vec3(tempfog(0),tempfog(1),tempfog(2)),
                                                          vec3(tempnoisecolor(0),tempnoisecolor(1),tempnoisecolor(2)),
                                                          noisephase(i), noiseintensity(i)));
        }
      } else {
        list[i] = entry;
      }
    }
  }
  return(new bvh_node(list, n, shutteropen, shutterclose));
}

// [[Rcpp::export]]
List generate_initial(int nx, int ny, int ns, float fov, bool ambient_light,
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
                      LogicalVector& isvolume, List& fogcolor, NumericVector& voldensity,
                      bool parallel) {
  NumericMatrix routput(nx,ny);
  NumericMatrix goutput(nx,ny);
  NumericMatrix boutput(nx,ny);
  vec3 lookfrom(lookfromvec[0],lookfromvec[1],lookfromvec[2]);
  vec3 lookat(lookatvec[0],lookatvec[1],lookatvec[2]);
  vec3 backgroundhigh(bghigh[0],bghigh[1],bghigh[2]);
  vec3 backgroundlow(bglow[0],bglow[1],bglow[2]);
  float dist_to_focus = focus_distance;
  camera cam(lookfrom, lookat, vec3(camera_up(0),camera_up(1),camera_up(2)), fov, float(nx)/float(ny), 
             aperture, dist_to_focus,
             shutteropen, shutterclose);
  hitable *world = build_scene(type, radius, shape, x, y, z, 
                                  properties, velocity, moving,
                                  n,shutteropen,shutterclose,
                                  ischeckered, checkercolors, 
                                  noise, isnoise,noisephase,noiseintensity, noisecolorlist,
                                  angle, 
                                  isimage, filelocation,
                                  islight, lightintensity,
                                  isflipped,
                                  isvolume, fogcolor, voldensity);
  if(!parallel) {
    for(int j = ny - 1; j >= 0; j--) {
      for(int i = 0; i < nx; i++) {
        vec3 col(0,0,0);
        for(int s = 0; s < ns; s++) {
          float u = float(i + drand48()) / float(nx);
          float v = float(j + drand48()) / float(ny);
          ray r = cam.get_ray(u,v);
          if(ambient_light) {
            col += color_amb(r, world, 0, backgroundhigh, backgroundlow);
          } else {
            col += color(r, world, 0);
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
                      cam, backgroundhigh, backgroundlow, world);
    RcppParallel::parallelFor(0, ny, color_worker);
  }
  return(List::create(_["r"] = routput, _["g"] = goutput, _["b"] = boutput));
}

