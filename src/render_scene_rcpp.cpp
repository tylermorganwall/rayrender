#define STB_IMAGE_IMPLEMENTATION 

#include "vec3.h"
#include "camera.h"
#include "float.h"
#include "buildscene.h"
#include "RProgress.h"
#include "rng.h"
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
#include "RcppThread.h"

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
  if(world->hit(r, 0.001, FLT_MAX, hrec, rng)) { //generated hit record, world space
    scatter_record srec;
    vec3 emitted = hrec.mat_ptr->emitted(r, hrec, hrec.u, hrec.v, hrec.p);
    float pdf_val;
    if(depth < 50 && hrec.mat_ptr->scatter(r, hrec, srec, rng)) { //generates scatter record, world space
      if(srec.is_specular) { //returns specular ray
        return(srec.attenuation * 
               color(srec.specular_ray, world, 
                     hlist, depth + 1, rng));
      }
      hitable_pdf p_imp(hlist, hrec.p); //creates pdf of all objects to be sampled
      mixture_pdf p(&p_imp, srec.pdf_ptr); //creates mixture pdf of surface intersected at hrec.p and all sampled objects
      
      //Generates a scatter direction (with origin hrec.p) from the mixture 
      //and saves surface normal from light to use in pdf_value calculation
      //(along with the scatter direction)
      //Translates the world space point into object space point, generates ray assuring intersection, and then translates 
      //ray back into world space
      ray scattered = ray(hrec.p, p.generate(rng), r.time()); //scatters a ray from hit point to direction
      pdf_val = p.value(scattered.direction(), rng); //generates a pdf value based the intersection point and the mixture pdf
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
                      LogicalVector& tri_normal_bools, LogicalVector& is_tri_color, List& tri_color_vert,
                      CharacterVector& fileinfo,
                      bool progress_bar, int numbercores) {
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
                                  group_angle, group_order_rotation, 
                                  tri_normal_bools, is_tri_color, tri_color_vert, fileinfo, rng);
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
  
  RProgress::RProgress pb("Raytracing [:bar] ETA: :eta");
  
  if(progress_bar) {
    pb.set_total(ny);
  }
  if(!parallel) {
    for(int j = ny - 1; j >= 0; j--) {
      Rcpp::checkUserInterrupt();
      if(progress_bar) {
        pb.tick();
      }
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
    RcppThread::ThreadPool pool(numbercores);
    auto worker = [&routput, &goutput, &boutput,
                   ambient_light, nx, ny, ns,
                   &cam, backgroundhigh, backgroundlow, &world, &hlist,
                   numbertosample, clampval, progress_bar, numbercores] (int j) {
    // auto worker = [nx, ns] (int j) {
      if(progress_bar && j % numbercores == 0) {
        RcppThread::Rcout << "Progress (" << numbercores << " cores): ";
        RcppThread::Rcout << (int)((1-(double)j/double(ny)) * 100) << "%\r";
      }
      random_gen rng;
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
              col += clamp(de_nan(color_amb_uniform(r, world, 0, backgroundhigh, backgroundlow, rng)),clampval);
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
    };
    for(int j = ny - 1; j >= 0; j--) {
      pool.push(worker,j);
    }
    pool.join();
  }
  return(List::create(_["r"] = routput, _["g"] = goutput, _["b"] = boutput));
}

