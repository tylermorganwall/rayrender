#include <Rcpp.h>
#include "sphere.h"
#include "hitablelist.h"
#include "camera.h"
#include "float.h"
#include "bvh_node.h"
#include "perlin.h"
using namespace Rcpp;

vec3 color(const ray& r, hitable *world, int depth, 
           const vec3& backgroundhigh, const vec3& backgroundlow) {
  hit_record rec;
  if(world->hit(r, 0.001, MAXFLOAT, rec)) {
    ray scattered;
    vec3 attenuation;
    if(depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
      return(attenuation * color(scattered, world, depth + 1, backgroundhigh, backgroundlow));
    } else {
      return(vec3(0,0,0));
    }
  } else {
    vec3 unit_direction = unit_vector(r.direction());
    float t = 0.5 * (unit_direction.y() + 1.0);
    return (1.0 - t) * backgroundlow + t * backgroundhigh;
  }
}

hitable *specific_scene(IntegerVector& type, 
                        NumericVector& radius,
                        NumericVector& x, NumericVector& y, NumericVector& z,
                        List& properties, List& velocity, LogicalVector& moving,
                        int n, float shutteropen, float shutterclose,
                        LogicalVector& ischeckered, List& checkercolors, 
                        NumericVector& noise, LogicalVector& isnoise) {
  hitable **list = new hitable*[n+1];
  NumericVector tempvector;
  NumericVector tempchecker;
  NumericVector tempvel;
  
  List templist;
  vec3 center(x(0), y(0), z(0));
  vec3 vel(x(0), y(0), z(0));
  for(int i = 0; i < n; i++) {
    tempvector = as<NumericVector>(properties(i));
    tempchecker = as<NumericVector>(checkercolors(i));
    tempvel = as<NumericVector>(velocity(i));
    center =  vec3(x(i), y(i), z(i));
    vel = vec3(tempvel(0),tempvel(1),tempvel(2));
    if (type(i) == 1) {
      if(ischeckered(i)) {
        texture *checker = new checker_texture(new constant_texture(vec3(tempchecker(0),tempchecker(1),tempchecker(2))),
                                               new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2)))); 
        if(!moving(i)) {
          list[i] = new sphere(center + vel * shutteropen, radius(i), 
                               new lambertian(checker));
        } else {
          list[i] = new moving_sphere(center + vel * shutteropen, center + vel*shutterclose, shutteropen, shutterclose, radius(i),
                                      new lambertian(checker));
        }
      } else if(isnoise(i)) {
        texture *perlinn = new noise_texture(noise(i),vec3(tempvector(0),tempvector(1),tempvector(2))); 
        if(!moving(i)) {
          list[i] = new sphere(center + vel * shutteropen, radius(i), 
                               new lambertian(perlinn));
        } else {
          list[i] = new moving_sphere(center + vel * shutteropen, center + vel*shutterclose, shutteropen, shutterclose, radius(i),
                                      new lambertian(perlinn));
        }
      } else {
        if(!moving(i)) {
          list[i] = new sphere(center + vel * shutteropen, radius(i), 
                               new lambertian(new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2)))));
        } else {
          list[i] = new moving_sphere(center + vel * shutteropen, center + vel*shutterclose, shutteropen, shutterclose, radius(i),
                                      new lambertian(new constant_texture(vec3(tempvector(0),tempvector(1),tempvector(2)))));
        }
      }
    } else if (type(i)  == 2) {
      if(!moving(i)) {
        list[i] = new sphere(center + vel * shutteropen, radius(i), 
                               new metal(vec3(tempvector(0),tempvector(1),tempvector(2)),tempvector(3)));
      } else {
        list[i] = new moving_sphere(center + vel * shutteropen, center + vel*shutterclose, shutteropen, shutterclose, radius(i),
                             new metal(vec3(tempvector(0),tempvector(1),tempvector(2)),tempvector(3)));
      }
    } else if (type(i)  == 3) {
      if(!moving(i)) {
        list[i] = new sphere(center + vel * shutteropen, radius(i), new dielectric(tempvector(0))); 
      } else {
        list[i] = new moving_sphere(center + vel * shutteropen, center + vel*shutterclose, shutteropen, shutterclose, radius(i),
                                    new dielectric(tempvector(0)));
      }
    }
  }
  return(new bvh_node(list, n, shutteropen, shutterclose));
}

// [[Rcpp::export]]
List generate_initial(int nx, int ny, int ns, float fov, 
                      NumericVector lookfromvec, NumericVector lookatvec, float aperture, 
                      IntegerVector type, 
                      NumericVector radius,
                      NumericVector x, NumericVector y, NumericVector z,
                      List properties, List velocity, LogicalVector moving,
                      int n,
                      NumericVector& bghigh, NumericVector& bglow,
                      float shutteropen, float shutterclose,
                      LogicalVector ischeckered, List checkercolors, 
                      NumericVector noise, LogicalVector isnoise) {
  NumericMatrix routput(nx,ny);
  NumericMatrix goutput(nx,ny);
  NumericMatrix boutput(nx,ny);
  vec3 lookfrom(lookfromvec[0],lookfromvec[1],lookfromvec[2]);
  vec3 lookat(lookatvec[0],lookatvec[1],lookatvec[2]);
  vec3 backgroundhigh(bghigh[0],bghigh[1],bghigh[2]);
  vec3 backgroundlow(bglow[0],bglow[1],bglow[2]);
  float dist_to_focus = (lookfrom-lookat).length();
  camera cam(lookfrom, lookat, vec3(0,1,0), fov, float(nx)/float(ny), 
             aperture, dist_to_focus,
             shutteropen, shutterclose);
  hitable *world = specific_scene(type, radius, x, y, z, 
                                  properties, velocity, moving,
                                  n,shutteropen,shutterclose,
                                  ischeckered, checkercolors, noise, isnoise);
  for(int j = ny - 1; j >= 0; j--) {
    for(int i = 0; i < nx; i++) {
      vec3 col(0,0,0);
      for(int s = 0; s < ns; s++) {
        float u = float(i + drand48()) / float(nx);
        float v = float(j + drand48()) / float(ny);
        ray r = cam.get_ray(u,v);
        col += color(r, world, 0, backgroundhigh, backgroundlow);
      }
      col /= float(ns);
      routput(i,j) = sqrt(col[0]);
      goutput(i,j) = sqrt(col[1]);
      boutput(i,j) = sqrt(col[2]);
    }
  }
  return(List::create(_["r"] = routput, _["g"] = goutput, _["b"] = boutput));
}

