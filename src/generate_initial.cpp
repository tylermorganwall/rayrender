#include <Rcpp.h>
#include "sphere.h"
#include "hitablelist.h"
#include "camera.h"
#include "float.h"
using namespace Rcpp;

vec3 color(const ray& r, hitable *world, int depth) {
  hit_record rec;
  if(world->hit(r, 0.001, MAXFLOAT, rec)) {
    ray scattered;
    vec3 attenuation;
    if(depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
      return(attenuation * color(scattered, world, depth + 1));
    } else {
      return(vec3(0,0,0));
    }
  } else {
    vec3 unit_direction = unit_vector(r.direction());
    float t = 0.5 * (unit_direction.y() + 1.0);
    return (1.0 - t) * vec3(1.0, 1.0, 1.0) + t * vec3(0.5,0.7,1.0);
  }
}

// [[Rcpp::export]]
List generate_initial(int nx = 200, int ny = 100, int ns = 100, float fov = 90.0) {
  
  NumericMatrix routput(nx,ny);
  NumericMatrix goutput(nx,ny);
  NumericMatrix boutput(nx,ny);
  camera cam(vec3(-2,2,1), vec3(0,0,-1), vec3(0,1,0),fov,float(nx)/float(ny));
  hitable *list[5];
  list[0] = new sphere(vec3(0,0,-1),0.5, new lambertian(vec3(0.8,0.3,0.3)));
  list[1] = new sphere(vec3(0,-100.5,-1),100, new lambertian(vec3(0.8,0.8,0.0)));
  list[2] = new sphere(vec3(1,0,-1),0.5, new metal(vec3(0.8,0.6,0.2), 0.8));
  list[3] = new sphere(vec3(-1,0,-1), 0.5, new dielectric(1.5));
  list[4] = new sphere(vec3(-1,0,-1), -0.45, new dielectric(1.5));
  hitable *world = new hitable_list(list,5);
  for(int j = ny - 1; j >= 0; j--) {
    for(int i = 0; i < nx; i++) {
      vec3 col(0,0,0);
      for(int s = 0; s < ns; s++) {
        float u = float(i + drand48()) / float(nx);
        float v = float(j + drand48()) / float(ny);
        ray r = cam.get_ray(u,v);
        col += color(r, world, 0);
      }
      col /= float(ns);
      routput(i,j) = sqrt(col[0]);
      goutput(i,j) = sqrt(col[1]);
      boutput(i,j) = sqrt(col[2]);
    }
  }
  return(List::create(_["r"] = routput, _["g"] = goutput, _["b"] = boutput));
}
