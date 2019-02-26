#include <Rcpp.h>
#include "sphere.h"
#include "hitablelist.h"
#include "float.h"
using namespace Rcpp;

vec3 color(const ray& r, hitable *world) {
  hit_record rec;
  if(world->hit(r, 0.0, MAXFLOAT, rec)) {
    return(0.5 * vec3(rec.normal.x()+1,rec.normal.y()+1,rec.normal.z()+1));
  } else {
    vec3 unit_direction = unit_vector(r.direction());
    float t = 0.5 * (unit_direction.y() + 1.0);
    return (1.0 - t) * vec3(1.0, 1.0, 1.0) + t * vec3(0.5,0.7,1.0);
  }
}

// [[Rcpp::export]]
List generate_initial() {
  int nx = 200;
  int ny = 100;
  NumericMatrix routput(nx,ny);
  NumericMatrix goutput(nx,ny);
  NumericMatrix boutput(nx,ny);
  vec3 lower_left_corner(-2.0, -1.0, -1.0);
  vec3 horizontal(4.0, 0.0, 0.0);
  vec3 vertical(0.0, 2.0, 0.0);
  vec3 origin(0.0, 0.0, 0.0);
  hitable *list[2];
  list[0] = new sphere(vec3(0,0,-1),0.5);
  list[1] = new sphere(vec3(0,-100.5,-1),100);
  hitable *world = new hitable_list(list,2);
  for(int j = ny - 1; j >= 0; j--) {
    for(int i = 0; i < nx; i++) {
      float u = float(i) / float(nx);
      float v = float(j) / float(ny);
      ray r(origin, lower_left_corner + u * horizontal + v * vertical);
      
      vec3 p = r.point_at_parameter(2.0);
      vec3 col = color(r, world);
      routput(i,j) = col[0];
      goutput(i,j) = col[1];
      boutput(i,j) = col[2];
    }
  }
  return(List::create(_["r"] = routput, _["g"] = goutput, _["b"] = boutput));
}
