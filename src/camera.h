#ifndef CAMERAH
#define CAMERAH

#include "ray.h"
#include "rng.h"
#include "Rcpp.h"
#include "onbh.h"

class camera {
  public:
    camera(vec3 lookfrom, vec3 lookat, vec3 vup, Float vfov, Float aspect, Float aperture, Float focus_dist,
           Float t0, Float t1);
    ray get_ray(Float s, Float t, vec3 u3, Float u1);
    
    vec3 origin;
    vec3 lower_left_corner;
    vec3 horizontal;
    vec3 vertical;
    vec3 u, v, w;
    Float time0, time1;
    Float lens_radius;
};

class ortho_camera {
public:
  ortho_camera(vec3 lookfrom, vec3 lookat, vec3 vup, 
               Float cam_width, Float cam_height, 
               Float t0, Float t1);
  ray get_ray(Float s, Float t, Float u);
  
  vec3 origin;
  vec3 lower_left_corner;
  vec3 horizontal;
  vec3 vertical;
  vec3 u, v, w;
  Float time0, time1;
};


class environment_camera {
  public:
    environment_camera(vec3 lookfrom, vec3 lookat, vec3 vup, 
                       Float t0, Float t1);
    ray get_ray(Float s, Float t, Float u1);
    
    vec3 origin;
    vec3 u, v, w;
    Float nx, ny;
    Float time0, time1;
    onb uvw;
};

  
#endif
