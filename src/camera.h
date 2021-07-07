#ifndef CAMERAH
#define CAMERAH

#include "ray.h"
#include "rng.h"
#include "Rcpp.h"
#include "onbh.h"

class camera {
  public:
    camera(vec3f lookfrom, vec3f lookat, vec3f vup, Float vfov, Float aspect, Float aperture, Float focus_dist,
           Float t0, Float t1);
    ray get_ray(Float s, Float t, vec3f u3, Float u1);
    
    vec3f origin;
    vec3f lower_left_corner;
    vec3f horizontal;
    vec3f vertical;
    vec3f u, v, w;
    Float time0, time1;
    Float lens_radius;
};

class ortho_camera {
public:
  ortho_camera(vec3f lookfrom, vec3f lookat, vec3f vup, 
               Float cam_width, Float cam_height, 
               Float t0, Float t1);
  ray get_ray(Float s, Float t, Float u);
  
  vec3f origin;
  vec3f lower_left_corner;
  vec3f horizontal;
  vec3f vertical;
  vec3f u, v, w;
  Float time0, time1;
};


class environment_camera {
  public:
    environment_camera(vec3f lookfrom, vec3f lookat, vec3f vup, 
                       Float t0, Float t1);
    ray get_ray(Float s, Float t, Float u1);
    
    vec3f origin;
    vec3f u, v, w;
    Float nx, ny;
    Float time0, time1;
    onb uvw;
};

  
#endif
