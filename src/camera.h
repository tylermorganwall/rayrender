#ifndef CAMERAH
#define CAMERAH

#include "ray.h"
#include "rng.h"

class camera {
  public:
    camera(vec3 lookfrom, vec3 lookat, vec3 vup, Float vfov, Float aspect, Float aperture, Float focus_dist,
           Float t0, Float t1, random_gen& rng) {
      rng2 = rng;
      time0 = t0;
      time1 = t1;
      lens_radius = aperture / 2;
      Float theta = vfov * M_PI/180;
      Float half_height = tan(theta/2);
      Float half_width = aspect * half_height;
      origin = lookfrom;
      w = unit_vector(lookfrom - lookat);
      u = unit_vector(cross(vup, w));
      v = cross(w, u);
      lower_left_corner = origin - half_width * focus_dist *  u - half_height * focus_dist * v - focus_dist * w;
      horizontal = 2 * half_width * focus_dist * u;
      vertical = 2 * half_height * focus_dist * v;
    }
    ray get_ray(Float s, Float t) {
      vec3 rd = lens_radius * rng2.random_in_unit_disk();
      vec3 offset = u * rd.x() + v * rd.y();
      Float time = time0 + rng2.unif_rand() * (time1 - time0);
      return(ray(origin + offset, lower_left_corner + s * horizontal + t * vertical - origin - offset, time)); 
    }
    
    vec3 origin;
    vec3 lower_left_corner;
    vec3 horizontal;
    vec3 vertical;
    vec3 u, v, w;
    Float time0, time1;
    Float lens_radius;
    random_gen rng2;
};
  
#endif
