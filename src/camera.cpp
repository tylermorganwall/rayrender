
#include "camera.h"

camera::camera(vec3 lookfrom, vec3 lookat, vec3 vup, Float vfov, Float aspect, Float aperture, Float focus_dist,
       Float t0, Float t1) {
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
  horizontal = 2.0f * half_width * focus_dist * u;
  vertical = 2.0f * half_height * focus_dist * v;
}

ray camera::get_ray(Float s, Float t, vec3 u3, Float u1) {
  vec3 rd = lens_radius * u3;
  vec3 offset = u * rd.x() + v * rd.y();
  Float time = time0 + u1 * (time1 - time0);
  return(ray(origin + offset, lower_left_corner + s * horizontal + t * vertical - origin - offset, time)); 
}

ortho_camera::ortho_camera(vec3 lookfrom, vec3 lookat, vec3 vup, 
             Float cam_width, Float cam_height, 
             Float t0, Float t1) {
  time0 = t0;
  time1 = t1;
  origin = lookfrom;
  w = unit_vector(lookfrom - lookat);
  u = unit_vector(cross(vup, w));
  v = cross(w, u);
  lower_left_corner = origin - cam_width/2 *  u - cam_height/2 * v;
  horizontal = cam_width * u;
  vertical = cam_height * v;
}

ray ortho_camera::get_ray(Float s, Float t, Float u) {
  Float time = time0 + u * (time1 - time0);
  return(ray(lower_left_corner + s * horizontal + t * vertical, -w, time)); 
}

environment_camera::environment_camera(vec3 lookfrom, vec3 lookat, vec3 vup, 
                   Float t0, Float t1) {
  time0 = t0;
  time1 = t1;
  origin = lookfrom;
  w = unit_vector(lookfrom - lookat);
  v = unit_vector(-cross(vup, w));
  u = cross(w, v);
  uvw = onb(w,v,u);
}

ray environment_camera::get_ray(Float s, Float t, Float u1) {
  Float time = time0 + u1 * (time1 - time0);
  Float theta = M_PI * t;
  Float phi = 2 * M_PI * s;
  vec3 dir(std::sin(theta) * std::cos(phi), 
           std::sin(theta) * std::sin(phi),
           std::cos(theta));
  dir = uvw.local_to_world(dir);
  return(ray(origin, dir, time)); 
}
