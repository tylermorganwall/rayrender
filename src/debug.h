#ifndef DEBUGH
#define DEBUGH

#include "vec3.h"
#include "vec2.h"
#include "mathinline.h"
#include "camera.h"
#include "hitable.h"
#include "material.h"
#include "RcppThread.h"
#include "rng.h"


#ifdef DEBUGBVH
inline Float debug_bvh(const ray& r, hitable *world, random_gen &rng) {
  hit_record hrec;
  hrec.bvh_nodes = 0.0;
  if(world->hit(r, 0.001, FLT_MAX, hrec, rng)) {
    return(hrec.bvh_nodes);
  } else {
    return(0);
  }
}
#endif


inline Float calculate_depth(const ray& r, hitable *world, random_gen &rng) {
  hit_record hrec;
  if(world->hit(r, 0.001, FLT_MAX, hrec, rng)) {
    return((r.origin()-hrec.p).length());
  } else {
    return(INFINITY);
  }
}

inline vec3f calculate_normals(const ray& r, hitable *world, size_t max_depth, random_gen &rng) {
  point3f final_color(0,0,0);
  ray r1 = r;
  ray r2 = r;
  bool diffuse_bounce = false;
  
  for(size_t i = 0; i < max_depth; i++) {
    bool is_invisible = false;
    hit_record hrec;
    if(world->hit(r2, 0.001, FLT_MAX, hrec, rng)) { //generated hit record, world space
      scatter_record srec;
      //Some lights can be invisible until after diffuse bounce
      //If so, generate new ray with intersection point and continue ray
      if(is_invisible && !diffuse_bounce) {
        r2.A = hrec.p;
        continue;
      }
      if(hrec.mat_ptr->scatter(r2, hrec, srec, rng)) { //generates scatter record, world space
        if(srec.is_specular) { //returns specular ray
          if(i == max_depth-1) {
            hrec.normal.make_unit_vector();
            return((vec3f(1,1,1) + vec3f(hrec.normal.x(),hrec.normal.y(),hrec.normal.z()))/2);
          }
          r2 = srec.specular_ray;
          continue;
        }
        // hitable_pdf p_imp(hlist, hrec.p); //creates pdf of all objects to be sampled
        // mixture_pdf p(&p_imp, srec.pdf_ptr); //creates mixture pdf of surface intersected at hrec.p and all sampled objects/lights
        
        point3f offset_p = offset_ray(hrec.p-r2.A, hrec.normal) + r2.A;
        
        r1 = r2;
        vec3f dir;
        dir = srec.pdf_ptr->generate(rng, diffuse_bounce, r2.time());
        if((dir.x() == 0 && dir.y() == 0 && dir.z() == 0)) {
          hrec.normal.make_unit_vector();
          return((vec3f(1,1,1) + vec3f(hrec.normal.x(),hrec.normal.y(),hrec.normal.z()))/2);
        }
        r2 = ray(OffsetRayOrigin(offset_p, hrec.pError, hrec.normal, dir), dir, r2.pri_stack, r2.time());
        
        if(i == max_depth-1) {
          hrec.normal.make_unit_vector();
          return((vec3f(1,1,1) + vec3f(hrec.normal.x(),hrec.normal.y(),hrec.normal.z()))/2);
        }
      } else {
        hrec.normal.make_unit_vector();
        return((vec3f(1,1,1) + vec3f(hrec.normal.x(),hrec.normal.y(),hrec.normal.z()))/2);
      }
    } else {
      return(vec3f(0,0,0));
    }
  }
  return(vec3f(0,0,0));
}

inline vec3f calculate_uv(const ray& r, hitable *world, random_gen &rng) {
  hit_record hrec;
  if(world->hit(r, 0.001, FLT_MAX, hrec, rng)) {
    return(vec3f(hrec.u,hrec.v,1-hrec.u-hrec.v));
  } else {
    return(vec3f(0,0,0));
  }
}

inline vec3f calculate_dpduv(const ray& r, hitable *world, random_gen &rng, bool u) {
  hit_record hrec;
  if(world->hit(r, 0.001, FLT_MAX, hrec, rng)) {
    if(u) {
      return((unit_vector(hrec.dpdu) + 1)/2);
    } else {
      return((unit_vector(hrec.dpdv) + 1)/2);
    }
  } else {
    return(vec3f(0,0,0));
  }
}

inline point3f calculate_color(const ray& r, hitable *world, random_gen &rng) {
  hit_record hrec;
  scatter_record srec;
  ray r2 = r;
  bool invisible = false;
  if(world->hit(r2, 0.001, FLT_MAX, hrec, rng)) {
    point3f emit = hrec.mat_ptr->emitted(r2, hrec, hrec.u, hrec.v, hrec.p, invisible);
    if(emit.x() != 0 || emit.y() != 0 || emit.z() != 0) {
      return(emit);
    }
    if(hrec.mat_ptr->scatter(r2, hrec, srec, rng)) { //generates scatter record, world space
      if(srec.is_specular) { 
        return(point3f(1,1,1));
      }
      return(hrec.mat_ptr->get_albedo(r2, hrec));
    } else {
      return(point3f(0,0,0));
    }
  } else {
    return(point3f(0,0,0));
  }
}

inline point3f quick_render(const ray& r, hitable *world, random_gen &rng, vec3f lightdir, Float n) {
  hit_record hrec;
  scatter_record srec;
  ray r2 = r;
  bool invisible = false;
  if(world->hit(r2, 0.001, FLT_MAX, hrec, rng)) {
    point3f emit = hrec.mat_ptr->emitted(r2, hrec, hrec.u, hrec.v, hrec.p, invisible);
    if(emit.x() != 0 || emit.y() != 0 || emit.z() != 0) {
      return(emit);
    }
    if(hrec.mat_ptr->scatter(r2, hrec, srec, rng)) { //generates scatter record, world space
      if(srec.is_specular) { 
        return(point3f(1,1,1));
      }
      normal3f normal = hrec.has_bump ? hrec.bump_normal : hrec.normal;
      vec3f R = Reflect(lightdir, hrec.normal); 
      return(hrec.mat_ptr->get_albedo(r2, hrec) * (dot(normal, lightdir)+1)/2 + 
             std::pow(std::fmax(0.f, dot(R, -unit_vector(r.direction()))), n));
    } else {
      return(point3f(0,0,0));
    }
  } else {
    return(point3f(0,0,0));
  }
}
// //Does not take into account moving objects
// void calculate_inside(const ray& r_in, hitable *world, random_gen rng) {
//   hit_record hrec;
//   if(world->hit(r_in, 0.001, FLT_MAX, hrec, rng)) {
//     if(hrec.mat_ptr->is_dielectric()) {
//       bool encountered = false;
//       int current_index = -1;
//       for(size_t i = 0; i < r_in.pri_stack->size(); i++) {
//         if(r_in.pri_stack->at(i) == hrec.mat_ptr) {
//           encountered = true;
//           current_index = i;
//         } 
//       }
//       if(dot(r_in.direction(),hrec.normal) > 0) {
//         if(encountered) {
//           r_in.pri_stack->erase(r_in.pri_stack->begin() + current_index);
//         } 
//       } else {
//         r_in.pri_stack->push_back((dielectric*)hrec.mat_ptr);
//       }
//     }
//     ray r = ray(hrec.p, r_in.direction(),  r_in.pri_stack);
//     calculate_inside(r, world, rng);
//   } 
// }

inline point3f calculate_position(const ray& r, hitable *world, hitable_list *hlist,
                                  size_t max_depth, random_gen &rng) {
  point3f final_color(0,0,0);
  ray r1 = r;
  ray r2 = r;
  bool diffuse_bounce = false;
  for(size_t i = 0; i < max_depth; i++) {
    bool is_invisible = false;
    hit_record hrec;
    if(world->hit(r2, 0.001, FLT_MAX, hrec, rng)) { //generated hit record, world space
      scatter_record srec;
      //Some lights can be invisible until after diffuse bounce
      //If so, generate new ray with intersection point and continue ray
      if(is_invisible && !diffuse_bounce) {
        r2.A = hrec.p;
        continue;
      }
      if(hrec.mat_ptr->scatter(r2, hrec, srec, rng)) { //generates scatter record, world space
        if(srec.is_specular) { //returns specular ray
          if(i == max_depth-1) {
            return(hrec.p);
          }
          r2 = srec.specular_ray;
          continue;
        }
        hitable_pdf p_imp(hlist, hrec.p); //creates pdf of all objects to be sampled
        mixture_pdf p(&p_imp, srec.pdf_ptr); //creates mixture pdf of surface intersected at hrec.p and all sampled objects/lights
        
        point3f offset_p = offset_ray(hrec.p-r2.A, hrec.normal) + r2.A;
        
        r1 = r2;
        vec3f dir;
        dir = p.generate(rng, diffuse_bounce, r2.time());
        if((dir.x() == 0 && dir.y() == 0 && dir.z() == 0)) {
          return(hrec.p);
        }
        r2 = ray(OffsetRayOrigin(offset_p, hrec.pError, hrec.normal, dir), dir, r2.pri_stack, r2.time());
        
        if(i == max_depth-1) {
          return(hrec.p);
        }
      } else {
        return(hrec.p);
      }
    } else {
      return(vec3f(0,0,0));
    }
  }
  return(vec3f(0,0,0));
}

inline point3f calculate_bounce_dir(const ray& r, hitable *world, hitable_list *hlist,
              size_t max_depth, random_gen& rng) {
  point3f final_color(0,0,0);
  ray r1 = r;
  ray r2 = r;
  bool diffuse_bounce = false;
  for(size_t i = 0; i < max_depth; i++) {
    bool is_invisible = false;
    hit_record hrec;
    if(world->hit(r2, 0.001, FLT_MAX, hrec, rng)) { //generated hit record, world space
      scatter_record srec;
      //Some lights can be invisible until after diffuse bounce
      //If so, generate new ray with intersection point and continue ray
      if(is_invisible && !diffuse_bounce) {
        r2.A = hrec.p;
        continue;
      }
      if(hrec.mat_ptr->scatter(r2, hrec, srec, rng)) { //generates scatter record, world space
        if(srec.is_specular) { //returns specular ray
          if(i == max_depth-1) {
            return(unit_vector(srec.specular_ray.direction()));
          }
          r2 = srec.specular_ray;
          continue;
        }
        hitable_pdf p_imp(hlist, hrec.p); //creates pdf of all objects to be sampled
        mixture_pdf p(&p_imp, srec.pdf_ptr); //creates mixture pdf of surface intersected at hrec.p and all sampled objects/lights
        
        point3f offset_p = offset_ray(hrec.p-r2.A, hrec.normal) + r2.A;
        
        r1 = r2;
        vec3f dir;
        dir = p.generate(rng, diffuse_bounce, r2.time());
        if((dir.x() == 0 && dir.y() == 0 && dir.z() == 0)) {
          return(dir);
        }
        r2 = ray(OffsetRayOrigin(offset_p, hrec.pError, hrec.normal, dir), dir, r2.pri_stack, r2.time());
        
        if(i == max_depth-1) {
          return(unit_vector(dir));
        }
      } else {
        return(vec3f(0,0,0));
      }
    } else {
      return(vec3f(0,0,0));
    }
  }
  return(vec3f(0,0,0));
}

inline Float calculate_time(const ray& r, hitable *world, hitable_list *hlist,
                                    size_t max_depth, random_gen& rng) {
  hit_record hrec;
  if(world->hit(r, 0.001, FLT_MAX, hrec, rng)) {
    return(hrec.t);
  } else {
    return(INFINITY);
  }
}

inline uint32_t hash32(uint32_t x) {
    x ^= x >> 16;
    x *= 0x45d9f3bU;
    x ^= x >> 16;
    x *= 0x45d9f3bU;
    x ^= x >> 16;
    return x;
}

inline point3f calculate_shape(const ray& r, hitable *world, hitable_list *hlist,
                            size_t max_depth, random_gen& rng) {
  hit_record hrec;
  if(world->hit(r, 0.001, FLT_MAX, hrec, rng)) {
    uint32_t r =  reinterpret_type<material*, uint32_t>(hrec.mat_ptr) % 65536;
    uint32_t g = (reinterpret_type<material*, uint32_t>(hrec.mat_ptr) % 65536) + 1;
    uint32_t b = (reinterpret_type<material*, uint32_t>(hrec.mat_ptr) % 65536) + 2;
    r = hash32(r) % 128;
    g = hash32(g) % 128;
    b = hash32(b) % 128;
    Float r2 = (Float)r/128;
    Float g2 = (Float)g/128;
    Float b2 = (Float)b/128;
    return(vec3f(r2,g2,b2));
  } else {
    return(vec3f(0,0,0));
  }
}

inline Float calculate_pdf(const ray& r, hitable *world, hitable_list *hlist,
                                    size_t max_depth, random_gen& rng) {
  point3f final_color(0,0,0);
  ray r1 = r;
  ray r2 = r;
  bool diffuse_bounce = false;
  for(size_t i = 0; i < max_depth; i++) {
    bool is_invisible = false;
    hit_record hrec;
    if(world->hit(r2, 0.001, FLT_MAX, hrec, rng)) { //generated hit record, world space
      scatter_record srec;
      //Some lights can be invisible until after diffuse bounce
      //If so, generate new ray with intersection point and continue ray
      if(is_invisible && !diffuse_bounce) {
        r2.A = OffsetRayOrigin(hrec.p, hrec.pError, hrec.normal, r2.direction());
        continue;
      }
      if(hrec.mat_ptr->scatter(r2, hrec, srec, rng)) { //generates scatter record, world space
        if(srec.is_specular) { //returns specular ray
          if(i == max_depth-1) {
            return(0);
          }
          r2 = srec.specular_ray;
          continue;
        }
        hitable_pdf p_imp(hlist, hrec.p); //creates pdf of all objects to be sampled
        mixture_pdf p(&p_imp, srec.pdf_ptr); //creates mixture pdf of surface intersected at hrec.p and all sampled objects/lights
        point3f offset_p = offset_ray(hrec.p-r2.A, hrec.normal) + r2.A;
        
        r1 = r2;
        vec3f dir;
        dir = p.generate(rng, diffuse_bounce, r2.time());
        Float pdf_val = p.value(dir, rng, r2.time()); //generates a pdf value based the intersection point and the mixture pdf
        if(i == max_depth-1) {
          return(pdf_val);
        }
        if((dir.x() == 0 && dir.y() == 0 && dir.z() == 0)) {
          return(0);
        }
        r2 = ray(OffsetRayOrigin(offset_p, hrec.pError, hrec.normal, dir), dir, r2.pri_stack, r2.time());
      } else {
        return(0);
      }
    } else {
      return(0);
    }
  }
  return(0);
}

inline Float calculate_error(const ray& r, hitable *world, hitable_list *hlist,
                           size_t max_depth, random_gen& rng) {
  point3f final_color(0,0,0);
  ray r1 = r;
  ray r2 = r;
  bool diffuse_bounce = false;
  for(size_t i = 0; i < 2; i++) {
    bool is_invisible = false;
    hit_record hrec;
    if(world->hit(r2, 0.001, FLT_MAX, hrec, rng)) { //generated hit record, world space
      scatter_record srec;
      //Some lights can be invisible until after diffuse bounce
      //If so, generate new ray with intersection point and continue ray
      if(is_invisible && !diffuse_bounce) {
        r2.A = hrec.p;
        continue;
      }
      if(hrec.mat_ptr->scatter(r2, hrec, srec, rng)) { //generates scatter record, world space
        if(srec.is_specular) { //returns specular ray
          return(0);
        }
        hitable_pdf p_imp(hlist, hrec.p); //creates pdf of all objects to be sampled
        mixture_pdf p(&p_imp, srec.pdf_ptr); //creates mixture pdf of surface intersected at hrec.p and all sampled objects/lights
        
        r1 = r2;
        vec3f dir;
        dir = p.generate(rng, diffuse_bounce, r2.time());
        Float offset = (OffsetRayOrigin(hrec.p, hrec.pError, hrec.normal, dir) - hrec.p).length();
        return(offset);
      } else {
        return(0);
      }
    } else {
      return(0);
    }
  }
  return(0);
}

inline Float calculate_bounces(const ray& r, hitable *world, hitable_list *hlist,
                             size_t max_depth, random_gen& rng) {
  point3f final_color(0,0,0);
  point3f emit_color(0,0,0);
  
  point3f throughput(1,1,1);
  ray r1 = r;
  ray r2 = r;
  bool diffuse_bounce = false;
  for(size_t i = 0; i < max_depth; i++) {
    // RcppThread::Rcout << "Ray origin: " << r2.A << "\n";
    bool is_invisible = false;
    hit_record hrec;
    if(world->hit(r2, 0.001, FLT_MAX, hrec, rng)) { //generated hit record, world space
      scatter_record srec;
      emit_color = throughput * hrec.mat_ptr->emitted(r2, hrec, hrec.u, hrec.v, hrec.p, is_invisible);
      //Some lights can be invisible until after diffuse bounce
      //If so, generate new ray with intersection point and continue ray
      if(is_invisible && !diffuse_bounce) {
        r2.A = OffsetRayOrigin(hrec.p, hrec.pError, hrec.normal, r2.direction());
        continue;
      }
      final_color += emit_color;
      if(throughput.x() == 0 && throughput.y() == 0 && throughput.z() == 0) {
        return(i);
      }
      float pdf_val;
      //generates scatter record and sends out new ray, otherwise exits out with accumulated color
      if(hrec.mat_ptr->scatter(r2, hrec, srec, rng)) { 
        if(srec.is_specular) { //returns specular ray
          r2 = srec.specular_ray;
          throughput *= srec.attenuation;
          continue;
        }
        hitable_pdf p_imp(hlist, hrec.p); //creates pdf of all objects to be sampled
        mixture_pdf p(&p_imp, srec.pdf_ptr); //creates mixture pdf of surface intersected at hrec.p and all sampled objects/lights
        
        //Generates a scatter direction (with origin hrec.p) from the mixture 
        //and saves surface normal from light to use in pdf_value calculation
        //(along with the scatter direction)
        
        //Translates the world space point into object space point, generates ray assuring intersection, and then translates 
        //ray back into world space
        //Remove this when error analysis fully implemented
        // point3f offset_p = offset_ray(hrec.p-r2.A, hrec.normal) + r2.A;
        
        r1 = r2;
        vec3f dir;
        dir = p.generate(rng, diffuse_bounce, r2.time()); //scatters a ray from hit point to random direction
        r2 = ray(OffsetRayOrigin(hrec.p, hrec.pError, hrec.normal, dir), dir, r2.pri_stack, r2.time());
        
        pdf_val = p.value(dir, rng, r2.time()); //generates a pdf value based the intersection point and the mixture pdf
        throughput *= hrec.mat_ptr->f(r1, hrec, r2) / pdf_val;
      } else {
        return(i);
      }
    } else {
      return(i);
    }
  }
  return(max_depth);
}

void debug_scene(size_t numbercores, size_t nx, size_t ny, size_t ns, int debug_channel,
                 Float min_variance, size_t min_adaptive_size, 
                 Rcpp::NumericMatrix& routput, Rcpp::NumericMatrix& goutput, Rcpp::NumericMatrix& boutput,
                 bool progress_bar, int sample_method, Rcpp::NumericVector& stratified_dim,
                 bool verbose, ortho_camera& ocam, camera &cam, environment_camera &ecam, Float fov,
                 hitable_list& world, hitable_list& hlist,
                 Float clampval, size_t max_depth, size_t roulette_active,
                 Rcpp::NumericVector& light_direction, random_gen& rng);



#endif
