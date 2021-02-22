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
inline double debug_bvh(const ray& r, hitable *world, random_gen &rng) {
  hit_record hrec;
  hrec.bvh_nodes = 0.0;
  if(world->hit(r, 0.001, FLT_MAX, hrec, rng)) {
    return(hrec.bvh_nodes);
  } else {
    return(0.0);
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

inline vec3 calculate_normals(const ray& r, hitable *world, random_gen &rng) {
  hit_record hrec;
  if(world->hit(r, 0.001, FLT_MAX, hrec, rng)) {
    hrec.normal.make_unit_vector();
    return((vec3(1,1,1) + hrec.normal)/2);
  } else {
    return(vec3(0,0,0));
  }
}

inline vec3 calculate_uv(const ray& r, hitable *world, random_gen &rng) {
  hit_record hrec;
  if(world->hit(r, 0.001, FLT_MAX, hrec, rng)) {
    return(vec3(hrec.u,hrec.v,1-hrec.u-hrec.v));
  } else {
    return(vec3(0,0,0));
  }
}

inline vec3 calculate_dpduv(const ray& r, hitable *world, random_gen &rng, bool u) {
  hit_record hrec;
  if(world->hit(r, 0.001, FLT_MAX, hrec, rng)) {
    if(u) {
      return((unit_vector(hrec.dpdu) + 1)/2);
    } else {
      return((unit_vector(hrec.dpdv) + 1)/2);
    }
  } else {
    return(vec3(0,0,0));
  }
}

inline vec3 calculate_color(const ray& r, hitable *world, random_gen &rng) {
  hit_record hrec;
  scatter_record srec;
  ray r2 = r;
  if(world->hit(r2, 0.001, FLT_MAX, hrec, rng)) {
    vec3 emit = hrec.mat_ptr->emitted(r2, hrec, hrec.u, hrec.v, hrec.p);
    if(emit.x() != 0 || emit.y() != 0 || emit.z() != 0) {
      return(emit);
    }
    if(hrec.mat_ptr->scatter(r2, hrec, srec, rng)) { //generates scatter record, world space
      if(srec.is_specular) { 
        return(vec3(1,1,1));
      }
      return(hrec.mat_ptr->get_albedo(r2, hrec));
    } else {
      return(vec3(0,0,0));
    }
  } else {
    return(vec3(0,0,0));
  }
}

inline vec3 quick_render(const ray& r, hitable *world, random_gen &rng, vec3 lightdir, Float n) {
  hit_record hrec;
  scatter_record srec;
  ray r2 = r;
  if(world->hit(r2, 0.001, FLT_MAX, hrec, rng)) {
    vec3 emit = hrec.mat_ptr->emitted(r2, hrec, hrec.u, hrec.v, hrec.p);
    if(emit.x() != 0 || emit.y() != 0 || emit.z() != 0) {
      return(emit);
    }
    if(hrec.mat_ptr->scatter(r2, hrec, srec, rng)) { //generates scatter record, world space
      if(srec.is_specular) { 
        return(vec3(1,1,1));
      }
      vec3 normal = hrec.has_bump ? hrec.bump_normal : hrec.normal;
      vec3 R = Reflect(lightdir, hrec.normal); 
      return(hrec.mat_ptr->get_albedo(r2, hrec) * (dot(normal, lightdir)+1)/2 + 
             std::pow(std::max(0.f, dot(R, -unit_vector(r.direction()))), n));
    } else {
      return(vec3(0,0,0));
    }
  } else {
    return(vec3(0,0,0));
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


void debug_scene(size_t numbercores, size_t nx, size_t ny, size_t ns, int debug_channel,
                 Float min_variance, size_t min_adaptive_size, 
                 Rcpp::NumericMatrix& routput, Rcpp::NumericMatrix& goutput, Rcpp::NumericMatrix& boutput,
                 bool progress_bar, int sample_method, Rcpp::NumericVector& stratified_dim,
                 bool verbose, ortho_camera& ocam, camera &cam, Float fov,
                 hitable_list& world, hitable_list& hlist,
                 Float clampval, size_t max_depth, size_t roulette_active,
                 Rcpp::NumericVector& light_direction, random_gen& rng);



#endif