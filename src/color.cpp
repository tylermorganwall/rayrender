#include "color.h"
#include "RcppThread.h"
#include "mathinline.h"

// #include "fstream"
// #define DEBUG

point3f color(const ray& r, hitable *world, hitable_list *hlist,
           size_t max_depth, size_t roulette_activate, random_gen& rng, Sampler* sampler) {
#ifdef DEBUG
  std::ofstream myfile;
  myfile.open("rays.txt", std::ios::app | std::ios::out);
#endif
  point3f final_color(0,0,0);
  point3f emit_color(0,0,0);
  
  point3f throughput(1,1,1);
  float prev_t = 1;
  ray r1 = r;
  ray r2 = r;
  bool diffuse_bounce = false;
  for(size_t i = 0; i < max_depth; i++) {
    bool is_invisible = false;
    hit_record hrec;
    if(world->hit(r2, 0.001, FLT_MAX, hrec, rng)) { //generated hit record, world space
      scatter_record srec;
      if(hrec.alpha_miss) {
        r2.A = hrec.p;
        continue;
      }
      emit_color = throughput * hrec.mat_ptr->emitted(r2, hrec, hrec.u, hrec.v, hrec.p, is_invisible);
      //Some lights can be invisible until after diffuse bounce
      //If so, generate new ray with intersection point and continue ray
      if(is_invisible && !diffuse_bounce) {
        r2.A = OffsetRayOrigin(hrec.p, hrec.pError, hrec.normal, r2.direction());
        continue;
      }
      final_color += emit_color;
      if(throughput.x() == 0 && throughput.y() == 0 && throughput.z() == 0) {
        return(point3f(0,0,0));
      }
      if(i > roulette_activate) {
        float t = std::fmax(throughput.x(), std::fmax(throughput.y(), throughput.z()));
        //From Szecsi, Szirmay-Kalos, and Kelemen
        float prob_continue = std::min(1.0f, std::sqrt(t/prev_t));
        prev_t = t;
        if(rng.unif_rand() > prob_continue) {
          return(final_color);
        }
        throughput *= 1 / prob_continue;
      }
#ifdef DEBUG
        myfile << i << ", " << r2.A << " ";
        myfile << ", " << hrec.p << ", " << hrec.normal << ", " << r2.direction() << ", " << throughput << "\n ";
#endif
      float pdf_val;
      //generates scatter record and sends out new ray, otherwise exits out with accumulated color
      if(hrec.mat_ptr->scatter(r2, hrec, srec, sampler)) { 
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
        r1 = r2;
        vec3f dir;
        if(!diffuse_bounce) {
          //`diffuse_bounce` switched by generate()
          dir = p.generate(sampler, diffuse_bounce, r2.time()); //scatters a ray from hit point to stratified direction
        } else {
          dir = p.generate(rng, diffuse_bounce, r2.time()); //scatters a ray from hit point to random direction
        }
        
        r2 = ray(OffsetRayOrigin(hrec.p, hrec.pError, hrec.normal, dir), dir, r2.pri_stack, r2.time());
        
        pdf_val = p.value(dir, rng, r2.time()); //generates a pdf value based the intersection point and the mixture pdf

        if(pdf_val == 0) {
          break;
        }
        
        if((dir.x() == 0 && dir.y() == 0 && dir.z() == 0)) {
          break;
        }
        
        throughput *= hrec.mat_ptr->f(r1, hrec, r2) / pdf_val;

      } else {
        return(final_color);
      }
    } else {
      return(final_color);
    }
  }
#ifdef DEBUG
  myfile.close();
#endif
  return(final_color);
}
