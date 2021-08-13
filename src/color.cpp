#include "color.h"
#include "RcppThread.h"

point3f color(const ray& r, hitable *world, hitable_list *hlist,
           size_t max_depth, size_t roulette_activate, random_gen& rng, Sampler* sampler) {
#ifdef DEBUG
  ofstream myfile;
  myfile.open("rays.txt", ios::app | ios::out);
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
#ifdef DEBUG
      myfile << i << ", " << r2.A << " ";
      myfile << ", " << hrec.p << "\n ";
#endif
      scatter_record srec;
      emit_color = throughput * hrec.mat_ptr->emitted(r2, hrec, hrec.u, hrec.v, hrec.p, is_invisible);
      //Some lights can be invisible until after diffuse bounce
      //If so, generate new ray with intersection point and continue ray
      if(is_invisible && !diffuse_bounce) {
        r2.A = hrec.p;
        continue;
      }
      final_color += emit_color;
      if(throughput.x() == 0 && throughput.y() == 0 && throughput.z() == 0) {
        return(point3f(0,0,0));
      }
      if(i > roulette_activate) {
        float t = std::max(throughput.x(), std::max(throughput.y(), throughput.z()));
        //From Szecsi, Szirmay-Kalos, and Kelemen
        float prob_continue = std::min(1.0f, std::sqrt(t/prev_t));
        prev_t = t;
        if(rng.unif_rand() > prob_continue) {
          return(final_color);
        }
        throughput *= 1 / prob_continue;
      }
      float pdf_val;
      if(hrec.mat_ptr->scatter(r2, hrec, srec, sampler)) { //generates scatter record, world space
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
        point3f offset_p = offset_ray(hrec.p-r2.A, hrec.normal) + r2.A;
        
        r1 = r2;
        if(!diffuse_bounce) {
          r2 = ray(offset_p, p.generate(sampler, diffuse_bounce, r2.time()), r2.pri_stack, r2.time()); //scatters a ray from hit point to stratified direction
        } else {
          r2 = ray(offset_p, p.generate(rng, diffuse_bounce, r2.time()), r2.pri_stack, r2.time()); //scatters a ray from hit point to direction
        }
        pdf_val = p.value(r2.direction(), sampler, r2.time()); //generates a pdf value based the intersection point and the mixture pdf
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
