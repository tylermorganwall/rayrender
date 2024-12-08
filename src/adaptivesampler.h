#ifndef ADAPTIVESAMPLERH
#define ADAPTIVESAMPLERH

#include "vectypes.h"
#include "RayMatrix.h"

using namespace Rcpp;

struct pixel_block {
  size_t startx, starty;
  size_t endx, endy;
  size_t split_axis;
  size_t split_pos;
  bool erase;
  bool split;
  float error;
};

class adaptive_sampler {
public:
  adaptive_sampler(size_t _numbercores, size_t nx, size_t ny, size_t ns, int debug_channel,
                   float min_variance, size_t min_adaptive_size, 
                   RayMatrix& rgb,  
                   RayMatrix& rgb2, 
                   RayMatrix& normalOutput, 
                   RayMatrix& albedoOutput,
                   RayMatrix& alpha, 
                   RayMatrix& draw_rgb_output,
                   bool adaptive_on);
  void reset();
  ~adaptive_sampler() {}
  void test_for_convergence(size_t k, size_t s,
                            size_t nx_end, size_t nx_begin,
                            size_t ny_end, size_t ny_begin);
  
  void split_remove_chunks(size_t s);
  void write_final_pixels();
  void add_color_main(size_t i, size_t j, point3f color);
  
  void add_color_sec(size_t i, size_t j, point3f color);
  //For use when s = 1 in small image preview
  void set_color_main(size_t i, size_t j, point3f color);
  void add_alpha_count(size_t i, size_t j);
  void add_albedo(size_t i, size_t j, point3f albedo);
  void add_normal(size_t i, size_t j, normal3f normal);

  size_t size() {return(pixel_chunks.size());}

  size_t numbercores;
  size_t nx, ny, ns;
  size_t max_s;
  int debug_channel;
  float min_variance;
  size_t min_adaptive_size;
  RayMatrix &rgb, &rgb2, &normalOutput, &albedoOutput, &draw_rgb_output;
  RayMatrix &a;
  std::vector<pixel_block> pixel_chunks;
  std::vector<bool> finalized;
  std::vector<bool> just_finalized;
  bool adaptive_on;
  
};

#endif
