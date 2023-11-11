#ifndef ADAPTIVESAMPLERH
#define ADAPTIVESAMPLERH

#include "vec3.h"
#include "point3.h"
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
                   RayMatrix& r, RayMatrix& g, RayMatrix& b,
                   RayMatrix& r2, RayMatrix& g2, RayMatrix& b2,
                   RayMatrix& alpha);
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
  size_t size() {return(pixel_chunks.size());}

  size_t numbercores;
  size_t nx, ny, ns;
  size_t max_s;
  int debug_channel;
  float min_variance;
  size_t min_adaptive_size;
  RayMatrix &r, &g, &b, &r2, &g2, &b2;
  RayMatrix &a;
  std::vector<pixel_block> pixel_chunks;
  std::vector<bool> finalized;
  std::vector<bool> just_finalized;
  
};

#endif
