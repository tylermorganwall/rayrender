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
                   RayMatrix& r2, RayMatrix& g2, RayMatrix& b2) : 
                   numbercores(_numbercores), nx(nx), ny(ny), ns(ns), max_s(0), debug_channel(debug_channel), 
                   min_variance(min_variance), min_adaptive_size(min_adaptive_size),
                   r(r), g(g), b(b), r2(r2), g2(g2), b2(b2) {
    size_t nx_chunk = nx / numbercores;
    size_t ny_chunk = ny / numbercores;
    size_t bonus_x = nx - nx_chunk * numbercores;
    size_t bonus_y = ny - ny_chunk * numbercores;
    finalized.resize(nx*ny, false);
    just_finalized.resize(nx*ny, true);
    for(size_t i = 0; i < numbercores; i++) {
      for(size_t j = 0; j < numbercores ; j++) {
        size_t extra_x = i == numbercores - 1 ? bonus_x : 0;
        size_t extra_y = j == numbercores - 1 ? bonus_y : 0;
        pixel_block chunk = {i*nx_chunk, j*ny_chunk,
                             (i+1)*nx_chunk + extra_x, (j+1)*ny_chunk  + extra_y,
                             0, 0, false, false, 0};
        pixel_chunks.push_back(chunk);
      }
    }
  }
  void reset() {
    pixel_chunks.clear();
    size_t nx_chunk = nx / numbercores;
    size_t ny_chunk = ny / numbercores;
    size_t bonus_x = nx - nx_chunk * numbercores;
    size_t bonus_y = ny - ny_chunk * numbercores;
    finalized.resize(nx*ny, false);
    just_finalized.resize(nx*ny, true);
    
    for(size_t i = 0; i < numbercores; i++) {
      for(size_t j = 0; j < numbercores ; j++) {
        size_t extra_x = i == numbercores - 1 ? bonus_x : 0;
        size_t extra_y = j == numbercores - 1 ? bonus_y : 0;
        pixel_block chunk = {i*nx_chunk, j*ny_chunk,
                             (i+1)*nx_chunk + extra_x, (j+1)*ny_chunk  + extra_y,
                             0, 0, false, false, 0};
        pixel_chunks.push_back(chunk);
      }
    }
    std::fill(finalized.begin(), finalized.end(), false);
    std::fill(just_finalized.begin(), just_finalized.end(), true);
    std::fill(r.begin(), r.end(), 0);
    std::fill(g.begin(), g.end(), 0);
    std::fill(b.begin(), b.end(), 0);
    std::fill(r2.begin(), r2.end(), 0);
    std::fill(g2.begin(), g2.end(), 0);
    std::fill(b2.begin(), b2.end(), 0);
  }
  ~adaptive_sampler() {}
  void test_for_convergence(size_t k, size_t s,
                            size_t nx_end, size_t nx_begin,
                            size_t ny_end, size_t ny_begin) {
    Float error_block = 0.0;
    size_t nx_block = (nx_end-nx_begin);
    size_t ny_block = (ny_end-ny_begin);
    Float N = (Float)nx_block * (Float)ny_block;
    Float r_b = std::sqrt(N / ((Float)nx * (Float)ny));
    std::vector<Float> error_sum(nx_block * ny_block, 0);
    for(size_t i = nx_begin; i < nx_end; i++) {
      for(size_t j = ny_begin; j < ny_end; j++) {
        error_sum[(i-nx_begin) + (j-ny_begin) * nx_block] = fabs(r(i,j) - 2 * r2(i,j)) +
          fabs(g(i,j) - 2 * g2(i,j)) +
          fabs(b(i,j) - 2 * b2(i,j));
        error_sum[(i-nx_begin) + (j-ny_begin) * nx_block] *= r_b / (s*N);
        Float normalize = sqrt(r(i,j) + g(i,j)  + b(i,j));
        if(normalize != 0) {
          error_sum[(i-nx_begin) + (j-ny_begin) * nx_block] /= normalize;
        }
        error_block += error_sum[(i-nx_begin) + (j-ny_begin) * nx_block];
      }
    }
    pixel_chunks[k].error = error_block;
    if(error_block < min_variance) {
      pixel_chunks[k].erase = true;
    } else if(error_block < min_variance*256) {
      pixel_chunks[k].split = true;
      Float error_half = 0.0f;
      if((nx_end-nx_begin) >= (ny_end-ny_begin)) {
        pixel_chunks[k].split_axis = 0;
        for(size_t i = nx_begin; i < nx_end; i++) {
          for(size_t j = ny_begin; j < ny_end; j++) {
            error_half += error_sum[(i-nx_begin) + (j-ny_begin) * nx_block];
          }
          if(error_half >= error_block/2) {
            pixel_chunks[k].split_pos = i;
            break;
          }
        }
      } else {
        pixel_chunks[k].split_axis = 1;
        for(size_t j = ny_begin; j < ny_end; j++) {
          for(size_t i = nx_begin; i < nx_end; i++) {
            error_half += error_sum[(i-nx_begin) + (j-ny_begin) * nx_block];
          }
          if(error_half >= error_block/2) {
            pixel_chunks[k].split_pos = j;
            break;
          }
        }
      }
    }
  }
  
  void split_remove_chunks(size_t s) {
    auto it = pixel_chunks.begin();
    std::vector<pixel_block> temppixels;
    while(it != pixel_chunks.end()) {
      if(it->erase) {
        for(size_t i = it->startx; i < it->endx; i++) {
          for(size_t j = it->starty; j < it->endy; j++) {
            r(i,j) /= (float)(s+1);
            g(i,j) /= (float)(s+1);
            b(i,j) /= (float)(s+1);
            if(debug_channel == 5) {
              r(i,j) = (float)(s+1)/(float)ns;
              g(i,j) = (float)(s+1)/(float)ns;
              b(i,j) = (float)(s+1)/(float)ns;
            }
            finalized[i + nx*j] = true;
          }
        }
      } else if(it->split &&
        (it->endx - it->startx) > min_adaptive_size &&
        (it->endy - it->starty) > min_adaptive_size) {
        if(it->split_axis == 1) {
          pixel_block b1 = {it->startx, it->starty,
                            it->endx, it->split_pos,
                            0, 0, false, false, 0};
          pixel_block b2 = {it->startx, it->split_pos,
                            it->endx, it->endy,
                            0, 0, false, false, 0};
          temppixels.push_back(b1);
          temppixels.push_back(b2);
        } else if(it->split_axis == 0) {
          pixel_block b1 = {it->startx, it->starty,
                            it->split_pos, it->endy,
                            0, 0, false, false, 0};
          pixel_block b2 = {it->split_pos, it->starty,
                            it->endx, it->endy,
                            0, 0, false, false, 0};
          temppixels.push_back(b1);
          temppixels.push_back(b2);
        } 
      } else {
        temppixels.push_back(*it);
      }
      it++;
    }
    pixel_chunks = temppixels;
  }
  void write_final_pixels() {
    auto it = pixel_chunks.begin();
    while(it != pixel_chunks.end()) {
      for(size_t i = it->startx; i < it->endx; i++) {
        for(size_t j = it->starty; j < it->endy; j++) {
          r(i,j) /= (float)ns;
          g(i,j) /= (float)ns;
          b(i,j) /= (float)ns;
          if(debug_channel == 5) {
            r(i,j) = (float)max_s/(float)ns;
            g(i,j) = (float)max_s/(float)ns;
            b(i,j) = (float)max_s/(float)ns;
          }
        }
      }
      it++;
    }
  }
  void add_color_main(size_t i, size_t j, point3f color) {
    r(i,j) += color.r();
    g(i,j) += color.g();
    b(i,j) += color.b();
  }
  
  void add_color_sec(size_t i, size_t j, point3f color) {
    r2(i,j) += color.r();
    g2(i,j) += color.g();
    b2(i,j) += color.b();
  }
  //For use when s = 1 in small image preview
  void set_color_main(size_t i, size_t j, point3f color) {
    r(i,j) = color.r();
    g(i,j) = color.g();
    b(i,j) = color.b();
  }
  size_t size() {return(pixel_chunks.size());}

  size_t numbercores;
  size_t nx, ny, ns;
  size_t max_s;
  int debug_channel;
  float min_variance;
  size_t min_adaptive_size;
  RayMatrix &r, &g, &b, &r2, &g2, &b2;
  std::vector<pixel_block> pixel_chunks;
  std::vector<bool> finalized;
  std::vector<bool> just_finalized;
  
};

#endif
