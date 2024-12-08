#include "adaptivesampler.h"
#include "RayMatrix.h"
#include "raylog.h"

adaptive_sampler::adaptive_sampler(size_t _numbercores, size_t nx, size_t ny, size_t ns, int debug_channel,
                 float min_variance, size_t min_adaptive_size, 
                 RayMatrix& rgb, 
                 RayMatrix& rgb2, 
                 RayMatrix& normalOutput, 
                 RayMatrix& albedoOutput,
                 RayMatrix& alpha, 
                 RayMatrix& draw_rgb_output, bool adaptive_on) : 
    numbercores(_numbercores), nx(nx), ny(ny), ns(ns), max_s(0), debug_channel(debug_channel), 
    min_variance(min_variance), min_adaptive_size(min_adaptive_size),
    rgb(rgb), rgb2(rgb2), normalOutput(normalOutput), albedoOutput(albedoOutput),
    draw_rgb_output(draw_rgb_output),
    a(alpha), adaptive_on(adaptive_on) {
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
void adaptive_sampler::reset() {
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
  rgb.reset();
  rgb2.reset();
  normalOutput.reset();
  albedoOutput.reset();
  draw_rgb_output.reset();
  a.reset();
}

void adaptive_sampler::test_for_convergence(size_t k, size_t s,
                          size_t nx_end, size_t nx_begin,
                          size_t ny_end, size_t ny_begin) {
  SCOPED_CONTEXT("Overall");
  SCOPED_TIMER_COUNTER("Convergence Test");
  if (!adaptive_on) {
    return; // Skip convergence test if adaptive sampling is off
  }
  Float error_block = 0.0;
  size_t nx_block = (nx_end-nx_begin);
  size_t ny_block = (ny_end-ny_begin);
  Float N = (Float)nx_block * (Float)ny_block;
  Float r_b = std::sqrt(N / ((Float)nx * (Float)ny));
  std::vector<Float> error_sum(nx_block * ny_block, 0);
  for(size_t i = nx_begin; i < nx_end; i++) {
    for(size_t j = ny_begin; j < ny_end; j++) {
      error_sum[(i-nx_begin) + (j-ny_begin) * nx_block] = 
        fabs(rgb(i,j,0) - 2 * rgb2(i,j,0)) +
        fabs(rgb(i,j,1) - 2 * rgb2(i,j,1)) +
        fabs(rgb(i,j,2) - 2 * rgb2(i,j,2));
      error_sum[(i-nx_begin) + (j-ny_begin) * nx_block] *= r_b / (s*N);
      Float normalize = sqrt(rgb(i,j,0) + rgb(i,j,1)  + rgb(i,j,2));
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
  
void adaptive_sampler::split_remove_chunks(size_t s) {
  SCOPED_CONTEXT("Overall");
  SCOPED_TIMER_COUNTER("Split/Remove Chunks");
  if (!adaptive_on) {
    return; // Skip convergence test if adaptive sampling is off
  }
  auto it = pixel_chunks.begin();
  std::vector<pixel_block> temppixels;
  while(it != pixel_chunks.end()) {
    if(it->erase) {
      for(size_t i = it->startx; i < it->endx; i++) {
        for(size_t j = it->starty; j < it->endy; j++) {
          rgb(i,j,0) /= (float)(s+1);
          rgb(i,j,1) /= (float)(s+1);
          rgb(i,j,2) /= (float)(s+1);
          a(i,j,0) = 1.f - a(i,j,0)/(float)(s+1);
          if(debug_channel == 5) {
            rgb(i,j,0) = (float)(s+1)/(float)ns;
            rgb(i,j,1) = (float)(s+1)/(float)ns;
            rgb(i,j,2) = (float)(s+1)/(float)ns;
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
void adaptive_sampler::write_final_pixels() {
  auto it = pixel_chunks.begin();
  while(it != pixel_chunks.end()) {
    for(size_t i = it->startx; i < it->endx; i++) {
      for(size_t j = it->starty; j < it->endy; j++) {
        rgb(i,j,0) /= (float)ns;
        rgb(i,j,1) /= (float)ns;
        rgb(i,j,2) /= (float)ns;
        normalOutput(i,j,0) /= (float)ns;
        normalOutput(i,j,1) /= (float)ns;
        normalOutput(i,j,2) /= (float)ns;
        albedoOutput(i,j,0) /= (float)ns;
        albedoOutput(i,j,1) /= (float)ns;
        albedoOutput(i,j,2) /= (float)ns;
        a(i,j,0)  = 1 - a(i,j,0)/(float)ns;
        if(debug_channel == 5) {
          rgb(i,j,0) = (float)max_s/(float)ns;
          rgb(i,j,1) = (float)max_s/(float)ns;
          rgb(i,j,2) = (float)max_s/(float)ns;
        }
      }
    }
    it++;
  }
}
void adaptive_sampler::add_color_main(size_t i, size_t j, point3f color) {
  rgb(i,j,0) += color.r();
  rgb(i,j,1) += color.g();
  rgb(i,j,2) += color.b();
}
  
void adaptive_sampler::add_color_sec(size_t i, size_t j, point3f color) {
  rgb2(i,j,0) += color.r();
  rgb2(i,j,1) += color.g();
  rgb2(i,j,2) += color.b();
}
//For use when s = 1 in small image preview
void adaptive_sampler::set_color_main(size_t i, size_t j, point3f color) {
  rgb(i,j,0) = color.r();
  rgb(i,j,1) = color.g();
  rgb(i,j,2) = color.b();
}

void adaptive_sampler::add_albedo(size_t i, size_t j, point3f albedo) {
  albedoOutput(i,j,0) += albedo.r();
  albedoOutput(i,j,1) += albedo.g();
  albedoOutput(i,j,2) += albedo.b();
}

void adaptive_sampler::add_normal(size_t i, size_t j, normal3f normal) {
  normalOutput(i,j,0) += normal.r();
  normalOutput(i,j,1) += normal.g();
  normalOutput(i,j,2) += normal.b();
}

void adaptive_sampler::add_alpha_count(size_t i, size_t j) {
  a.add_one(i,j,0); 
}
