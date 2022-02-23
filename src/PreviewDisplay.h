#ifndef PREVIEWDISPLAYH
#define PREVIEWDISPLAYH

#include <memory>
#include "Rcpp.h"
#include "RProgress.h"
#include "adaptivesampler.h"
#include "camera.h"
#include "hitable.h"

#ifdef RAY_HAS_X11

#include <X11/Xlib.h>

class PreviewDisplay {
public: 
  PreviewDisplay(unsigned int _width, unsigned int _height, bool preview, bool _interactive,
                 Float initial_lookat_distance);
  ~PreviewDisplay();
  void DrawImage(adaptive_sampler& adaptive_pixel_sampler, size_t &ns,
                 RProgress::RProgress &pb, bool progress,
                 RayCamera* cam,Float percent_done,
                 hitable* world, random_gen& rng);
  Display *d;
  XImage *img;
  std::unique_ptr<char[]> data;
  Window w;
  XEvent e;
  unsigned int width;
  unsigned int height;
  int s;
  bool terminate;
  Float speed;
  bool interactive;
  bool orbit;
  Float base_step;
};

#else

class PreviewDisplay {
public: 
  PreviewDisplay(unsigned int _width, unsigned int _height, bool preview, bool _interactive,
                 Float initial_lookat_distance) {};
  void DrawImage(adaptive_sampler& adaptive_pixel_sampler, size_t &ns,
                 RProgress::RProgress &pb, bool progress,
                 RayCamera* cam,Float percent_done,
                 hitable* world, random_gen& rng) {};
  bool terminate;
  bool interactive;
};

#endif

#endif