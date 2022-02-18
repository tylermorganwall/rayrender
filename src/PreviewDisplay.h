#ifndef PREVIEWDISPLAYH
#define PREVIEWDISPLAYH

#ifdef RAY_HAS_X11

#include <X11/Xlib.h>
#include <memory>
#include "Rcpp.h"
#include "RProgress.h"
#include "adaptivesampler.h"
#include "camera.h"
#include "hitable.h"

class PreviewDisplay {
public: 
  PreviewDisplay(unsigned int _width, unsigned int _height, bool preview, bool _interactive);
  ~PreviewDisplay();
  void DrawImage(Rcpp::NumericMatrix& r, Rcpp::NumericMatrix& g, Rcpp::NumericMatrix& b, size_t &ns,
                 std::vector<bool>& finalized, RProgress::RProgress &pb, bool progress,
                 RayCamera* cam, adaptive_sampler& adaptive_pixel_sampler, Float percent_done,
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
};

#else

class PreviewDisplay {
public: 
  PreviewDisplay(unsigned int _width, unsigned int _height);
  ~PreviewDisplay();
  void DrawImage(Rcpp::NumericMatrix& r, Rcpp::NumericMatrix& g, Rcpp::NumericMatrix& b, size_t &ns,
                 std::vector<bool>& finalized, RProgress::RProgress &pb, bool progress,
                 RayCamera* cam, adaptive_sampler& adaptive_pixel_sampler, Float percent_done,
                 hitable* world) {};
  bool terminate;
  bool interactive;
};

#endif

#endif