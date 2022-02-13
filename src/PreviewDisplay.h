#ifndef PREVIEWDISPLAYH
#define PREVIEWDISPLAYH

#ifdef RAY_HAS_X11

#include <X11/Xlib.h>
#include <memory>
#include "Rcpp.h"

class PreviewDisplay {
public: 
  PreviewDisplay(unsigned int _width, unsigned int _height);
  ~PreviewDisplay();
  void DrawImage(Rcpp::NumericMatrix& r, Rcpp::NumericMatrix& g, Rcpp::NumericMatrix& b, int ns,
                 std::vector<bool>& finalized);
  Display *d;
  XImage *img;
  std::unique_ptr<char[]> data;
  Window w;
  XEvent e;
  unsigned int width;
  unsigned int height;
  int s;
};

#else

class PreviewDisplay {
public: 
  PreviewDisplay(unsigned int _width, unsigned int _height);
  ~PreviewDisplay();
  void DrawImage(Rcpp::NumericMatrix& r, Rcpp::NumericMatrix& g, Rcpp::NumericMatrix& b, int ns,
                 std::vector<bool>& finalized) {};
};

#endif

#endif