#ifndef PREVIEWDISPLAYH
#define PREVIEWDISPLAYH

#ifdef RAY_HAS_X11

#include <X11/Xlib.h>
#include <memory>
#include "Rcpp.h"
#include "RProgress.h"

class PreviewDisplay {
public: 
  PreviewDisplay(unsigned int _width, unsigned int _height);
  ~PreviewDisplay();
  void DrawImage(Rcpp::NumericMatrix& r, Rcpp::NumericMatrix& g, Rcpp::NumericMatrix& b, size_t &ns,
                 std::vector<bool>& finalized, RProgress::RProgress &pb);
  Display *d;
  XImage *img;
  std::unique_ptr<char[]> data;
  Window w;
  XEvent e;
  unsigned int width;
  unsigned int height;
  int s;
  bool terminate;
};

#else

class PreviewDisplay {
public: 
  PreviewDisplay(unsigned int _width, unsigned int _height);
  ~PreviewDisplay();
  void DrawImage(Rcpp::NumericMatrix& r, Rcpp::NumericMatrix& g, Rcpp::NumericMatrix& b, size_t &ns,
                 std::vector<bool>& finalized, RProgress::RProgress &pb) {};
  bool terminate;
};

#endif

#endif