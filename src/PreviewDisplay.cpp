#ifdef RAY_HAS_X11
#include "PreviewDisplay.h"
#include <string.h>
#include "mathinline.h"
#include "Rcpp.h"
#include <X11/Xutil.h>
#include "X11/keysym.h"

void PreviewDisplay::DrawImage(Rcpp::NumericMatrix& r, Rcpp::NumericMatrix& g, Rcpp::NumericMatrix& b,
                               size_t &ns, std::vector<bool>& finalized, RProgress::RProgress &pb) {
  if (d) {
    for(unsigned int i = 0; i < 4*width; i += 4 ) {
      for(unsigned int j = 0; j < height; j++) {
        int samples;
        if(finalized[i/4 + width * (height-1-j)]) {
          samples = 1;
        } else {
          samples = ns+1;
        }
        data[i + 4*width*j]   = 255.999*sqrt(clamp(b(i/4,height-1-j)/samples,0,1));
        data[i + 4*width*j+1] = 255.999*sqrt(clamp(g(i/4,height-1-j)/samples,0,1));
        data[i + 4*width*j+2] = 255.999*sqrt(clamp(r(i/4,height-1-j)/samples,0,1));
      }
    }
    KeySym esc = XKeysymToKeycode(d, XK_Escape);
    // while (1) {
      XPutImage(d,w,DefaultGC(d,s),
                img,0,0,0,0,width,height);
      if(XCheckWindowEvent(d, w, KeyPressMask, &e)) {
        if (e.type == KeyPress) {
          if (e.xkey.keycode == esc ) {
            terminate = true;
          }
          ns = 0;
          std::fill(finalized.begin(), finalized.end(), false);
          std::fill(r.begin(), r.end(), 0);
          std::fill(g.begin(), g.end(), 0);
          std::fill(b.begin(), b.end(), 0);
          pb.update(0);
          while(XPending(d)) {
            XNextEvent(d, &e);
            if (e.xkey.keycode == esc ) {
              terminate = true;
            }
          }
        }
      }
    // }
  }
}

PreviewDisplay::PreviewDisplay(unsigned int _width, unsigned int _height) {
  terminate = false;
  d = XOpenDisplay(NULL);
  if (d) {
    s = DefaultScreen(d);
    Visual *visual = DefaultVisual(d,s);
    
    width = _width;
    height = _height;
    
    data = std::unique_ptr<char[]>(new char[width*height*4]);
    for(unsigned int i = 0; i < 4*width; i += 4 ) {
      for(unsigned int j = 0; j < height; j++) {
        data[i + 4*width*j]   = 0;
        data[i + 4*width*j+1] = 0;
        data[i + 4*width*j+2] = 0;
      }
    }
    
    img = XCreateImage(d,visual,
                       DefaultDepth(d, s),
                       ZPixmap,
                       0,data.get(),width,height,32,0);
    
    w = XCreateSimpleWindow(d, RootWindow(d, s), 100, 100, width, height, 1,
                            BlackPixel(d, s),  BlackPixel(d, s));
    XSelectInput(d, w, ExposureMask | KeyPressMask);
    XMapWindow(d, w);
    XFlush(d);
  }
}

PreviewDisplay::~PreviewDisplay() {
  if (d) {
    XDestroyWindow(d, w);
    XCloseDisplay(d);
    // XDestroyImage(img);
  }
}

#endif