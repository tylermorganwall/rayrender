#include "PreviewDisplay.h"
#include <string.h>
#include <X11/Xutil.h>
#include "mathinline.h"
#include "Rcpp.h"


#ifdef RAY_HAS_X11

void PreviewDisplay::DrawImage(Rcpp::NumericMatrix& r, Rcpp::NumericMatrix& g, Rcpp::NumericMatrix& b,
                               int ns, std::vector<bool>& finalized) {
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
    // while (1) {
      XFlush(d);
      XPutImage(d,w,DefaultGC(d,s),
                img,0,0,0,0,width,height);
      // if(XCheckWindowEvent(d, w, KeyPressMask, &e)) {
      //   if (e.type == KeyPress) break;
      // }
    // }
  }
}

PreviewDisplay::PreviewDisplay(unsigned int _width, unsigned int _height) {
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
    
    // XFillRectangle(d, w, DefaultGC(d, s),
    //                0, 0, width, height);
    // XPutImage(d,w,DefaultGC(d,s),
    //           img,0,0,0,0,width,height);
    // while (1) {
    //   XNextEvent(d, &e);
    //   if (e.type == Expose) {
    //     XPutImage(d,w,DefaultGC(d,s),
    //               img,0,0,0,0,width,height);
    //     // XFillRectangle(d, w, DefaultGC(d, s),
    //     //                20, 20, 10, 10);
    //     // XDrawString(d, w, DefaultGC(d, s),
    //     //             10, 50, msg, strlen(msg));
    //   }
    //   if (e.type == KeyPress)
    //     break;
    // }
    
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