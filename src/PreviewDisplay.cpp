#ifdef RAY_HAS_X11
#include "PreviewDisplay.h"
#include <string.h>
#include "mathinline.h"
#include "Rcpp.h"
#include <X11/Xutil.h>
#include "X11/keysym.h"

void PreviewDisplay::DrawImage(Rcpp::NumericMatrix& r, Rcpp::NumericMatrix& g, Rcpp::NumericMatrix& b,
                               size_t &ns, std::vector<bool>& finalized, RProgress::RProgress &pb,
                               camera& cam, adaptive_sampler& adaptive_pixel_sampler) {
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
    KeyCode esc = XKeysymToKeycode(d, XK_Escape);
    //Movement
    KeyCode W_key = XKeysymToKeycode(d,XStringToKeysym("w"));
    KeyCode A_key = XKeysymToKeycode(d,XStringToKeysym("a"));
    KeyCode S_key = XKeysymToKeycode(d,XStringToKeysym("s"));
    KeyCode D_key = XKeysymToKeycode(d,XStringToKeysym("d"));
    KeyCode Q_key = XKeysymToKeycode(d,XStringToKeysym("q"));
    KeyCode Z_key = XKeysymToKeycode(d,XStringToKeysym("z"));
    
    //Speed control
    KeyCode E_key = XKeysymToKeycode(d,XStringToKeysym("e"));
    KeyCode C_key = XKeysymToKeycode(d,XStringToKeysym("c"));
    
    //Reset
    KeyCode R_key = XKeysymToKeycode(d,XStringToKeysym("r"));
    
    XPutImage(d,w,DefaultGC(d,s),
              img,0,0,0,0,width,height);
    if(XCheckWindowEvent(d, w, KeyPressMask, &e)) {
      if (e.type == KeyPress) {
        if (e.xkey.keycode == esc ) {
          terminate = true;
        }
        
        if(interactive) {
          if (e.xkey.keycode == W_key ) {
            cam.update_position(-speed * cam.w);
          }
          if (e.xkey.keycode == A_key ) {
            cam.update_position(-speed * cam.u);
          }
          if (e.xkey.keycode == S_key ) {
            cam.update_position(speed * cam.w);
          }
          if (e.xkey.keycode == D_key ) {
            cam.update_position(speed * cam.u);
          }
          if (e.xkey.keycode == Q_key ) {
            cam.update_position(speed * cam.v);
          }
          if (e.xkey.keycode == Z_key ) {
            cam.update_position(-speed * cam.v);
          }
          if (e.xkey.keycode == E_key ) {
            speed = 2 * speed;
          }
          if (e.xkey.keycode == C_key ) {
            speed = 0.5 * speed;
          }
          if (e.xkey.keycode == R_key ) {
            cam.reset();
          }
          ns = 0;
          std::fill(finalized.begin(), finalized.end(), false);
          std::fill(r.begin(), r.end(), 0);
          std::fill(g.begin(), g.end(), 0);
          std::fill(b.begin(), b.end(), 0);
          adaptive_pixel_sampler.reset();
          pb.update(0);
          while(XPending(d)) {
            XNextEvent(d, &e);
            if (e.xkey.keycode == esc ) {
              terminate = true;
            }
            if (e.xkey.keycode == W_key ) {
              cam.update_position(-speed * cam.w);
            }
            if (e.xkey.keycode == A_key ) {
              cam.update_position(-speed * cam.u);
            }
            if (e.xkey.keycode == S_key ) {
              cam.update_position(speed * cam.w);
            }
            if (e.xkey.keycode == D_key ) {
              cam.update_position(speed * cam.u);
            }
            if (e.xkey.keycode == Q_key ) {
              cam.update_position(speed * cam.v);
            }
            if (e.xkey.keycode == Z_key ) {
              cam.update_position(-speed * cam.v);
            }
            if (e.xkey.keycode == E_key ) {
              speed = 2 * speed;
            }
            if (e.xkey.keycode == C_key ) {
              speed = 0.5 * speed;
            }
            if (e.xkey.keycode == R_key ) {
              cam.reset();
            }
          }
        }
      }
    }
  }
}

PreviewDisplay::PreviewDisplay(unsigned int _width, unsigned int _height, 
                               bool preview, bool _interactive) {
  speed = 1.f;
  interactive = _interactive;
  terminate = false;
  if(preview) {
    d = XOpenDisplay(NULL);
  } else {
    d = nullptr;
  }
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