#ifdef RAY_HAS_X11
#include "PreviewDisplay.h"
#include <string.h>
#include "mathinline.h"
#include "Rcpp.h"
#include <X11/Xutil.h>
#include "X11/keysym.h"


void PreviewDisplay::DrawImage(Rcpp::NumericMatrix& r, Rcpp::NumericMatrix& g, Rcpp::NumericMatrix& b,
                               size_t &ns, std::vector<bool>& finalized, RProgress::RProgress &pb, bool progress,
                               RayCamera* cam, adaptive_sampler& adaptive_pixel_sampler, Float percent_done) {
  if (d) {
    for(unsigned int i = 0; i < 4*width; i += 4 ) {
      for(unsigned int j = 0; j < height; j++) {
        int samples;
        if(finalized[i/4 + width * (height-1-j)]) {
          samples = 1;
        } else {
          samples = ns+1;
        }
        data[i + 4*width*j]   = 255*sqrt(clamp(b(i/4,height-1-j)*cam->get_iso()/samples,0,1));
        data[i + 4*width*j+1] = 255*sqrt(clamp(g(i/4,height-1-j)*cam->get_iso()/samples,0,1));
        data[i + 4*width*j+2] = 255*sqrt(clamp(r(i/4,height-1-j)*cam->get_iso()/samples,0,1));
      }
    }
    if(progress) {
      for(unsigned int i = 0; i < 4*width*percent_done; i += 4 ) {
        for(unsigned int j = 0; j < 2; j++) {
          data[i + 4*width*j]   = 0;
          data[i + 4*width*j+1] = 0;
          data[i + 4*width*j+2] = 255;
        }
      }
    }
    KeyCode tab = XKeysymToKeycode(d, XK_Tab);
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
    
    //Fov control
    KeyCode Up_key = XKeysymToKeycode(d,XK_Up);
    KeyCode Down_key = XKeysymToKeycode(d,XK_Down);
    
    //Aperture control
    KeyCode Left_key = XKeysymToKeycode(d,XK_Left);
    KeyCode Right_key = XKeysymToKeycode(d,XK_Right);
    
    
    //Focus Distance control
    KeyCode One_key = XKeysymToKeycode(d,XStringToKeysym("1"));
    KeyCode Two_key = XKeysymToKeycode(d,XStringToKeysym("2"));
    
    //Reset
    KeyCode R_key = XKeysymToKeycode(d,XStringToKeysym("r"));
    
    //Print Position
    KeyCode P_key = XKeysymToKeycode(d,XStringToKeysym("p"));
    
    
    XPutImage(d,w,DefaultGC(d,s),
              img,0,0,0,0,width,height);
    if(XCheckWindowEvent(d, w, KeyPressMask, &e)) {
      if (e.type == KeyPress) {
        if (e.xkey.keycode == esc ) {
          terminate = true;
        }
        bool one_orbit = false;
        if (e.xkey.keycode == tab ) {
          orbit = !orbit;
          one_orbit = true;
        }
        
        vec3f w = orbit ? cam->get_w() : vec3f(0,0,1);
        vec3f u = orbit ? cam->get_u() : vec3f(1,0,0);
        vec3f v = orbit ? cam->get_v() : vec3f(0,1,0);
        bool blanked = false;
        if(interactive) {
          if (e.xkey.keycode == W_key ) {
            cam->update_position(-speed * w);
          }
          if (e.xkey.keycode == A_key ) {
            cam->update_position(-speed * u);
          }
          if (e.xkey.keycode == S_key ) {
            cam->update_position(speed * w);
          }
          if (e.xkey.keycode == D_key ) {
            cam->update_position(speed * u);
          }
          if (e.xkey.keycode == Q_key ) {
            cam->update_position(speed * v);
          }
          if (e.xkey.keycode == Z_key ) {
            cam->update_position(-speed * v);
          }
          if (e.xkey.keycode == E_key ) {
            speed = 2 * speed;
            speed = std::fmin(speed,128);
          }
          if (e.xkey.keycode == C_key ) {
            speed = 0.5 * speed;
          }
          if (e.xkey.keycode == Up_key ) {
            cam->update_fov(speed*1.f);
          }
          if (e.xkey.keycode == Down_key ) {
            cam->update_fov(speed*-1.f);
          }
          if (e.xkey.keycode == Left_key ) {
            cam->update_aperture(speed*-0.1f);
          }
          if (e.xkey.keycode == Right_key ) {
            cam->update_aperture(speed*0.1f);
          }
          if (e.xkey.keycode == One_key ) {
            cam->update_focal_distance(speed*-1.f);
          }
          if (e.xkey.keycode == Two_key ) {
            cam->update_focal_distance(speed*1.f);
          }
          if (e.xkey.keycode == R_key ) {
            cam->reset();
            speed = 1;
          }
          if (e.xkey.keycode == P_key ) {
            point3f origin = cam->get_origin();
            Rprintf("\nCamera Position: c(%.2f, %.2f, %.2f) FOV: %.1f Aperture: %0.1f Focal Dist: %.1f Step Multiplier: %.2f",
                    origin.x(), origin.y(), origin.z(), cam->get_fov(), 
                    cam->get_aperture(), cam->get_focal_distance(), speed);
          } else {
            if(!blanked) {
              blanked = true;
              ns = 0;
              std::fill(finalized.begin(), finalized.end(), false);
              std::fill(r.begin(), r.end(), 0);
              std::fill(g.begin(), g.end(), 0);
              std::fill(b.begin(), b.end(), 0);
              adaptive_pixel_sampler.reset();
              if(progress && !interactive) {
                pb.update(0);
              }
            }
          }
          while(XPending(d)) {
            XNextEvent(d, &e);
            if (e.xkey.keycode == esc ) {
              terminate = true;
            }
            if (e.xkey.keycode == tab && !one_orbit) {
              orbit = !orbit;
              one_orbit = true;
            }
            
            w = orbit ? cam->get_w() : vec3f(0,0,1);
            u = orbit ? cam->get_u() : vec3f(1,0,0);
            v = orbit ? cam->get_v() : vec3f(0,1,0);
            
            if (e.xkey.keycode == W_key ) {
              cam->update_position(-speed * w);
            }
            if (e.xkey.keycode == A_key ) {
              cam->update_position(-speed * u);
            }
            if (e.xkey.keycode == S_key ) {
              cam->update_position(speed * w);
            }
            if (e.xkey.keycode == D_key ) {
              cam->update_position(speed * u);
            }
            if (e.xkey.keycode == Q_key ) {
              cam->update_position(speed * v);
            }
            if (e.xkey.keycode == Z_key ) {
              cam->update_position(-speed * v);
            }
            if (e.xkey.keycode == E_key ) {
              speed = 2 * speed;
            }
            if (e.xkey.keycode == C_key ) {
              speed = 0.5 * speed;
              speed = std::fmin(speed,128);
            }
            if (e.xkey.keycode == Up_key ) {
              cam->update_fov(speed*1.f);
            }
            if (e.xkey.keycode == Down_key ) {
              cam->update_fov(speed*-1.f);
            }
            if (e.xkey.keycode == Left_key ) {
              cam->update_aperture(speed*-0.1f);
            }
            if (e.xkey.keycode == Right_key ) {
              cam->update_aperture(speed*0.1f);
            }
            if (e.xkey.keycode == One_key ) {
              cam->update_focal_distance(speed*-1.f);
            }
            if (e.xkey.keycode == Two_key ) {
              cam->update_focal_distance(speed*1.f);
            }
            if (e.xkey.keycode == R_key ) {
              cam->reset();
              speed = 1;
            }
            if (e.xkey.keycode == P_key ) {
              point3f origin = cam->get_origin();
              Rprintf("\nCamera Position: c(%.2f, %.2f, %.2f) FOV: %.1f Aperture: %0.1f Focal Dist: %.1f Step Multiplier: %.2f",
                      origin.x(), origin.y(), origin.z(), cam->get_fov(), 
                      cam->get_aperture(), cam->get_focal_distance(), speed);
            } else {
              if(!blanked) {
                blanked = true;
                ns = 0;
                std::fill(finalized.begin(), finalized.end(), false);
                std::fill(r.begin(), r.end(), 0);
                std::fill(g.begin(), g.end(), 0);
                std::fill(b.begin(), b.end(), 0);
                adaptive_pixel_sampler.reset();
                if(progress && !interactive) {
                  pb.update(0);
                }
              }
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
  orbit = true;
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