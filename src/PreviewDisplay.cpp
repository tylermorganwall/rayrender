#ifdef RAY_HAS_X11
#include "PreviewDisplay.h"
#include <string.h>
#include "mathinline.h"
#include "Rcpp.h"
#include <X11/Xutil.h>
#include "X11/keysym.h"
#include "X11/Xatom.h"


void PreviewDisplay::DrawImage(adaptive_sampler& adaptive_pixel_sampler,
                               size_t &ns, RProgress::RProgress &pb, bool progress,
                               RayCamera* cam,  Float percent_done,
                               hitable *world, random_gen& rng) {
  if (d) {
    Rcpp::NumericMatrix &r  = adaptive_pixel_sampler.r;
    Rcpp::NumericMatrix &g  = adaptive_pixel_sampler.g;
    Rcpp::NumericMatrix &b  = adaptive_pixel_sampler.b;
    std::vector<bool>& finalized = adaptive_pixel_sampler.finalized;
    std::vector<bool>& just_finalized = adaptive_pixel_sampler.just_finalized;
    
    for(unsigned int i = 0; i < 4*width; i += 4 ) {
      for(unsigned int j = 0; j < height; j++) {
        int samples;
        Float r_col,g_col,b_col;
        if(finalized[i/4 + width * (height-1-j)]) {
          samples = 1;
          if(just_finalized[i/4 + width * (height-1-j)] && !interactive ) {
            r_col = 0;
            g_col = 1;
            b_col = 0;
          } else {
            r_col = interactive ? std::sqrt((r(i/4,height-1-j))) : std::sqrt((r(i/4,height-1-j))/4);
            g_col = interactive ? std::sqrt((g(i/4,height-1-j))) : std::sqrt((g(i/4,height-1-j))/4);
            b_col = interactive ? std::sqrt((b(i/4,height-1-j))) : std::sqrt((b(i/4,height-1-j))/4);
          }
        } else {
          samples = ns+1;
          r_col = std::sqrt((r(i/4,height-1-j))/samples);
          g_col = std::sqrt((g(i/4,height-1-j))/samples);
          b_col = std::sqrt((b(i/4,height-1-j))/samples);
        }

        data[i + 4*width*j]   = (unsigned char)(255*clamp(b_col,0,1));
        data[i + 4*width*j+1] = (unsigned char)(255*clamp(g_col,0,1));
        data[i + 4*width*j+2] = (unsigned char)(255*clamp(r_col,0,1));

        if(finalized[i/4 + width * (height-1-j)]) {
          just_finalized[i/4 + width * (height-1-j)] = false;
        }
      }
    }
    if(progress) {
      for(unsigned int i = 0; i < 4*width*percent_done; i += 4 ) {
        for(unsigned int j = 0; j < 3; j++) {
          data[i + 4*width*j]   = (unsigned char)0;
          data[i + 4*width*j+1] = (unsigned char)0;
          data[i + 4*width*j+2] = (unsigned char)255;
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
    while (XPending(d)) {
      XNextEvent(d, &e);
      if (e.type == KeyPress) {
        if (e.xkey.keycode == esc ) {
          terminate = true;
        }
        bool one_orbit = false;
        if (e.xkey.keycode == tab ) {
          orbit = !orbit;
          one_orbit = true;
        }
        
        vec3f w = cam->get_w();
        vec3f u = cam->get_u();
        vec3f v = cam->get_v();
        
        bool blanked = false;
        if(interactive) {
          if (e.xkey.keycode == W_key ) {
            cam->update_position(-speed * w * base_step, orbit);
          }
          if (e.xkey.keycode == A_key ) {
            cam->update_position(-speed * u * base_step, orbit);
          }
          if (e.xkey.keycode == S_key ) {
            cam->update_position(speed * w * base_step, orbit);
          }
          if (e.xkey.keycode == D_key ) {
            cam->update_position(speed * u * base_step, orbit);
          }
          if (e.xkey.keycode == Q_key ) {
            cam->update_position(speed * v * base_step, orbit);
          }
          if (e.xkey.keycode == Z_key ) {
            cam->update_position(-speed * v * base_step, orbit);
          }
          if (e.xkey.keycode == E_key ) {
            speed = 2 * speed;
            speed = std::fmin(speed,128);
          }
          if (e.xkey.keycode == C_key ) {
            speed = 0.5 * speed;
          }
          if (e.xkey.keycode == Down_key ) {
            cam->update_fov(speed*1.f);
          }
          if (e.xkey.keycode == Up_key ) {
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
            Float fov =  cam->get_fov();
            point3f cam_direction = -cam->get_w();
            Float fd = cam->get_focal_distance();
            
            if(fov >= 0) {
              Rprintf("\nLookfrom: c(%.2f, %.2f, %.2f) LookAt: c(%.2f, %.2f, %.2f) FOV: %.1f Aperture: %0.3f Focal Dist: %.1f Step Multiplier: %.2f",
                      origin.x(), origin.y(), origin.z(), 
                      origin.x()+cam_direction.x()*fd, origin.y()+cam_direction.y()*fd, origin.z()+cam_direction.z()*fd, 
                      fov, 
                      cam->get_aperture(), fd, speed);
            }  else {
              Rprintf("\nLookfrom: c(%.2f, %.2f, %.2f) LookAt: c(%.2f, %.2f, %.2f)  Focal Dist: %0.3f Step Multiplier: %.2f",
                      origin.x(), origin.y(), origin.z(), 
                      origin.x()+cam_direction.x()*fd, origin.y()+cam_direction.y()*fd, origin.z()+cam_direction.z()*fd, 
                      fd, speed);
            }
          } else {
            if(!blanked && !terminate) {
              blanked = true;
              ns = 0;
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
            
            w = cam->get_w();
            u = cam->get_u();
            v = cam->get_v();
            
            if (e.xkey.keycode == W_key ) {
              cam->update_position(-speed * w * base_step, orbit);
            }
            if (e.xkey.keycode == A_key ) {
              cam->update_position(-speed * u * base_step, orbit);
            }
            if (e.xkey.keycode == S_key ) {
              cam->update_position(speed * w * base_step, orbit);
            }
            if (e.xkey.keycode == D_key ) {
              cam->update_position(speed * u * base_step, orbit);
            }
            if (e.xkey.keycode == Q_key ) {
              cam->update_position(speed * v * base_step, orbit);
            }
            if (e.xkey.keycode == Z_key ) {
              cam->update_position(-speed * v * base_step, orbit);
            }
            if (e.xkey.keycode == E_key ) {
              speed = 2 * speed;
            }
            if (e.xkey.keycode == C_key ) {
              speed = 0.5 * speed;
              speed = std::fmin(speed,128);
            }
            if (e.xkey.keycode == Down_key ) {
              cam->update_fov(speed*1.f);
            }
            if (e.xkey.keycode == Up_key ) {
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
              Float fov =  cam->get_fov();
              point3f cam_direction = -cam->get_w();
              Float fd = cam->get_focal_distance();
              
              if(fov >= 0) {
                Rprintf("\nLookfrom: c(%.2f, %.2f, %.2f) LookAt: c(%.2f, %.2f, %.2f) FOV: %.1f Aperture: %0.3f Focal Dist: %.1f Step Multiplier: %.2f",
                        origin.x(), origin.y(), origin.z(), 
                        origin.x()+cam_direction.x()*fd, origin.y()+cam_direction.y()*fd, origin.z()+cam_direction.z()*fd, 
                        fov, 
                        cam->get_aperture(), fd, speed);
              }  else {
                Rprintf("\nLookfrom: c(%.2f, %.2f, %.2f) LookAt: c(%.2f, %.2f, %.2f)  Focal Dist: %0.3f Step Multiplier: %.2f",
                        origin.x(), origin.y(), origin.z(), 
                        origin.x()+cam_direction.x()*fd, origin.y()+cam_direction.y()*fd, origin.z()+cam_direction.z()*fd, 
                        fd, speed);
              }
            } else {
              if(!blanked && !terminate) {
                blanked = true;
                ns = 0;
                adaptive_pixel_sampler.reset();
                if(progress && !interactive) {
                  pb.update(0);
                }
              }
            }
          }
        }
      } else if (e.type == ButtonPress) {
        bool left = e.xbutton.button == Button1;
        bool right = e.xbutton.button == Button3;
        if(!left && !right) {
          return;
        }
        Float x = e.xbutton.x;
        Float y = e.xbutton.y;
        Float fov = cam->get_fov();
        Float u = (Float(x)) / Float(width);
        Float v = (Float(y)) / Float(height);
        vec3f dir;
        hit_record hrec;
        if(fov < 0) {
          CameraSample samp({1-u,v},point2f(0.5,0.5), 0.5);
          ray r2;
          cam->GenerateRay(samp,&r2);
          if(world->hit(r2, 0.001, FLT_MAX, hrec, rng)) {
            if( hrec.shape->GetName() != "EnvironmentLight") {
              dir = -hrec.p;
            }  else {
              return;
            }
          }
        } else if (fov > 0) {
          ray r2 = cam->get_ray(u,1-v, point3f(0.5),
                               0.5f);
          if(left) {
            world->hit(r2, 0.001, FLT_MAX, hrec, rng);
          }
          dir = r2.direction();
        } else {
          ray r2 = cam->get_ray(u,1-v, point3f(0.5),
                               0.5f);
          if(world->hit(r2, 0.001, FLT_MAX, hrec, rng)) {
            if( hrec.shape->GetName() != "EnvironmentLight") {
              dir = -(cam->get_origin()-hrec.p);
            } else {
              dir = -cam->get_w();
            }
          } else {
            dir = -cam->get_w();
          }
        }
        if(left && !right) {
          Float current_fd = cam->get_focal_distance();
          Float new_fd = (hrec.p-cam->get_origin()).length();
          cam->update_focal_distance(new_fd- current_fd);
          cam->update_lookat(hrec.p);
        }
        cam->update_look_direction(-dir);
        ns = 0;
        adaptive_pixel_sampler.reset();
        if(progress && !interactive) {
          pb.update(0);
        }
      } else if (e.type == ClientMessage) {
        terminate = true;
      } 
    }
  }
}

PreviewDisplay::PreviewDisplay(unsigned int _width, unsigned int _height, 
                               bool preview, bool _interactive,
                               Float initial_lookat_distance) {
  speed = 1.f;
  interactive = _interactive;
  orbit = true;
  terminate = false;
  base_step = initial_lookat_distance/20;
  if(preview) {
    d = XOpenDisplay(NULL);
  } else {
    d = nullptr;
  }
  if (d) {
    s = DefaultScreen(d);
    XVisualInfo vinfo;
    if (!XMatchVisualInfo(d, s, 24, TrueColor, &vinfo)) {
      Rprintf("No X11 `visual` object found matching display requirements (24 bit depth and True Color)");
      d = nullptr;
      XCloseDisplay(d);
      return;
    }
    Visual *visual = vinfo.visual;
    
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
    XSelectInput(d, w, ExposureMask | KeyPressMask | ButtonPress);
    XMapWindow(d, w);
    Atom WM_DELETE_WINDOW = XInternAtom(d, "WM_DELETE_WINDOW", False); 
    XSetWMProtocols(d, w, &WM_DELETE_WINDOW, 1);
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