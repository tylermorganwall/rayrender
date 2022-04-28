#include "PreviewDisplay.h"
#include "mathinline.h"
#include "Rcpp.h"


#ifdef RAY_HAS_X11

//Had to undefine Xlib.h's Status because RcppThread also defines Status 
#define Status int
#include <string.h>
#include <X11/Xutil.h>
#include "X11/keysym.h"
#include "X11/Xatom.h"
static Float env_y_angle;;
#undef Status

#endif

#ifdef RAY_WINDOWS

#include <windows.h>
#include <winuser.h>
#include "float.h"
#include <wingdi.h>
#include <windowsx.h>
#include "RProgress.h"

static unsigned int width;
static unsigned int height;
static bool term;

static std::vector<Float> rgb;
static RayCamera* cam_w;
static Float speed;
static bool preview;
static bool orbit;
static Float base_step;
static bool blanked;
static bool interactive_w;
static adaptive_sampler* aps;
static adaptive_sampler* aps_small;
static size_t* ns_w;
static hitable* world_w;
static random_gen* rng_w;
static RProgress::RProgress* pb_w;
static bool progress_w;
static Float env_y_angle;;
static Transform* EnvWorldToObject_w;
static Transform* EnvObjectToWorld_w;
static Transform Start_EnvWorldToObject_w;
static Transform Start_EnvObjectToWorld_w;
static std::vector<Rcpp::List>* Keyframes_w;
static bool* write_fast_output_w;



LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam);

#endif


void PreviewDisplay::DrawImage(adaptive_sampler& adaptive_pixel_sampler,
                               adaptive_sampler& adaptive_pixel_sampler_small,
                               size_t &ns, RProgress::RProgress &pb, bool progress,
                               Float percent_done,
                               hitable *world, random_gen& rng) {
#ifdef RAY_HAS_X11
  if (d) {
    RayMatrix &r  = adaptive_pixel_sampler.r;
    RayMatrix &g  = adaptive_pixel_sampler.g;
    RayMatrix &b  = adaptive_pixel_sampler.b;
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
    
    
    //Environment Rotate control
    KeyCode Three_key = XKeysymToKeycode(d,XStringToKeysym("3"));
    KeyCode Four_key = XKeysymToKeycode(d,XStringToKeysym("4"));
    
    //Reset
    KeyCode R_key = XKeysymToKeycode(d,XStringToKeysym("r"));
    
    //Print Position
    KeyCode P_key = XKeysymToKeycode(d,XStringToKeysym("p"));
    
    //Save Keyframe
    KeyCode K_key = XKeysymToKeycode(d,XStringToKeysym("k"));
    
    //Move to Last Keyframe
    KeyCode L_key = XKeysymToKeycode(d,XStringToKeysym("l"));
    
    //Fast Movement Key
    KeyCode F_key = XKeysymToKeycode(d,XStringToKeysym("f"));
    
    
    XPutImage(d,w,DefaultGC(d,s),
              img,0,0,0,0,width,height);
    while (XPending(d)) {
      XNextEvent(d, &e);
      if (e.type == KeyPress) {
        if (e.xkey.keycode == esc ) {
          terminate = true;
          break;
        }
        if(interactive) {
          vec3f w = cam->get_w();
          vec3f u = cam->get_u();
          vec3f v = cam->get_v();
        
          bool blanked = false;
          bool one_orbit = false;
          bool one_fast = false;
          
          if (e.xkey.keycode == tab ) {
            orbit = !orbit;
            one_orbit = true;
          }
          if (e.xkey.keycode == F_key ) {
            write_fast_output = !write_fast_output;
            one_fast  = true;
          }
          if (e.xkey.keycode == W_key ) {
            vec3f step = -speed * w * base_step;
            if(orbit) {
              Float dist_to_orbit = (cam->get_origin() - cam->get_lookat()).length();
              if(dist_to_orbit <= base_step * speed) {
                Rprintf("Moving forward will overstep orbit point, stopping (decrease step size to move closer).\n");
                step = vec3f(0);
              }
            } 
            cam->update_position(step, orbit, false);
          }
          if (e.xkey.keycode == A_key ) {
            cam->update_position(-speed * u * base_step, orbit);
          }
          if (e.xkey.keycode == S_key ) {
            cam->update_position(speed * w * base_step, orbit, false);
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
            Rprintf("Step Multiplier: %.3f\n", speed);
          }
          if (e.xkey.keycode == C_key ) {
            speed = 0.5 * speed;
            Rprintf("Step Multiplier: %.3f\n", speed);
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
          if (e.xkey.keycode == Three_key ) {
            (*EnvObjectToWorld) =  RotateY(speed*8) * (*EnvObjectToWorld);
            (*EnvWorldToObject) =  RotateY(-speed*8) * (*EnvWorldToObject);
            env_y_angle -= speed*8;
          }
          if (e.xkey.keycode == Four_key ) {
            (*EnvObjectToWorld) =  RotateY(-speed*8) * (*EnvObjectToWorld);
            (*EnvWorldToObject) =  RotateY(speed*8) * (*EnvWorldToObject);
            env_y_angle += speed*8;
          }
          if (e.xkey.keycode == R_key ) {
            cam->reset();
            speed = 1;
            (*EnvWorldToObject) = Start_EnvWorldToObject;
            (*EnvObjectToWorld) = Start_EnvObjectToWorld;
          }
          if (e.xkey.keycode == L_key) {
            if(Keyframes.size() > 0) {
              Rcpp::List LastKeyframe = Keyframes.at(Keyframes.size()-1);
              point3f last_pos = point3f(Rcpp::as<Float>(LastKeyframe["x"]),
                                       Rcpp::as<Float>(LastKeyframe["y"]),
                                       Rcpp::as<Float>(LastKeyframe["z"]));
              point3f last_lookat = point3f(Rcpp::as<Float>(LastKeyframe["dx"]),
                                            Rcpp::as<Float>(LastKeyframe["dy"]),
                                            Rcpp::as<Float>(LastKeyframe["dz"]));
              Float last_aperture = Rcpp::as<Float>(LastKeyframe["aperture"]);
              Float last_fov = Rcpp::as<Float>(LastKeyframe["fov"]);
              Float last_focal = Rcpp::as<Float>(LastKeyframe["focal"]);
              vec2f last_ortho = vec2f(Rcpp::as<Float>(LastKeyframe["orthox"]),
                                       Rcpp::as<Float>(LastKeyframe["orthoy"]));
              cam->update_focal_absolute(last_focal);
              cam->update_position_absolute(last_pos);
              cam->update_lookat(last_lookat);
              cam->update_aperture_absolute(last_aperture);
              cam->update_fov_absolute(last_fov);
              cam->update_ortho_absolute(last_ortho);
            } else {
              Rprintf("Can't reset to last keyframe: No keyframes have been saved. Use the R key to reset camera.");
            }
          }
          if (e.xkey.keycode == P_key || e.xkey.keycode == K_key ) {
            point3f origin = cam->get_origin();
            Float fov =  cam->get_fov();
            point3f cam_direction = -cam->get_w();
            Float cam_aperture = cam->get_aperture();
            Float fd = cam->get_focal_distance();
            vec3f cam_up = cam->get_up();
            point2f ortho = cam->get_ortho();
            point3f cam_lookat = cam->get_lookat();
            
            if(e.xkey.keycode == K_key) {
              if(fov > 0) {
                Keyframes.push_back(Rcpp::List::create(Named("x")  = origin.x(),
                                                       Named("y")  = origin.y(),
                                                       Named("z")  = origin.z(),
                                                       Named("dx") = origin.x()+cam_direction.x()*fd,
                                                       Named("dy") = origin.y()+cam_direction.y()*fd,
                                                       Named("dz") = origin.z()+cam_direction.z()*fd,
                                                       Named("aperture") = cam_aperture,
                                                       Named("fov") = fov,
                                                       Named("focal") = fd,
                                                       Named("orthox") = ortho.x(),
                                                       Named("orthoy") = ortho.y(),
                                                       Named("upx")  = cam_up.x(),
                                                       Named("upy")  = cam_up.y(),
                                                       Named("upz")  = cam_up.z()));
              } else if (fov < 0) {
                Keyframes.push_back(Rcpp::List::create(Named("x")  = origin.x(),
                                                       Named("y")  = origin.y(),
                                                       Named("z")  = origin.z(),
                                                       Named("dx") = origin.x()+cam_direction.x()*fd,
                                                       Named("dy") = origin.y()+cam_direction.y()*fd,
                                                       Named("dz") = origin.z()+cam_direction.z()*fd,
                                                       Named("aperture") = 0,
                                                       Named("fov") = 0,
                                                       Named("focal") = fd,
                                                       Named("orthox") = ortho.x(),
                                                       Named("orthoy") = ortho.y(),
                                                       Named("upx")  = cam_up.x(),
                                                       Named("upy")  = cam_up.y(),
                                                       Named("upz")  = cam_up.z()));
              } else {
                Keyframes.push_back(Rcpp::List::create(Named("x")  = origin.x(),
                                                       Named("y")  = origin.y(),
                                                       Named("z")  = origin.z(),
                                                       Named("dx") = cam_lookat.x(),
                                                       Named("dy") = cam_lookat.y(),
                                                       Named("dz") = cam_lookat.z(),
                                                       Named("aperture") = 0,
                                                       Named("fov") = 0,
                                                       Named("focal") = fd,
                                                       Named("orthox") = ortho.x(),
                                                       Named("orthoy") = ortho.y(),
                                                       Named("upx")  = cam_up.x(),
                                                       Named("upy")  = cam_up.y(),
                                                       Named("upz")  = cam_up.z()));
              }
            }
            if(fov > 0) {
              Rprintf("Lookfrom: c(%.2f, %.2f, %.2f) LookAt: c(%.2f, %.2f, %.2f) FOV: %.1f Aperture: %0.3f Focal Dist: %0.3f Env Rotation: %.2f\n",
                      origin.x(), origin.y(), origin.z(), 
                      origin.x()+cam_direction.x()*fd, origin.y()+cam_direction.y()*fd, origin.z()+cam_direction.z()*fd, 
                      fov, 
                      cam_aperture, fd, env_y_angle);
            } else if (fov < 0) {
              Rprintf("Lookfrom: c(%.2f, %.2f, %.2f) LookAt: c(%.2f, %.2f, %.2f) Focal Dist: %0.3f Env Rotation: %.2f\n",
                      origin.x(), origin.y(), origin.z(), 
                      origin.x()+cam_direction.x()*fd, origin.y()+cam_direction.y()*fd, origin.z()+cam_direction.z()*fd, 
                      fd, env_y_angle);
            } else {
              Rprintf("Lookfrom: c(%.2f, %.2f, %.2f) LookAt: c(%.2f, %.2f, %.2f) Focal Dist: %0.3f Env Rotation: %.2f\n",
                      origin.x(), origin.y(), origin.z(), 
                      cam_lookat.x(), cam_lookat.y(), cam_lookat.z(), 
                      fd, env_y_angle);
            }
            
          } else {
            if(!blanked && !terminate && e.xkey.keycode != C_key && e.xkey.keycode != E_key && e.xkey.keycode != tab) {
              blanked = true;
              ns = 0;
              adaptive_pixel_sampler.reset();
              adaptive_pixel_sampler_small.reset();
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
            
            if (e.xkey.keycode == F_key && !one_fast) {
              write_fast_output = !write_fast_output;
              one_fast  = true;
            }
            
            w = cam->get_w();
            u = cam->get_u();
            v = cam->get_v();
            
            if (e.xkey.keycode == W_key ) {
              vec3f step = -speed * w * base_step;
              if(orbit) {
                Float dist_to_orbit = (cam->get_origin() - cam->get_lookat()).length();
                if(dist_to_orbit <= base_step * speed) {
                  Rprintf("Moving forward will overstep orbit point, stopping (decrease step size to move closer).\n");
                  step = vec3f(0);
                }
              } 
              cam->update_position(step, orbit, false);
            }
            if (e.xkey.keycode == A_key ) {
              cam->update_position(-speed * u * base_step, orbit);
            }
            if (e.xkey.keycode == S_key ) {
              cam->update_position(speed * w * base_step, orbit, false);
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
              Rprintf("Step Multiplier: %.3f\n", speed);
              
            }
            if (e.xkey.keycode == C_key ) {
              speed = 0.5 * speed;
              speed = std::fmin(speed,128);
              Rprintf("Step Multiplier: %.3f\n", speed);
              
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
            if (e.xkey.keycode == Three_key ) {
              (*EnvObjectToWorld) =  RotateY(speed*8) * (*EnvObjectToWorld);
              (*EnvWorldToObject) =  RotateY(-speed*8) * (*EnvWorldToObject);
              env_y_angle -= speed*8;
              
            }
            if (e.xkey.keycode == Four_key ) {
              (*EnvObjectToWorld) =  RotateY(-speed*8) * (*EnvObjectToWorld);
              (*EnvWorldToObject) =  RotateY(speed*8) * (*EnvWorldToObject);
              env_y_angle += speed*8;
              
            }
            if (e.xkey.keycode == R_key ) {
              cam->reset();
              speed = 1;
              (*EnvWorldToObject) = Start_EnvWorldToObject;
              (*EnvObjectToWorld) = Start_EnvObjectToWorld;
            }
            if (e.xkey.keycode == P_key ) {
              point3f origin = cam->get_origin();
              Float fov =  cam->get_fov();
              point3f cam_direction = -cam->get_w();
              Float fd = cam->get_focal_distance();
              point3f cam_lookat = cam->get_lookat();
              Float cam_aperture = cam->get_aperture();
              if(fov > 0) {
                Rprintf("Lookfrom: c(%.2f, %.2f, %.2f) LookAt: c(%.2f, %.2f, %.2f) FOV: %.1f Aperture: %0.3f Focal Dist: %0.3f Env Rotation: %.2f\n",
                        origin.x(), origin.y(), origin.z(), 
                        origin.x()+cam_direction.x()*fd, origin.y()+cam_direction.y()*fd, origin.z()+cam_direction.z()*fd, 
                        fov, 
                        cam_aperture, fd, env_y_angle);
              } else if (fov < 0) {
                Rprintf("Lookfrom: c(%.2f, %.2f, %.2f) LookAt: c(%.2f, %.2f, %.2f) Focal Dist: %0.3f Env Rotation: %.2f\n",
                        origin.x(), origin.y(), origin.z(), 
                        origin.x()+cam_direction.x()*fd, origin.y()+cam_direction.y()*fd, origin.z()+cam_direction.z()*fd, 
                        fd, env_y_angle);
              } else {
                Rprintf("Lookfrom: c(%.2f, %.2f, %.2f) LookAt: c(%.2f, %.2f, %.2f) Focal Dist: %0.3f Env Rotation: %.2f\n",
                        origin.x(), origin.y(), origin.z(), 
                        cam_lookat.x(), cam_lookat.y(), cam_lookat.z(), 
                        fd, env_y_angle);
              }
            } else {
              if(!blanked && !terminate && e.xkey.keycode != C_key && e.xkey.keycode != E_key && e.xkey.keycode != tab) {
                blanked = true;
                ns = 0;
                adaptive_pixel_sampler.reset();
                adaptive_pixel_sampler_small.reset();
                if(progress && !interactive) {
                  pb.update(0);
                }
              }
            }
          }
        }
      } else if (e.type == ButtonPress) {
        if(interactive) {
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
                left = false;
                right = true;
                dir = -hrec.p;
              }
            }
          } else if (fov > 0) {
            ray r2 = cam->get_ray(u,1-v, point3f(0),
                                 0.5f);
            if(left) {
              world->hit(r2, 0.001, FLT_MAX, hrec, rng);
              if( hrec.shape->GetName() == "EnvironmentLight") {
                right = true;
              }
            }
            dir = r2.direction();
          } else {
            ray r2 = cam->get_ray(u,1-v, point3f(0.5),
                                 0.5f);
            if(world->hit(r2, 0.001, FLT_MAX, hrec, rng)) {
              if( hrec.shape->GetName() != "EnvironmentLight") {
                dir = -(cam->get_origin()-hrec.p);
              } else {
                Rprintf("Clicking on the environment light while using a orthographic camera does not change the view.\n");
                dir = -cam->get_w();
              }
            } else {
              dir = -cam->get_w();
            }
          }
          if(left && !right) {
            if(fov != 0 && fov != 360) {
              Float current_fd = cam->get_focal_distance();
              Float new_fd = (hrec.p-cam->get_origin()).length();
              cam->update_focal_distance(new_fd- current_fd);
            }
            cam->update_lookat(hrec.p);
          }
          cam->update_look_direction(-dir);
          ns = 0;
          adaptive_pixel_sampler.reset();
          adaptive_pixel_sampler_small.reset();
          
          if(progress && !interactive) {
            pb.update(0);
          }
        }
      } else if (e.type == ClientMessage) {
        terminate = true;
      } 
    }
  }
#endif
#ifdef RAY_WINDOWS
  if(hwnd) {
    aps = &adaptive_pixel_sampler;
    aps_small = &adaptive_pixel_sampler_small;
    ns_w = &ns;
    pb_w = &pb;
    progress_w = progress;
    interactive_w = interactive;
    RayMatrix &r  = adaptive_pixel_sampler.r;
    RayMatrix &g  = adaptive_pixel_sampler.g;
    RayMatrix &b  = adaptive_pixel_sampler.b;
    EnvWorldToObject_w = EnvWorldToObject;
    EnvObjectToWorld_w = EnvObjectToWorld;
    Start_EnvWorldToObject_w = Start_EnvWorldToObject;
    Start_EnvObjectToWorld_w = Start_EnvObjectToWorld;
    std::vector<bool>& finalized = adaptive_pixel_sampler.finalized;
    std::vector<bool>& just_finalized = adaptive_pixel_sampler.just_finalized;
    write_fast_output_w = &write_fast_output;
    world_w = world;
    rng_w = &rng;
    height = (unsigned int)r.cols();
    width = (unsigned int)r.rows();
    rgb.resize(width*height*3);
    for(unsigned int i = 0; i < width*3; i += 3) {
      for(unsigned int j = 0; j < height; j++) {
        Float samples;
        Float r_col,g_col,b_col;
        if(finalized[i/3 + width * (height-1-j)]) {
          if(just_finalized[i/3 + width * (height-1-j)] && !interactive ) {
            r_col = 0.f;
            g_col = 1.f;
            b_col = 0.f;
          } else {
            r_col = interactive ? std::sqrt((r(i/3,height-1-j))) : std::sqrt((r(i/3,height-1-j))/4.f);
            g_col = interactive ? std::sqrt((g(i/3,height-1-j))) : std::sqrt((g(i/3,height-1-j))/4.f);
            b_col = interactive ? std::sqrt((b(i/3,height-1-j))) : std::sqrt((b(i/3,height-1-j))/4.f);
          }
        } else {
          samples = (Float)ns+1.f;
          r_col = std::sqrt((r(i/3,height-1-j))/samples);
          g_col = std::sqrt((g(i/3,height-1-j))/samples);
          b_col = std::sqrt((b(i/3,height-1-j))/samples);
        }
        rgb[i+3*width*j]   = clamp(r_col,0.f,1.f);
        rgb[i+3*width*j+1] = clamp(g_col,0.f,1.f);
        rgb[i+3*width*j+2] = clamp(b_col,0.f,1.f);
        
        if(finalized[i/3 + width * (height-1-j)]) {
          just_finalized[i/3 + width * (height-1-j)] = false;
        }
      }
    }
    blanked = false;
    
    if(progress) {
      for(unsigned int i = 0; i < 3*width*percent_done; i += 3 ) {
        for(unsigned int j = 0; j < 3; j++) {
          rgb[i + 3*width*j]   = 1.f;
          rgb[i + 3*width*j+1] = 0.f;
          rgb[i + 3*width*j+2] = 0.f;
        }
      }
    }
    
    InvalidateRect(hwnd, NULL, 0);
    while (PeekMessage (&msg, NULL, 0, 0, PM_REMOVE) > 0) {
      TranslateMessage(&msg);
      DispatchMessage(&msg); 
    }
    terminate = term;
  }
#endif
}

PreviewDisplay::PreviewDisplay(unsigned int _width, unsigned int _height, 
                               bool preview, bool _interactive,
                               Float initial_lookat_distance, RayCamera* _cam,
                               Transform* _EnvObjectToWorld, Transform* _EnvWorldToObject) :
  preview(preview), EnvObjectToWorld(_EnvObjectToWorld), EnvWorldToObject(_EnvWorldToObject),
  Start_EnvObjectToWorld(*_EnvObjectToWorld), Start_EnvWorldToObject(*_EnvWorldToObject) {
  Keyframes.clear();
  write_fast_output = false;
#ifdef RAY_HAS_X11
  speed = 1.f;
  interactive = _interactive;
  orbit = true;
  terminate = false;
  env_y_angle = 0;
  base_step = initial_lookat_distance/20;
  cam = _cam;
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
#endif
#ifdef RAY_WINDOWS
  speed = 1.f;
  interactive = _interactive;
  orbit = true;
  terminate = false;
  term = false;
  env_y_angle = 0;
  Keyframes_w = &Keyframes;
  if(preview) {
    width = _width;
    height = _height;
    base_step = initial_lookat_distance/20;
    rgb.resize(width*height*3);
    cam_w = _cam;
    hInstance = (HINSTANCE)GetModuleHandle(NULL);
    // Register the window class.
    const wchar_t CLASS_NAME[]  = L"Rayrender";
  
    wc = { };
  
    wc.lpfnWndProc   = WindowProc;
    wc.hInstance     = hInstance;
    wc.lpszClassName = CLASS_NAME;
  
    RegisterClass(&wc);
  
    // Create the window.
    RECT rect = {0, 0, (long int)width, (long int)height};
    AdjustWindowRect(&rect, WS_THICKFRAME | WS_VISIBLE | WS_SYSMENU, true);
    
    hwnd = CreateWindowEx(
      0,                              // Optional window styles.
      CLASS_NAME,                     // Window class
      L"Rayrender",    // Window text
      WS_THICKFRAME | WS_VISIBLE | WS_SYSMENU ,            // Window style
  
      // Size and position
      0, 0, rect.right - rect.left, rect.bottom - rect.top, 
  
      NULL,       // Parent window
      NULL,       // Menu
      hInstance,  // Instance handle
      NULL        // Additional application data
    );
  
  
    if (hwnd == NULL) {
      throw std::runtime_error("Can't open window");
    }
    ShowWindow(hwnd, SW_SHOW);
    // SetForegroundWindow(hwnd)
    // BringWindowToTop(hwnd);
  } else {
    hwnd = nullptr;
  }
#endif
}

PreviewDisplay::~PreviewDisplay() {
#ifdef RAY_HAS_X11
  if (d) {
    XDestroyWindow(d, w);
    XCloseDisplay(d);
  }
#endif
#ifdef RAY_WINDOWS
  if (hwnd != NULL) {
    DestroyWindow(hwnd);
  }
  rgb.resize(0);
#endif
}

#ifdef RAY_WINDOWS

#define VK_KEY_0 48
#define VK_KEY_1 49
#define VK_KEY_2 50
#define VK_KEY_3 51
#define VK_KEY_4 52
#define VK_KEY_5 53
#define VK_KEY_6 54
#define VK_KEY_7 55
#define VK_KEY_8 56
#define VK_KEY_9 57
//These are lower case
#define VK_KEY_A 65
#define VK_KEY_B 66
#define VK_KEY_C 67
#define VK_KEY_D 68
#define VK_KEY_E 69
#define VK_KEY_F 70
#define VK_KEY_G 71
#define VK_KEY_H 72
#define VK_KEY_I 73
#define VK_KEY_J 74
#define VK_KEY_K 75
#define VK_KEY_L 76
#define VK_KEY_M 77
#define VK_KEY_N 78
#define VK_KEY_O 79
#define VK_KEY_P 80
#define VK_KEY_Q 81
#define VK_KEY_R 82
#define VK_KEY_S 83
#define VK_KEY_T 84
#define VK_KEY_U 85
#define VK_KEY_V 86
#define VK_KEY_W 87
#define VK_KEY_X 88
#define VK_KEY_Y 89
#define VK_KEY_Z 90

LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam) {
  switch (uMsg) {
    case WM_DESTROY: {
      PostQuitMessage(0);
      term = true;
      DestroyWindow(hwnd);
      return 0;
    }
    case WM_KEYDOWN: {
      vec3f w(1,0,0);
      vec3f u(0,1,0);
      vec3f v(0,0,1);
      
      if(interactive_w) {
        w = cam_w->get_w();
        u = cam_w->get_u();
        v = cam_w->get_v();
      }

      switch (wParam) {
        case VK_ESCAPE: {
          PostQuitMessage(0);
          term = true;
          DestroyWindow(hwnd);
          return 0;
        }
        case VK_TAB: {
          orbit = !orbit;
          break;
        }
        case VK_KEY_F: {
          if(interactive_w) {
            (*write_fast_output_w) = !(*write_fast_output_w);
          }
          break;
        }
        case VK_KEY_W: {
          if(interactive_w) {
            vec3f step = -speed * w * base_step;
            if(orbit) {
              Float dist_to_orbit = (cam_w->get_origin() - cam_w->get_lookat()).length();
              if(dist_to_orbit <= base_step * speed) {
                Rprintf("Moving forward will overstep orbit point, stopping (decrease step size to move closer).\n");
                step = vec3f(0);
              }
            } 
            cam_w->update_position(step, orbit, false);
          }
          break;
        }

        case VK_KEY_A: {
          if(interactive_w) {
          cam_w->update_position(-speed * u * base_step, orbit);
        }
          break;
        }
        case VK_KEY_S: {
          if(interactive_w) {
          cam_w->update_position(speed * w * base_step, orbit, false);
        }
          break;
        }
        case VK_KEY_D: {
          if(interactive_w) {
          cam_w->update_position(speed * u * base_step, orbit);
        }
          break;
        }
        case VK_KEY_Q: { 
          if(interactive_w) {
            cam_w->update_position(speed * v * base_step, orbit);
        }
          break;
          }
        case VK_KEY_Z: { 
          if(interactive_w) {
            cam_w->update_position(-speed * v * base_step, orbit);
        }
          break;
          }
        case VK_KEY_E: { 
          if(interactive_w) {
            speed = 2 * speed;
            speed = std::fmin(speed,128);
            Rprintf("Step Multiplier: %.3f\n", speed);
          }
          break;
          }
        case VK_KEY_C: { 
          if(interactive_w) {
            speed = 0.5 * speed;
            Rprintf("Step Multiplier: %.3f\n", speed);
          }
          break;
          }
        case VK_DOWN: {
          if(interactive_w) {
            cam_w->update_fov(speed*1.f);
        }
          break;
          }
        case VK_UP: {
          if(interactive_w) {
            cam_w->update_fov(speed*-1.f);
        }
          break;
          }
        case VK_LEFT: {
          if(interactive_w) {
            cam_w->update_aperture(speed*-0.1f);
        }
          break;
          }
        case VK_RIGHT: {
          if(interactive_w) {
            cam_w->update_aperture(speed*0.1f);
        }
          break;
          }
        case VK_KEY_1: {
          if(interactive_w) {
            cam_w->update_focal_distance(speed*-1.f);
          }
          break;
          }
        case VK_KEY_2: {
          if(interactive_w) {
            cam_w->update_focal_distance(speed*1.f);
          }
          break;
          }
        case VK_KEY_3: {
          if(interactive_w) {
            (*EnvObjectToWorld_w) =  RotateY(speed*8) * (*EnvObjectToWorld_w);
            (*EnvWorldToObject_w) =  RotateY(-speed*8) * (*EnvWorldToObject_w);
          }
          break;
        }
        case VK_KEY_4: {
          if(interactive_w) {
          (*EnvObjectToWorld_w) =  RotateY(-speed*8) * (*EnvObjectToWorld_w);
          (*EnvWorldToObject_w) =  RotateY(speed*8) * (*EnvWorldToObject_w);
          }
          break;
        }
        case VK_KEY_R: {
          if(interactive_w) {
            cam_w->reset();
            speed = 1;
            (*EnvWorldToObject_w) = Start_EnvWorldToObject_w;
            (*EnvObjectToWorld_w) = Start_EnvObjectToWorld_w;
        }
          break;
          }
        case VK_KEY_P: {
          point3f origin = cam_w->get_origin();
          Float fov =  cam_w->get_fov();
          point3f cam_direction = -cam_w->get_w();
          Float cam_aperture = cam_w->get_aperture();
          Float fd = cam_w->get_focal_distance();
          point2f ortho = cam_w->get_ortho();
          point3f cam_lookat = cam_w->get_lookat();
          if(fov > 0) {
            Rprintf("Lookfrom: c(%.2f, %.2f, %.2f) LookAt: c(%.2f, %.2f, %.2f) FOV: %.1f Aperture: %0.3f Focal Dist: %0.3f Env Rotation: %.2f\n",
                    origin.x(), origin.y(), origin.z(), 
                    origin.x()+cam_direction.x()*fd, origin.y()+cam_direction.y()*fd, origin.z()+cam_direction.z()*fd, 
                    fov, 
                    cam_aperture, fd, env_y_angle);
          } else if (fov < 0) {
            Rprintf("Lookfrom: c(%.2f, %.2f, %.2f) LookAt: c(%.2f, %.2f, %.2f) Focal Dist: %0.3f Env Rotation: %.2f\n",
                    origin.x(), origin.y(), origin.z(), 
                    origin.x()+cam_direction.x()*fd, origin.y()+cam_direction.y()*fd, origin.z()+cam_direction.z()*fd, 
                    fd, env_y_angle);
          } else {
            Rprintf("Lookfrom: c(%.2f, %.2f, %.2f) LookAt: c(%.2f, %.2f, %.2f) Focal Dist: %0.3f Env Rotation: %.2f\n",
                    origin.x(), origin.y(), origin.z(), 
                    cam_lookat.x(), cam_lookat.y(), cam_lookat.z(), 
                    fd, env_y_angle);
          }
            break;
          } 
          case VK_KEY_L: {
            if(Keyframes_w->size() > 0) {
              Rcpp::List LastKeyframe = Keyframes_w->at(Keyframes_w->size()-1);
              point3f last_pos = point3f(Rcpp::as<Float>(LastKeyframe["x"]),
                                         Rcpp::as<Float>(LastKeyframe["y"]),
                                         Rcpp::as<Float>(LastKeyframe["z"]));
              point3f last_lookat = point3f(Rcpp::as<Float>(LastKeyframe["dx"]),
                                            Rcpp::as<Float>(LastKeyframe["dy"]),
                                            Rcpp::as<Float>(LastKeyframe["dz"]));
              Float last_aperture = Rcpp::as<Float>(LastKeyframe["aperture"]);
              Float last_fov = Rcpp::as<Float>(LastKeyframe["fov"]);
              Float last_focal = Rcpp::as<Float>(LastKeyframe["focal"]);
              vec2f last_ortho = vec2f(Rcpp::as<Float>(LastKeyframe["orthox"]),
                                       Rcpp::as<Float>(LastKeyframe["orthoy"]));
              cam_w->update_focal_absolute(last_focal);
              cam_w->update_position_absolute(last_pos);
              cam_w->update_lookat(last_lookat);
              cam_w->update_aperture_absolute(last_aperture);
              cam_w->update_fov_absolute(last_fov);
              cam_w->update_ortho_absolute(last_ortho);
            } else {
              Rprintf("Can't reset to last keyframe: No keyframes have been saved. Use the R key to reset camera.");
            }
            break;
          }
          case VK_KEY_K: {
            point3f origin = cam_w->get_origin();
            Float fov =  cam_w->get_fov();
            point3f cam_direction = -cam_w->get_w();
            Float cam_aperture = cam_w->get_aperture();
            Float fd = cam_w->get_focal_distance();
            vec3f cam_up = cam_w->get_up();
            point2f ortho = cam_w->get_ortho();
            point3f cam_lookat = cam_w->get_lookat();
            
            if(fov > 0) {
              Keyframes_w->push_back(Rcpp::List::create(Named("x")  = origin.x(),
                                                     Named("y")  = origin.y(),
                                                     Named("z")  = origin.z(),
                                                     Named("dx") = origin.x()+cam_direction.x()*fd,
                                                     Named("dy") = origin.y()+cam_direction.y()*fd,
                                                     Named("dz") = origin.z()+cam_direction.z()*fd,
                                                     Named("aperture") = cam_aperture,
                                                     Named("fov") = fov,
                                                     Named("focal") = fd,
                                                     Named("orthox") = ortho.x(),
                                                     Named("orthoy") = ortho.y(),
                                                     Named("upx")  = cam_up.x(),
                                                     Named("upy")  = cam_up.y(),
                                                     Named("upz")  = cam_up.z()));
            } else if (fov < 0) {
              Keyframes_w->push_back(Rcpp::List::create(Named("x")  = origin.x(),
                                                     Named("y")  = origin.y(),
                                                     Named("z")  = origin.z(),
                                                     Named("dx") = origin.x()+cam_direction.x()*fd,
                                                     Named("dy") = origin.y()+cam_direction.y()*fd,
                                                     Named("dz") = origin.z()+cam_direction.z()*fd,
                                                     Named("aperture") = 0,
                                                     Named("fov") = 0,
                                                     Named("focal") = fd,
                                                     Named("orthox") = ortho.x(),
                                                     Named("orthoy") = ortho.y(),
                                                     Named("upx")  = cam_up.x(),
                                                     Named("upy")  = cam_up.y(),
                                                     Named("upz")  = cam_up.z()));
            } else {
              Keyframes_w->push_back(Rcpp::List::create(Named("x")  = origin.x(),
                                                     Named("y")  = origin.y(),
                                                     Named("z")  = origin.z(),
                                                     Named("dx") = cam_lookat.x(),
                                                     Named("dy") = cam_lookat.y(),
                                                     Named("dz") = cam_lookat.z(),
                                                     Named("aperture") = 0,
                                                     Named("fov") = 0,
                                                     Named("focal") = fd,
                                                     Named("orthox") = ortho.x(),
                                                     Named("orthoy") = ortho.y(),
                                                     Named("upx")  = cam_up.x(),
                                                     Named("upy")  = cam_up.y(),
                                                     Named("upz")  = cam_up.z()));
            }
            if(fov > 0) {
              Rprintf("Lookfrom: c(%.2f, %.2f, %.2f) LookAt: c(%.2f, %.2f, %.2f) FOV: %.1f Aperture: %0.3f Focal Dist: %0.3f Env Rotation: %.2f\n",
                      origin.x(), origin.y(), origin.z(), 
                      origin.x()+cam_direction.x()*fd, origin.y()+cam_direction.y()*fd, origin.z()+cam_direction.z()*fd, 
                      fov, 
                      cam_aperture, fd, env_y_angle);
            } else if (fov < 0) {
              Rprintf("Lookfrom: c(%.2f, %.2f, %.2f) LookAt: c(%.2f, %.2f, %.2f) Focal Dist: %0.3f Env Rotation: %.2f\n",
                      origin.x(), origin.y(), origin.z(), 
                      origin.x()+cam_direction.x()*fd, origin.y()+cam_direction.y()*fd, origin.z()+cam_direction.z()*fd, 
                      fd, env_y_angle);
            } else {
              Rprintf("Lookfrom: c(%.2f, %.2f, %.2f) LookAt: c(%.2f, %.2f, %.2f) Focal Dist: %0.3f Env Rotation: %.2f\n",
                      origin.x(), origin.y(), origin.z(), 
                      cam_lookat.x(), cam_lookat.y(), cam_lookat.z(), 
                      fd, env_y_angle);
            }
            break;
          } 
          
          default: 
            break;
      }
      if(interactive_w) {
        if(wParam != VK_KEY_P && wParam != VK_KEY_K && wParam != VK_TAB) {
          if(!blanked && !term) {
            blanked = true;
            *ns_w = 0;
            aps->reset();
            aps_small->reset();
            if(progress_w && !interactive_w) {
              pb_w->update(0);
            }
          }
        }
      }
      break;
    }
  case WM_LBUTTONDOWN: {
    if(interactive_w) {
        Float x = GET_X_LPARAM(lParam);
        Float y = GET_Y_LPARAM(lParam);
        Float fov = cam_w->get_fov();
        Float u = (Float(x)) / Float(width);
        Float v = (Float(y)) / Float(height);
        vec3f dir;
        bool just_direction = false;
        hit_record hrec;
        if(fov < 0) {
          CameraSample samp({1-u,v},point2f(0.5,0.5), 0.5);
          ray r2;
          cam_w->GenerateRay(samp,&r2);
          if(world_w->hit(r2, 0.001, FLT_MAX, hrec, *rng_w)) {
            if( hrec.shape->GetName() != "EnvironmentLight") {
              dir = -hrec.p;
            }  else {
              dir = -hrec.p;
              just_direction = true;
            }
          }
        } else if (fov > 0) {
          ray r2 = cam_w->get_ray(u,1-v, point3f(0),
                                0.5f);
          world_w->hit(r2, 0.001, FLT_MAX, hrec, *rng_w);
          dir = r2.direction();
          if( hrec.shape->GetName() == "EnvironmentLight") {
            just_direction = true;
          } 
        } else {
          ray r2 = cam_w->get_ray(u,1-v, point3f(0),
                                0.5f);
          if(world_w->hit(r2, 0.001, FLT_MAX, hrec, *rng_w)) {
            if( hrec.shape->GetName() != "EnvironmentLight") {
              dir = -(cam_w->get_origin()-hrec.p);
            } else {
              dir = -cam_w->get_w();
              just_direction = true;
            }
          } else {
            dir = -cam_w->get_w();
          }
        }
        if(!just_direction) {
          if(fov != 0 && fov != 360) {
            Float current_fd = cam_w->get_focal_distance();
            Float new_fd = (hrec.p-cam_w->get_origin()).length();
            cam_w->update_focal_distance(new_fd- current_fd);
          }
          cam_w->update_lookat(hrec.p);
        } else {
          cam_w->update_look_direction(-dir);
        }
        *ns_w = 0;
        aps->reset();
        aps_small->reset();
        if(progress_w && !interactive_w) {
          pb_w->update(0);
        }
      }
    break;
  }
  case WM_RBUTTONDOWN: {
    if(interactive_w) {
      Float x = GET_X_LPARAM(lParam);
      Float y = GET_Y_LPARAM(lParam);
      Float fov = cam_w->get_fov();
      Float u = (Float(x)) / Float(width);
      Float v = (Float(y)) / Float(height);
      vec3f dir;
      hit_record hrec;
      if(fov < 0) {
        CameraSample samp({1-u,v},point2f(0.5,0.5), 0.5);
        ray r2;
        cam_w->GenerateRay(samp,&r2);
        if(world_w->hit(r2, 0.001, FLT_MAX, hrec, *rng_w)) {
          if( hrec.shape->GetName() != "EnvironmentLight") {
            dir = -hrec.p;
          }  else {
            dir = -hrec.p;
          }
        }
      } else if (fov > 0) {
        ray r2 = cam_w->get_ray(u,1-v, point3f(0),
                                0.5f);
        dir = r2.direction();
      } else {
        ray r2 = cam_w->get_ray(u,1-v, point3f(0),
                              0.5f);
        if(world_w->hit(r2, 0.001, FLT_MAX, hrec, *rng_w)) {
          if( hrec.shape->GetName() != "EnvironmentLight") {
            dir = -(cam_w->get_origin()-hrec.p);
          } else {
            Rprintf("Clicking on the environment light while using a orthographic camera does not change the view.\n");
            dir = -cam_w->get_w();
          }
        } else {
          dir = -cam_w->get_w();
        }
      }
      cam_w->update_look_direction(-dir);
      *ns_w = 0;
      aps->reset();
      aps_small->reset();
      if(progress_w && !interactive_w) {
        pb_w->update(0);
      }
    }
    break;
  }
    case WM_PAINT: {
      PAINTSTRUCT ps;
      HDC hdc = BeginPaint(hwnd, &ps);

      COLORREF *arr = (COLORREF*) calloc(width*height, sizeof(COLORREF));

      for(unsigned int i = 0; i < width*height*3; i += 3) {
        arr[i/3] = ((unsigned int)(255*rgb[i]) << 16) | ((unsigned int)(255*rgb[i+1]) << 8) | (unsigned int)(255*rgb[i+2]);
      }

      HBITMAP map = CreateBitmap(width,
                                 height,
                                 1,
                                 8*4,
                                 (void*) arr);

      HDC src = CreateCompatibleDC(hdc);
      SelectObject(src, map);

      // Copy image from temp HDC to window
      BitBlt(hdc,
             0,
             0,
             width,
             height,
             src,
             0,
             0,
             SRCCOPY); // Defined DWORD to just copy pixels.
      DeleteDC(src);
      DeleteObject(map);
      free(arr);

      EndPaint(hwnd, &ps);
    }
  return 0;
  }
  return DefWindowProc(hwnd, uMsg, wParam, lParam);
}
#endif

