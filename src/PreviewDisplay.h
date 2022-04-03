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
#endif

#ifdef RAY_WINDOWS
#ifndef UNICODE
#define UNICODE
#endif 


#include <windows.h>
#include <winuser.h>
#include "float.h"
#include <wingdi.h>

#endif

class PreviewDisplay {
public: 
  PreviewDisplay(unsigned int _width, unsigned int _height, bool preview, bool _interactive,
                 Float initial_lookat_distance, RayCamera* _cam,
                 Transform* _EnvObjectToWorld, Transform* _EnvWorldToObject);
  ~PreviewDisplay();
  void DrawImage(adaptive_sampler& adaptive_pixel_sampler, 
                 adaptive_sampler& adaptive_pixel_sampler_small,
                 size_t &ns,
                 RProgress::RProgress &pb, bool progress,
                 Float percent_done,
                 hitable* world, random_gen& rng);
  std::vector<Rcpp::List> GetKeyframes() {return(Keyframes);}
#ifdef RAY_HAS_X11
  Display *d;
  XImage *img;
  std::unique_ptr<char[]> data;
  Window w;
  XEvent e;
  unsigned int width;
  unsigned int height;
  Float speed;
  bool preview;
  bool orbit;
  Float base_step;
  int s;
#endif
#ifdef RAY_WINDOWS
  HWND hwnd;
  HINSTANCE hInstance;
  MSG msg;
  WNDCLASS wc;
#endif
  bool write_fast_output;
  bool interactive;
  bool terminate;
  RayCamera* cam;
  Transform* EnvObjectToWorld;
  Transform* EnvWorldToObject;
  Transform Start_EnvObjectToWorld;
  Transform Start_EnvWorldToObject;
  std::vector<Rcpp::List> Keyframes;
  
};

#endif