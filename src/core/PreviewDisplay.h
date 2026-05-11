#ifndef PREVIEWDISPLAYH
#define PREVIEWDISPLAYH

#include <memory>
#include <vector>
#include "Rcpp.h"
#include "RProgress.h"
#include "../core/adaptivesampler.h"
#include "../core/camera.h"
#include "../hitables/hitable.h"

struct PreviewTextOverlay {
  point3f anchor;
  int x_offset;
  int y_offset;
  Float hjust;
  Float vjust;
  bool clip;
  bool occlusion;
  bool partial_occlusion;
  Float occlusion_tolerance;
  unsigned int width;
  unsigned int height;
  std::vector<unsigned char> rgba;
};

struct PreviewLineOverlay {
  point3f start;
  point3f end;
  int x_offset;
  int y_offset;
  int xend_offset;
  int yend_offset;
  Float width;
  Float red;
  Float green;
  Float blue;
  Float alpha;
  int lineend;
  bool clip;
  bool occlusion;
  bool partial_occlusion;
  Float occlusion_tolerance;
};

#ifdef RAY_HAS_X11
#include <X11/Xlib.h>
#undef Status

#endif

#ifdef HAS_OIDN
#undef None
#include <OpenImageDenoise/oidn.hpp>
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
#ifdef HAS_OIDN
  PreviewDisplay(unsigned int _width, unsigned int _height, bool preview, bool _interactive,
                 bool _deferred_render, Float initial_lookat_distance, RayCamera* _cam,
                 Transform* _EnvObjectToWorld, Transform* _EnvWorldToObject, oidn::FilterRef& _filter,
                 bool denoise, bool _auto_exposure);
#else
  PreviewDisplay(unsigned int _width, unsigned int _height, bool preview, bool _interactive,
                 bool _deferred_render, Float initial_lookat_distance, RayCamera* _cam,
                 Transform* _EnvObjectToWorld, Transform* _EnvWorldToObject,
                 bool _auto_exposure);
#endif
  ~PreviewDisplay();
  void DrawImage(adaptive_sampler& adaptive_pixel_sampler, 
                 adaptive_sampler& adaptive_pixel_sampler_small,
                 size_t &ns,
                 RProgress::RProgress &pb, bool progress,
                 Float percent_done,
                 hitable* world, random_gen& rng);
  std::vector<Rcpp::List> GetKeyframes() {return(Keyframes);}
  void CalibratePreviewExposure(adaptive_sampler& adaptive_pixel_sampler,
                                RayMatrix& rgb,
                                size_t ns);
  void ResetPreviewExposure();
  Float ApplyPreviewExposure(Float value, Float sample_count) const;
  void IncreasePreviewExposure();
  void DecreasePreviewExposure();
  void SetTextOverlays(const std::vector<PreviewTextOverlay>& overlays);
  void SetLineOverlays(const std::vector<PreviewLineOverlay>& overlays);
  bool ProjectTextAnchor(const PreviewTextOverlay& overlay,
                         Float& screen_x,
                         Float& screen_y) const;
  bool ProjectWorldPoint(const point3f& point,
                         bool clip,
                         Float& screen_x,
                         Float& screen_y,
                         Float& depth) const;
  bool IsTextAnchorOccluded(const PreviewTextOverlay& overlay,
                            hitable* world,
                            random_gen& rng) const;
  bool IsTextPixelOccluded(const PreviewTextOverlay& overlay,
                           Float screen_x,
                           Float screen_y,
                           hitable* world,
                           random_gen& rng) const;
  bool IsLineAnchorOccluded(const PreviewLineOverlay& overlay,
                            hitable* world,
                            random_gen& rng) const;
  bool IsLinePixelOccluded(const PreviewLineOverlay& overlay,
                           Float screen_x,
                           Float screen_y,
                           Float line_depth,
                           hitable* world,
                           random_gen& rng) const;
#ifdef RAY_HAS_X11
  void CompositeTextOverlaysToX11Buffer(hitable* world, random_gen& rng);
  void CompositeLineOverlaysToX11Buffer(hitable* world, random_gen& rng);
#endif
#ifdef RAY_WINDOWS
  void CompositeTextOverlaysToFloatBuffer(std::vector<Float>& rgb,
                                          hitable* world,
                                          random_gen& rng);
  void CompositeLineOverlaysToFloatBuffer(std::vector<Float>& rgb,
                                          hitable* world,
                                          random_gen& rng);
#endif
#ifdef RAY_HAS_X11
  Display *d;
  XImage *img;
  std::unique_ptr<char[]> data;
  Window w;
  XEvent e;
  unsigned int width;
  unsigned int height;
  Float speed;
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
  bool preview;
  bool auto_exposure;
  bool preview_exposure_calibrated;
  Float preview_exposure_scale;
  Float preview_exposure_adjustment;
  bool write_fast_output;
  bool interactive;
  bool deferred_render;
  bool render_requested;
  bool terminate;
  RayCamera* cam;
  Transform* EnvObjectToWorld;
  Transform* EnvWorldToObject;
  Transform Start_EnvObjectToWorld;
  Transform Start_EnvWorldToObject;
  #ifdef HAS_OIDN
  oidn::FilterRef& filter;
  bool denoise;
  #endif
  std::vector<Rcpp::List> Keyframes;
  std::vector<PreviewTextOverlay> text_overlays;
  std::vector<PreviewLineOverlay> line_overlays;
};

#endif
