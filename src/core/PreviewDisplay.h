#ifndef PREVIEWDISPLAYH
#define PREVIEWDISPLAYH

#include <memory>
#include <string>
#include <vector>
#include "Rcpp.h"
#include "RProgress.h"
#include "../core/adaptivesampler.h"
#include "../core/camera.h"
#include "../hitables/hitable.h"

class VolumeScene;

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
#include "../core/oidn_denoiser.h"
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
                 Transform* _EnvObjectToWorld, Transform* _EnvWorldToObject,
                 RayOidnDenoiser* _denoiser,
                 RayMatrix* _oidn_albedo_output,
                 RayMatrix* _oidn_normal_output,
                 bool denoise, bool _auto_exposure);
#else
  PreviewDisplay(unsigned int _width, unsigned int _height, bool preview, bool _interactive,
                 bool _deferred_render, Float initial_lookat_distance, RayCamera* _cam,
                 Transform* _EnvObjectToWorld, Transform* _EnvWorldToObject,
                 bool _auto_exposure);
#endif
  ~PreviewDisplay();
  void SetCamera(RayCamera* _cam);
  // Coordinates match the renderer's film samples, independent of window API.
  bool PickCameraTarget(Float u, Float v, bool update_focus, hitable* world);
  std::shared_ptr<VolumeScene> volume_scene;
  void SetSnapshotFilename(const std::string& filename);
  void SavePreviewSnapshot() const;
  void CaptureVolumeSnapshot(adaptive_sampler&, RayMatrix&, size_t samples, hitable*, random_gen&);
  bool transparent_volume_background = false;
  std::vector<Float> snapshot_alpha;
  bool PollCloseEvent();
#ifdef HAS_OIDN
  void SetDenoiser(RayOidnDenoiser* _denoiser,
                   RayMatrix* _oidn_albedo_output,
                   RayMatrix* _oidn_normal_output,
                   bool _denoise);
  void InvalidateOidnAux();
  void MarkOidnAuxClean(bool fast_preview);
  void MarkDenoisedPreviewReady(size_t sample_count);
  bool HasDenoisedPreview() const { return has_denoised_preview; }
  size_t DenoisedPreviewSampleCount() const { return denoised_preview_sample_count; }
#endif
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
  void SetShutterSpeed(Float value);
  Float GetShutterSpeed() const;
  void AdjustShutterSpeedStops(Float stops);
  void ApplyShutterSpeedToCameras();
  void PrintShutterSpeed() const;
  Rcpp::List CreateCurrentKeyframe(Float env_rotation) const;
  void SaveCurrentKeyframe(Float env_rotation);
  bool ApplyCameraState(const Rcpp::List& state, Float* env_rotation);
  bool ApplyKeyframe(int index, Float* env_rotation);
  bool JumpKeyframe(int step, Float* env_rotation);
  bool DeleteCurrentKeyframe(Float* env_rotation);
  void PrintCameraInfo(Float env_rotation) const;
  std::string PreviewStatusText(Float env_rotation) const;
  void SetKeyframeMotionArgs(const Rcpp::List& args);
  bool ToggleKeyframeMotionClosed();
  bool KeyframeMotionClosed() const { return keyframe_motion_closed; }
  Rcpp::DataFrame KeyframesDataFrame() const;
  bool StartPreviewMotion(Float env_rotation);
  bool CancelPreviewMotion(Float* env_rotation);
  bool AdvancePreviewMotion(Float* env_rotation);
  bool IsPreviewMotionActive() const { return preview_motion_active; }
  bool ToggleCameraMotionBlur();
  void SetCameraMotionBlur(bool enabled);
  bool CameraMotionBlurEnabled() const { return camera_motion_blur_enabled; }
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
  void DrawStatusBarX11(Float env_rotation);
  void CompositeTextOverlaysToX11Buffer(hitable* world, random_gen& rng);
  void CompositeLineOverlaysToX11Buffer(hitable* world, random_gen& rng);
#endif
#ifdef RAY_WINDOWS
  void DrawStatusBarWindows(HDC hdc, Float env_rotation) const;
#endif
  void CompositeTextOverlaysToFloatBuffer(std::vector<Float>& rgb,
                                          hitable* world,
                                          random_gen& rng, std::vector<Float>* coverage = nullptr);
  void CompositeLineOverlaysToFloatBuffer(std::vector<Float>& rgb,
                                          hitable* world,
                                          random_gen& rng, std::vector<Float>* coverage = nullptr);
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
  RayOidnDenoiser* denoiser;
  RayMatrix* oidn_albedo_output;
  RayMatrix* oidn_normal_output;
  bool denoise;
  bool oidn_aux_dirty;
  bool oidn_fast_aux_dirty;
  bool has_denoised_preview;
  size_t denoised_preview_sample_count;
  #endif
  std::vector<Rcpp::List> Keyframes;
  int current_keyframe;
  Rcpp::List keyframe_motion_args;
  bool keyframe_motion_closed;
  Rcpp::DataFrame preview_motion;
  Rcpp::List preview_motion_restore_state;
  int preview_motion_frame;
  int preview_motion_restore_keyframe;
  bool preview_motion_active;
  bool camera_motion_blur_enabled;
  Float shutter_speed;
  std::string snapshot_filename;
  std::vector<unsigned char> snapshot_pixels;
  unsigned int snapshot_width;
  unsigned int snapshot_height;
  std::vector<PreviewTextOverlay> text_overlays;
  std::vector<PreviewLineOverlay> line_overlays;
};

#endif
