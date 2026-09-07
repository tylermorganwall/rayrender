#include "Rcpp.h"
#include "../core/PreviewDisplay.h"
#include "../math/mathinline.h"
#include "../utils/raylog.h"
#include "../volumes/picking.h"
#include "RcppThread.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

#ifdef NOT_CRAN
#include "../hitables/sphere.h"
#include "../materials/material.h"
#include <testthat.h>
#endif

static const unsigned int PREVIEW_STATUS_MIN_WIDTH = 640;
static const unsigned int PREVIEW_STATUS_BAR_HEIGHT = 24;
static const Float PREVIEW_SHUTTER_SPEED_MIN = static_cast<Float>(1);
static const Float PREVIEW_SHUTTER_SPEED_MAX = static_cast<Float>(4096);

static point3f PreviewCameraLookat(RayCamera* cam) {
  point3f origin = cam->get_origin(), pivot = cam->get_lookat();
  vec3f forward = unit_vector(cam->get_w()), to_pivot = pivot - origin;
  if(to_pivot.length() > 0 && dot(unit_vector(to_pivot), forward) > .999999f) {
    return pivot;
  }
  // Free-flight translation can leave the stored orbit point off the viewing
  // axis. Capture that view without replacing an aligned, explicitly picked pivot.
  return origin + forward * std::max(cam->get_focal_distance(), Float(.001));
}

struct PreviewCameraState {
  point3f origin;
  point3f lookat;
  vec3f up;
  Float focal;
};

static PreviewCameraState CapturePreviewCameraState(RayCamera* cam) {
  PreviewCameraState state;
  state.origin = cam->get_origin();
  state.lookat = PreviewCameraLookat(cam);
  state.up = cam->get_up();
  state.focal = cam->get_focal_distance();
  return state;
}

static void ApplyPreviewCameraMotionRange(RayCamera* cam,
                                          const PreviewCameraState& start,
                                          const PreviewCameraState& end) {
  cam->set_camera_motion_blur_range(start.origin,
                                    start.lookat,
                                    start.up,
                                    start.focal,
                                    end.origin,
                                    end.lookat,
                                    end.up,
                                    end.focal);
}

static void ApplyStaticPreviewCameraMotionRange(RayCamera* cam) {
  if(cam == nullptr) {
    return;
  }
  PreviewCameraState state = CapturePreviewCameraState(cam);
  ApplyPreviewCameraMotionRange(cam, state, state);
}

bool PreviewDisplay::PickCameraTarget(Float u, Float v, bool update_focus, hitable* world) {
  if(!cam || !world || u < 0 || u > 1 || v < 0 || v > 1) {
    return false;
  }
  Ray ray;
  Float fov = cam->get_fov();
  if(fov < 0) {
    // Match GenerateRay's film convention, using the center of the lens and
    // a fixed shutter sample rather than the current rendering sample.
    CameraSample sample(point2f(1 - u, 1 - v), point2f(.5f, .5f), .5f);
    if(!(cam->GenerateRay(sample, &ray) > 0)) {
      return false;
    }
  } else {
    ray = cam->get_ray(u, v, point3f(0), .5f);
  }
  bool cancelled = false;
  size_t polls = 0;
  auto cancel = [&] {
    if(!cancelled && polls++ % 64 == 0) {
      cancelled = PollCloseEvent() || RcppThread::isInterrupted();
    }
    return cancelled;
  };
  try {
    auto target = PickRay(ray, world, volume_scene.get(), .15, cancel);
    if(!target || (target->p - cam->get_origin()).length() <= 0) {
      return false;
    }
    if(update_focus && fov != 0 && fov != 360) {
      cam->update_focal_distance((target->p - cam->get_origin()).length() -
                                cam->get_focal_distance());
    }
    // This updates both the stored orbit point and the camera frame. Updating
    // just the direction would snap back to the previous pivot on the next orbit.
    cam->update_lookat(target->p);
    ApplyStaticPreviewCameraMotionRange(cam);
    return true;
  } catch(const std::exception& error) {
    Rprintf("Unable to pick preview target: %s\n", error.what());
    return false;
  }
}

static bool IsKeyframeSuppliedMotionArg(const std::string& name) {
  return name == "positions" ||
         name == "lookats" ||
         name == "apertures" ||
         name == "fovs" ||
         name == "focal_distances" ||
         name == "ortho_dims" ||
         name == "camera_ups";
}


#ifdef RAY_HAS_X11

//Had to undefine Xlib.h's Status because RcppThread also defines Status 
#define Status int
#include <string.h>
#include <X11/Xutil.h>
#include "X11/keysym.h"
#include "X11/Xatom.h"
static Float env_y_angle;;
#undef Status

static bool PreviewDisplayHasKeyboardFocus(Display* display, Window window) {
  Window focus_window;
  int revert_to;
  XGetInputFocus(display, &focus_window, &revert_to);
  return focus_window == window;
}

static bool IsX11RenderInvalidatingKey(Display* display, KeyCode keycode) {
  return keycode == XKeysymToKeycode(display, XStringToKeysym("w")) ||
         keycode == XKeysymToKeycode(display, XStringToKeysym("a")) ||
         keycode == XKeysymToKeycode(display, XStringToKeysym("s")) ||
         keycode == XKeysymToKeycode(display, XStringToKeysym("d")) ||
         keycode == XKeysymToKeycode(display, XStringToKeysym("q")) ||
         keycode == XKeysymToKeycode(display, XStringToKeysym("z")) ||
         keycode == XKeysymToKeycode(display, XK_Up) ||
         keycode == XKeysymToKeycode(display, XK_Down) ||
         keycode == XKeysymToKeycode(display, XK_Left) ||
         keycode == XKeysymToKeycode(display, XK_Right) ||
         keycode == XKeysymToKeycode(display, XStringToKeysym("1")) ||
         keycode == XKeysymToKeycode(display, XStringToKeysym("2")) ||
         keycode == XKeysymToKeycode(display, XStringToKeysym("3")) ||
         keycode == XKeysymToKeycode(display, XStringToKeysym("4")) ||
         keycode == XKeysymToKeycode(display, XStringToKeysym("f")) ||
         keycode == XKeysymToKeycode(display, XStringToKeysym("b")) ||
         keycode == XKeysymToKeycode(display, XStringToKeysym("l")) ||
         keycode == XKeysymToKeycode(display, XK_less) ||
         keycode == XKeysymToKeycode(display, XK_greater) ||
         keycode == XKeysymToKeycode(display, XK_slash) ||
         keycode == XKeysymToKeycode(display, XStringToKeysym("r"));
}

static KeySym X11EventKeysym(XKeyEvent event) {
  return XLookupKeysym(&event, (event.state & ShiftMask) != 0 ? 1 : 0);
}

static bool IsX11ShiftedLeftBracket(Display* display, XKeyEvent event) {
  bool shift_pressed = (event.state & ShiftMask) != 0;
  KeySym key_symbol = X11EventKeysym(event);
  return key_symbol == XK_braceleft ||
         (shift_pressed && event.keycode == XKeysymToKeycode(display, XK_bracketleft));
}

static bool IsX11ShiftedRightBracket(Display* display, XKeyEvent event) {
  bool shift_pressed = (event.state & ShiftMask) != 0;
  KeySym key_symbol = X11EventKeysym(event);
  return key_symbol == XK_braceright ||
         (shift_pressed && event.keycode == XKeysymToKeycode(display, XK_bracketright));
}

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
static bool deferred_render_w;
static bool* render_requested_w;
static Float* preview_exposure_adjustment_w;
static PreviewDisplay* preview_display_w;



LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam);

#endif


void PreviewDisplay::CalibratePreviewExposure(adaptive_sampler& adaptive_pixel_sampler,
                                              RayMatrix& rgb,
                                              size_t ns) {
  if(!auto_exposure || preview_exposure_calibrated) {
    return;
  }
  std::vector<bool>& finalized = adaptive_pixel_sampler.finalized;
  std::vector<Float> luminance;
  luminance.reserve(rgb.size());

  for(unsigned int x = 0; x < rgb.rows(); x++) {
    for(unsigned int y = 0; y < rgb.cols(); y++) {
      Float sample_count = (Float)ns + 1.f;
      if(finalized[x + rgb.rows() * y]) {
        sample_count = interactive ? 1.f : 4.f;
      }
      Float r = std::fmax((Float)0, rgb(x, y, 0) / sample_count);
      Float g = std::fmax((Float)0, rgb(x, y, 1) / sample_count);
      Float b = std::fmax((Float)0, rgb(x, y, 2) / sample_count);
      Float l = 0.2126f * r + 0.7152f * g + 0.0722f * b;
      if(std::isfinite(l)) {
        luminance.push_back(l);
      }
    }
  }

  preview_exposure_scale = 1.f;
  if(!luminance.empty()) {
    size_t quantile_index = static_cast<size_t>(std::ceil(0.9 * luminance.size())) - 1;
    quantile_index = std::min(quantile_index, luminance.size() - 1);
    std::nth_element(luminance.begin(),
                     luminance.begin() + quantile_index,
                     luminance.end());
    Float quantile = luminance[quantile_index];
    if(std::isfinite(quantile) && quantile > 0.f) {
      preview_exposure_scale = 1.f / quantile;
    }
  }
  preview_exposure_calibrated = true;
}

Float PreviewDisplay::ApplyPreviewExposure(Float value, Float sample_count) const {
  return std::sqrt(std::fmax((Float)0,
                            value * preview_exposure_scale *
                            preview_exposure_adjustment / sample_count));
}

void PreviewDisplay::ResetPreviewExposure() {
  preview_exposure_calibrated = false;
  preview_exposure_scale = 1.f;
}

void PreviewDisplay::IncreasePreviewExposure() {
  preview_exposure_adjustment *= 2.f;
  Rprintf("Preview Exposure: %.3f\n", preview_exposure_adjustment);
}

void PreviewDisplay::DecreasePreviewExposure() {
  preview_exposure_adjustment *= 0.5f;
  Rprintf("Preview Exposure: %.3f\n", preview_exposure_adjustment);
}

void PreviewDisplay::SetShutterSpeed(Float value) {
  if(std::isnan(value) ||
     value < static_cast<Float>(1) ||
     (std::isinf(value) && value < static_cast<Float>(0))) {
    throw std::runtime_error("shutter_speed must be greater than or equal to 1, or Inf.");
  }
  shutter_speed = value;
  ApplyShutterSpeedToCameras();
}

Float PreviewDisplay::GetShutterSpeed() const {
  return shutter_speed;
}

void PreviewDisplay::AdjustShutterSpeedStops(Float stops) {
  Float next = shutter_speed;
  if(std::isinf(shutter_speed)) {
    if(stops < static_cast<Float>(0)) {
      next = PREVIEW_SHUTTER_SPEED_MAX;
    }
  } else {
    next = shutter_speed * std::pow(static_cast<Float>(2), stops);
    next = clamp(next, PREVIEW_SHUTTER_SPEED_MIN, PREVIEW_SHUTTER_SPEED_MAX);
  }
  SetShutterSpeed(next);
  PrintShutterSpeed();
}

void PreviewDisplay::ApplyShutterSpeedToCameras() {
  if(cam != nullptr) {
    cam->set_shutter_speed(shutter_speed);
  }
#ifdef RAY_WINDOWS
  if(cam_w != nullptr) {
    cam_w->set_shutter_speed(shutter_speed);
  }
#endif
}

void PreviewDisplay::PrintShutterSpeed() const {
  if(std::isinf(shutter_speed)) {
    Rprintf("Shutter speed: Inf; motion blur interval: 0\n");
    return;
  }
  Float interval = static_cast<Float>(1) / shutter_speed;
  Float angle = static_cast<Float>(360) * interval;
  Rprintf("Shutter speed: %.3f; interval: %.3f frame; angle: %.0f deg\n",
          shutter_speed, interval, angle);
}

void PreviewDisplay::SetCameraMotionBlur(bool enabled) {
  camera_motion_blur_enabled = enabled;
  if(cam != nullptr) {
    cam->set_camera_motion_blur(enabled);
    if(enabled && interactive && !preview_motion_active) {
      ApplyStaticPreviewCameraMotionRange(cam);
    }
  }
#ifdef RAY_WINDOWS
  if(cam_w != nullptr) {
    cam_w->set_camera_motion_blur(enabled);
    if(enabled && interactive && !preview_motion_active) {
      ApplyStaticPreviewCameraMotionRange(cam_w);
    }
  }
#endif
}

bool PreviewDisplay::ToggleCameraMotionBlur() {
  SetCameraMotionBlur(!camera_motion_blur_enabled);
  Rprintf("Camera Motion Blur: %s\n", camera_motion_blur_enabled ? "ON" : "OFF");
  return camera_motion_blur_enabled;
}

Rcpp::List PreviewDisplay::CreateCurrentKeyframe(Float env_rotation) const {
  point3f origin = cam->get_origin();
  Float fov = cam->get_fov();
  Float cam_aperture = cam->get_aperture();
  Float fd = cam->get_focal_distance();
  vec3f cam_up = cam->get_up();
  point2f ortho = cam->get_ortho();
  point3f key_lookat = PreviewCameraLookat(cam);
  Float key_aperture = fov > 0 ? cam_aperture : 0;
  Float key_fov = fov > 0 ? fov : 0;

  return Rcpp::List::create(Named("x") = origin.xyz.x,
                            Named("y") = origin.xyz.y,
                            Named("z") = origin.xyz.z,
                            Named("dx") = key_lookat.xyz.x,
                            Named("dy") = key_lookat.xyz.y,
                            Named("dz") = key_lookat.xyz.z,
                            Named("aperture") = key_aperture,
                            Named("fov") = key_fov,
                            Named("focal") = fd,
                            Named("exposure") = preview_exposure_adjustment,
                            Named("env_rotation") = env_rotation,
                            Named("orthox") = ortho.xy.x,
                            Named("orthoy") = ortho.xy.y,
                            Named("upx") = cam_up.xyz.x,
                            Named("upy") = cam_up.xyz.y,
                            Named("upz") = cam_up.xyz.z);
}

void PreviewDisplay::SaveCurrentKeyframe(Float env_rotation) {
  Keyframes.push_back(CreateCurrentKeyframe(env_rotation));
  current_keyframe = static_cast<int>(Keyframes.size()) - 1;
}

bool PreviewDisplay::ApplyCameraState(const Rcpp::List& state,
                                      Float* env_rotation) {
  point3f key_pos = point3f(Rcpp::as<Float>(state["x"]),
                            Rcpp::as<Float>(state["y"]),
                            Rcpp::as<Float>(state["z"]));
  point3f key_lookat = point3f(Rcpp::as<Float>(state["dx"]),
                               Rcpp::as<Float>(state["dy"]),
                               Rcpp::as<Float>(state["dz"]));
  Float key_aperture = Rcpp::as<Float>(state["aperture"]);
  Float key_fov = Rcpp::as<Float>(state["fov"]);
  Float key_focal = Rcpp::as<Float>(state["focal"]);
  vec3f key_up = cam->get_up();
  if(state.containsElementNamed("upx") &&
     state.containsElementNamed("upy") &&
     state.containsElementNamed("upz")) {
    key_up = vec3f(Rcpp::as<Float>(state["upx"]),
                   Rcpp::as<Float>(state["upy"]),
                   Rcpp::as<Float>(state["upz"]));
  }
  if(state.containsElementNamed("exposure")) {
    preview_exposure_adjustment = Rcpp::as<Float>(state["exposure"]);
  }
  if(env_rotation != nullptr && state.containsElementNamed("env_rotation")) {
    *env_rotation = Rcpp::as<Float>(state["env_rotation"]);
    (*EnvObjectToWorld) = RotateY(-*env_rotation) * Start_EnvObjectToWorld;
    (*EnvWorldToObject) = RotateY(*env_rotation) * Start_EnvWorldToObject;
  }
  vec2f key_ortho = vec2f(Rcpp::as<Float>(state["orthox"]),
                          Rcpp::as<Float>(state["orthoy"]));
  cam->update_focal_absolute(key_focal);
  cam->update_position_absolute(key_pos);
  cam->update_lookat(key_lookat);
  cam->update_up(key_up);
  cam->update_aperture_absolute(key_aperture);
  cam->update_fov_absolute(key_fov);
  cam->update_ortho_absolute(key_ortho);
  return true;
}

bool PreviewDisplay::ApplyKeyframe(int index, Float* env_rotation) {
  if(Keyframes.empty()) {
    current_keyframe = -1;
    Rprintf("Can't jump to keyframe: No keyframes have been saved. Use K to save a keyframe.\n");
    return false;
  }
  if(index < 0 || index >= static_cast<int>(Keyframes.size())) {
    return false;
  }

  Rcpp::List keyframe = Keyframes.at(index);
  ApplyCameraState(keyframe, env_rotation);
  current_keyframe = index;
  return true;
}

bool PreviewDisplay::JumpKeyframe(int step, Float* env_rotation) {
  if(Keyframes.empty()) {
    current_keyframe = -1;
    Rprintf("Can't jump to keyframe: No keyframes have been saved. Use K to save a keyframe.\n");
    return false;
  }

  int keyframe_count = static_cast<int>(Keyframes.size());
  if(current_keyframe < 0 || current_keyframe >= keyframe_count) {
    current_keyframe = step >= 0 ? 0 : keyframe_count - 1;
  } else {
    current_keyframe = (current_keyframe + step + keyframe_count) % keyframe_count;
  }
  return ApplyKeyframe(current_keyframe, env_rotation);
}

bool PreviewDisplay::DeleteCurrentKeyframe(Float* env_rotation) {
  if(Keyframes.empty()) {
    current_keyframe = -1;
    Rprintf("Can't delete keyframe: No keyframes have been saved. Use K to save a keyframe.\n");
    return false;
  }

  int keyframe_count = static_cast<int>(Keyframes.size());
  if(current_keyframe < 0 || current_keyframe >= keyframe_count) {
    current_keyframe = keyframe_count - 1;
  }
  int deleted_keyframe = current_keyframe;
  Keyframes.erase(Keyframes.begin() + current_keyframe);

  if(Keyframes.empty()) {
    current_keyframe = -1;
    Rprintf("Deleted keyframe %d. 0 keyframes saved.\n", deleted_keyframe + 1);
    return true;
  }

  if(current_keyframe >= static_cast<int>(Keyframes.size())) {
    current_keyframe = static_cast<int>(Keyframes.size()) - 1;
  }
  Rprintf("Deleted keyframe %d. %zu keyframes saved.\n",
          deleted_keyframe + 1,
          Keyframes.size());
  return ApplyKeyframe(current_keyframe, env_rotation);
}

void PreviewDisplay::SetKeyframeMotionArgs(const Rcpp::List& args) {
  keyframe_motion_args = Rcpp::clone(args);
  keyframe_motion_closed = args.containsElementNamed("closed") ?
    Rcpp::as<bool>(args["closed"]) : false;
}

bool PreviewDisplay::ToggleKeyframeMotionClosed() {
  keyframe_motion_closed = !keyframe_motion_closed;
  Rprintf("Keyframe motion loop: %s.\n",
          keyframe_motion_closed ? "CLOSED" : "OPEN");
  return keyframe_motion_closed;
}

Rcpp::DataFrame PreviewDisplay::KeyframesDataFrame() const {
  size_t keyframe_count = Keyframes.size();
  Rcpp::NumericVector x(keyframe_count);
  Rcpp::NumericVector y(keyframe_count);
  Rcpp::NumericVector z(keyframe_count);
  Rcpp::NumericVector dx(keyframe_count);
  Rcpp::NumericVector dy(keyframe_count);
  Rcpp::NumericVector dz(keyframe_count);
  Rcpp::NumericVector aperture(keyframe_count);
  Rcpp::NumericVector fov(keyframe_count);
  Rcpp::NumericVector focal(keyframe_count);
  Rcpp::NumericVector exposure(keyframe_count);
  Rcpp::NumericVector env_rotation(keyframe_count);
  Rcpp::NumericVector orthox(keyframe_count);
  Rcpp::NumericVector orthoy(keyframe_count);
  Rcpp::NumericVector upx(keyframe_count);
  Rcpp::NumericVector upy(keyframe_count);
  Rcpp::NumericVector upz(keyframe_count);

  for(size_t i = 0; i < keyframe_count; i++) {
    Rcpp::List keyframe = Keyframes.at(i);
    x[i] = Rcpp::as<Float>(keyframe["x"]);
    y[i] = Rcpp::as<Float>(keyframe["y"]);
    z[i] = Rcpp::as<Float>(keyframe["z"]);
    dx[i] = Rcpp::as<Float>(keyframe["dx"]);
    dy[i] = Rcpp::as<Float>(keyframe["dy"]);
    dz[i] = Rcpp::as<Float>(keyframe["dz"]);
    aperture[i] = Rcpp::as<Float>(keyframe["aperture"]);
    fov[i] = Rcpp::as<Float>(keyframe["fov"]);
    focal[i] = Rcpp::as<Float>(keyframe["focal"]);
    exposure[i] = keyframe.containsElementNamed("exposure") ?
      Rcpp::as<Float>(keyframe["exposure"]) : preview_exposure_adjustment;
    env_rotation[i] = keyframe.containsElementNamed("env_rotation") ?
      Rcpp::as<Float>(keyframe["env_rotation"]) : 0;
    orthox[i] = Rcpp::as<Float>(keyframe["orthox"]);
    orthoy[i] = Rcpp::as<Float>(keyframe["orthoy"]);
    upx[i] = Rcpp::as<Float>(keyframe["upx"]);
    upy[i] = Rcpp::as<Float>(keyframe["upy"]);
    upz[i] = Rcpp::as<Float>(keyframe["upz"]);
  }

  return Rcpp::DataFrame::create(Named("x") = x,
                                 Named("y") = y,
                                 Named("z") = z,
                                 Named("dx") = dx,
                                 Named("dy") = dy,
                                 Named("dz") = dz,
                                 Named("aperture") = aperture,
                                 Named("fov") = fov,
                                 Named("focal") = focal,
                                 Named("exposure") = exposure,
                                 Named("env_rotation") = env_rotation,
                                 Named("orthox") = orthox,
                                 Named("orthoy") = orthoy,
                                 Named("upx") = upx,
                                 Named("upy") = upy,
                                 Named("upz") = upz);
}

bool PreviewDisplay::StartPreviewMotion(Float env_rotation) {
  if(preview_motion_active) {
    return false;
  }
  if(Keyframes.size() < 2) {
    Rprintf("Can't preview keyframe motion: Save at least two keyframes with K.\n");
    return false;
  }

  try {
    Rcpp::DataFrame keyframes = KeyframesDataFrame();
    Rcpp::List call_args;
    call_args.push_back(keyframes, "positions");

    bool frames_supplied = false;
    Rcpp::CharacterVector arg_names = keyframe_motion_args.names();
    for(int i = 0; i < keyframe_motion_args.size(); i++) {
      std::string arg_name = Rcpp::as<std::string>(arg_names[i]);
      if(arg_name == "frames") {
        frames_supplied = true;
      }
      if(arg_name == "closed") {
        continue;
      }
      if(!IsKeyframeSuppliedMotionArg(arg_name)) {
        call_args.push_back(keyframe_motion_args[i], arg_name);
      }
    }
    if(!frames_supplied) {
      call_args.push_back(static_cast<int>(Keyframes.size() * 30), "frames");
    }
    call_args.push_back(keyframe_motion_closed, "closed");

    Rcpp::Environment pkg = Rcpp::Environment::namespace_env("rayrender");
    Rcpp::Function generate_camera_motion = pkg["generate_camera_motion"];
    Rcpp::Environment base = Rcpp::Environment::base_env();
    Rcpp::Function do_call = base["do.call"];
    preview_motion = Rcpp::as<Rcpp::DataFrame>(
      do_call(generate_camera_motion, call_args)
    );
  } catch(std::exception& ex) {
    Rprintf("Can't preview keyframe motion: %s\n", ex.what());
    return false;
  }

  if(preview_motion.nrows() < 1) {
    Rprintf("Can't preview keyframe motion: generate_camera_motion() returned no frames.\n");
    return false;
  }

  preview_motion_restore_state = CreateCurrentKeyframe(env_rotation);
  preview_motion_restore_keyframe = current_keyframe;
  preview_motion_frame = 0;
  current_keyframe = 0;
  preview_motion_active = true;
  Rprintf("Previewing keyframe motion (%d frames). Press M to cancel.\n",
          static_cast<int>(preview_motion.nrows()));
  return true;
}

bool PreviewDisplay::CancelPreviewMotion(Float* env_rotation) {
  if(!preview_motion_active) {
    return false;
  }
  ApplyCameraState(preview_motion_restore_state, env_rotation);
  current_keyframe = preview_motion_restore_keyframe;
  preview_motion_active = false;
  preview_motion_frame = 0;
  ApplyStaticPreviewCameraMotionRange(cam);
  Rprintf("Cancelled keyframe motion preview. Restored original camera.\n");
  return true;
}

bool PreviewDisplay::AdvancePreviewMotion(Float* env_rotation) {
  if(!preview_motion_active) {
    return false;
  }

  if(preview_motion_frame >= preview_motion.nrows()) {
    ApplyCameraState(preview_motion_restore_state, env_rotation);
    current_keyframe = preview_motion_restore_keyframe;
    preview_motion_active = false;
    preview_motion_frame = 0;
    ApplyStaticPreviewCameraMotionRange(cam);
    Rprintf("Finished keyframe motion preview. Restored original camera.\n");
    return true;
  }

  PreviewCameraState camera_state_before = CapturePreviewCameraState(cam);
  Rcpp::NumericVector x = preview_motion["x"];
  Rcpp::NumericVector y = preview_motion["y"];
  Rcpp::NumericVector z = preview_motion["z"];
  Rcpp::NumericVector dx = preview_motion["dx"];
  Rcpp::NumericVector dy = preview_motion["dy"];
  Rcpp::NumericVector dz = preview_motion["dz"];
  Rcpp::NumericVector aperture = preview_motion["aperture"];
  Rcpp::NumericVector fov = preview_motion["fov"];
  Rcpp::NumericVector focal = preview_motion["focal"];
  Rcpp::NumericVector orthox = preview_motion["orthox"];
  Rcpp::NumericVector orthoy = preview_motion["orthoy"];
  Rcpp::NumericVector upx = preview_motion["upx"];
  Rcpp::NumericVector upy = preview_motion["upy"];
  Rcpp::NumericVector upz = preview_motion["upz"];
  int keyframe_count = static_cast<int>(Keyframes.size());
  if(keyframe_count > 0) {
    double frames_per_keyframe =
      static_cast<double>(preview_motion.nrows()) / static_cast<double>(keyframe_count);
    int playback_keyframe = frames_per_keyframe > 0 ?
      static_cast<int>(std::floor(preview_motion_frame / frames_per_keyframe)) : 0;
    current_keyframe = std::max(0, std::min(playback_keyframe, keyframe_count - 1));
  }

  Rcpp::List state = Rcpp::List::create(
    Named("x") = x[preview_motion_frame],
    Named("y") = y[preview_motion_frame],
    Named("z") = z[preview_motion_frame],
    Named("dx") = dx[preview_motion_frame],
    Named("dy") = dy[preview_motion_frame],
    Named("dz") = dz[preview_motion_frame],
    Named("aperture") = aperture[preview_motion_frame],
    Named("fov") = fov[preview_motion_frame],
    Named("focal") = focal[preview_motion_frame],
    Named("orthox") = orthox[preview_motion_frame],
    Named("orthoy") = orthoy[preview_motion_frame],
    Named("upx") = upx[preview_motion_frame],
    Named("upy") = upy[preview_motion_frame],
    Named("upz") = upz[preview_motion_frame]
  );
  ApplyCameraState(state, env_rotation);
  ApplyPreviewCameraMotionRange(cam,
                                camera_state_before,
                                CapturePreviewCameraState(cam));
  preview_motion_frame++;
  return true;
}

void PreviewDisplay::PrintCameraInfo(Float env_rotation) const {
  point3f origin = cam->get_origin();
  Float fov = cam->get_fov();
  Float cam_aperture = cam->get_aperture();
  Float fd = cam->get_focal_distance();
  point3f key_lookat = PreviewCameraLookat(cam);
  const char* shutter_label = std::isinf(shutter_speed) ? "Inf" : "";

  if(fov > 0) {
    if(std::isinf(shutter_speed)) {
      Rprintf("Lookfrom: c(%.2f, %.2f, %.2f) LookAt: c(%.2f, %.2f, %.2f) FOV: %.1f Aperture: %0.3f Focal Dist: %0.3f Env Rotation: %.2f Exposure: %.3f Camera Motion Blur: %s Shutter Speed: %s\n",
              origin.xyz.x, origin.xyz.y, origin.xyz.z,
              key_lookat.xyz.x, key_lookat.xyz.y, key_lookat.xyz.z,
              fov,
              cam_aperture, fd, env_rotation, preview_exposure_adjustment,
              camera_motion_blur_enabled ? "ON" : "OFF",
              shutter_label);
    } else {
      Rprintf("Lookfrom: c(%.2f, %.2f, %.2f) LookAt: c(%.2f, %.2f, %.2f) FOV: %.1f Aperture: %0.3f Focal Dist: %0.3f Env Rotation: %.2f Exposure: %.3f Camera Motion Blur: %s Shutter Speed: %.3f\n",
              origin.xyz.x, origin.xyz.y, origin.xyz.z,
              key_lookat.xyz.x, key_lookat.xyz.y, key_lookat.xyz.z,
              fov,
              cam_aperture, fd, env_rotation, preview_exposure_adjustment,
              camera_motion_blur_enabled ? "ON" : "OFF",
              shutter_speed);
    }
  } else {
    if(std::isinf(shutter_speed)) {
      Rprintf("Lookfrom: c(%.2f, %.2f, %.2f) LookAt: c(%.2f, %.2f, %.2f) Focal Dist: %0.3f Env Rotation: %.2f Exposure: %.3f Camera Motion Blur: %s Shutter Speed: %s\n",
              origin.xyz.x, origin.xyz.y, origin.xyz.z,
              key_lookat.xyz.x, key_lookat.xyz.y, key_lookat.xyz.z,
              fd, env_rotation, preview_exposure_adjustment,
              camera_motion_blur_enabled ? "ON" : "OFF",
              shutter_label);
    } else {
      Rprintf("Lookfrom: c(%.2f, %.2f, %.2f) LookAt: c(%.2f, %.2f, %.2f) Focal Dist: %0.3f Env Rotation: %.2f Exposure: %.3f Camera Motion Blur: %s Shutter Speed: %.3f\n",
              origin.xyz.x, origin.xyz.y, origin.xyz.z,
              key_lookat.xyz.x, key_lookat.xyz.y, key_lookat.xyz.z,
              fd, env_rotation, preview_exposure_adjustment,
              camera_motion_blur_enabled ? "ON" : "OFF",
              shutter_speed);
    }
  }
}

std::string PreviewDisplay::PreviewStatusText(Float env_rotation) const {
  if(cam == nullptr) {
    return "";
  }

  point3f origin = cam->get_origin();
  point3f key_lookat = PreviewCameraLookat(cam);
  int keyframe_number = 0;
  if(current_keyframe >= 0 &&
     current_keyframe < static_cast<int>(Keyframes.size())) {
    keyframe_number = current_keyframe + 1;
  }

  char buffer[320];
  if(std::isinf(shutter_speed)) {
    std::snprintf(buffer,
                  sizeof(buffer),
                  "Loop %s | Cam x/y/z %.2f %.2f %.2f | Look %.2f %.2f %.2f | Exp %.3f | Shutter Inf | Env %.1f | Blur %s | Key %d/%zu",
                  keyframe_motion_closed ? "CLOSED" : "OPEN",
                  origin.xyz.x,
                  origin.xyz.y,
                  origin.xyz.z,
                  key_lookat.xyz.x,
                  key_lookat.xyz.y,
                  key_lookat.xyz.z,
                  preview_exposure_adjustment,
                  env_rotation,
                  camera_motion_blur_enabled ? "ON" : "OFF",
                  keyframe_number,
                  Keyframes.size());
  } else {
    std::snprintf(buffer,
                  sizeof(buffer),
                  "Loop %s | Cam x/y/z %.2f %.2f %.2f | Look %.2f %.2f %.2f | Exp %.3f | Shutter %.3f | Env %.1f | Blur %s | Key %d/%zu",
                  keyframe_motion_closed ? "CLOSED" : "OPEN",
                  origin.xyz.x,
                  origin.xyz.y,
                  origin.xyz.z,
                  key_lookat.xyz.x,
                  key_lookat.xyz.y,
                  key_lookat.xyz.z,
                  preview_exposure_adjustment,
                  shutter_speed,
                  env_rotation,
                  camera_motion_blur_enabled ? "ON" : "OFF",
                  keyframe_number,
                  Keyframes.size());
  }
  return std::string(buffer);
}

void PreviewDisplay::SetTextOverlays(const std::vector<PreviewTextOverlay>& overlays) {
  text_overlays = overlays;
}

void PreviewDisplay::SetLineOverlays(const std::vector<PreviewLineOverlay>& overlays) {
  line_overlays = overlays;
}

bool PreviewDisplay::ProjectWorldPoint(const point3f& point,
                                       bool clip,
                                       Float& screen_x,
                                       Float& screen_y,
                                       Float& depth) const {
  if(!cam) {
    return false;
  }
  unsigned int display_width = 0;
  unsigned int display_height = 0;
#ifdef RAY_HAS_X11
  display_width = width;
  display_height = height;
#endif
#ifdef RAY_WINDOWS
  display_width = ::width;
  display_height = ::height;
#endif
  if(display_width == 0 || display_height == 0) {
    return false;
  }
  Float fov = cam->get_fov();
  if(fov < 0) {
    return false;
  }
  point3f origin = cam->get_origin();
  vec3f relative = point - origin;
  Float rel_len2 = relative.squared_length();
  if(rel_len2 <= 0) {
    return false;
  }
  vec3f right = cam->get_u();
  vec3f up = cam->get_v();
  vec3f forward = cam->get_w();
  Float x_camera = dot(relative, right);
  Float y_camera = dot(relative, up);
  Float z_camera = dot(relative, forward);
  Float s = 0.5f;
  Float t = 0.5f;
  bool in_front = true;
  const Float pi_val = static_cast<Float>(3.14159265358979323846);

  if(fov == 0) {
    point2f ortho = cam->get_ortho();
    s = 0.5f + x_camera / ortho.xy.x;
    t = 0.5f + y_camera / ortho.xy.y;
    in_front = z_camera >= 0;
  } else if(fov == 360) {
    vec3f direction = relative / std::sqrt(rel_len2);
    Float local_x = dot(direction, right);
    Float local_y = dot(direction, up);
    Float local_z = dot(direction, forward);
    Float theta = std::acos(clamp(local_y, -1.f, 1.f));
    Float phi = std::atan2(local_x, local_z);
    s = std::fmod((phi - pi_val) /
                    (2.f * pi_val) + 1.f,
                  1.f);
    t = 1.f - theta / pi_val;
  } else {
    if(z_camera <= 0) {
      return false;
    }
    Float aspect = static_cast<Float>(display_width) /
      static_cast<Float>(display_height);
    Float half_height = std::tan(fov * pi_val / 360.f);
    Float half_width = aspect * half_height;
    s = 0.5f + x_camera / (2.f * z_camera * half_width);
    t = 0.5f + y_camera / (2.f * z_camera * half_height);
  }
  if(!in_front) {
    return false;
  }
  if(clip && (s < 0 || s > 1 || t < 0 || t > 1)) {
    return false;
  }
  screen_x = (1.f - s) * static_cast<Float>(display_width - 1);
  screen_y = (1.f - t) * static_cast<Float>(display_height - 1);
  depth = z_camera;
  return std::isfinite(screen_x) && std::isfinite(screen_y);
}

bool PreviewDisplay::ProjectTextAnchor(const PreviewTextOverlay& overlay,
                                       Float& screen_x,
                                       Float& screen_y) const {
  if(!cam) {
    return false;
  }
  unsigned int display_width = 0;
  unsigned int display_height = 0;
#ifdef RAY_HAS_X11
  display_width = width;
  display_height = height;
#endif
#ifdef RAY_WINDOWS
  display_width = ::width;
  display_height = ::height;
#endif
  if(display_width == 0 || display_height == 0) {
    return false;
  }
  Float fov = cam->get_fov();
  if(fov < 0) {
    return false;
  }
  point3f origin = cam->get_origin();
  vec3f relative = overlay.anchor - origin;
  Float rel_len2 = relative.squared_length();
  if(rel_len2 <= 0) {
    return false;
  }
  vec3f right = cam->get_u();
  vec3f up = cam->get_v();
  vec3f forward = cam->get_w();
  Float x_camera = dot(relative, right);
  Float y_camera = dot(relative, up);
  Float z_camera = dot(relative, forward);
  Float s = 0.5f;
  Float t = 0.5f;
  bool in_front = true;
  const Float pi_val = static_cast<Float>(3.14159265358979323846);

  if(fov == 0) {
    point2f ortho = cam->get_ortho();
    s = 0.5f + x_camera / ortho.xy.x;
    t = 0.5f + y_camera / ortho.xy.y;
    in_front = z_camera >= 0;
  } else if(fov == 360) {
    vec3f direction = relative / std::sqrt(rel_len2);
    Float local_x = dot(direction, right);
    Float local_y = dot(direction, up);
    Float local_z = dot(direction, forward);
    Float theta = std::acos(clamp(local_y, -1.f, 1.f));
    Float phi = std::atan2(local_x, local_z);
    s = std::fmod((phi - pi_val) /
                    (2.f * pi_val) + 1.f,
                  1.f);
    t = 1.f - theta / pi_val;
  } else {
    if(z_camera <= 0) {
      return false;
    }
    Float aspect = static_cast<Float>(display_width) /
      static_cast<Float>(display_height);
    Float half_height = std::tan(fov * pi_val / 360.f);
    Float half_width = aspect * half_height;
    s = 0.5f + x_camera / (2.f * z_camera * half_width);
    t = 0.5f + y_camera / (2.f * z_camera * half_height);
  }
  if(!in_front) {
    return false;
  }
  if(overlay.clip && (s < 0 || s > 1 || t < 0 || t > 1)) {
    return false;
  }
  screen_x = (1.f - s) * static_cast<Float>(display_width - 1);
  screen_y = (1.f - t) * static_cast<Float>(display_height - 1);
  return std::isfinite(screen_x) && std::isfinite(screen_y);
}

bool PreviewDisplay::IsTextAnchorOccluded(const PreviewTextOverlay& overlay,
                                          hitable* world,
                                          random_gen& rng) const {
  if(!overlay.occlusion || !cam || !world) {
    return false;
  }
  point3f origin = cam->get_origin();
  vec3f relative = overlay.anchor - origin;
  Float distance_to_anchor = relative.length();
  Float occlusion_tolerance = overlay.occlusion_tolerance;
  if(occlusion_tolerance < 0.f) {
    occlusion_tolerance = 0.f;
  }
  Float endpoint_tolerance = occlusion_tolerance < 1.f ?
    distance_to_anchor * occlusion_tolerance : occlusion_tolerance;
  if(distance_to_anchor <= endpoint_tolerance) {
    return false;
  }
  vec3f direction = relative / distance_to_anchor;
  Float t_max = distance_to_anchor - endpoint_tolerance;

  if(cam->get_fov() == 0) {
    vec3f right = cam->get_u();
    vec3f up = cam->get_v();
    vec3f forward = cam->get_w();
    Float x_camera = dot(relative, right);
    Float y_camera = dot(relative, up);
    Float z_camera = dot(relative, forward);
    endpoint_tolerance = occlusion_tolerance < 1.f ?
      z_camera * occlusion_tolerance : occlusion_tolerance;
    if(z_camera <= endpoint_tolerance) {
      return false;
    }
    origin = origin + x_camera * right + y_camera * up;
    direction = forward;
    t_max = z_camera - endpoint_tolerance;
  }

  hit_record hrec;
  Ray visibility_ray(origin, direction, 0.5f);
  return world->hit(visibility_ray, 0.001f, t_max, hrec, rng) &&
    hrec.shape->GetName() != "EnvironmentLight";
}

bool PreviewDisplay::IsTextPixelOccluded(const PreviewTextOverlay& overlay,
                                         Float screen_x,
                                         Float screen_y,
                                         hitable* world,
                                         random_gen& rng) const {
  if(!overlay.partial_occlusion || !cam || !world) {
    return false;
  }
  unsigned int display_width = 0;
  unsigned int display_height = 0;
#ifdef RAY_HAS_X11
  display_width = width;
  display_height = height;
#endif
#ifdef RAY_WINDOWS
  display_width = ::width;
  display_height = ::height;
#endif
  if(display_width == 0 || display_height == 0 || cam->get_fov() < 0) {
    return false;
  }

  vec3f forward = cam->get_w();
  Float label_depth = dot(overlay.anchor - cam->get_origin(), forward);
  Float occlusion_tolerance = overlay.occlusion_tolerance;
  if(occlusion_tolerance < 0.f) {
    occlusion_tolerance = 0.f;
  }
  Float endpoint_tolerance = occlusion_tolerance < 1.f ?
    label_depth * occlusion_tolerance : occlusion_tolerance;
  if(label_depth <= endpoint_tolerance) {
    return false;
  }

  Float denom_x = static_cast<Float>(std::max(1u, display_width - 1));
  Float denom_y = static_cast<Float>(std::max(1u, display_height - 1));
  Float s = 1.f - screen_x / denom_x;
  Float t = 1.f - screen_y / denom_y;
  Ray visibility_ray = cam->get_ray(s, t, point3f(0.f, 0.f, 0.f), 0.5f);

  hit_record hrec;
  if(!world->hit(visibility_ray, 0.001f, MaxT, hrec, rng) ||
     hrec.shape->GetName() == "EnvironmentLight") {
    return false;
  }
  Float scene_depth = dot(hrec.p - cam->get_origin(), forward);
  return scene_depth < label_depth - endpoint_tolerance;
}

bool PreviewDisplay::IsLineAnchorOccluded(const PreviewLineOverlay& overlay,
                                          hitable* world,
                                          random_gen& rng) const {
  if(!overlay.occlusion || !cam || !world) {
    return false;
  }
  point3f anchor = overlay.start + (overlay.end - overlay.start) * 0.5f;
  point3f origin = cam->get_origin();
  vec3f relative = anchor - origin;
  Float distance_to_anchor = relative.length();
  Float occlusion_tolerance = overlay.occlusion_tolerance;
  if(occlusion_tolerance < 0.f) {
    occlusion_tolerance = 0.f;
  }
  Float endpoint_tolerance = occlusion_tolerance < 1.f ?
    distance_to_anchor * occlusion_tolerance : occlusion_tolerance;
  if(distance_to_anchor <= endpoint_tolerance) {
    return false;
  }
  vec3f direction = relative / distance_to_anchor;
  Float t_max = distance_to_anchor - endpoint_tolerance;

  if(cam->get_fov() == 0) {
    vec3f right = cam->get_u();
    vec3f up = cam->get_v();
    vec3f forward = cam->get_w();
    Float x_camera = dot(relative, right);
    Float y_camera = dot(relative, up);
    Float z_camera = dot(relative, forward);
    endpoint_tolerance = occlusion_tolerance < 1.f ?
      z_camera * occlusion_tolerance : occlusion_tolerance;
    if(z_camera <= endpoint_tolerance) {
      return false;
    }
    origin = origin + x_camera * right + y_camera * up;
    direction = forward;
    t_max = z_camera - endpoint_tolerance;
  }

  hit_record hrec;
  Ray visibility_ray(origin, direction, 0.5f);
  return world->hit(visibility_ray, 0.001f, t_max, hrec, rng) &&
    hrec.shape->GetName() != "EnvironmentLight";
}

bool PreviewDisplay::IsLinePixelOccluded(const PreviewLineOverlay& overlay,
                                         Float screen_x,
                                         Float screen_y,
                                         Float line_depth,
                                         hitable* world,
                                         random_gen& rng) const {
  if(!overlay.partial_occlusion || !cam || !world) {
    return false;
  }
  unsigned int display_width = 0;
  unsigned int display_height = 0;
#ifdef RAY_HAS_X11
  display_width = width;
  display_height = height;
#endif
#ifdef RAY_WINDOWS
  display_width = ::width;
  display_height = ::height;
#endif
  if(display_width == 0 || display_height == 0 || cam->get_fov() < 0) {
    return false;
  }
  Float occlusion_tolerance = overlay.occlusion_tolerance;
  if(occlusion_tolerance < 0.f) {
    occlusion_tolerance = 0.f;
  }
  Float endpoint_tolerance = occlusion_tolerance < 1.f ?
    line_depth * occlusion_tolerance : occlusion_tolerance;
  if(line_depth <= endpoint_tolerance) {
    return false;
  }

  Float denom_x = static_cast<Float>(std::max(1u, display_width - 1));
  Float denom_y = static_cast<Float>(std::max(1u, display_height - 1));
  Float s = 1.f - screen_x / denom_x;
  Float t = 1.f - screen_y / denom_y;
  Ray visibility_ray = cam->get_ray(s, t, point3f(0.f, 0.f, 0.f), 0.5f);

  hit_record hrec;
  if(!world->hit(visibility_ray, 0.001f, MaxT, hrec, rng) ||
     hrec.shape->GetName() == "EnvironmentLight") {
    return false;
  }
  Float scene_depth = dot(hrec.p - cam->get_origin(), cam->get_w());
  return scene_depth < line_depth - endpoint_tolerance;
}

static Float LineCoverageAndDepth(Float x0,
                                  Float y0,
                                  Float depth0,
                                  Float x1,
                                  Float y1,
                                  Float depth1,
                                  Float px,
                                  Float py,
                                  Float width,
                                  int lineend,
                                  Float& line_depth) {
  Float radius = width / 2.f;
  Float dx = x1 - x0;
  Float dy = y1 - y0;
  Float len2 = dx * dx + dy * dy;
  if(len2 <= static_cast<Float>(1e-12)) {
    Float dist = std::sqrt((px - x0) * (px - x0) + (py - y0) * (py - y0));
    line_depth = (depth0 + depth1) * 0.5f;
    return clamp(radius + 0.5f - dist, 0.f, 1.f);
  }

  Float t_depth = ((px - x0) * dx + (py - y0) * dy) / len2;
  line_depth = depth0 + clamp(t_depth, 0.f, 1.f) * (depth1 - depth0);
  Float cx0 = x0;
  Float cy0 = y0;
  Float cx1 = x1;
  Float cy1 = y1;

  if(lineend == 2) {
    Float len = std::sqrt(len2);
    Float ux = dx / len;
    Float uy = dy / len;
    cx0 -= ux * radius;
    cy0 -= uy * radius;
    cx1 += ux * radius;
    cy1 += uy * radius;
    dx = cx1 - cx0;
    dy = cy1 - cy0;
    len2 = dx * dx + dy * dy;
  }

  Float t = ((px - cx0) * dx + (py - cy0) * dy) / len2;
  if(lineend == 1 && (t < 0.f || t > 1.f)) {
    return 0.f;
  }
  t = clamp(t, 0.f, 1.f);
  Float closest_x = cx0 + t * dx;
  Float closest_y = cy0 + t * dy;
  Float dist = std::sqrt((px - closest_x) * (px - closest_x) +
                         (py - closest_y) * (py - closest_y));
  return clamp(radius + 0.5f - dist, 0.f, 1.f);
}

#ifdef RAY_HAS_X11
void PreviewDisplay::CompositeTextOverlaysToX11Buffer(hitable* world, random_gen& rng) {
  if(text_overlays.empty() || !data) {
    return;
  }
  for(const PreviewTextOverlay& overlay : text_overlays) {
    Float screen_x;
    Float screen_y;
    if(!ProjectTextAnchor(overlay, screen_x, screen_y)) {
      continue;
    }
    if(IsTextAnchorOccluded(overlay, world, rng)) {
      continue;
    }
    int overlay_width = static_cast<int>(overlay.width);
    int overlay_height = static_cast<int>(overlay.height);
    int left = static_cast<int>(std::round(screen_x + overlay.x_offset -
                                           overlay.hjust * overlay_width));
    int top = static_cast<int>(std::round(screen_y + overlay.y_offset -
                                          overlay.vjust * overlay_height));
    int right_bound = std::min<int>(left + overlay_width, width);
    int bottom_bound = std::min<int>(top + overlay_height, height);
    int x_start = std::max(0, left);
    int y_start = std::max(0, top);
    if(x_start >= right_bound || y_start >= bottom_bound) {
      continue;
    }
    for(int y = y_start; y < bottom_bound; y++) {
      unsigned int source_y = static_cast<unsigned int>(y - top);
      for(int x = x_start; x < right_bound; x++) {
        unsigned int source_x = static_cast<unsigned int>(x - left);
        size_t src_idx = 4 * (source_x + overlay.width * source_y);
        Float alpha = overlay.rgba[src_idx + 3] / 255.f;
        if(alpha <= 0) {
          continue;
        }
        if(IsTextPixelOccluded(overlay,
                               static_cast<Float>(x),
                               static_cast<Float>(y),
                               world,
                               rng)) {
          continue;
        }
        size_t dst_idx = 4 * (x + width * y);
        Float src_r = overlay.rgba[src_idx] / 255.f;
        Float src_g = overlay.rgba[src_idx + 1] / 255.f;
        Float src_b = overlay.rgba[src_idx + 2] / 255.f;
        Float dst_b = static_cast<unsigned char>(data[dst_idx]) / 255.f;
        Float dst_g = static_cast<unsigned char>(data[dst_idx + 1]) / 255.f;
        Float dst_r = static_cast<unsigned char>(data[dst_idx + 2]) / 255.f;
        data[dst_idx] = static_cast<unsigned char>(
          255.f * clamp(src_b * alpha + dst_b * (1.f - alpha), 0.f, 1.f)
        );
        data[dst_idx + 1] = static_cast<unsigned char>(
          255.f * clamp(src_g * alpha + dst_g * (1.f - alpha), 0.f, 1.f)
        );
        data[dst_idx + 2] = static_cast<unsigned char>(
          255.f * clamp(src_r * alpha + dst_r * (1.f - alpha), 0.f, 1.f)
        );
      }
    }
  }
}

void PreviewDisplay::CompositeLineOverlaysToX11Buffer(hitable* world, random_gen& rng) {
  if(line_overlays.empty() || !data) {
    return;
  }
  for(const PreviewLineOverlay& overlay : line_overlays) {
    Float x0;
    Float y0;
    Float depth0;
    Float x1;
    Float y1;
    Float depth1;
    if(!ProjectWorldPoint(overlay.start, overlay.clip, x0, y0, depth0) ||
       !ProjectWorldPoint(overlay.end, overlay.clip, x1, y1, depth1)) {
      continue;
    }
    x0 += overlay.x_offset;
    y0 += overlay.y_offset;
    x1 += overlay.xend_offset;
    y1 += overlay.yend_offset;
    if(IsLineAnchorOccluded(overlay, world, rng)) {
      continue;
    }
    Float pad = overlay.width / 2.f + 1.f;
    int left = std::max(0, static_cast<int>(std::floor(std::min(x0, x1) - pad)));
    int right = std::min<int>(width - 1, static_cast<int>(std::ceil(std::max(x0, x1) + pad)));
    int top = std::max(0, static_cast<int>(std::floor(std::min(y0, y1) - pad)));
    int bottom = std::min<int>(height - 1, static_cast<int>(std::ceil(std::max(y0, y1) + pad)));
    if(left > right || top > bottom) {
      continue;
    }
    for(int y = top; y <= bottom; y++) {
      for(int x = left; x <= right; x++) {
        Float line_depth;
        Float alpha = overlay.alpha * LineCoverageAndDepth(
          x0, y0, depth0, x1, y1, depth1,
          static_cast<Float>(x), static_cast<Float>(y),
          overlay.width, overlay.lineend, line_depth
        );
        if(alpha <= 0.f) {
          continue;
        }
        if(IsLinePixelOccluded(overlay,
                               static_cast<Float>(x),
                               static_cast<Float>(y),
                               line_depth,
                               world,
                               rng)) {
          continue;
        }
        size_t dst_idx = 4 * (x + width * y);
        Float dst_b = static_cast<unsigned char>(data[dst_idx]) / 255.f;
        Float dst_g = static_cast<unsigned char>(data[dst_idx + 1]) / 255.f;
        Float dst_r = static_cast<unsigned char>(data[dst_idx + 2]) / 255.f;
        data[dst_idx] = static_cast<unsigned char>(
          255.f * clamp(overlay.blue * alpha + dst_b * (1.f - alpha), 0.f, 1.f)
        );
        data[dst_idx + 1] = static_cast<unsigned char>(
          255.f * clamp(overlay.green * alpha + dst_g * (1.f - alpha), 0.f, 1.f)
        );
        data[dst_idx + 2] = static_cast<unsigned char>(
          255.f * clamp(overlay.red * alpha + dst_r * (1.f - alpha), 0.f, 1.f)
        );
      }
    }
  }
}

void PreviewDisplay::DrawStatusBarX11(Float env_rotation) {
  if(!interactive ||
     width < PREVIEW_STATUS_MIN_WIDTH ||
     height <= PREVIEW_STATUS_BAR_HEIGHT) {
    return;
  }

  std::string status_text = PreviewStatusText(env_rotation);
  if(status_text.empty()) {
    return;
  }

  GC gc = DefaultGC(d, s);
  XSetForeground(d, gc, BlackPixel(d, s));
  XFillRectangle(d,
                 w,
                 gc,
                 0,
                 height - PREVIEW_STATUS_BAR_HEIGHT,
                 width,
                 PREVIEW_STATUS_BAR_HEIGHT);
  XSetForeground(d, gc, WhitePixel(d, s));
  XDrawString(d,
              w,
              gc,
              8,
              height - 7,
              status_text.c_str(),
              static_cast<int>(status_text.size()));
  XFlush(d);
}
#endif

void PreviewDisplay::CompositeTextOverlaysToFloatBuffer(std::vector<Float>& rgb,
                                                        hitable* world,
                                                        random_gen& rng, std::vector<Float>* coverage) {
  if(text_overlays.empty() || rgb.empty()) {
    return;
  }
  for(const PreviewTextOverlay& overlay : text_overlays) {
    Float screen_x;
    Float screen_y;
    if(!ProjectTextAnchor(overlay, screen_x, screen_y)) {
      continue;
    }
    if(IsTextAnchorOccluded(overlay, world, rng)) {
      continue;
    }
    int overlay_width = static_cast<int>(overlay.width);
    int overlay_height = static_cast<int>(overlay.height);
    int left = static_cast<int>(std::round(screen_x + overlay.x_offset -
                                           overlay.hjust * overlay_width));
    int top = static_cast<int>(std::round(screen_y + overlay.y_offset -
                                          overlay.vjust * overlay_height));
    int right_bound = std::min<int>(left + overlay_width, width);
    int bottom_bound = std::min<int>(top + overlay_height, height);
    int x_start = std::max(0, left);
    int y_start = std::max(0, top);
    if(x_start >= right_bound || y_start >= bottom_bound) {
      continue;
    }
    for(int y = y_start; y < bottom_bound; y++) {
      unsigned int source_y = static_cast<unsigned int>(y - top);
      for(int x = x_start; x < right_bound; x++) {
        unsigned int source_x = static_cast<unsigned int>(x - left);
        size_t src_idx = 4 * (source_x + overlay.width * source_y);
        Float alpha = overlay.rgba[src_idx + 3] / 255.f;
        if(alpha <= 0) {
          continue;
        }
        if(IsTextPixelOccluded(overlay,
                               static_cast<Float>(x),
                               static_cast<Float>(y),
                               world,
                               rng)) {
          continue;
        }
        size_t dst_idx = 3 * (x + width * y);
        Float src_r = overlay.rgba[src_idx] / 255.f;
        Float src_g = overlay.rgba[src_idx + 1] / 255.f;
        Float src_b = overlay.rgba[src_idx + 2] / 255.f;
        Float dest_alpha=coverage ? (*coverage)[dst_idx/3] : 1;
        Float out_alpha=alpha+dest_alpha*(1-alpha);
        Float src_weight=out_alpha>0 ? alpha/out_alpha : 0;
        Float dst_weight=out_alpha>0 ? dest_alpha*(1-alpha)/out_alpha : 0;
        if(coverage) (*coverage)[dst_idx/3]=out_alpha;
        rgb[dst_idx] = clamp(src_r * src_weight + rgb[dst_idx] * dst_weight, 0.f, 1.f);
        rgb[dst_idx + 1] = clamp(src_g * src_weight + rgb[dst_idx + 1] * dst_weight, 0.f, 1.f);
        rgb[dst_idx + 2] = clamp(src_b * src_weight + rgb[dst_idx + 2] * dst_weight, 0.f, 1.f);
      }
    }
  }
}

void PreviewDisplay::CompositeLineOverlaysToFloatBuffer(std::vector<Float>& rgb,
                                                        hitable* world,
                                                        random_gen& rng, std::vector<Float>* coverage) {
  if(line_overlays.empty() || rgb.empty()) {
    return;
  }
  for(const PreviewLineOverlay& overlay : line_overlays) {
    Float x0;
    Float y0;
    Float depth0;
    Float x1;
    Float y1;
    Float depth1;
    if(!ProjectWorldPoint(overlay.start, overlay.clip, x0, y0, depth0) ||
       !ProjectWorldPoint(overlay.end, overlay.clip, x1, y1, depth1)) {
      continue;
    }
    x0 += overlay.x_offset;
    y0 += overlay.y_offset;
    x1 += overlay.xend_offset;
    y1 += overlay.yend_offset;
    if(IsLineAnchorOccluded(overlay, world, rng)) {
      continue;
    }
    Float pad = overlay.width / 2.f + 1.f;
    int left = std::max(0, static_cast<int>(std::floor(std::min(x0, x1) - pad)));
    int right = std::min<int>(width - 1, static_cast<int>(std::ceil(std::max(x0, x1) + pad)));
    int top = std::max(0, static_cast<int>(std::floor(std::min(y0, y1) - pad)));
    int bottom = std::min<int>(height - 1, static_cast<int>(std::ceil(std::max(y0, y1) + pad)));
    if(left > right || top > bottom) {
      continue;
    }
    for(int y = top; y <= bottom; y++) {
      for(int x = left; x <= right; x++) {
        Float line_depth;
        Float alpha = overlay.alpha * LineCoverageAndDepth(
          x0, y0, depth0, x1, y1, depth1,
          static_cast<Float>(x), static_cast<Float>(y),
          overlay.width, overlay.lineend, line_depth
        );
        if(alpha <= 0.f) {
          continue;
        }
        if(IsLinePixelOccluded(overlay,
                               static_cast<Float>(x),
                               static_cast<Float>(y),
                               line_depth,
                               world,
                               rng)) {
          continue;
        }
        size_t dst_idx = 3 * (x + width * y);
        Float dest_alpha=coverage ? (*coverage)[dst_idx/3] : 1;
        Float out_alpha=alpha+dest_alpha*(1-alpha);
        Float src_weight=out_alpha>0 ? alpha/out_alpha : 0;
        Float dst_weight=out_alpha>0 ? dest_alpha*(1-alpha)/out_alpha : 0;
        if(coverage) (*coverage)[dst_idx/3]=out_alpha;
        rgb[dst_idx] = clamp(overlay.red * src_weight + rgb[dst_idx] * dst_weight, 0.f, 1.f);
        rgb[dst_idx + 1] = clamp(overlay.green * src_weight + rgb[dst_idx + 1] * dst_weight, 0.f, 1.f);
        rgb[dst_idx + 2] = clamp(overlay.blue * src_weight + rgb[dst_idx + 2] * dst_weight, 0.f, 1.f);
      }
    }
  }
}

#ifdef RAY_WINDOWS
void PreviewDisplay::DrawStatusBarWindows(HDC hdc, Float env_rotation) const {
  if(!interactive ||
     width < PREVIEW_STATUS_MIN_WIDTH ||
     height <= PREVIEW_STATUS_BAR_HEIGHT) {
    return;
  }

  std::string status_text = PreviewStatusText(env_rotation);
  if(status_text.empty()) {
    return;
  }

  RECT bar_rect = {0,
                   static_cast<LONG>(height - PREVIEW_STATUS_BAR_HEIGHT),
                   static_cast<LONG>(width),
                   static_cast<LONG>(height)};
  HBRUSH brush = CreateSolidBrush(RGB(0, 0, 0));
  FillRect(hdc, &bar_rect, brush);
  DeleteObject(brush);

  RECT text_rect = bar_rect;
  text_rect.left += 8;
  SetBkMode(hdc, TRANSPARENT);
  SetTextColor(hdc, RGB(255, 255, 255));
  DrawTextA(hdc,
            status_text.c_str(),
            -1,
            &text_rect,
            DT_SINGLELINE | DT_VCENTER | DT_LEFT | DT_END_ELLIPSIS);
}
#endif

void PreviewDisplay::DrawImage(adaptive_sampler& adaptive_pixel_sampler,
                               adaptive_sampler& adaptive_pixel_sampler_small,
                               size_t &ns, RProgress::RProgress &pb, bool progress,
                               Float percent_done,
                               hitable *world, random_gen& rng) {
  SCOPED_CONTEXT("Overall");
  SCOPED_TIMER_COUNTER("Draw Image");
  auto reset_preview_render = [&]() {
    ns = 0;
    adaptive_pixel_sampler.reset();
    adaptive_pixel_sampler_small.reset();
    ResetPreviewExposure();
#ifdef HAS_OIDN
    InvalidateOidnAux();
#endif
    if(progress && !interactive) {
      pb.update(0);
    }
  };
  if(PollCloseEvent()) {
    return;
  }
#ifdef RAY_HAS_X11
  if (d) {
#ifdef HAS_OIDN
    bool use_denoised_preview = denoise &&
      denoiser != nullptr &&
      denoiser->Ready();
    if(use_denoised_preview && !write_fast_output) {
      if(PollCloseEvent()) {
        return;
      }
      denoiser->Execute();
      if(!denoiser->ReportError()) {
        MarkDenoisedPreviewReady(ns + 1);
      }
    }
    RayMatrix& rgb = use_denoised_preview ? adaptive_pixel_sampler.draw_rgb_output : adaptive_pixel_sampler.rgb;
#else
    RayMatrix &rgb  = adaptive_pixel_sampler.rgb;
#endif
    CalibratePreviewExposure(adaptive_pixel_sampler, rgb, ns);
    std::vector<bool>& finalized = adaptive_pixel_sampler.finalized;
    std::vector<bool>& just_finalized = adaptive_pixel_sampler.just_finalized;
    
    for(unsigned int i = 0; i < 4*width; i += 4 ) {
      for(unsigned int j = 0; j < height; j++) {
        int samples;
        Float r_col,g_col,b_col;
        if(finalized[(width - 1 - i/4) + width * (height-1-j)]) {
          samples = 1;
          if (just_finalized[(width - 1 - i / 4) + width * (height - 1 - j)] &&
              !interactive) {
            r_col = 0;
            g_col = 1;
            b_col = 0;
          } else {
            Float sample_count = interactive ? 1.f : 4.f;
            r_col = ApplyPreviewExposure(rgb((width - 1 - i/4),height-1-j,0), sample_count);
            g_col = ApplyPreviewExposure(rgb((width - 1 - i/4),height-1-j,1), sample_count);
            b_col = ApplyPreviewExposure(rgb((width - 1 - i/4),height-1-j,2), sample_count);
          }
        } else {
          samples = ns+1;
          r_col = ApplyPreviewExposure(rgb((width - 1 - i/4),height-1-j,0), samples);
          g_col = ApplyPreviewExposure(rgb((width - 1 - i/4),height-1-j,1), samples);
          b_col = ApplyPreviewExposure(rgb((width - 1 - i/4),height-1-j,2), samples);
        }

        data[i + 4*width*j]   = (unsigned char)(255*clamp(b_col,0,1));
        data[i + 4*width*j+1] = (unsigned char)(255*clamp(g_col,0,1));
        data[i + 4*width*j+2] = (unsigned char)(255*clamp(r_col,0,1));

        if(finalized[(width - 1 - i/4) + width * (height-1-j)]) {
          just_finalized[(width - 1 - i/4) + width * (height-1-j)] = false;
        }
      }
    }
    CompositeTextOverlaysToX11Buffer(world, rng);
    CompositeLineOverlaysToX11Buffer(world, rng);
    snapshot_width = width;
    snapshot_height = height;
    snapshot_pixels.resize(static_cast<size_t>(width) *
                           static_cast<size_t>(height) * 3);
    for(size_t pixel = 0;
        pixel < static_cast<size_t>(width) * static_cast<size_t>(height);
        pixel++) {
      snapshot_pixels[3 * pixel] =
        static_cast<unsigned char>(data[4 * pixel + 2]);
      snapshot_pixels[3 * pixel + 1] =
        static_cast<unsigned char>(data[4 * pixel + 1]);
      snapshot_pixels[3 * pixel + 2] =
        static_cast<unsigned char>(data[4 * pixel]);
    }
    CaptureVolumeSnapshot(adaptive_pixel_sampler, rgb, ns+1, world, rng);
    // Paint display-only UI after populating the snapshot buffer.
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

    //Preview exposure control
    KeyCode RightBracket_key = XKeysymToKeycode(d, XK_bracketright);
    KeyCode LeftBracket_key = XKeysymToKeycode(d, XK_bracketleft);
    
    //Reset
    KeyCode R_key = XKeysymToKeycode(d,XStringToKeysym("r"));

    //Toggle deferred/final render
    KeyCode Return_key = XKeysymToKeycode(d,XK_Return);
    KeyCode KeypadEnter_key = XKeysymToKeycode(d,XK_KP_Enter);
    
    //Print Position
    KeyCode P_key = XKeysymToKeycode(d,XStringToKeysym("p"));
    
    //Save Keyframe
    KeyCode K_key = XKeysymToKeycode(d,XStringToKeysym("k"));
    
    //Move to Last Keyframe
    KeyCode L_key = XKeysymToKeycode(d,XStringToKeysym("l"));

    //Keyframe navigation
    KeyCode PreviousKeyframe_key = XKeysymToKeycode(d, XK_less);
    KeyCode NextKeyframe_key = XKeysymToKeycode(d, XK_greater);
    KeyCode DeleteKeyframe_key = XKeysymToKeycode(d, XK_slash);
    KeyCode MotionPreview_key = XKeysymToKeycode(d,XStringToKeysym("m"));
    KeyCode CameraMotionBlur_key = XKeysymToKeycode(d,XStringToKeysym("b"));
    
    //Fast Movement Key
    KeyCode F_key = XKeysymToKeycode(d,XStringToKeysym("f"));
    
    
    XPutImage(d,w,DefaultGC(d,s),
              img,0,0,0,0,width,height);
    DrawStatusBarX11(env_y_angle);
    while (XPending(d)) {
      XNextEvent(d, &e);
      if (e.type == KeyPress) {
        if(!PreviewDisplayHasKeyboardFocus(d, this->w)) {
          continue;
        }
        if (e.xkey.keycode == esc ) {
          terminate = true;
          break;
        }
        if (e.xkey.keycode == CameraMotionBlur_key ) {
          ToggleCameraMotionBlur();
          reset_preview_render();
          continue;
        }
        if(interactive &&
           (e.xkey.keycode == Return_key ||
            e.xkey.keycode == KeypadEnter_key) &&
           (e.xkey.state & ShiftMask) != 0) {
          SavePreviewSnapshot();
          continue;
        }
        if(interactive && IsPreviewMotionActive()) {
          if(e.xkey.keycode == MotionPreview_key) {
            CancelPreviewMotion(&env_y_angle);
            reset_preview_render();
          }
          continue;
        }
        if(interactive) {
          vec3f w = cam->get_w();
          vec3f u = cam->get_u();
          vec3f v = cam->get_v();
        
          bool blanked = false;
          bool one_orbit = false;
          bool one_fast = false;
          bool shift_pressed = (e.xkey.state & ShiftMask) != 0;
          if(IsX11ShiftedLeftBracket(d, e.xkey)) {
            AdjustShutterSpeedStops(static_cast<Float>(-1) / static_cast<Float>(3));
            reset_preview_render();
            continue;
          }
          if(IsX11ShiftedRightBracket(d, e.xkey)) {
            AdjustShutterSpeedStops(static_cast<Float>(1) / static_cast<Float>(3));
            reset_preview_render();
            continue;
          }
          
          if (e.xkey.keycode == tab ) {
            orbit = !orbit;
            one_orbit = true;
          }
          if (e.xkey.keycode == F_key ) {
            write_fast_output = !write_fast_output;
            one_fast  = true;
          }
          if (e.xkey.keycode == W_key ) {
            if(shift_pressed) {
              cam->rotate_forward(speed * 1.f);
            } else {
              vec3f step = speed * w * base_step;
              if(orbit) {
                Float dist_to_orbit = (cam->get_origin() - cam->get_lookat()).length();
                if(dist_to_orbit <= base_step * speed) {
                  Rprintf("Moving forward will overstep orbit point, stopping (decrease step size to move closer).\n");
                  step = vec3f(0);
                }
              }
              cam->update_position(step, orbit, false);
            }
          }
          if (e.xkey.keycode == A_key ) {
            if(shift_pressed) {
              cam->rotate_up(speed * -1.f);
            } else {
              cam->update_position(speed * u * base_step, orbit);
            }
          }
          if (e.xkey.keycode == S_key ) {
            if(shift_pressed) {
              cam->rotate_forward(speed * -1.f);
            } else {
              cam->update_position(-speed * w * base_step, orbit, false);
            }
          }
          if (e.xkey.keycode == D_key ) {
            if(shift_pressed) {
              cam->rotate_up(speed * 1.f);
            } else {
              cam->update_position(-speed * u * base_step, orbit);
            }
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
          if (e.xkey.keycode == RightBracket_key ) {
            IncreasePreviewExposure();
          }
          if (e.xkey.keycode == LeftBracket_key ) {
            DecreasePreviewExposure();
          }
          if (e.xkey.keycode == Return_key || e.xkey.keycode == KeypadEnter_key ) {
            if(deferred_render) {
              render_requested = !render_requested;
              if(render_requested) {
                Rprintf("Starting final render...\n");
              } else {
                Rprintf("Returning to deferred render...\n");
              }
            }
          }
          if (e.xkey.keycode == R_key ) {
            cam->reset();
            speed = 1;
            (*EnvWorldToObject) = Start_EnvWorldToObject;
            (*EnvObjectToWorld) = Start_EnvObjectToWorld;
            env_y_angle = 0;
          }
          if (e.xkey.keycode == L_key) {
            if(shift_pressed) {
              ToggleKeyframeMotionClosed();
              DrawStatusBarX11(env_y_angle);
              continue;
            } else if(Keyframes.size() > 0) {
              ApplyKeyframe(static_cast<int>(Keyframes.size()) - 1, &env_y_angle);
            } else {
              Rprintf("Can't reset to last keyframe: No keyframes have been saved. Use the R key to reset camera.");
            }
          }
          if (e.xkey.keycode == PreviousKeyframe_key ) {
            JumpKeyframe(-1, &env_y_angle);
          }
          if (e.xkey.keycode == NextKeyframe_key ) {
            JumpKeyframe(1, &env_y_angle);
          }
          if (e.xkey.keycode == DeleteKeyframe_key ) {
            DeleteCurrentKeyframe(&env_y_angle);
          }
          if (e.xkey.keycode == MotionPreview_key ) {
            if(StartPreviewMotion(env_y_angle)) {
              blanked = true;
              reset_preview_render();
            }
          }
          if (e.xkey.keycode == P_key || e.xkey.keycode == K_key ) {
            if(e.xkey.keycode == K_key) {
              SaveCurrentKeyframe(env_y_angle);
            }
            PrintCameraInfo(env_y_angle);
            
          } else {
            if(!blanked && !terminate && IsX11RenderInvalidatingKey(d, e.xkey.keycode)) {
              ApplyStaticPreviewCameraMotionRange(cam);
              blanked = true;
              ns = 0;
              adaptive_pixel_sampler.reset();
              adaptive_pixel_sampler_small.reset();
              ResetPreviewExposure();
#ifdef HAS_OIDN
              InvalidateOidnAux();
#endif
              if(progress && !interactive) {
                pb.update(0);
              }
            }
          }
          while(XPending(d)) {
            XNextEvent(d, &e);
            if (e.type != KeyPress) {
              XPutBackEvent(d, &e);
              break;
            }
            if(!PreviewDisplayHasKeyboardFocus(d, this->w)) {
              continue;
            }
            if (e.xkey.keycode == esc ) {
              terminate = true;
            }
            if (e.xkey.keycode == CameraMotionBlur_key ) {
              ToggleCameraMotionBlur();
              reset_preview_render();
              continue;
            }
            if(interactive &&
               (e.xkey.keycode == Return_key ||
                e.xkey.keycode == KeypadEnter_key) &&
               (e.xkey.state & ShiftMask) != 0) {
              SavePreviewSnapshot();
              continue;
            }
            if(interactive && IsPreviewMotionActive()) {
              if(e.xkey.keycode == MotionPreview_key) {
                CancelPreviewMotion(&env_y_angle);
                reset_preview_render();
              }
              continue;
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
            shift_pressed = (e.xkey.state & ShiftMask) != 0;
            if(IsX11ShiftedLeftBracket(d, e.xkey)) {
              AdjustShutterSpeedStops(static_cast<Float>(-1) / static_cast<Float>(3));
              reset_preview_render();
              continue;
            }
            if(IsX11ShiftedRightBracket(d, e.xkey)) {
              AdjustShutterSpeedStops(static_cast<Float>(1) / static_cast<Float>(3));
              reset_preview_render();
              continue;
            }
            
            if (e.xkey.keycode == W_key ) {
              if(shift_pressed) {
                cam->rotate_forward(speed * 1.f);
              } else {
                vec3f step = speed * w * base_step;
                if(orbit) {
                  Float dist_to_orbit = (cam->get_origin() - cam->get_lookat()).length();
                  if(dist_to_orbit <= base_step * speed) {
                    Rprintf("Moving forward will overstep orbit point, stopping (decrease step size to move closer).\n");
                    step = vec3f(0);
                  }
                }
                cam->update_position(step, orbit, false);
              }
            }
            if (e.xkey.keycode == A_key ) {
              if(shift_pressed) {
                cam->rotate_up(speed * -1.f);
              } else {
                cam->update_position(speed * u * base_step, orbit);
              }
            }
            if (e.xkey.keycode == S_key ) {
              if(shift_pressed) {
                cam->rotate_forward(speed * -1.f);
              } else {
                cam->update_position(-speed * w * base_step, orbit, false);
              }
            }
            if (e.xkey.keycode == D_key ) {
              if(shift_pressed) {
                cam->rotate_up(speed * 1.f);
              } else {
                cam->update_position(-speed * u * base_step, orbit);
              }
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
            if (e.xkey.keycode == RightBracket_key ) {
              IncreasePreviewExposure();
            }
            if (e.xkey.keycode == LeftBracket_key ) {
              DecreasePreviewExposure();
            }
            if (e.xkey.keycode == Return_key || e.xkey.keycode == KeypadEnter_key ) {
              if(deferred_render) {
                render_requested = !render_requested;
                if(render_requested) {
                  Rprintf("Starting final render...\n");
                } else {
                  Rprintf("Returning to deferred render...\n");
                }
              }
            }
            if (e.xkey.keycode == R_key ) {
              cam->reset();
              speed = 1;
              (*EnvWorldToObject) = Start_EnvWorldToObject;
              (*EnvObjectToWorld) = Start_EnvObjectToWorld;
              env_y_angle = 0;
            }
            if (e.xkey.keycode == L_key) {
              if(shift_pressed) {
                ToggleKeyframeMotionClosed();
                DrawStatusBarX11(env_y_angle);
                continue;
              } else if(Keyframes.size() > 0) {
                ApplyKeyframe(static_cast<int>(Keyframes.size()) - 1, &env_y_angle);
              } else {
                Rprintf("Can't reset to last keyframe: No keyframes have been saved. Use the R key to reset camera.");
              }
            }
            if (e.xkey.keycode == PreviousKeyframe_key ) {
              JumpKeyframe(-1, &env_y_angle);
            }
            if (e.xkey.keycode == NextKeyframe_key ) {
              JumpKeyframe(1, &env_y_angle);
            }
            if (e.xkey.keycode == DeleteKeyframe_key ) {
              DeleteCurrentKeyframe(&env_y_angle);
            }
            if (e.xkey.keycode == MotionPreview_key ) {
              if(StartPreviewMotion(env_y_angle)) {
                blanked = true;
                reset_preview_render();
              }
            }
            if (e.xkey.keycode == P_key || e.xkey.keycode == K_key ) {
              if(e.xkey.keycode == K_key) {
                SaveCurrentKeyframe(env_y_angle);
              }
              PrintCameraInfo(env_y_angle);
            } else {
              if(!blanked && !terminate && IsX11RenderInvalidatingKey(d, e.xkey.keycode)) {
                ApplyStaticPreviewCameraMotionRange(cam);
                blanked = true;
                ns = 0;
                adaptive_pixel_sampler.reset();
                adaptive_pixel_sampler_small.reset();
                ResetPreviewExposure();
#ifdef HAS_OIDN
                InvalidateOidnAux();
#endif
                if(progress && !interactive) {
                  pb.update(0);
                }
              }
            }
          }
        }
      } else if (e.type == ButtonPress) {
        if(interactive && IsPreviewMotionActive()) {
          continue;
        }
        if(interactive) {
          bool left = e.xbutton.button == Button1;
          bool right = e.xbutton.button == Button3;
          if(!left && !right) {
            return;
          }
          Float x = e.xbutton.x;
          Float y = e.xbutton.y;
          Float u = 1 - (x + .5f) / Float(width);
          Float v = 1 - (y + .5f) / Float(height);
          if(!PickCameraTarget(u, v, left, world)) {
            continue;
          }
          ns = 0;
          adaptive_pixel_sampler.reset();
          adaptive_pixel_sampler_small.reset();
          ResetPreviewExposure();
#ifdef HAS_OIDN
          InvalidateOidnAux();
#endif
          
          if(progress && !interactive) {
            pb.update(0);
          }
        }
      } else if (e.type == ClientMessage) {
        terminate = true;
        break;
      } 
    }
    if(AdvancePreviewMotion(&env_y_angle)) {
      reset_preview_render();
    }
  }
#endif
#ifdef RAY_WINDOWS
  if(hwnd) {
    preview_display_w = this;
    aps = &adaptive_pixel_sampler;
    aps_small = &adaptive_pixel_sampler_small;
    ns_w = &ns;
    pb_w = &pb;
    progress_w = progress;
    interactive_w = interactive;
#ifdef HAS_OIDN
    bool use_denoised_preview = denoise &&
      denoiser != nullptr &&
      denoiser->Ready();
    if(use_denoised_preview && !write_fast_output) {
      if(PollCloseEvent()) {
        return;
      }
      denoiser->Execute();
      if(!denoiser->ReportError()) {
        MarkDenoisedPreviewReady(ns + 1);
      }
    }
    RayMatrix& rgb_s = use_denoised_preview ? adaptive_pixel_sampler.draw_rgb_output : adaptive_pixel_sampler.rgb;
#else
    RayMatrix &rgb_s  = adaptive_pixel_sampler.rgb;
#endif
    CalibratePreviewExposure(adaptive_pixel_sampler, rgb_s, ns);
    EnvWorldToObject_w = EnvWorldToObject;
    EnvObjectToWorld_w = EnvObjectToWorld;
    Start_EnvWorldToObject_w = Start_EnvWorldToObject;
    Start_EnvObjectToWorld_w = Start_EnvObjectToWorld;
    deferred_render_w = deferred_render;
    render_requested_w = &render_requested;
    preview_exposure_adjustment_w = &preview_exposure_adjustment;
    std::vector<bool>& finalized = adaptive_pixel_sampler.finalized;
    std::vector<bool>& just_finalized = adaptive_pixel_sampler.just_finalized;
    write_fast_output_w = &write_fast_output;
    world_w = world;
    rng_w = &rng;
    height = (unsigned int)rgb_s.cols();
    width = (unsigned int)rgb_s.rows();
    rgb.resize(width*height*3);
    for(unsigned int i = 0; i < width*3; i += 3) {
      for(unsigned int j = 0; j < height; j++) {
        Float samples;
        Float r_col,g_col,b_col;
        if(finalized[(width - 1 - i/3) + width * (height-1-j)]) {
          if(just_finalized[(width - 1 - i/3) + width * (height-1-j)] && !interactive ) {
            r_col = 0.f;
            g_col = 1.f;
            b_col = 0.f;
          } else {
            Float sample_count = interactive ? 1.f : 4.f;
            r_col = ApplyPreviewExposure(rgb_s((width - 1 - i/3),height-1-j,0), sample_count);
            g_col = ApplyPreviewExposure(rgb_s((width - 1 - i/3),height-1-j,1), sample_count);
            b_col = ApplyPreviewExposure(rgb_s((width - 1 - i/3),height-1-j,2), sample_count);
          }
        } else {
          samples = (Float)ns+1.f;
          r_col = ApplyPreviewExposure(rgb_s((width - 1 - i/3),height-1-j,0), samples);
          g_col = ApplyPreviewExposure(rgb_s((width - 1 - i/3),height-1-j,1), samples);
          b_col = ApplyPreviewExposure(rgb_s((width - 1 - i/3),height-1-j,2), samples);
        }
        rgb[i+3*width*j]   = clamp(r_col,0.f,1.f);
        rgb[i+3*width*j+1] = clamp(g_col,0.f,1.f);
        rgb[i+3*width*j+2] = clamp(b_col,0.f,1.f);
        
        if(finalized[(width - 1 - i/3) + width * (height-1-j)]) {
          just_finalized[(width - 1 - i/3) + width * (height-1-j)] = false;
        }
      }
    }
    blanked = false;
    
    CompositeTextOverlaysToFloatBuffer(rgb, world, rng);
    CompositeLineOverlaysToFloatBuffer(rgb, world, rng);
    snapshot_width = width;
    snapshot_height = height;
    snapshot_pixels.resize(static_cast<size_t>(width) *
                           static_cast<size_t>(height) * 3);
    for(size_t pixel = 0;
        pixel < static_cast<size_t>(width) * static_cast<size_t>(height);
        pixel++) {
      snapshot_pixels[3 * pixel] = static_cast<unsigned char>(
        255.f * clamp(rgb[3 * pixel], 0.f, 1.f));
      snapshot_pixels[3 * pixel + 1] = static_cast<unsigned char>(
        255.f * clamp(rgb[3 * pixel + 1], 0.f, 1.f));
      snapshot_pixels[3 * pixel + 2] = static_cast<unsigned char>(
        255.f * clamp(rgb[3 * pixel + 2], 0.f, 1.f));
    }
    CaptureVolumeSnapshot(adaptive_pixel_sampler, rgb_s, ns+1, world, rng);
    // Paint display-only UI after populating the snapshot buffer.
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
    if(AdvancePreviewMotion(&env_y_angle)) {
      blanked = true;
      reset_preview_render();
    }
    terminate = term;
  }
#endif
}
#ifdef HAS_OIDN
PreviewDisplay::PreviewDisplay(unsigned int _width, unsigned int _height, 
                               bool preview, bool _interactive,
                               bool _deferred_render, Float initial_lookat_distance, RayCamera* _cam,
                               Transform* _EnvObjectToWorld, Transform* _EnvWorldToObject, 
                               RayOidnDenoiser* _denoiser,
                               RayMatrix* _oidn_albedo_output,
                               RayMatrix* _oidn_normal_output,
                               bool denoise, bool _auto_exposure) :
  preview(preview), auto_exposure(_auto_exposure), preview_exposure_calibrated(false),
  preview_exposure_scale(1.f), preview_exposure_adjustment(1.f),
  EnvObjectToWorld(_EnvObjectToWorld), EnvWorldToObject(_EnvWorldToObject),
  Start_EnvObjectToWorld(*_EnvObjectToWorld), Start_EnvWorldToObject(*_EnvWorldToObject),
  denoiser(_denoiser), oidn_albedo_output(_oidn_albedo_output),
  oidn_normal_output(_oidn_normal_output),
  denoise(denoise && _denoiser != nullptr &&
          _oidn_albedo_output != nullptr &&
          _oidn_normal_output != nullptr),
  oidn_aux_dirty(true), oidn_fast_aux_dirty(true),
  has_denoised_preview(false), denoised_preview_sample_count(0) {
#else
PreviewDisplay::PreviewDisplay(unsigned int _width, unsigned int _height, 
                               bool preview, bool _interactive,
                               bool _deferred_render, Float initial_lookat_distance, RayCamera* _cam,
                               Transform* _EnvObjectToWorld, Transform* _EnvWorldToObject,
                               bool _auto_exposure) :
  preview(preview), auto_exposure(_auto_exposure), preview_exposure_calibrated(false),
  preview_exposure_scale(1.f), preview_exposure_adjustment(1.f),
  EnvObjectToWorld(_EnvObjectToWorld), EnvWorldToObject(_EnvWorldToObject),
  Start_EnvObjectToWorld(*_EnvObjectToWorld), Start_EnvWorldToObject(*_EnvWorldToObject) {
#endif
  width = _width; height = _height;
  Keyframes.clear();
  current_keyframe = -1;
  keyframe_motion_args = Rcpp::List::create(
    Named("type") = "spline",
    Named("smooth_orientation") = true,
    Named("damp_motion") = true
  );
  keyframe_motion_closed = false;
  preview_motion = Rcpp::DataFrame::create();
  preview_motion_restore_state = Rcpp::List::create();
  preview_motion_frame = 0;
  preview_motion_restore_keyframe = -1;
  preview_motion_active = false;
  snapshot_width = 0;
  snapshot_height = 0;
  write_fast_output = false;
  terminate = false;
  deferred_render = _deferred_render && preview && _interactive;
  render_requested = !deferred_render;
  cam = _cam;
  camera_motion_blur_enabled = cam != nullptr && cam->get_camera_motion_blur();
  shutter_speed = cam != nullptr ? cam->get_shutter_speed() : static_cast<Float>(2);
  ApplyShutterSpeedToCameras();
#ifdef RAY_HAS_X11
  speed = 1.f;
  interactive = _interactive;
  env_y_angle = 0;
  orbit = true;
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
  term = false;
  env_y_angle = 0;
  Keyframes_w = &Keyframes;
  preview_display_w = this;
  if(preview) {
    width = _width;
    height = _height;
    base_step = initial_lookat_distance/20;
    rgb.resize(width*height*3);
    cam_w = _cam;
    ApplyShutterSpeedToCameras();
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

void PreviewDisplay::SetCamera(RayCamera* _cam) {
  cam = _cam;
  if(cam != nullptr) {
    cam->set_camera_motion_blur(camera_motion_blur_enabled);
    cam->set_shutter_speed(shutter_speed);
  }
#ifdef RAY_WINDOWS
  if(hwnd != NULL) {
    cam_w = _cam;
    if(cam_w != nullptr) {
      cam_w->set_camera_motion_blur(camera_motion_blur_enabled);
      cam_w->set_shutter_speed(shutter_speed);
    }
  }
#endif
}

void PreviewDisplay::SetSnapshotFilename(const std::string& filename) {
  snapshot_filename = filename;
}

void PreviewDisplay::CaptureVolumeSnapshot(adaptive_sampler& sampler, RayMatrix& color,
                                           size_t samples, hitable* world, random_gen& rng) {
  snapshot_alpha.clear();
  if(!transparent_volume_background) return;
  snapshot_width=width; snapshot_height=height;
  snapshot_alpha.resize(size_t(width)*height);
  std::vector<Float> rgb(size_t(width)*height*3);
  for(size_t y=0;y<height;++y) for(size_t x=0;x<width;++x) {
    size_t sx=width-1-x, sy=height-1-y, pixel=x+width*y;
    bool final=sampler.finalized[sx+width*sy];
    Float count=final ? 1 : std::max(size_t(1),samples);
    Float alpha=clamp(final ? sampler.a(sx,sy,0) : 1-sampler.a(sx,sy,0)/count,0.f,1.f);
    snapshot_alpha[pixel]=alpha;
    for(int channel=0;channel<3;++channel)
      rgb[3*pixel+channel]=alpha>0 ? ApplyPreviewExposure(color(sx,sy,channel)/alpha,count) : 0;
  }
  CompositeTextOverlaysToFloatBuffer(rgb,world,rng,&snapshot_alpha);
  CompositeLineOverlaysToFloatBuffer(rgb,world,rng,&snapshot_alpha);
  snapshot_pixels.resize(rgb.size());
  for(size_t i=0;i<rgb.size();++i) snapshot_pixels[i]=static_cast<unsigned char>(255*clamp(rgb[i],0.f,1.f));
}

void PreviewDisplay::SavePreviewSnapshot() const {
  if(snapshot_width == 0 || snapshot_height == 0 || snapshot_pixels.empty()) {
    Rprintf("Unable to save preview snapshot: no preview image is available.\n");
    return;
  }

  try {
    size_t channels=snapshot_alpha.empty() ? 3 : 4;
    Rcpp::NumericVector image(
      static_cast<R_xlen_t>(snapshot_width) *
      static_cast<R_xlen_t>(snapshot_height) * channels
    );
    size_t channel_size = static_cast<size_t>(snapshot_width) *
      static_cast<size_t>(snapshot_height);
    for(unsigned int y = 0; y < snapshot_height; y++) {
      for(unsigned int x = 0; x < snapshot_width; x++) {
        size_t source_pixel = static_cast<size_t>(x) +
          static_cast<size_t>(snapshot_width) * y;
        size_t target_pixel = static_cast<size_t>(y) +
          static_cast<size_t>(snapshot_height) * x;
        if(channels==4) image[target_pixel + 3*channel_size]=snapshot_alpha[source_pixel];
        for(size_t channel = 0; channel < 3; channel++) {
          image[target_pixel + channel_size * channel] =
            static_cast<Float>(snapshot_pixels[3 * source_pixel + channel]) /
            255.f;
        }
      }
    }
    image.attr("dim") = Rcpp::IntegerVector::create(
      snapshot_height,
      snapshot_width,
      channels
    );

    Rcpp::CharacterVector source_filename(1);
    if(snapshot_filename.empty()) {
      source_filename[0] = NA_STRING;
    } else {
      source_filename[0] = snapshot_filename;
    }
    Rcpp::Environment pkg = Rcpp::Environment::namespace_env("rayrender");
    Rcpp::Function save_snapshot = pkg["save_preview_snapshot"];
    save_snapshot(image, source_filename);
  } catch(const std::exception& error) {
    Rprintf("Unable to save preview snapshot: %s\n", error.what());
  } catch(...) {
    Rprintf("Unable to save preview snapshot.\n");
  }
}

bool PreviewDisplay::PollCloseEvent() {
#ifdef RAY_HAS_X11
  if(d != nullptr && !terminate) {
    std::vector<XEvent> deferred_events;
    KeyCode esc = XKeysymToKeycode(d, XK_Escape);
    while(XPending(d)) {
      XEvent poll_event;
      XNextEvent(d, &poll_event);
      if(poll_event.type == ClientMessage) {
        terminate = true;
      } else if(poll_event.type == KeyPress &&
                poll_event.xkey.keycode == esc &&
                PreviewDisplayHasKeyboardFocus(d, this->w)) {
        terminate = true;
      } else {
        deferred_events.push_back(poll_event);
      }
    }
    if(!terminate) {
      for(auto event_iter = deferred_events.rbegin();
          event_iter != deferred_events.rend();
          ++event_iter) {
        XPutBackEvent(d, &(*event_iter));
      }
    }
  }
#endif
#ifdef RAY_WINDOWS
  if(hwnd != NULL && !terminate) {
    if(GetForegroundWindow() == hwnd && (GetAsyncKeyState(VK_ESCAPE) & 0x8000)) {
      term = true;
      terminate = true;
      PostMessage(hwnd, WM_CLOSE, 0, 0);
    }
    MSG close_msg;
    while(PeekMessage(&close_msg, hwnd, WM_SYSCOMMAND, WM_SYSCOMMAND, PM_REMOVE) > 0) {
      TranslateMessage(&close_msg);
      DispatchMessage(&close_msg);
    }
    while(PeekMessage(&close_msg, hwnd, WM_CLOSE, WM_CLOSE, PM_REMOVE) > 0) {
      TranslateMessage(&close_msg);
      DispatchMessage(&close_msg);
    }
    while(PeekMessage(&close_msg, NULL, WM_QUIT, WM_QUIT, PM_REMOVE) > 0) {
      term = true;
    }
    terminate = term;
  }
#endif
  return terminate;
}

#ifdef HAS_OIDN
void PreviewDisplay::SetDenoiser(RayOidnDenoiser* _denoiser,
                                 RayMatrix* _oidn_albedo_output,
                                 RayMatrix* _oidn_normal_output,
                                 bool _denoise) {
  denoiser = _denoiser;
  oidn_albedo_output = _oidn_albedo_output;
  oidn_normal_output = _oidn_normal_output;
  denoise = _denoise &&
    denoiser != nullptr &&
    oidn_albedo_output != nullptr &&
    oidn_normal_output != nullptr;
  InvalidateOidnAux();
}

void PreviewDisplay::InvalidateOidnAux() {
  oidn_aux_dirty = true;
  oidn_fast_aux_dirty = true;
  has_denoised_preview = false;
  denoised_preview_sample_count = 0;
}

void PreviewDisplay::MarkOidnAuxClean(bool fast_preview) {
  if(fast_preview) {
    oidn_fast_aux_dirty = false;
  } else {
    oidn_aux_dirty = false;
  }
}

void PreviewDisplay::MarkDenoisedPreviewReady(size_t sample_count) {
  has_denoised_preview = true;
  denoised_preview_sample_count = std::max<size_t>(sample_count, 1);
}
#endif

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
  if (preview_display_w == this) {
    preview_display_w = nullptr;
  }
  rgb.resize(0);
#endif
}

#ifdef RAY_WINDOWS

#define VK_KEY_1 49
#define VK_KEY_2 50
#define VK_KEY_3 51
#define VK_KEY_4 52
#define VK_KEY_5 53
#define VK_KEY_6 54
#define VK_KEY_7 55
#define VK_KEY_8 56
#define VK_KEY_COMMA 188
#define VK_KEY_PERIOD 190
#define VK_KEY_SLASH 191
#define VK_KEY_LEFT_BRACKET 219
#define VK_KEY_RIGHT_BRACKET 221
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

static bool PreviewWindowHasKeyboardFocus(HWND hwnd) {
  return GetForegroundWindow() == hwnd && GetFocus() == hwnd;
}

static bool IsWindowsRenderInvalidatingKey(WPARAM key) {
  return key == VK_KEY_W ||
         key == VK_KEY_A ||
         key == VK_KEY_S ||
         key == VK_KEY_D ||
         key == VK_KEY_Q ||
         key == VK_KEY_Z ||
         key == VK_UP ||
         key == VK_DOWN ||
         key == VK_LEFT ||
         key == VK_RIGHT ||
         key == VK_KEY_1 ||
         key == VK_KEY_2 ||
         key == VK_KEY_3 ||
         key == VK_KEY_4 ||
         key == VK_KEY_F ||
         key == VK_KEY_B ||
         key == VK_KEY_L ||
         key == VK_KEY_M ||
         key == VK_KEY_COMMA ||
         key == VK_KEY_PERIOD ||
         key == VK_KEY_SLASH ||
         key == VK_KEY_R;
}

static void ResetWindowsPreviewRenderState(bool update_camera_motion_range = true) {
  if(!blanked && !term) {
    if(update_camera_motion_range) {
      ApplyStaticPreviewCameraMotionRange(cam_w);
    }
    blanked = true;
    if(ns_w != nullptr) {
      *ns_w = 0;
    }
    if(aps != nullptr) {
      aps->reset();
    }
    if(aps_small != nullptr) {
      aps_small->reset();
    }
    if(preview_display_w != nullptr) {
      preview_display_w->ResetPreviewExposure();
#ifdef HAS_OIDN
      preview_display_w->InvalidateOidnAux();
#endif
    }
    if(progress_w && !interactive_w && pb_w != nullptr) {
      pb_w->update(0);
    }
  }
}

LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam) {
  switch (uMsg) {
    case WM_DESTROY: {
      PostQuitMessage(0);
      term = true;
      return 0;
    }
    case WM_KEYDOWN: {
      if(!PreviewWindowHasKeyboardFocus(hwnd)) {
        return 0;
      }
      bool shift_pressed = (GetKeyState(VK_SHIFT) & 0x8000) != 0;
      if(interactive_w &&
         shift_pressed &&
         wParam == VK_RETURN &&
         preview_display_w != nullptr) {
        preview_display_w->SavePreviewSnapshot();
        return 0;
      }
      if(interactive_w &&
         preview_display_w != nullptr &&
         preview_display_w->IsPreviewMotionActive() &&
         wParam != VK_ESCAPE &&
         wParam != VK_KEY_M) {
        return 0;
      }
      if(interactive_w &&
         preview_display_w != nullptr &&
         shift_pressed &&
         wParam == VK_KEY_L) {
        preview_display_w->ToggleKeyframeMotionClosed();
        InvalidateRect(hwnd, nullptr, FALSE);
        return 0;
      }
      vec3f w(1,0,0);
      vec3f u(0,1,0);
      vec3f v(0,0,1);
      
      if(interactive_w) {
        w = cam_w->get_w();
        u = cam_w->get_u();
        v = cam_w->get_v();
      }
      if(interactive_w &&
         preview_display_w != nullptr &&
         shift_pressed &&
         wParam == VK_KEY_LEFT_BRACKET) {
        preview_display_w->AdjustShutterSpeedStops(static_cast<Float>(-1) / static_cast<Float>(3));
        ResetWindowsPreviewRenderState(false);
        return 0;
      }
      if(interactive_w &&
         preview_display_w != nullptr &&
         shift_pressed &&
         wParam == VK_KEY_RIGHT_BRACKET) {
        preview_display_w->AdjustShutterSpeedStops(static_cast<Float>(1) / static_cast<Float>(3));
        ResetWindowsPreviewRenderState(false);
        return 0;
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
        case VK_KEY_B: {
          if(preview_display_w != nullptr) {
            preview_display_w->ToggleCameraMotionBlur();
          }
          break;
        }
        case VK_KEY_W: {
          if(interactive_w) {
            if(shift_pressed) {
              cam_w->rotate_forward(speed * 1.f);
            } else {
              vec3f step = speed * w * base_step;
              if(orbit) {
                Float dist_to_orbit = (cam_w->get_origin() - cam_w->get_lookat()).length();
                if(dist_to_orbit <= base_step * speed) {
                  Rprintf("Moving forward will overstep orbit point, stopping (decrease step size to move closer).\n");
                  step = vec3f(0);
                }
              }
              cam_w->update_position(step, orbit, false);
            }
          }
          break;
        }

        case VK_KEY_A: {
          if(interactive_w) {
            if(shift_pressed) {
              cam_w->rotate_up(speed * -1.f);
            } else {
              cam_w->update_position(-speed * u * base_step, orbit);
            }
          }
          break;
        }
        case VK_KEY_S: {
          if(interactive_w) {
            if(shift_pressed) {
              cam_w->rotate_forward(speed * -1.f);
            } else {
              cam_w->update_position(-speed * w * base_step, orbit, false);
            }
          }
          break;
        }
        case VK_KEY_D: {
          if(interactive_w) {
            if(shift_pressed) {
              cam_w->rotate_up(speed * 1.f);
            } else {
              cam_w->update_position(speed * u * base_step, orbit);
            }
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
            env_y_angle -= speed*8;
          }
          break;
        }
        case VK_KEY_4: {
          if(interactive_w) {
          (*EnvObjectToWorld_w) =  RotateY(-speed*8) * (*EnvObjectToWorld_w);
          (*EnvWorldToObject_w) =  RotateY(speed*8) * (*EnvWorldToObject_w);
          env_y_angle += speed*8;
          }
          break;
        }
        case VK_KEY_RIGHT_BRACKET: {
          if(interactive_w) {
            (*preview_exposure_adjustment_w) *= 2.f;
            Rprintf("Preview Exposure: %.3f\n", *preview_exposure_adjustment_w);
          }
          break;
        }
        case VK_KEY_LEFT_BRACKET: {
          if(interactive_w) {
            (*preview_exposure_adjustment_w) *= 0.5f;
            Rprintf("Preview Exposure: %.3f\n", *preview_exposure_adjustment_w);
          }
          break;
        }
        case VK_RETURN: {
          if(interactive_w) {
            if(deferred_render_w) {
              (*render_requested_w) = !(*render_requested_w);
              if(*render_requested_w) {
                Rprintf("Starting final render...\n");
              } else {
                Rprintf("Returning to deferred render...\n");
              }
            }
          }
          break;
        }
        case VK_KEY_R: {
          if(interactive_w) {
            cam_w->reset();
            speed = 1;
            (*EnvWorldToObject_w) = Start_EnvWorldToObject_w;
            (*EnvObjectToWorld_w) = Start_EnvObjectToWorld_w;
            env_y_angle = 0;
          }
          break;
        }
        case VK_KEY_P: {
          if(preview_display_w != nullptr) {
            preview_display_w->PrintCameraInfo(env_y_angle);
          }
          break;
        }
        case VK_KEY_L: {
          if(preview_display_w != nullptr) {
            if(preview_display_w->Keyframes.size() > 0) {
              preview_display_w->ApplyKeyframe(
                static_cast<int>(preview_display_w->Keyframes.size()) - 1,
                &env_y_angle
              );
            } else {
              Rprintf("Can't reset to last keyframe: No keyframes have been saved. Use the R key to reset camera.");
            }
          }
          break;
        }
        case VK_KEY_M: {
          if(preview_display_w != nullptr) {
            if(preview_display_w->IsPreviewMotionActive()) {
              preview_display_w->CancelPreviewMotion(&env_y_angle);
            } else {
              preview_display_w->StartPreviewMotion(env_y_angle);
            }
          }
          break;
        }
        case VK_KEY_COMMA: {
          if(preview_display_w != nullptr) {
            preview_display_w->JumpKeyframe(-1, &env_y_angle);
          }
          break;
        }
        case VK_KEY_PERIOD: {
          if(preview_display_w != nullptr) {
            preview_display_w->JumpKeyframe(1, &env_y_angle);
          }
          break;
        }
        case VK_KEY_SLASH: {
          if(preview_display_w != nullptr) {
            preview_display_w->DeleteCurrentKeyframe(&env_y_angle);
          }
          break;
        }
        case VK_KEY_K: {
          if(preview_display_w != nullptr) {
            preview_display_w->SaveCurrentKeyframe(env_y_angle);
            preview_display_w->PrintCameraInfo(env_y_angle);
          }
          break;
        }
          
          default: 
            break;
      }
      if(interactive_w) {
        if(IsWindowsRenderInvalidatingKey(wParam)) {
          ResetWindowsPreviewRenderState();
        }
      }
      break;
    }
  case WM_LBUTTONDOWN:
  case WM_RBUTTONDOWN: {
    if(interactive_w && preview_display_w != nullptr &&
       !preview_display_w->IsPreviewMotionActive()) {
      Float u = 1 - (Float(GET_X_LPARAM(lParam)) + .5f) / Float(width);
      Float v = 1 - (Float(GET_Y_LPARAM(lParam)) + .5f) / Float(height);
      if(preview_display_w->PickCameraTarget(u, v, uMsg == WM_LBUTTONDOWN, world_w)) {
        ResetWindowsPreviewRenderState();
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
      if(preview_display_w != nullptr) {
        preview_display_w->DrawStatusBarWindows(hdc, env_y_angle);
      }
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

#ifdef NOT_CRAN
namespace {
std::unique_ptr<PreviewDisplay> MakeTestPreviewDisplay(RayCamera& cam,
                                                       Transform& env_transform) {
#ifdef HAS_OIDN
  return std::unique_ptr<PreviewDisplay>(new PreviewDisplay(
    4, 4, false, true, false, 10.f, &cam,
    &env_transform, &env_transform,
    nullptr, nullptr, nullptr, false, false));
#else
  return std::unique_ptr<PreviewDisplay>(new PreviewDisplay(
    4, 4, false, true, false, 10.f, &cam,
    &env_transform, &env_transform, false));
#endif
}
}

context("Preview picking and orbit targets") {
  test_that("volume clicks update the pivot and sparse misses leave the camera untouched") {
    Transform identity;
    camera cam(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0),
               60.f, 1.f, 1.f, 10.f, 0.f, 1.f, 1.f);
    auto display = MakeTestPreviewDisplay(cam, identity);
    auto scene = std::make_shared<VolumeScene>();
    Rcpp::Function describe = Rcpp::Environment::namespace_env("rayrender")["homogeneous_medium"];
    auto medium = LoadMedium(describe());
    auto mat = std::make_shared<lambertian>(std::make_shared<constant_texture>(point3f(1)));
    auto sphere_geometry = std::make_shared<sphere>(2, mat, nullptr, nullptr,
                                                  &identity, &identity, false);
    auto boundary = std::make_shared<MediumBoundary>(sphere_geometry, medium, identity, false,
                                                    scene->NextBoundaryId());
    scene->boundaries.add(boundary);
    scene->Finish(0, 1);
    display->volume_scene = scene;
    hitable_list world;
    world.add(boundary);
    expect_true(display->PickCameraTarget(.5, .5, true, &world));
    point3f pivot = cam.get_lookat();
    expect_true(std::abs(pivot[2] - (-2 - std::log(.85))) < 1e-5);
    expect_true(std::abs(cam.get_focal_distance() - (pivot - cam.get_origin()).length()) < 1e-5);
    cam.update_position(cam.get_u(), true);
    expect_true((cam.get_lookat() - pivot).length() == 0);
    expect_true(dot(unit_vector(pivot - cam.get_origin()), cam.get_w()) > .99999);
    point3f origin = cam.get_origin();
    vec3f direction = cam.get_w();
    Float focal = cam.get_focal_distance();
    expect_false(display->PickCameraTarget(0, 0, true, &world));
    expect_true((cam.get_origin() - origin).length() == 0);
    expect_true((cam.get_lookat() - pivot).length() == 0);
    expect_true((cam.get_w() - direction).length() == 0);
    expect_true(cam.get_focal_distance() == focal);
  }
  test_that("right clicks preserve focus while orbiting about the selected surface") {
    Transform identity;
    auto mat = std::make_shared<lambertian>(std::make_shared<constant_texture>(point3f(1)));
    hitable_list world;
    world.add(std::make_shared<sphere>(2, mat, nullptr, nullptr, &identity, &identity, false));
    camera cam(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0),
               60.f, 1.f, 0.f, 10.f, 0.f, 1.f, 1.f);
    auto display = MakeTestPreviewDisplay(cam, identity);
    expect_true(display->PickCameraTarget(.5, .5, false, &world));
    point3f pivot = cam.get_lookat();
    Float distance = (pivot - cam.get_origin()).length();
    expect_true(pivot[2] == Approx(-2));
    expect_true(cam.get_focal_distance() == 10);
    cam.update_position(cam.get_u(), true);
    expect_true((cam.get_lookat() - pivot).length() == 0);
    expect_true((cam.get_origin() - pivot).length() == Approx(distance));
    expect_true(cam.get_focal_distance() == 10);
    Rcpp::List keyframe = display->CreateCurrentKeyframe(0);
    expect_true(Rcpp::as<Float>(keyframe["dz"]) == pivot[2]);
    expect_true(Rcpp::as<Float>(keyframe["focal"]) == 10);
  }
  test_that("free-flight keyframes still capture the current viewing direction") {
    Transform identity;
    camera cam(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0),
               60.f, 1.f, 0.f, 10.f, 0.f, 1.f, 1.f);
    auto display = MakeTestPreviewDisplay(cam, identity);
    cam.update_position(vec3f(2, 0, 0), false);
    Rcpp::List keyframe = display->CreateCurrentKeyframe(0);
    point3f target(Rcpp::as<Float>(keyframe["dx"]), Rcpp::as<Float>(keyframe["dy"]),
                    Rcpp::as<Float>(keyframe["dz"]));
    expect_true(dot(unit_vector(target - cam.get_origin()), cam.get_w()) > .999999);
  }
}

context("Preview keyframe loop controls") {
  test_that("preview keyframe loop state toggles and appears in status text") {
    Matrix4x4 identity_matrix;
    Transform env_transform(identity_matrix);
    camera cam(point3f(0, 0, -10), point3f(0, 0, 0), vec3f(0, 1, 0),
               60.f, 1.f, 0.f, 10.f, 0.f, 1.f, 1.f);
    auto display = MakeTestPreviewDisplay(cam, env_transform);

    expect_false(display->KeyframeMotionClosed());
    expect_true(display->PreviewStatusText(0.f).find("Loop OPEN") !=
                std::string::npos);

    display->SetKeyframeMotionArgs(
      Rcpp::List::create(Rcpp::Named("closed") = true)
    );
    expect_true(display->KeyframeMotionClosed());
    expect_true(display->PreviewStatusText(0.f).find("Loop CLOSED") !=
                std::string::npos);

    expect_false(display->ToggleKeyframeMotionClosed());
    expect_true(display->PreviewStatusText(0.f).find("Loop OPEN") !=
                std::string::npos);
  }
}

context("Preview shutter speed controls") {
  test_that("preview shutter speed adjustments use one-third stop increments") {
    Matrix4x4 identity_matrix;
    Transform env_transform(identity_matrix);
    camera cam(point3f(0, 0, -10), point3f(0, 0, 0), vec3f(0, 1, 0),
               60.f, 1.f, 0.f, 10.f, 0.f, 1.f, 1.f);
    auto display = MakeTestPreviewDisplay(cam, env_transform);

    display->SetShutterSpeed(2.f);
    display->AdjustShutterSpeedStops(static_cast<Float>(1) / static_cast<Float>(3));
    expect_true(display->GetShutterSpeed() == Approx(2.f * std::pow(2.f, 1.f / 3.f)));
    expect_true(cam.get_shutter_speed() == Approx(display->GetShutterSpeed()));

    display->AdjustShutterSpeedStops(static_cast<Float>(-1) / static_cast<Float>(3));
    expect_true(display->GetShutterSpeed() == Approx(2.f));
  }

  test_that("preview shutter speed clamps finite keyboard adjustments") {
    Matrix4x4 identity_matrix;
    Transform env_transform(identity_matrix);
    camera cam(point3f(0, 0, -10), point3f(0, 0, 0), vec3f(0, 1, 0),
               60.f, 1.f, 0.f, 10.f, 0.f, 1.f, 1.f);
    auto display = MakeTestPreviewDisplay(cam, env_transform);

    display->SetShutterSpeed(1.f);
    display->AdjustShutterSpeedStops(static_cast<Float>(-1) / static_cast<Float>(3));
    expect_true(display->GetShutterSpeed() == Approx(1.f));

    display->SetShutterSpeed(4096.f);
    display->AdjustShutterSpeedStops(static_cast<Float>(1) / static_cast<Float>(3));
    expect_true(display->GetShutterSpeed() == Approx(4096.f));
  }

  test_that("preview shutter speed exits Inf through slower-shutter adjustment") {
    Matrix4x4 identity_matrix;
    Transform env_transform(identity_matrix);
    camera cam(point3f(0, 0, -10), point3f(0, 0, 0), vec3f(0, 1, 0),
               60.f, 1.f, 0.f, 10.f, 0.f, 1.f, 1.f);
    auto display = MakeTestPreviewDisplay(cam, env_transform);

    display->SetShutterSpeed(Infinity);
    display->AdjustShutterSpeedStops(static_cast<Float>(1) / static_cast<Float>(3));
    expect_true(std::isinf(display->GetShutterSpeed()));

    display->AdjustShutterSpeedStops(static_cast<Float>(-1) / static_cast<Float>(3));
    expect_true(display->GetShutterSpeed() == Approx(4096.f));
  }
}
#endif
