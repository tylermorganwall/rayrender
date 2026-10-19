#include "PreviewDisplay.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>

void PreviewDisplay::ConfigureTimelineCamera(std::shared_ptr<const animation::CameraTrack> track,
                                             std::unique_ptr<RayCamera>* owner,
                                             const RealisticCamera *optics) {
  if (!timeline) throw std::runtime_error("Camera playback requires a scene timeline.");
  if (owner) {
    if (owner->get() != cam) throw std::runtime_error("Timeline camera owner does not match the preview.");
    timeline_camera_owner = owner;
  }
  std::unique_ptr<const RealisticCamera> physical;
  if (auto *lens = optics ? optics : dynamic_cast<RealisticCamera *>(cam)) {
    physical = std::make_unique<const RealisticCamera>(*lens);
    if (!timeline_camera_owner)
      throw std::runtime_error("Physical-camera playback requires an owned preview camera.");
  }
  // Validate a replacement before changing the active track. The owner is
  // replaced only by UpdateTimeline at the render handoff.
  if (timeline_camera_owner)
    track->MakeCamera(track->At(timeline->Requested().time +
                        (timeline_open == timeline_close ? timeline_open : 0)),Float(width)/Float(height),
                      cam->get_iso(),timeline_open,timeline_close,cam->get_camera_flip_x(),physical.get());
  else track->Apply(track->At(timeline->Requested().time), *cam);
  timeline_physical_camera = std::move(physical);
  timeline_camera = std::move(track);
  cam->scene_camera_track = timeline_camera;
  cam->scene_camera_time = timeline_scene_time;
  cam->scene_camera_aspect = Float(width)/Float(height);
  timeline_camera_override = false;
  timeline->Invalidate(animation::Camera);
}

void PreviewDisplay::TimelineCameraEdited() {
  if (!timeline || !timeline_camera) return;
  timeline->Pause(true);
  // Hold this manual view until an explicit scene-time command resumes the
  // immutable camera track. Navigation never writes back into authored keys.
  timeline_camera_override = true;
  cam->scene_camera_track.reset();
  timeline->Invalidate(animation::Camera);
}

void PreviewDisplay::ConfigureTimeline(double begin, double end, double initial, double fps,
                                      std::shared_ptr<SceneTimeContext> scene_time,
                                      Float shutter_open, Float shutter_close) {
  timeline = std::make_unique<animation::TimelineController>(begin, end, initial, fps);
  timeline_scene_time = std::move(scene_time);
  timeline_epoch = initial;
  timeline_open = shutter_open;
  timeline_close = shutter_close;
  timeline->Publish(timeline->Request());
  if (preview && interactive) {
#if defined(RAY_HAS_X11)
    if (!d) throw std::runtime_error("Interactive timeline requires an available X11 display.");
#elif defined(RAY_WINDOWS)
    if (!hwnd) throw std::runtime_error("Interactive timeline requires an available preview window.");
#else
    throw std::runtime_error("This build does not have an interactive preview backend.");
#endif
  }
#ifdef RAY_HAS_X11
  if (d) XResizeWindow(d, w, width, height + StatusBarHeight());
#endif
#ifdef RAY_WINDOWS
  if (hwnd) {
    RECT rect = {0, 0, LONG(width), LONG(height + StatusBarHeight())};
    AdjustWindowRect(&rect, WS_THICKFRAME | WS_VISIBLE | WS_SYSMENU, FALSE);
    SetWindowPos(hwnd, NULL, 0, 0, rect.right-rect.left, rect.bottom-rect.top, SWP_NOMOVE | SWP_NOZORDER);
  }
#endif
}

void PreviewDisplay::TimelineAction(animation::TimelineCommand command) noexcept {
  if (!timeline) return;
  try {
    if (timeline_camera_override) {
      timeline_camera_override = false;
      timeline->Invalidate(animation::Camera);
    }
    timeline->Command(command);
    timeline_error.clear();
  } catch (const std::exception &error) {
    timeline_error = error.what();
    Rprintf("Timeline command failed: %s\n", error.what());
  } catch (...) { Rprintf("Timeline command failed.\n"); }
}

bool PreviewDisplay::TimelinePointer(int x, int y, bool press, bool release) noexcept {
  if (!timeline) return false;
  try {
    bool dragging = timeline->Requested().dragging;
    if (release && dragging) {
      timeline->EndDrag();
      return true;
    }
    if (dragging || (press && y >= int(height) + 96 && y < int(height) + 120)) {
      timeline_camera_override = false;
      double fraction = std::clamp((x - 12.) / std::max(1., double(width) - 24), 0., 1.);
      timeline->Drag(timeline->Begin() + fraction * (timeline->End() - timeline->Begin()));
      return true;
    }
    if (press && y >= int(height) + 72 && y < int(height) + 96) {
      static const animation::TimelineCommand commands[] = {
        animation::TimelineCommand::Begin, animation::TimelineCommand::StepBack,
        animation::TimelineCommand::PlayPause, animation::TimelineCommand::Reverse,
        animation::TimelineCommand::StepForward, animation::TimelineCommand::End,
        animation::TimelineCommand::Slower, animation::TimelineCommand::Faster,
        animation::TimelineCommand::Loop};
      TimelineAction(commands[std::clamp(int(9. * x / std::max(1u, width)), 0, 8)]);
      return true;
    }
  } catch (const std::exception &error) {
    timeline_error = error.what();
    Rprintf("Timeline input failed: %s\n", error.what());
  } catch (...) { Rprintf("Timeline input failed.\n"); }
  return false;
}

void PreviewDisplay::TimelineFocusLost() noexcept {
  if (!timeline) return;
  try {
    timeline->EndDrag();
    timeline->Pause(true);
  } catch (...) { Rprintf("Unable to pause timeline after focus loss.\n"); }
}

bool PreviewDisplay::UpdateTimeline() noexcept {
  if (!timeline || !timeline_scene_time || !cam) return false;
  try {
    auto request = timeline->Request();
    if (timeline->Pending()) {
      // Validate the complete replacement before touching the live camera or
      // epoch. No worker or device batch is active at this handoff.
      double offset = request.time - timeline_epoch;
      Float anchor = Float(offset);
      double open = offset + (request.dragging ? 0 : timeline_open);
      double close = offset + (request.dragging ? 0 : timeline_close);
      if (!std::isfinite(offset) || std::abs(double(anchor) - offset) > 1e-5 ||
          !std::isfinite(Float(open)) || !std::isfinite(Float(close)) ||
          std::abs(double(Float(open)) - open) > 1e-5 || std::abs(double(Float(close)) - close) > 1e-5 ||
          (open < close && Float(open) >= Float(close)))
        throw std::runtime_error("Requested timeline exposure exceeds float time precision.");
      std::unique_ptr<RayCamera> replacement;
      if (timeline_camera && !timeline_camera_override && timeline_camera_owner) {
        const double camera_time = request.time +
          (!request.dragging && timeline_open == timeline_close ? timeline_open : 0);
        replacement = timeline_camera->MakeCamera(timeline_camera->At(camera_time),
          Float(width)/Float(height),cam->get_iso(),open,close,cam->get_camera_flip_x(),
          timeline_physical_camera.get());
      }
      RayCamera* next = replacement ? replacement.get() : cam;
      auto perspective = dynamic_cast<camera *>(next);
      auto orthographic = dynamic_cast<ortho_camera *>(next);
      auto panoramic = dynamic_cast<environment_camera *>(next);
      auto physical = dynamic_cast<RealisticCamera *>(next);
      if (!perspective && !orthographic && !panoramic && !physical)
        throw std::runtime_error("Timeline does not support this camera type.");
      if (timeline_camera && !timeline_camera_override) {
        if (!replacement) timeline_camera->Apply(timeline_camera->At(request.time), *next);
        next->scene_camera_track = timeline_camera;
        next->scene_camera_time = timeline_scene_time;
        next->scene_camera_aspect = Float(width)/Float(height);
        if (timeline_saved_camera && !Keyframes.empty()) {
          double row = (request.time-timeline->Begin()) * timeline->Fps();
          double key = row * Keyframes.size() / std::max<R_xlen_t>(1, preview_motion.nrows());
          current_keyframe = int(std::clamp(key, 0., double(Keyframes.size()-1)));
        }
      }
      timeline_scene_time->UpdateGeometry(request.time, anchor, Float(open), Float(close));
      if (!timeline->Publish(request)) return false;
      if (replacement) {
        *timeline_camera_owner = std::move(replacement);
        SetCamera(timeline_camera_owner->get());
      }
      navigation.SetDistance((cam->get_lookat() - cam->get_origin()).length());
      if (perspective) { perspective->time0 = open; perspective->time1 = close; }
      if (orthographic) { orthographic->time0 = open; orthographic->time1 = close; }
      if (panoramic) { panoramic->time0 = open; panoramic->time1 = close; }
      if (physical) physical->SetShutterInterval(open, close);
      timeline_scene_time->epoch = request.time;
      timeline_scene_time->ray_anchor = anchor;
      cam->set_shutter_speed(1);
      cam->set_camera_motion_blur(false);
      timeline_changed = true;
      timeline_error.clear();
    }
    bool moving = timeline->Moving() || timeline->Requested().dragging;
    if (moving != timeline_fast) {
      if (moving) timeline_restore_fast = write_fast_output;
      write_fast_output = moving ? true : timeline_restore_fast;
      timeline_fast = moving;
      timeline_changed = true;
    }
    return timeline_changed;
  } catch (const std::exception &error) {
    timeline_error = error.what();
    timeline->Reject(timeline->Requested().revision);
    Rprintf("Timeline update rejected: %s\n", error.what());
  } catch (...) {
    timeline->Reject(timeline->Requested().revision);
    Rprintf("Timeline update rejected.\n");
  }
  return false;
}

bool PreviewDisplay::ConsumeTimelineChange() {
  bool changed = timeline_changed;
  timeline_changed = false;
  return changed;
}

std::string PreviewDisplay::TimelineStatusText() const {
  if (!timeline) return "";
  char text[256];
  std::snprintf(text, sizeof(text), "Time %.3f / %.3f | %s %.2fx | Scene loop %s",
    timeline->Displayed().time, timeline->End(), timeline->Paused() ? "Paused" : "Play",
    timeline->Rate(), timeline->Looping() ? "ON" : "OFF");
  return timeline_error.empty() ? text : timeline_error;
}
