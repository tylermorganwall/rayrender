#include "PreviewDisplay.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>

void PreviewDisplay::ConfigureTimelineCamera(std::shared_ptr<const animation::CameraTrack> track) {
  if (!timeline) throw std::runtime_error("Camera playback requires a scene timeline.");
  if (timeline_open != 0 || timeline_close != 0)
    throw std::runtime_error("Animated camera timeline requires a zero-width exposure until exact camera exposure sampling is implemented.");
  // Probe a copy before replacing the published camera or the active track.
  track->Apply(track->At(timeline->Requested().time), *cam);
  timeline_camera = std::move(track);
  timeline_camera_override = false;
  timeline->Invalidate(animation::Camera);
}

void PreviewDisplay::TimelineCameraEdited() {
  if (!timeline || !timeline_camera) return;
  timeline->Pause(true);
  // Hold this manual view until an explicit scene-time command resumes the
  // immutable camera track. Navigation never writes back into authored keys.
  timeline_camera_override = true;
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
      auto perspective = dynamic_cast<camera *>(cam);
      auto orthographic = dynamic_cast<ortho_camera *>(cam);
      auto panoramic = dynamic_cast<environment_camera *>(cam);
      auto physical = dynamic_cast<RealisticCamera *>(cam);
      if (!perspective && !orthographic && !panoramic && !physical)
        throw std::runtime_error("Timeline does not support this camera type.");
      if (timeline_camera && !timeline_camera_override) {
        timeline_camera->Apply(timeline_camera->At(request.time), *cam);
        navigation.SetDistance((cam->get_lookat() - cam->get_origin()).length());
        if (timeline_saved_camera && !Keyframes.empty()) {
          double row = (request.time-timeline->Begin()) * timeline->Fps();
          double key = row * Keyframes.size() / std::max<R_xlen_t>(1, preview_motion.nrows());
          current_keyframe = int(std::clamp(key, 0., double(Keyframes.size()-1)));
        }
      }
      if (!timeline->Publish(request)) return false;
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
