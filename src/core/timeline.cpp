#include "timeline.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace animation {
TimelineController::TimelineController(double begin, double end, double initial, double fps, Clock clock)
  : begin(begin), end(end), fps(fps), clock(std::move(clock)) {
  if (!std::isfinite(begin) || !std::isfinite(end) || begin > end ||
      !std::isfinite(end - begin) || !std::isfinite(initial) || initial < begin || initial > end ||
      !std::isfinite(fps) || fps <= 0)
    throw std::runtime_error("Timeline requires a finite ordered range, contained initial time, and positive fps.");
  if (!this->clock) this->clock = [] {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
  };
  requested.time = initial;
  published = displayed = requested;
  published.revision = displayed.revision = 0;
  Anchor(this->clock());
}

void TimelineController::Anchor(double now) {
  if (!std::isfinite(now)) throw std::runtime_error("Timeline wall clock must be finite.");
  anchor_wall = now;
  anchor_scene = requested.time;
}

void TimelineController::Change(double time, uint32_t dirty) {
  requested.time = time;
  requested.dirty |= dirty;
  ++requested.revision;
}

void TimelineController::Advance(double now) {
  if (!std::isfinite(now) || now < anchor_wall)
    throw std::runtime_error("Timeline requires a monotonic finite wall clock.");
  if (!Moving() || requested.dragging) return;
  double time = anchor_scene + rate * (now - anchor_wall);
  if (!std::isfinite(time)) throw std::runtime_error("Timeline playback exceeds finite time precision.");
  if (loop) {
    double duration = end - begin;
    double phase = std::fmod(time - begin, duration);
    if (phase < 0) phase += duration;
    if (phase >= duration) phase = 0;
    time = begin + phase;
  } else {
    time = std::clamp(time, begin, end);
    if ((rate > 0 && time == end) || (rate < 0 && time == begin)) paused = true;
  }
  if (time != requested.time) Change(time, TimePlacement | Camera);
}

void TimelineController::Seek(double time) {
  if (!std::isfinite(time)) throw std::runtime_error("Timeline seek time must be finite.");
  time = std::clamp(time, begin, end);
  if (time != requested.time) Change(time, TimePlacement | Camera);
  Anchor(clock());
}

void TimelineController::Pause(bool value) {
  double now = clock();
  Advance(now);
  paused = value || begin == end;
  Anchor(now);
}

void TimelineController::SetRate(double value) {
  if (!std::isfinite(value)) throw std::runtime_error("Timeline rate must be finite.");
  double now = clock();
  Advance(now);
  rate = value;
  Anchor(now);
}

void TimelineController::SetLoop(bool value) {
  double now = clock();
  Advance(now);
  loop = value;
  Anchor(now);
}

void TimelineController::Drag(double time) {
  if (begin == end) return;
  Pause(true);
  if (!requested.dragging) {
    requested.dragging = true;
    Change(requested.time, TimePlacement | Camera);
  }
  Seek(time);
}

void TimelineController::EndDrag() {
  if (!requested.dragging) return;
  requested.dragging = false;
  Change(requested.time, TimePlacement | Camera);
  Anchor(clock());
}

void TimelineController::Command(TimelineCommand command) {
  switch (command) {
    case TimelineCommand::PlayPause: Pause(!paused); break;
    case TimelineCommand::Reverse: SetRate(rate == 0 ? -1 : -rate); Pause(false); break;
    case TimelineCommand::StepBack:
    case TimelineCommand::StepForward:
      Pause(true);
      Seek(requested.time + (command == TimelineCommand::StepForward ? 1 : -1) / fps);
      break;
    case TimelineCommand::Begin: Pause(true); Seek(begin); break;
    case TimelineCommand::End: Pause(true); Seek(end); break;
    case TimelineCommand::Slower: SetRate(rate / 2); break;
    case TimelineCommand::Faster: SetRate(rate == 0 ? 1 : rate * 2); break;
    case TimelineCommand::Loop: SetLoop(!loop); break;
  }
}

TimelineRequest TimelineController::Request() {
  Advance(clock());
  return requested;
}

bool TimelineController::Pending() const {
  return requested.revision != published.revision && requested.revision != rejected;
}

bool TimelineController::Publish(const TimelineRequest &prepared) {
  if (prepared.revision != requested.revision || prepared.revision == rejected) return false;
  published = prepared;
  requested.dirty = 0;
  return true;
}

bool TimelineController::Display(uint64_t revision) {
  if (revision != published.revision) return false;
  displayed = published;
  return true;
}

void TimelineController::Reject(uint64_t revision) {
  if (revision != requested.revision) return;
  rejected = revision;
  paused = true;
  Anchor(clock());
}

void TimelineController::Invalidate(uint32_t dirty) {
  Change(requested.time, dirty);
}
} // namespace animation
