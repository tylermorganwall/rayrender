#ifndef RAYRENDER_TIMELINE_H
#define RAYRENDER_TIMELINE_H

#include <cstdint>
#include <functional>

namespace animation {
enum class TimelineCommand { PlayPause, Reverse, StepBack, StepForward, Begin, End, Slower, Faster, Loop };
enum TimelineDirty : uint32_t { TimePlacement = 1, Pose = 2, Topology = 4, Camera = 8, Shading = 16 };

struct TimelineRequest {
  double time = 0;
  uint64_t revision = 1;
  uint32_t dirty = TimePlacement;
  bool dragging = false;
};

// The owning render thread polls this controller between completed batches.
// Tests inject a clock: there are no sleeps, windows, frame counters, or R APIs.
// Requested, published, and displayed revisions are deliberately separate.
class TimelineController {
public:
  using Clock = std::function<double()>;
  TimelineController(double begin, double end, double initial, double fps, Clock clock = {});
  void Seek(double time);
  void Pause(bool paused);
  void SetRate(double rate);
  void SetLoop(bool loop);
  void Drag(double time);
  void EndDrag();
  void Command(TimelineCommand command);
  TimelineRequest Request();
  bool Pending() const;
  bool Publish(const TimelineRequest &prepared);
  bool Display(uint64_t revision);
  void Reject(uint64_t revision);
  void Invalidate(uint32_t dirty);

  const TimelineRequest &Requested() const { return requested; }
  const TimelineRequest &Published() const { return published; }
  const TimelineRequest &Displayed() const { return displayed; }
  double Begin() const { return begin; }
  double End() const { return end; }
  double Fps() const { return fps; }
  double Rate() const { return rate; }
  bool Paused() const { return paused; }
  bool Looping() const { return loop; }
  bool Moving() const { return !paused && rate != 0 && begin < end; }

private:
  void Advance(double now);
  void Anchor(double now);
  void Change(double time, uint32_t dirty);
  const double begin, end, fps;
  Clock clock;
  double anchor_wall = 0, anchor_scene = 0, rate = 1;
  bool paused = true, loop = false;
  TimelineRequest requested, published, displayed;
  uint64_t rejected = 0;
};
} // namespace animation
#endif
