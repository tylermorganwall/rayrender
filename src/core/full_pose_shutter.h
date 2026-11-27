#ifndef RAYRENDER_FULL_POSE_SHUTTER_H
#define RAYRENDER_FULL_POSE_SHUTTER_H
#include "camera.h"
#include <algorithm>
#include <cstdint>

// One exact pose serves all paths of an image sample. Film and pixel samplers
// advance normally; only geometry and existing bounds change between samples.
class FullPoseShutter {
public:
  static bool Enabled(const RayCamera &camera) {
    auto clock = camera.scene_camera_time;
    auto interval = camera.ShutterInterval();
    return clock && clock->full_pose_blur && clock->has_model_poses &&
      interval.first != interval.second;
  }

  // Digit-scrambled base-two sequence: every power-of-two prefix covers all
  // equal-width strata, even when preview rendering is interrupted early.
  static Float Unit(uint32_t sample, uint32_t seed) {
    sample = (sample << 16) | (sample >> 16);
    sample = ((sample & 0x00ff00ffu) << 8) | ((sample & 0xff00ff00u) >> 8);
    sample = ((sample & 0x0f0f0f0fu) << 4) | ((sample & 0xf0f0f0f0u) >> 4);
    sample = ((sample & 0x33333333u) << 2) | ((sample & 0xccccccccu) >> 2);
    sample = ((sample & 0x55555555u) << 1) | ((sample & 0xaaaaaaaau) >> 1);
    // Use exactly representable float cell centers. Rounding a 32-bit fraction
    // can otherwise cross a stratum boundary and duplicate the last stratum.
    return (Float((sample ^ seed) >> 9) + Float(.5)) * Float(1.0 / 8388608.0);
  }

  bool Prepare(RayCamera &camera, size_t sample, uint32_t seed) {
    if (!Enabled(camera)) return false;
    if (!clock) {
      clock = camera.scene_camera_time;
      interval = camera.ShutterInterval();
      epoch = clock->epoch;
      anchor = clock->ray_anchor;
    }
    Float unit = Unit(uint32_t(sample), seed);
    // Subtract the frame anchor before interpolation, so adding that anchor
    // cannot erase fractional bits from the sampled offset.
    Float offset = Float((double(interval.first) - double(anchor)) +
      (double(interval.second) - double(interval.first)) *
      double(unit / camera.get_shutter_speed()));
    clock->UpdateGeometry(epoch + double(offset), anchor, anchor, anchor);
    clock->shutter_unit = unit;
    ++clock->shutter_samples;
    if (clock->benchmark_timing) {
      clock->shutter_times.push_back(clock->epoch);
      clock->shutter_pose_seconds += clock->pose_seconds;
      clock->shutter_refit_seconds += clock->refit_seconds;
    }
    return true;
  }

  // Restore the requested frame before UI picking or returning a retained scene.
  void Restore() {
    if (!clock) return;
    clock->shutter_unit = -1;
    clock->UpdateGeometry(epoch, anchor, interval.first, interval.second);
    clock.reset();
  }

  ~FullPoseShutter() {
    // On an exception the render scene is destroyed by its owner. Do not attempt
    // a potentially throwing deformation update while unwinding that exception.
    if (clock) clock->shutter_unit = -1;
  }
private:
  std::shared_ptr<SceneTimeContext> clock;
  std::pair<Float, Float> interval;
  double epoch = 0;
  Float anchor = 0;
};
#endif
