#ifndef RAYRENDER_WAVEFRONT_MOTION_LAYOUT_H
#define RAYRENDER_WAVEFRONT_MOTION_LAYOUT_H
#include "wavefront.h"

// ABI 3: byte offsets address aligned, immutable storage in one owned buffer.
// Times are seconds relative to AnimatedTransform::TimeOrigin(), validated at
// export so float packing cannot collapse neighboring keys.
struct WFMotionHeader {
  uint32_t count, owners[5], nodes, bounds;
  uint32_t version = 3, byte_size = 0, reserved[2]{};
};
struct WFMotionNode {
  WFVector endpoints[2][3], translation[2], rotation[2], scale[2][3];
  WFVector before[3], reference_inverse[3];
  uint32_t parent;
  float start, end;
  uint32_t track = 0;
};
struct WFTransformTrack {
  uint32_t channels[3]{}, version = 2;
  WFVector pivot{};
  // Phase at the scene epoch, speed, duration, mode (0 = scene time).
  // Clocked channel times are relative to the clip's local beginning.
  WFVector clock{};
  WFVector frame{}; // x: float ray anchor subtracted before clock mapping.
};
struct WFPropertyTrack {
  uint32_t count, interpolation, property, keys;
};
struct WFTrackKey {
  WFVector value{}, incoming{}, outgoing{};
  float time;
  uint32_t reserved[3]{};
};
static_assert(sizeof(WFMotionHeader) == 48 && sizeof(WFMotionNode) == 368, "motion ABI 3 mismatch");
static_assert(sizeof(WFTransformTrack) == 64 && sizeof(WFPropertyTrack) == 16 && sizeof(WFTrackKey) == 64,
              "track ABI 2 mismatch");
#endif
