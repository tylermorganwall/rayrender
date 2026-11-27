#ifndef RAYRENDER_WAVEFRONT_MOTION_LAYOUT_H
#define RAYRENDER_WAVEFRONT_MOTION_LAYOUT_H
#include "wavefront.h"

// ABI 4: byte offsets address aligned, immutable storage in one owned buffer.
// Times are seconds relative to AnimatedTransform::TimeOrigin(), validated at
// export so float packing cannot collapse neighboring keys.
struct WFMotionHeader {
  uint32_t count, owners[5], nodes, bounds;
  uint32_t version = 4, byte_size = 0, reserved[2]{};
};
struct WFMotionNode {
  WFVector endpoints[2][3], translation[2], rotation[2], scale[2][3];
  WFVector before[3], reference_inverse[3];
  uint32_t parent;
  float start, end;
  uint32_t track = 0;
};
struct WFTransformTrack {
  uint32_t channels[3]{}, version = 3;
  WFVector pivot{};
  // Phase at the scene epoch, speed, duration, mode (0 = scene time).
  // Clocked channel times are relative to the clip's local beginning.
  WFVector clock{};
  WFVector frame{}; // x: float ray anchor subtracted before clock mapping.
  WFVector root_position{}; // xyz reference position; w translation bit mask.
  WFVector root_yaw{}; // y,w reference twist; x enables yaw removal.
};
struct WFPropertyTrack {
  uint32_t count, interpolation, property, keys;
};
struct WFTrackKey {
  WFVector value{}, incoming{}, outgoing{};
  float time;
  uint32_t reserved[3]{};
};
static_assert(sizeof(WFMotionHeader) == 48 && sizeof(WFMotionNode) == 368, "motion ABI 4 mismatch");
static_assert(sizeof(WFTransformTrack) == 96 && sizeof(WFPropertyTrack) == 16 && sizeof(WFTrackKey) == 64,
              "track ABI 3 mismatch");
#endif
