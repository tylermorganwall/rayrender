#include "wavefront_motion.h"
#include "wavefront_motion_layout.h"
#include <cstring>
#include <set>
#include <map>

namespace {
Transform unpack(const WFVector *rows) {
  Float m[4][4] = {{0}};
  m[3][3] = 1;
  for (int i = 0; i < 3; ++i) {
    m[i][0] = rows[i].x;
    m[i][1] = rows[i].y;
    m[i][2] = rows[i].z;
    m[i][3] = rows[i].w;
  }
  return Transform(m);
}
void pack(const Matrix4x4 &m, WFVector *rows) {
  for (int i = 0; i < 3; ++i)
    rows[i] = {float(m.m[i][0]), float(m.m[i][1]), float(m.m[i][2]), float(m.m[i][3])};
}
void extend(WFVector &lo, WFVector &hi, const WFVector &p) {
  lo.x = std::min(lo.x, p.x);
  lo.y = std::min(lo.y, p.y);
  lo.z = std::min(lo.z, p.z);
  hi.x = std::max(hi.x, p.x);
  hi.y = std::max(hi.y, p.y);
  hi.z = std::max(hi.z, p.z);
}
double clock_phase(const animation::ClipClock &clock, double epoch) {
  double duration = clock.end - clock.begin;
  if (duration == 0) return 0;
  double phase = clock.Phase(epoch);
  if (clock.loop != animation::Loop::Hold) {
    double period = clock.loop == animation::Loop::Repeat ? duration : 2 * duration;
    phase = std::fmod(phase, period);
    if (phase < 0) phase += period;
    if (phase >= period) phase = 0;
  }
  return phase;
}
} // namespace

void WavefrontMotion::ValidateTimes(double begin, double end) const {
  for (const auto &node : nodes) {
    const auto &track = node.animation->Track();
    if (!track || !track->clock) continue;
    const auto &clock = *track->clock;
    double duration = clock.end - clock.begin;
    if (duration == 0) continue;
    double phase = clock_phase(clock, node.animation->FrameEpoch());
    double anchor = node.animation->FrameAnchor();
    double extent = std::max(std::abs(begin - anchor), std::abs(end - anchor));
    double distance = std::abs(phase) + std::abs(clock.speed) * extent;
    // Bound float input, multiply/add, and period quantization error before a
    // batch is published. Excessive cycle counts must not silently lose phase.
    double error = std::abs(double(float(phase)) - phase) +
      extent * std::abs(double(float(clock.speed)) - clock.speed) +
      (2 + distance / duration) * std::abs(double(float(duration)) - duration) +
      4 * std::numeric_limits<float>::epsilon() * distance;
    double tolerance = std::min(1e-5, duration * 1e-5);
    for (const auto &channel : {track->position, track->rotation, track->scale}) {
      if (!channel) continue;
      for (size_t i = 1; i < channel->times.size(); ++i)
        tolerance = std::min(tolerance, (channel->times[i] - channel->times[i - 1]) * 1e-5);
    }
    if (!std::isfinite(error) || !std::isfinite(float(duration)) ||
        float(duration) <= 0 || error > tolerance)
      throw std::runtime_error("Device clip clock exceeds key-interval or ten-microsecond precision; use CPU or a shorter scene interval.");
    // A small time error can be a large pose error at a discontinuity. Refuse
    // snapshots whose packed phase would cross a repeat seam or a STEP key.
    float packed_duration = float(duration), local = float(phase);
    if (clock.loop == animation::Loop::Hold) local = std::clamp(local, 0.f, packed_duration);
    else {
      float period = clock.loop == animation::Loop::Repeat ? packed_duration : 2 * packed_duration;
      if (local >= period && phase < (clock.loop == animation::Loop::Repeat ? duration : 2 * duration))
        throw std::runtime_error("Device clock rounding crosses a loop boundary; use CPU for this exact time.");
      local = std::fmod(local, period);
      if (clock.loop == animation::Loop::PingPong && local > packed_duration) local = period - local;
    }
    double native_time = clock.At(node.animation->FrameEpoch());
    for (const auto &channel : {track->position, track->rotation, track->scale}) {
      if (!channel || channel->interpolation != animation::Interpolation::Step) continue;
      for (double key : channel->times)
        if ((native_time >= key) != (local >= float(key - clock.begin)))
          throw std::runtime_error("Device clock rounding crosses a STEP key; use CPU for this exact time.");
    }
  }
}

void WavefrontMotion::UpdateClocks() {
  for (size_t i = 0; i < nodes.size(); ++i) {
    const auto &animation = *nodes[i].animation;
    if (!animation.Track() || !animation.Track()->clock) continue;
    const auto &clock = *animation.Track()->clock;
    WFTransformTrack descriptor;
    std::memcpy(&descriptor, packed.data() + packed_transform_offsets[i], sizeof(descriptor));
    descriptor.clock.x = float(clock_phase(clock, animation.FrameEpoch()));
    descriptor.frame.x = animation.FrameAnchor();
    std::memcpy(packed.data() + packed_transform_offsets[i], &descriptor, sizeof(descriptor));
  }
}

void WavefrontMotion::CopyClockState(void *destination) const {
  for (size_t i = 0; i < nodes.size(); ++i) {
    if (!nodes[i].animation->Track() || !nodes[i].animation->Track()->clock) continue;
    size_t offset = packed_transform_offsets[i] + offsetof(WFTransformTrack, clock);
    std::memcpy(static_cast<unsigned char *>(destination) + offset, packed.data() + offset, 2 * sizeof(WFVector));
  }
}

std::vector<unsigned char> WavefrontMotion::Export(WavefrontScene &scene) {
  packed_transform_offsets.clear();
  WFMotionHeader header{};
  std::vector<unsigned char> bytes(sizeof(header));
  auto append = [&](const void *data, size_t size) {
    size_t offset = (bytes.size() + 15) & ~size_t(15);
    if (offset > UINT32_MAX || size > UINT32_MAX - offset)
      throw std::runtime_error("Metal motion metadata exceeds 32-bit capacity");
    bytes.resize(offset + size);
    if (size)
      std::memcpy(bytes.data() + offset, data, size);
    return uint32_t(offset);
  };
  auto counts = sizes(scene);
  for (size_t k = 0; k < 5; ++k) {
    owners[k].resize(counts[k]);
    header.owners[k] = append(owners[k].data(), owners[k].size() * sizeof(uint32_t));
  }
  // Deduplicate property storage by immutable identity and epoch. Distinct
  // instances can share keys while retaining independent owner transforms.
  std::map<std::pair<const animation::Track *, double>, uint32_t> track_offsets;
  auto export_track = [&](const std::shared_ptr<const animation::Track> &track, double epoch) {
    if (!track) return uint32_t(0);
    auto key = std::make_pair(track.get(), epoch);
    auto found = track_offsets.find(key);
    if (found != track_offsets.end()) return found->second;
    std::vector<WFTrackKey> keys(track->times.size());
    std::vector<double> quantized_times;
    std::vector<animation::Value> quantized_values, quantized_incoming, quantized_outgoing;
    auto vector = [](const animation::Value &v) { return WFVector{float(v[0]), float(v[1]), float(v[2]), float(v[3])}; };
    for (size_t i = 0; i < keys.size(); ++i) {
      double relative = track->times[i] - epoch;
      keys[i].time = float(relative);
      double interval = keys.size() == 1 ? 1. : i == 0 ? track->times[1] - track->times[0] :
        track->times[i] - track->times[i - 1];
      if (i + 1 < keys.size()) interval = std::min(interval, track->times[i + 1] - track->times[i]);
      if (!std::isfinite(keys[i].time) || std::abs(double(keys[i].time) - relative) > interval * 1e-5 ||
          (i && !(keys[i].time > keys[i - 1].time)))
        throw std::runtime_error("Device key times lose precision at this scene epoch; use CPU or a nearer time interval.");
      keys[i].value = vector(track->values[i]);
      quantized_times.push_back(keys[i].time);
      auto native = [](WFVector v) { return animation::Value{v.x, v.y, v.z, v.w}; };
      quantized_values.push_back(native(keys[i].value));
      if (!track->incoming.empty()) {
        keys[i].incoming = vector(track->incoming[i]); keys[i].outgoing = vector(track->outgoing[i]);
        quantized_incoming.push_back(native(keys[i].incoming));
        quantized_outgoing.push_back(native(keys[i].outgoing));
        for (int k = 0; k < 4; ++k)
          if (!std::isfinite(float(track->incoming[i][k])) || !std::isfinite(float(track->outgoing[i][k])))
            throw std::runtime_error("Device cubic tangents exceed float precision range.");
      }
    }
    auto method = track->interpolation == animation::Interpolation::Smooth ? animation::Interpolation::Cubic : track->interpolation;
    animation::Track quantized(track->property, method, std::move(quantized_times), std::move(quantized_values),
      std::move(quantized_incoming), std::move(quantized_outgoing));
    quantized.ValidateFloatPacking();
    WFPropertyTrack descriptor{uint32_t(keys.size()), uint32_t(track->interpolation), uint32_t(track->property),
      append(keys.data(), keys.size() * sizeof(WFTrackKey))};
    auto offset = append(&descriptor, sizeof(descriptor));
    track_offsets.emplace(key, offset);
    return offset;
  };
  std::vector<WFMotionNode> animations;
  for (const auto &node : nodes) {
    WFMotionNode a{};
    a.parent = node.parent;
    a.start = double(node.animation->StartTime()) - node.animation->TimeOrigin();
    a.end = double(node.animation->EndTime()) - node.animation->TimeOrigin();
    pack(node.before.GetMatrix(), a.before);
    pack(node.reference_inverse.GetMatrix(), a.reference_inverse);
    if (const auto &track = node.animation->Track()) {
      WFTransformTrack descriptor{};
      double key_epoch = node.animation->TimeOrigin();
      if (track->clock) {
        const auto &clock = *track->clock;
        key_epoch = clock.begin;
        descriptor.clock = {float(clock_phase(clock, node.animation->FrameEpoch())),
          clock.end == clock.begin ? 0.f : float(clock.speed), float(clock.end - clock.begin), float(clock.loop)};
        descriptor.frame.x = node.animation->FrameAnchor();
      }
      descriptor.channels[0] = export_track(track->position, key_epoch);
      descriptor.channels[1] = export_track(track->rotation, key_epoch);
      descriptor.channels[2] = export_track(track->scale, key_epoch);
      descriptor.pivot = {float(track->pivot[0]), float(track->pivot[1]), float(track->pivot[2]), 0};
      a.track = append(&descriptor, sizeof(descriptor));
    }
    Quaternion rotations[2];
    for (int k = 0; k < 2; ++k) {
      const auto &t = k ? node.animation->EndTransform() : node.animation->StartTransform();
      pack(t.GetMatrix(), a.endpoints[k]);
      vec3f translation;
      Matrix4x4 scale;
      AnimatedTransform::Decompose(t.GetMatrix(), &translation, &rotations[k], &scale);
      a.translation[k] = {float(translation[0]), float(translation[1]), float(translation[2]), 0};
      pack(scale, a.scale[k]);
    }
    if (dot(rotations[0], rotations[1]) < 0)
      rotations[1] = -rotations[1];
    for (int k = 0; k < 2; ++k)
      a.rotation[k] = {float(rotations[k].v[0]), float(rotations[k].v[1]), float(rotations[k].v[2]),
                       float(rotations[k].w)};
    animations.push_back(a);
    packed_transform_offsets.push_back(a.track);
  }
  header.count = nodes.size();
  header.nodes = append(animations.data(), animations.size() * sizeof(WFMotionNode));

  // Metal instance masks provide eight disjoint hardware traversal partitions.
  // Keep the existing owner-filtered path for larger hierarchies or animated
  // loose triangles, whose shared BLAS cannot be masked independently.
  const size_t loose_count = scene.meshes.empty() ? scene.triangles.size() : scene.loose_triangles;
  const bool static_loose = std::all_of(owners[0].begin(), owners[0].begin() + loose_count,
                                       [](uint32_t owner) { return owner == 0; });
  scene.texture_features &= ~4096u;
  if (!scene.meshes.empty() && nodes.size() < 8 && static_loose)
    scene.texture_features |= 4096;


  // Bounds in each owner's reference frame reject unrelated animation groups
  // before querying the shared immutable hardware BVH. Prototype bounds are
  // computed once, irrespective of the number of placements.
  std::vector<aabb> groups(nodes.size() + 1), meshes(scene.meshes.size());
  for (size_t i = 0; i < scene.meshes.size(); ++i) {
    const auto &mesh = scene.meshes[i];
    for (uint32_t j = 0; j < mesh.count; ++j)
      for (const auto &p : scene.triangles[mesh.start + j].p)
        meshes[i] = surrounding_box(meshes[i], aabb(point3f(p.x, p.y, p.z)));
  }
  size_t loose = scene.meshes.empty() ? scene.triangles.size() : scene.loose_triangles;
  for (size_t i = 0; i < loose; ++i)
    for (const auto &p : scene.triangles[i].p)
      groups[owners[0][i]] = surrounding_box(groups[owners[0][i]], aabb(point3f(p.x, p.y, p.z)));
  for (size_t i = 0; i < scene.quadrics.size(); ++i) {
    const auto &lo = scene.quadric_bounds[2 * i], &hi = scene.quadric_bounds[2 * i + 1];
    groups[owners[1][i]] =
        surrounding_box(groups[owners[1][i]], aabb(point3f(lo.x, lo.y, lo.z), point3f(hi.x, hi.y, hi.z)));
  }
  for (size_t i = 0; i < scene.instances.size(); ++i) {
    const auto &instance = scene.instances[i];
    if (instance.triangle_count)
      groups[owners[3][i]] = surrounding_box(groups[owners[3][i]], unpack(instance.forward)(meshes[instance.mesh]));
  }
  std::vector<WFVector> boxes;
  for (const auto &box : groups) {
    auto lo = box.min(), hi = box.max();
    if (lo[0] > hi[0]) {
      boxes.push_back({INFINITY, INFINITY, INFINITY, 0});
      boxes.push_back({-INFINITY, -INFINITY, -INFINITY, 0});
      continue;
    }
    Float pad = 64 * std::numeric_limits<float>::epsilon() *
                std::max(Float(1), std::max(MaxComponent(Abs(lo)), MaxComponent(Abs(hi))));
    boxes.push_back({float(lo[0] - pad), float(lo[1] - pad), float(lo[2] - pad), 0});
    boxes.push_back({float(hi[0] + pad), float(hi[1] + pad), float(hi[2] + pad), 0});
  }
  header.bounds = append(boxes.data(), boxes.size() * sizeof(WFVector));
  header.byte_size = bytes.size();
  std::memcpy(bytes.data(), &header, sizeof(header));

  // Membership and diffusion bounds enclose the full native motion, including
  // nested rotation. They are only broad-phase bounds; membership still tests
  // the exact geometry at the path's time.
  std::set<std::pair<uint32_t, uint32_t>> bodies;
  auto collect = [&](const auto &values, size_t kind) {
    for (size_t i = 0; i < values.size(); ++i)
      if (values[i].boundary && owners[kind][i])
        bodies.emplace(values[i].boundary, owners[kind][i]);
  };
  collect(scene.triangles, 0);
  collect(scene.quadrics, 1);
  collect(scene.instances, 3);
  const auto reference_bounds = scene.boundaries;
  for (const auto &body : bodies) {
    const auto &b = reference_bounds[body.first];
    if (b.lo.w != 0)
      continue; // Open interfaces have no interior bounds to sweep.
    aabb box(point3f(b.lo.x, b.lo.y, b.lo.z), point3f(b.hi.x, b.hi.y, b.hi.z));
    for (uint32_t owner = body.second; owner; owner = nodes[owner - 1].parent) {
      const auto &node = nodes[owner - 1];
      box = node.before(node.animation->MotionBounds(node.reference_inverse(box)));
    }
    auto lo = box.min(), hi = box.max();
    float pad = 64 * std::numeric_limits<float>::epsilon();
    for (int k = 0; k < 3; ++k)
      pad =
          std::max(pad, 64 * std::numeric_limits<float>::epsilon() * float(std::max(std::abs(lo[k]), std::abs(hi[k]))));
    WFVector l{float(lo[0]) - pad, float(lo[1]) - pad, float(lo[2]) - pad, 0};
    WFVector h{float(hi[0]) + pad, float(hi[1]) + pad, float(hi[2]) + pad, 0};
    extend(scene.boundaries[body.first].lo, scene.boundaries[body.first].hi, l);
    extend(scene.boundaries[body.first].lo, scene.boundaries[body.first].hi, h);
    extend(scene.boundary_lo, scene.boundary_hi, l);
    extend(scene.boundary_lo, scene.boundary_hi, h);
  }
  return bytes;
}
