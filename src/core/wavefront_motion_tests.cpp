#ifdef NOT_CRAN
#include "wavefront_motion.h"
#include "wavefront_motion_layout.h"
#include "hlbvh.h"
#include <testthat.h>
#include <cstring>

namespace {
template<class T> T Read(const std::vector<unsigned char> &bytes, size_t offset) {
  if (offset % 16 || offset > bytes.size() || sizeof(T) > bytes.size() - offset)
    throw std::runtime_error("Invalid packed animation offset.");
  T result;
  std::memcpy(&result, bytes.data() + offset, sizeof(T));
  return result;
}
animation::Value Unpack(WFVector v) { return {v.x, v.y, v.z, v.w}; }
std::shared_ptr<const animation::Track> Decode(const std::vector<unsigned char> &bytes, uint32_t offset) {
  if (!offset) return nullptr;
  auto descriptor = Read<WFPropertyTrack>(bytes, offset);
  std::vector<double> times;
  std::vector<animation::Value> values, incoming, outgoing;
  for (uint32_t i = 0; i < descriptor.count; ++i) {
    auto key = Read<WFTrackKey>(bytes, descriptor.keys + size_t(i) * sizeof(WFTrackKey));
    times.push_back(key.time); values.push_back(Unpack(key.value));
    if (descriptor.interpolation >= 2) { incoming.push_back(Unpack(key.incoming)); outgoing.push_back(Unpack(key.outgoing)); }
  }
  auto method = descriptor.interpolation == 3 ? animation::Interpolation::Cubic :
    static_cast<animation::Interpolation>(descriptor.interpolation);
  return std::make_shared<const animation::Track>(static_cast<animation::Property>(descriptor.property),
    method, times, values, incoming, outgoing);
}
}
context("Packed scene-time tracks") {
  test_that("frame clock publication changes phase without changing resident keys") {
    using namespace animation;
    auto track = std::make_shared<TransformTrack>();
    track->position = std::make_shared<Track>(Property::Position, Interpolation::Linear,
      std::vector<double>{0, 1, 2}, std::vector<Value>{{-1,0,0,0},{0,0,0,0},{1,0,0,0}});
    track->clock = std::make_shared<ClipClock>(0,2,0,0,1,Loop::Repeat);
    auto frame = std::make_shared<SceneTimeContext>();
    frame->epoch = -.001;
    Transform identity;
    AnimatedTransform animation(&identity,0,&identity,0,track,-.001,frame);
    WavefrontMotion motion;
    motion.Add(0,animation,identity);
    WavefrontScene scene;
    motion.Seal(scene);
    auto header = Read<WFMotionHeader>(motion.packed,0);
    auto node = Read<WFMotionNode>(motion.packed,header.nodes);
    auto original = motion.packed;
    for (double time : {4., 0., 2., -.001, 4.}) {
      frame->epoch = time;
      frame->ray_anchor = Float(time + .001);
      motion.ValidateTimes(frame->ray_anchor,frame->ray_anchor);
      motion.UpdateClocks();
      auto descriptor = Read<WFTransformTrack>(motion.packed,node.track);
      expect_true(descriptor.frame.x == frame->ray_anchor);
      expect_true(std::abs(descriptor.clock.x - track->clock->At(time)) < 1e-6);
      auto copy = original;
      motion.CopyClockState(copy.data());
      expect_true(copy == motion.packed);
      size_t clock_offset = node.track + offsetof(WFTransformTrack,clock);
      std::copy(original.begin()+clock_offset,original.begin()+clock_offset+2*sizeof(WFVector),copy.begin()+clock_offset);
      expect_true(copy == original);
    }
    expect_error(motion.ValidateTimes(-1e10,1e10));
  }

  test_that("versioned storage preserves asynchronous curves and transformed ray parameters") {
    using namespace animation;
    auto authored = std::make_shared<TransformTrack>();
    double epoch = 1e12;
    authored->position = std::make_shared<Track>(Property::Position, Interpolation::Cubic,
      std::vector<double>{epoch, epoch + .25, epoch + .75, epoch + 1},
      std::vector<Value>{{0,0,0,0}, {1,2,0,0}, {-1,1,1,0}, {0,0,0,0}});
    authored->rotation = std::make_shared<Track>(Property::Rotation, Interpolation::Smooth,
      std::vector<double>{epoch, epoch + .5, epoch + 1},
      std::vector<Value>{{0,0,0,1}, {0,.6,0,.8}, {0,0,0,1}});
    authored->scale = std::make_shared<Track>(Property::Scale, Interpolation::Linear,
      std::vector<double>{epoch, epoch + 1}, std::vector<Value>{{1,2,-1,0}, {2,1,-.5,0}});
    authored->pivot = {.3, .1, -.2, 0};
    Transform identity;
    AnimatedTransform animation(&identity, 0, &identity, 0, authored, epoch);
    WavefrontMotion motion;
    expect_true(motion.Add(0, animation, identity) == 1);
    expect_true(motion.Add(0, animation, identity) == 1);
    WavefrontScene scene;
    scene.quadrics.resize(1); scene.quadric_bounds = {{-1,-1,-1,0}, {1,1,1,0}};
    motion.owners[1] = {1};
    motion.Seal(scene);
    auto header = Read<WFMotionHeader>(motion.packed, 0);
    expect_true(header.version == 3);
    expect_true(header.byte_size == motion.packed.size());
    auto node = Read<WFMotionNode>(motion.packed, header.nodes);
    auto descriptor = Read<WFTransformTrack>(motion.packed, node.track);
    expect_true(descriptor.version == 2);
    TransformTrack decoded;
    decoded.position = Decode(motion.packed, descriptor.channels[0]);
    decoded.rotation = Decode(motion.packed, descriptor.channels[1]);
    decoded.scale = Decode(motion.packed, descriptor.channels[2]);
    decoded.pivot = Unpack(descriptor.pivot);
    for (double time : {-.5, 0., .125, .25, .5, .625, .75, 1., 2., .125}) {
      auto expected = authored->At(epoch + time), actual = decoded.At(time);
      for (int r = 0; r < 4; ++r) for (int c = 0; c < 4; ++c)
        expect_true(std::abs(expected.GetMatrix().m[r][c] - actual.GetMatrix().m[r][c]) < 3e-6);
      // A transformed ray must retain its parameter even under reflection and
      // nonuniform scale; normalizing its inverse-transformed direction fails.
      point3f local(.2, .3, .4), world = actual(local);
      vec3f direction(.3, -.1, .7);
      Ray ray(world - 2.5f * direction, direction, time);
      Ray inverse_ray = Inverse(actual)(ray);
      expect_true(inverse_ray.time() == ray.time());
      expect_true(inverse_ray.WithDirection(vec3f(1, 0, 0)).time() == ray.time());
      auto hit = inverse_ray(2.5f);
      for (int k = 0; k < 3; ++k) expect_true(std::abs(hit[k] - local[k]) < 2e-6);
    }
  }
  test_that("packing rejects collapsed nearby keys instead of changing the curve") {
    auto track = std::make_shared<animation::TransformTrack>();
    track->position = std::make_shared<animation::Track>(animation::Property::Position, animation::Interpolation::Linear,
      std::vector<double>{1e12, 1e12 + .25}, std::vector<animation::Value>{{0,0,0,0}, {1,0,0,0}});
    Transform identity;
    AnimatedTransform animation(&identity, 0, &identity, 0, track);
    WavefrontMotion motion;
    motion.Add(0, animation, identity);
    WavefrontScene scene;
    expect_error(motion.Seal(scene));
    WavefrontMotion invalid;
    expect_error(invalid.Add(1, animation, identity));
    invalid.owners[0] = {1};
    expect_error(invalid.Seal(scene));
  }
#ifdef RAY_HAS_METAL_BVH
  test_that("device keyframe intersections match native transforms at keys and midpoints") {
    if (!MetalBVHAvailable()) {
      Rcpp::Rcout << "Keyframe device probes: unavailable\n";
      return;
    }
    using namespace animation;
    for (int clock_case = 0; clock_case < 7; ++clock_case) {
      auto track = std::make_shared<TransformTrack>();
      track->position = std::make_shared<Track>(Property::Position, Interpolation::Cubic,
        std::vector<double>{0, .7, 2.3, 4}, std::vector<Value>{{0,0,0,0}, {1,2,0,0}, {-1,1,1,0}, {0,0,0,0}});
      track->rotation = std::make_shared<Track>(Property::Rotation, Interpolation::Smooth,
        std::vector<double>{0, 1.3, 4}, std::vector<Value>{{0,0,0,1}, {0,.6,0,.8}, {0,0,0,1}});
      track->scale = std::make_shared<Track>(Property::Scale, Interpolation::Linear,
        std::vector<double>{0,4}, std::vector<Value>{{1,1,1,0}, {1.5,.7,1,0}});
      if (clock_case) {
        Loop mode = clock_case == 6 ? Loop::Hold : clock_case % 2 ? Loop::Repeat : Loop::PingPong;
        double speed = clock_case == 5 ? 0 : clock_case >= 3 && clock_case <= 4 ? -1 : 1;
        track->clock = std::make_shared<ClipClock>(0, 4, 0, clock_case == 5 ? .5 : 0, speed, mode);
      }
      Transform identity;
      AnimatedTransform animation(&identity, 0, &identity, 0, track);
      WavefrontScene scene;
      scene.motion = std::make_shared<WavefrontMotion>();
      scene.motion->Add(0, animation, identity);
      scene.materials.resize(1); scene.surfaces.resize(1);
      auto &surface = scene.surfaces[0];
      for (int k = 0; k < 3; ++k) {
        WFVector row{}; (&row.x)[k] = 1;
        surface.forward[k] = surface.inverse[k] = surface.bump_forward[k] = surface.bump_inverse[k] = row;
      }
      WFTriangle triangle{};
      triangle.p[0] = {-2,-2,0,0}; triangle.p[1] = {2,-2,0,0}; triangle.p[2] = {0,2,0,0};
      auto reference = track->At(0);
      for (auto &p : triangle.p) {
        auto point = reference(point3f(p.x,p.y,p.z));
        p = {point[0],point[1],point[2],0};
      }
      for (auto &n : triangle.n) n = {0,0,1,0};
      triangle.geometric = {0,0,1,0}; triangle.uv01 = {0,0,1,0}; triangle.uv2 = {.5,1,0,0};
      scene.triangles.push_back(triangle);
      scene.motion->owners[0] = {1}; scene.motion->owners[2] = {1}; scene.motion->Seal(scene);
      std::vector<WFGeometryProbe> probes;
      std::vector<point3f> points;
      scene.motion->ValidateTimes(-4.001, 8.001);
      for (float time : {-4.001f, -4.f, -3.999f, -1.f, 0.f, .125f, .7f, 1.3f, 2.3f, 3.f, 3.999f, 4.f, 4.001f, 5.f, 8.f}) {
        auto pose = track->At(time);
        point3f target = pose(point3f(0,0,0));
        vec3f direction = unit_vector(pose(vec3f(0,0,-1)));
        point3f origin = target - 2.5f * direction;
        WFGeometryProbe probe;
        probe.origin = {origin[0],origin[1],origin[2],0};
        probe.direction = {direction[0],direction[1],direction[2],10};
        std::memcpy(&probe.padding[1], &time, sizeof(time));
        probes.push_back(probe); points.push_back(target);
      }
      auto hits = ProbeMetalGeometry(scene, probes);
      expect_true(hits.size() == probes.size());
      for (size_t i = 0; i < hits.size(); ++i) {
        expect_true(hits[i].normal.w != 0);
        expect_true(std::abs(hits[i].point.w - 2.5f) < 5e-5);
        for (int k = 0; k < 3; ++k) expect_true(std::abs((&hits[i].point.x)[k] - points[i][k]) < 5e-5);
      }
    }
    Rcpp::Rcout << "Keyframe device probes executed: 105 scene-time and loop-boundary hits\n";
  }
#endif
}
#endif
