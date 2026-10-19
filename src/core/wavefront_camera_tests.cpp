#if defined(NOT_CRAN) && defined(RAY_HAS_METAL_BVH)
#include "camera.h"
#include "camera_track.h"
#include "hlbvh.h"
#include "wavefront.h"
#include "wavefront_motion.h"
#include <testthat.h>

namespace {
void compare_camera(RayCamera &camera) {
  WFCamera exported;
  std::vector<WFVector> data;
  ExportWavefrontCamera(camera, exported, data);
  // Production rendering writes atmosphere bookkeeping after camera export.
  // Camera probes must remain independent of those shared runtime fields.
  exported.padding[0] = 1024;
  exported.padding[1] = 4;
  std::vector<WFCameraProbe> probes;
  for (int i = 0; i < 127; ++i)
    probes.push_back(
        {{float((i % 11 + .3) / 11), float((i % 13 + .7) / 13), float((i % 7 + .5) / 7), float((i % 17 + .5) / 17)},
         {float((i % 19) / 18.), 0, 0, 0}});
  auto results = ProbeMetalCamera(exported, data, probes);
  double max_position = 0, max_direction = 0, max_weight = 0;
  size_t acceptance = 0, valid = 0;
  for (size_t i = 0; i < probes.size(); ++i) {
    auto p = probes[i];
    Ray ray;
    Float weight = 1;
    if (camera.get_fov() < 0)
      weight = camera.GenerateRay(
          CameraSample(point2f(1 - p.film_lens.x, 1 - p.film_lens.y), point2f(p.film_lens.z, p.film_lens.w), p.time.x),
          &ray);
    else {
      float radius = std::sqrt(p.film_lens.z), phi = 2 * M_PI * p.film_lens.w;
      ray = camera.get_ray(p.film_lens.x, p.film_lens.y, point3f(radius * std::cos(phi), radius * std::sin(phi), 0),
                           p.time.x);
    }
    const auto &r = results[i];
    acceptance += (weight > 0) != (r.origin.w > 0);
    if (!(weight > 0) || !(r.origin.w > 0))
      continue;
    ++valid;
    max_weight = std::max(max_weight, std::abs(double(weight - r.origin.w)) / std::max(1e-8, double(weight)));
    auto d = unit_vector(ray.d);
    max_position = std::max(max_position, double((ray.o - point3f(r.origin.x, r.origin.y, r.origin.z)).length()));
    max_direction = std::max(max_direction, double((d - vec3f(r.direction.x, r.direction.y, r.direction.z)).length()));
  }
  expect_true(valid > 0);
  expect_true(acceptance == 0);
  expect_true(max_position < 1e-4);
  expect_true(max_direction < 1e-4);
  expect_true(max_weight < 1e-4);
}
} // namespace
context("Scene-wide motion identity") {
  test_that("projection cuts sample the correct camera on either side of the key") {
    if(!MetalBVHAvailable()) return;
    animation::CameraKey a{0,{0,0,-10, 0,0,0, .1,60,10, 4,3, 0,1,0}};
    auto b=a,c=a,d=a;
    b.time=1; b.values[7]=0; b.values[0]=2;
    c.time=2; c.values[7]=360;
    d.time=3;
    auto track=std::make_shared<const animation::CameraTrack>(std::vector<animation::CameraKey>{a,b,c,d});
    auto clock=std::make_shared<SceneTimeContext>();
    for (const auto &key : {a,b,c}) {
      auto cam=track->MakeCamera(key,1.5,1,0,3,false);
      cam->scene_camera_track=track; cam->scene_camera_time=clock;
      compare_camera(*cam);
      cam->set_camera_flip_x(true);
      compare_camera(*cam);
    }
    auto cam=track->MakeCamera(a,1.5,1,0,3,false);
    cam->scene_camera_track=track; cam->scene_camera_time=clock;
    clock->epoch=.1;
    WFCamera descriptor;
    std::vector<WFVector> data;
    expect_error(ExportWavefrontCamera(*cam,descriptor,data));
  }
  test_that("multi-key camera exposures agree with device rays") {
    if(!MetalBVHAvailable()) return;
    for (double projection : {60.,0.,360.}) {
      animation::CameraKey a{100,{0,0,-10, 0,0,0, .2,projection,10, 4,3, 0,1,0}};
      auto b=a,c=a,d=a;
      b.time=100.13; b.values[0]=2; b.values[4]=1; b.values[12]=.8;
      c.time=100.37; c.values[0]=-2; c.values[6]=.5; c.values[9]=6;
      d.time=100.6;
      auto track=std::make_shared<const animation::CameraTrack>(std::vector<animation::CameraKey>{a,b,c,d});
      auto clock=std::make_shared<SceneTimeContext>(); clock->epoch=100; clock->ray_anchor=5;
      std::unique_ptr<RayCamera> cam;
      if(projection==60) cam=std::make_unique<camera>(point3f(0,0,-10),point3f(0),vec3f(0,1,0),60,1.5,.2,10,5,5.6,1);
      else if(projection==0) cam=std::make_unique<ortho_camera>(point3f(0,0,-10),point3f(0),vec3f(0,1,0),4,3,5,5.6,1);
      else cam=std::make_unique<environment_camera>(point3f(0,0,-10),point3f(0),vec3f(0,1,0),5,5.6,1);
      cam->set_shutter_speed(1);
      cam->scene_camera_track=track; cam->scene_camera_time=clock;
      compare_camera(*cam);
      cam->set_camera_flip_x(true);
      compare_camera(*cam);
    }
  }
  test_that("nonadjacent repeats preserve exact transforms, times and parents") {
    Transform start, end = Translate(vec3f(1, 2, 3)) * RotateY(35) * Scale(.7, 1.3, .9);
    Transform other_end = Translate(vec3f(2, 1, -1));
    AnimatedTransform first(&start, 0, &end, 1), copy(&start, 0, &end, 1);
    AnimatedTransform other(&start, 0, &other_end, 1);
    WavefrontMotion motion;
    const auto first_id = motion.Add(0, first, Transform());
    const auto other_id = motion.Add(0, other, Transform());
    expect_true(first_id == 1);
    expect_true(other_id == 2);
    expect_true(motion.Add(0, copy, Transform()) == first_id);
    expect_true(motion.nodes.size() == 2);

    // Equal numerical values must share a hash even when their bits differ.
    Matrix4x4 signed_zero;
    signed_zero.m[0][1] = -Float(0);
    Transform zero_transform(signed_zero, signed_zero);
    AnimatedTransform zero_time(&start, -Float(0), &end, 1);
    expect_true(motion.Add(0, zero_time, zero_transform) == first_id);

    AnimatedTransform later_start(&start, .1, &end, 1);
    AnimatedTransform later_end(&start, 0, &end, 2);
    Transform shifted_start = Translate(vec3f(.01, 0, 0));
    AnimatedTransform different_start(&shifted_start, 0, &end, 1);
    const auto start_id = motion.Add(0, later_start, Transform());
    const auto end_id = motion.Add(0, later_end, Transform());
    const auto endpoint_id = motion.Add(0, different_start, Transform());
    const auto placement_id = motion.Add(0, first, Translate(vec3f(.01, 0, 0)));
    const auto parent_id = motion.Add(other_id, first, Transform());
    expect_true(start_id == 3);
    expect_true(end_id == 4);
    expect_true(endpoint_id == 5);
    expect_true(placement_id == 6);
    expect_true(parent_id == 7);
    expect_true(motion.Add(0, copy, Transform()) == first_id);
    expect_true(motion.Add(other_id, copy, Transform()) == parent_id);

    // Explicit inverses participate in Transform::operator== independently
    // from the forward matrix, even if a caller supplies an unusual pair.
    Matrix4x4 inverse = end.GetInverseMatrix();
    inverse.m[0][3] = std::nextafter(inverse.m[0][3], std::numeric_limits<Float>::infinity());
    Transform different_inverse(end.GetMatrix(), inverse);
    AnimatedTransform inverse_motion(&start, 0, &different_inverse, 1);
    const auto inverse_id = motion.Add(0, inverse_motion, Transform());
    expect_true(inverse_id == 8);
    Transform before_inverse(start.GetMatrix(), inverse);
    expect_true(motion.Add(0, first, before_inverse) == 9);
    expect_true(motion.Add(0, inverse_motion, Transform()) == inverse_id);

    for (Float time : {Float(-.2), Float(0), Float(.37), Float(1), Float(1.2)}) {
      Transform a, b;
      first.Interpolate(time, &a);
      other.Interpolate(time, &b);
      const auto expected = (b * a)(point3f(.2, -.5, 1));
      const auto actual = motion.At(parent_id, time)(point3f(.2, -.5, 1));
      expect_true((expected - actual).length() < 1e-5);
    }
    WavefrontScene scene;
    scene.triangles.resize(3);
    motion.Claim({0, 0, 0, 0, 0}, {1, 0, 0, 0, 0}, first_id);
    motion.Claim({1, 0, 0, 0, 0}, {2, 0, 0, 0, 0}, other_id);
    motion.Claim({2, 0, 0, 0, 0}, {3, 0, 0, 0, 0}, motion.Add(0, copy, Transform()));
    motion.Seal(scene);
    expect_true(motion.owners[0] == std::vector<uint32_t>({first_id, other_id, first_id}));
  }
  test_that("motion IDs survive lookup rehashing and node vector growth") {
    constexpr size_t count = 2048;
    Transform start;
    std::vector<Transform> ends;
    std::vector<AnimatedTransform> animations;
    ends.reserve(count);
    animations.reserve(count);
    WavefrontMotion motion;
    for (size_t i = 0; i < count; ++i) {
      ends.push_back(Translate(vec3f(Float(i + 1), .2, -.3)));
      animations.emplace_back(&start, 0, &ends.back(), 1);
      expect_true(motion.Add(0, animations.back(), Transform()) == i + 1);
    }
    for (size_t i = count; i-- > 0;)
      expect_true(motion.Add(0, animations[i], Transform()) == i + 1);
    expect_true(motion.nodes.size() == count);
  }
}

context("Metal hardware motion partitions") {
  test_that("masked traversal matches owner filtering for hits and boundaries") {
    if (!MetalBVHAvailable()) return;
    for (int group_count : {1, 7, 8}) {
      WavefrontScene scene;
      scene.motion = std::make_shared<WavefrontMotion>();
      auto &motion = *scene.motion;
      std::vector<Transform> starts(group_count), ends;
      std::vector<AnimatedTransform> animations;
      ends.reserve(group_count);
      animations.reserve(group_count);
      for (int i = 0; i < group_count; ++i) {
        ends.push_back(Translate(vec3f(.1 * (i + 1), -.1, .3)) * RotateY(10 + i));
        animations.emplace_back(&starts[i], 0, &ends[i], 1);
        motion.Add(0, animations.back(), Transform());
      }
      // One static loose triangle and one shared triangle prototype. Include
      // static and moving instances and analytic spheres in the same TLAS.
      WFTriangle triangle;
      triangle.p[0] = {-1, -1, 0, 0};
      triangle.p[1] = {1, -1, 0, 0};
      triangle.p[2] = {0, 1, 0, 0};
      for (auto &n : triangle.n) n = {0, 0, 1, 0};
      triangle.geometric = {0, 0, 1, 0};
      triangle.uv01 = {0, 0, 1, 0};
      triangle.uv2 = {.5, 1, 0, 0};
      scene.triangles = {triangle, triangle};
      for (auto &p : scene.triangles[0].p) p.z = -2;
      scene.loose_triangles = 1;
      scene.meshes.push_back({1, 1});
      scene.instances.resize(2); // Loose triangle and analytic BLAS prefixes.
      motion.owners[3] = {0, 0};
      scene.boundaries.resize(2);
      for (int axis = 0; axis < 3; ++axis) {
        WFVector row{};
        (&row.x)[axis] = 1;
        scene.surfaces[0].forward[axis] = scene.surfaces[0].inverse[axis] = row;
      }
      scene.surfaces[0].axes = {1, 1, 1, 1};
      for (int owner = 0; owner <= group_count; ++owner) {
        const float x = (owner - group_count / 2.f) * 1.1f;
        WFMeshInstance instance;
        instance.triangle_start = 1;
        instance.triangle_count = 1;
        instance.boundary = owner % 2;
        for (int axis = 0; axis < 3; ++axis)
          instance.forward[axis] = instance.inverse[axis] = scene.surfaces[0].forward[axis];
        instance.forward[0].w = x;
        instance.inverse[0].w = -x;
        scene.instances.push_back(instance);
        motion.owners[3].push_back(owner);
        WFQuadric sphere;
        for (int axis = 0; axis < 3; ++axis)
          sphere.forward[axis] = sphere.inverse[axis] = scene.surfaces[0].forward[axis];
        sphere.forward[0].w = x;
        sphere.inverse[0].w = -x;
        sphere.forward[1].w = -1.4f;
        sphere.inverse[1].w = 1.4f;
        sphere.measure = {1, 1, 1, 1};
        sphere.boundary = owner % 2;
        scene.quadrics.push_back(sphere);
        scene.quadric_bounds.push_back({x - 1, -2.4f, -1, 0});
        scene.quadric_bounds.push_back({x + 1, -.4f, 1, 0});
        motion.owners[1].push_back(owner);
      }
      motion.Seal(scene);
      expect_true(bool(scene.texture_features & 4096) == (group_count < 8));
      std::vector<WFGeometryProbe> probes;
      for (float time : {0.f, .23f, .71f, 1.f})
        for (uint32_t body : {0u, 1u})
          for (int i = 0; i < 1024; ++i) {
            WFGeometryProbe probe;
            probe.origin = {float(i % 64 - 32) * .17f, float(i / 64 - 8) * .25f, 5, 0};
            probe.direction = {0, 0, -1, 12};
            probe.body = body;
            std::memcpy(&probe.padding[1], &time, sizeof(time));
            probes.push_back(probe);
          }
      auto partitioned = ProbeMetalGeometry(scene, probes);
      scene.texture_features &= ~4096u;
      auto reference = ProbeMetalGeometry(scene, probes);
      size_t errors = 0, hits = 0;
      double maximum = 0;
      for (size_t i = 0; i < probes.size(); ++i) {
        const auto &a = partitioned[i], &b = reference[i];
        errors += (a.normal.w != b.normal.w) || a.primitive != b.primitive || a.boundary != b.boundary;
        if (!a.normal.w || !b.normal.w) continue;
        ++hits;
        for (int k = 0; k < 4; ++k) {
          maximum = std::max(maximum, std::abs(double((&a.point.x)[k]) - (&b.point.x)[k]));
          maximum = std::max(maximum, std::abs(double((&a.normal.x)[k]) - (&b.normal.x)[k]));
          maximum = std::max(maximum, std::abs(double((&a.uv.x)[k]) - (&b.uv.x)[k]));
        }
      }
      Rcpp::Rcout << "Motion mask groups=" << group_count << " hits=" << hits
                  << " errors=" << errors << " maximum difference=" << maximum << "\n";
      expect_true(hits > 0);
      expect_true(errors == 0);
      expect_true(maximum < 1e-6);
      // Animated loose triangles cannot be given an independent instance mask.
      // Re-sealing must also clear a previously enabled optimization flag.
      motion.owners[0][0] = 1;
      motion.Seal(scene);
      expect_true((scene.texture_features & 4096) == 0);
    }
  }
}

context("Metal cameras and density threshold intervals") {
  test_that("motion ownership coalesces identical animations and retains native affine parameters") {
    Transform start = Translate(vec3f(2, -1, 3));
    Transform end = Translate(vec3f(-1, 2, 4)) * RotateY(155) * Scale(.8, 1.2, 1.1);
    AnimatedTransform animation(&start, 0, &end, 1);
    AnimatedTransform duplicate(&start, 0, &end, 1);
    AnimatedTransform different_time(&start, 0, &end, 2);
    WavefrontMotion motion;
    expect_true(motion.Add(0, animation, Transform()) == 1);
    expect_true(motion.Add(0, duplicate, Transform()) == 1);
    WavefrontScene scene;
    scene.triangles.resize(1);
    scene.triangles[0].p[0] = {1, 2, 3, 0};
    motion.owners[0] = {1};
    motion.Seal(scene);
    expect_true((scene.texture_features & 1024) != 0);
    // Keep t in world-distance units under nonuniform scaling. A normalized
    // reference-frame direction would silently move intersections and shadows.
    double maximum_error = 0;
    for (int i = 0; i <= 100; ++i) {
      float time = i / 100.f;
      const auto inverse = motion.At(1, time).GetInverseMatrix();
      Transform pose;
      animation.Interpolate(time, &pose);
      Transform delta = pose * Inverse(start);
      point3f reference(1, 2, 3);
      point3f world = delta(reference);
      vec3f direction(.2, -.3, .7);
      point3f origin = world - Float(2.5) * direction;
      for (int a = 0; a < 3; ++a) {
        const auto *r = inverse.m[a];
        double local_origin = r[0] * origin[0] + r[1] * origin[1] + r[2] * origin[2] + r[3];
        double local_direction = r[0] * direction[0] + r[1] * direction[1] + r[2] * direction[2];
        maximum_error = std::max(maximum_error, std::abs(local_origin + 2.5 * local_direction - reference[a]));
      }
    }
    expect_true(maximum_error < 2e-5);
    scene.triangles.resize(2);
    motion.Seal(scene);
    expect_true(motion.owners[0][1] == 0);
    expect_true(motion.Add(0, different_time, Transform()) == 2);
    expect_true(motion.Add(1, different_time, Transform()) == 3);
    expect_true(motion.Add(1, different_time, Translate(vec3f(1, 0, 0))) == 4);
  }
  test_that("per-ray motion intersections match native nested transforms") {
    if (!MetalBVHAvailable())
      return;
    Transform start, end = Translate(vec3f(1, -.4, .7)) * RotateY(137) * Scale(.7, 1.3, .9);
    Float shear[4][4] = {{1, .2, 0, 0}, {0, 1, .15, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
    Transform child_end = RotateX(-71) * Transform(shear) * Scale(1.1, .8, 1.4);
    AnimatedTransform parent(&start, 0, &end, 1), child(&start, .1, &child_end, .9);
    double position_error = 0, normal_error = 0, distance_error = 0;
    size_t misses = 0;
    for (int kind = 0; kind < 3; ++kind) {
      WavefrontScene scene;
      scene.motion = std::make_shared<WavefrontMotion>();
      auto &motion = *scene.motion;
      motion.Add(0, parent, Transform());
      motion.Add(1, child, Translate(vec3f(.3, -.2, .1)));
      WFTriangle triangle{};
      triangle.p[0] = {-2, -2, 0, 0};
      triangle.p[1] = {2, -2, 0, 0};
      triangle.p[2] = {0, 2, 0, 0};
      for (auto &n : triangle.n)
        n = {0, 0, 1, 0};
      triangle.geometric = {0, 0, 1, 0};
      triangle.uv01 = {0, 0, 1, 0};
      triangle.uv2 = {.5, 1, 0, 0};
      auto &surface = scene.surfaces[0];
      for (int k = 0; k < 3; ++k) {
        WFVector row{};
        (&row.x)[k] = 1;
        surface.forward[k] = surface.inverse[k] = surface.bump_forward[k] = surface.bump_inverse[k] = row;
      }
      if (kind == 1) {
        WFQuadric q;
        for (int k = 0; k < 3; ++k)
          q.forward[k] = q.inverse[k] = surface.forward[k];
        q.measure = {1, 1, 1, 1};
        scene.quadrics.push_back(q);
        scene.quadric_bounds = {{-1, -1, -1, 0}, {1, 1, 1, 0}};
        surface.axes = {1, 1, 1, 1};
        motion.owners[1] = {2};
      } else {
        scene.triangles.push_back(triangle);
        if (kind == 2) {
          scene.meshes.push_back({0, 1});
          WFMeshInstance instance;
          for (int k = 0; k < 3; ++k)
            instance.forward[k] = instance.inverse[k] = surface.forward[k];
          instance.triangle_count = 1;
          scene.instances.push_back(instance);
          motion.owners[3] = {2};
        } else
          motion.owners[0] = {2};
      }
      motion.owners[2] = {2};
      motion.Seal(scene);
      std::vector<WFGeometryProbe> probes;
      std::vector<point3f> points;
      std::vector<normal3f> normals;
      for (float time : {-.2f, 0.f, .001f, .1f, .27f, .6f, .9f, 1.f, 1.2f}) {
        Transform pose = motion.At(2, time);
        point3f target = pose(point3f(0, 0, kind == 1 ? 1 : 0));
        vec3f d = unit_vector(pose(vec3f(0, 0, -1)));
        point3f origin = target - Float(2.5) * d;
        WFGeometryProbe probe;
        probe.origin = {float(origin[0]), float(origin[1]), float(origin[2]), 0};
        probe.direction = {float(d[0]), float(d[1]), float(d[2]), 10};
        std::memcpy(&probe.padding[1], &time, sizeof(time));
        probes.push_back(probe);
        points.push_back(target);
        normals.push_back(unit_vector(pose(normal3f(0, 0, 1))));
      }
      auto results = ProbeMetalGeometry(scene, probes);
      for (size_t i = 0; i < results.size(); ++i) {
        const auto &r = results[i];
        if (!r.normal.w) {
          ++misses;
          continue;
        }
        distance_error = std::max(distance_error, std::abs(double(r.point.w) - 2.5));
        for (int k = 0; k < 3; ++k) {
          position_error = std::max(position_error, std::abs(double((&r.point.x)[k]) - points[i][k]));
          normal_error = std::max(normal_error, std::abs(double((&r.normal.x)[k]) - normals[i][k]));
        }
      }
    }
    Rcpp::Rcout << "Per-ray nested motion: misses=" << misses << ", position=" << position_error
                << ", distance=" << distance_error << ", normal=" << normal_error << "\n";
    expect_true(misses == 0);
    expect_true(position_error < 5e-5);
    expect_true(distance_error < 5e-5);
    expect_true(normal_error < 5e-5);
  }
  test_that("camera projections and motion match native rays") {
    if (!MetalBVHAvailable())
      return;
    point3f origin(2, 3, 6), target(0, 1, 0);
    vec3f up(0, 1, 0);
    camera perspective(origin, target, up, 45, 1.3, .12, 6, 0, 1, 1);
    ortho_camera ortho(origin, target, up, 4, 3, 0, 1, 1);
    environment_camera panorama(origin, target, up, 0, 1, 1);
    for (RayCamera *cam : {static_cast<RayCamera *>(&perspective), static_cast<RayCamera *>(&ortho),
                           static_cast<RayCamera *>(&panorama)}) {
      compare_camera(*cam);
      cam->set_camera_flip_x(true);
      compare_camera(*cam);
      cam->set_camera_motion_blur_range(origin, target, up, 6, point3f(-3, 4, 2), point3f(1, 0, 0), vec3f(.2, 1, .1),
                                        4);
      cam->set_camera_motion_blur(true);
      cam->set_shutter_speed(1);
      compare_camera(*cam);
      cam->set_shutter_speed(INFINITY);
      compare_camera(*cam);
    }
  }
  test_that("real lens rejection and radiometric weights match native lens tracing") {
    if (!MetalBVHAvailable())
      return;
    Transform placement = Inverse(LookAt(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0)));
    AnimatedTransform animation(&placement, 0, &placement, 1);
    std::vector<Float> lens = {29.475f,  3.76f,  1.67f,  25.2f,  84.83f,  .12f,     1.f,      25.2f,  19.275f,
                               4.025f,   1.67f,  23.f,   40.77f, 3.275f,  1.699f,   23.f,     12.75f, 5.705f,
                               1.f,      18.f,   0.f,    4.5f,   0.f,     17.1f,    -14.495f, 1.18f,  1.603f,
                               17.f,     40.77f, 6.065f, 1.658f, 20.f,    -20.385f, .19f,     1.f,    20.f,
                               437.065f, 3.22f,  1.717f, 20.f,   -39.73f, 0.f,      1.f,      20.f};
    for (bool simple : {false, true}) {
      RealisticCamera cam(animation, 2, 6, 10, 10, 10, 10, simple, lens, .022, 1, 1, vec3f(0, 1, 0), placement,
                          point3f(0));
      compare_camera(cam);
      cam.set_camera_flip_x(true);
      compare_camera(cam);
      cam.set_camera_motion_blur_range(point3f(0, 0, -10), point3f(0), vec3f(0, 1, 0), 10, point3f(3, 1, -8),
                                       point3f(0), vec3f(0, 1, 0), 8);
      cam.set_camera_motion_blur(true);
      compare_camera(cam);
      cam.set_camera_motion_blur(false);
      cam.scene_camera_time = std::make_shared<SceneTimeContext>();
      cam.SetShutterInterval(0,0);
      compare_camera(cam);
      cam.SetShutterInterval(-.1,.2);
      compare_camera(cam);
      animation::CameraKey a{100,{0,0,-10, 0,0,0, 10,-1,10, 4,3, 0,1,0}};
      auto b=a,c=a;
      b.time=100.2; b.values[0]=3; b.values[4]=1; b.values[11]=.2;
      c.time=100.4;
      cam.scene_camera_track=std::make_shared<const animation::CameraTrack>(
        std::vector<animation::CameraKey>{a,b,c});
      auto clock=std::make_shared<SceneTimeContext>();
      clock->epoch=100; clock->ray_anchor=5;
      cam.scene_camera_time=clock;
      cam.SetShutterInterval(5,5.4);
      for (bool flip : {false,true}) {
        cam.set_camera_flip_x(flip);
        compare_camera(cam);
      }
      // Seek order and a large absolute epoch must not alter the lens model.
      for (double time : {100.2,100.4,100.,100.1,100.2}) {
        clock->epoch=time;
        cam.SetShutterInterval(5,5);
        compare_camera(cam);
      }
      b.values[6]=4; b.values[8]=6;
      auto variable=std::make_shared<const animation::CameraTrack>(
        std::vector<animation::CameraKey>{a,b,c});
      expect_error(variable->MakeCamera(b,1,1,5,5.4,false,&cam));
      for (double time : {100.2,100.,100.2}) {
        clock->epoch=time;
        auto snapshot=variable->MakeCamera(variable->At(time),1,1,5,5,false,&cam);
        snapshot->scene_camera_time=clock;
        snapshot->scene_camera_track=variable;
        compare_camera(*snapshot);
      }
    }
  }
  test_that("cubic cutoff keeps all three crossings and constant equality") {
    if (!MetalBVHAvailable())
      return;
    WFHazeProbe cubic = {{{.92f, 1.14f, 1.14f, .86f}, {1.14f, .86f, .86f, 1.08f}}, {}, {1, 1, 1, 1}};
    WFHazeProbe equal = {{{1, 1, 1, 1}, {1, 1, 1, 1}}, {}, {1, 1, 1, 1}};
    auto result = ProbeMetalHaze({cubic, equal});
    expect_true(result[0].direction.y == 5);
    expect_true(std::abs(result[0].origin.y - .2f) < 2e-5);
    expect_true(std::abs(result[0].origin.z - .5f) < 2e-5);
    expect_true(std::abs(result[0].origin.w - .8f) < 2e-5);
    expect_true((result[1].direction.y == 2 && result[1].origin.y == 1));
    expect_true(result[0].direction.z == 5); // Below the cutoff in intervals 0 and 2.
    expect_true(result[1].direction.z == 0); // Equality excludes haze.
  }
}
#endif
