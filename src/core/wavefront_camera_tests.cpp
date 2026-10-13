#if defined(NOT_CRAN) && defined(RAY_HAS_METAL_BVH)
#include "camera.h"
#include "hlbvh.h"
#include "wavefront.h"
#include "wavefront_motion.h"
#include <testthat.h>

namespace {
void compare_camera(RayCamera &camera) {
  WFCamera exported;
  std::vector<WFVector> data;
  ExportWavefrontCamera(camera, exported, data);
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
