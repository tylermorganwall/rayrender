#if defined(NOT_CRAN) && defined(RAY_HAS_METAL_BVH)
#include "../hitables/rectangle.h"
#include "../hitables/sphere.h"
#include "camera.h"
#include "hlbvh.h"
#include "wavefront.h"
#include <testthat.h>

context("Metal analytic intersections") {
  test_that("mixed traversal preserves both roots, transforms, normals, and clipped segments") {
    if (!MetalBVHAvailable())
      return;
    auto material =
        std::make_shared<diffuse_material>(std::make_shared<constant_texture>(point3f(.5)));
    Transform identity;
    camera cam({0, 0, 4}, {0, 0, 0}, {0, 1, 0}, 30, 1, 0, 4, 0, 1, 1);
    for (Transform transform :
         {Transform(), Translate(vec3f(.3, -.2, .4)) * RotateY(31) * Scale(-1.3, .7, .9),
          Translate(vec3f(0, -1001, 0)) * Scale(1000, 1000, 1000)}) {
      Transform inverse = Inverse(transform);
      hitable_list world, lights;
      world.add(
          std::make_shared<sphere>(1, material, nullptr, nullptr, &transform, &inverse, false));
      world.add(std::make_shared<xy_rect>(-3, 3, -3, 3, -3000, material, nullptr, nullptr,
                                          &identity, &identity, false));
      WavefrontScene scene;
      WavefrontReport report;
      auto session = PrepareWavefront(world, lights, cam, 1, 0, scene, report);
      expect_true(bool(session));
      if (!session)
        continue;
      expect_true(scene.quadrics.size() == 1);
      expect_true(scene.triangles.size() == 2);
      const auto &q = scene.quadrics[0];
      std::vector<WFGeometryProbe> probes;
      std::vector<double> distances;
      std::vector<vec3<double>> normals;
      // Transform rays that hit near/far roots, start inside, narrowly miss,
      // graze the silhouette, and clip away the near root with t_min.
      for (int i = 0; i < 183; ++i) {
        double x = (i % 31 - 15) / 14.5, y = (i % 7 - 3) / 13.;
        point3f origin = transform(point3f(x, y, i % 3 == 0 ? 0 : 3));
        vec3f direction = unit_vector(transform(vec3f(0, 0, -1)));
        // Near-pole hits on the enormous ground sphere expose cancellation
        // between its center and radius, which a mid-latitude probe misses.
        if (i >= 151) {
          origin = transform(point3f((i - 167) * .0001, 1.002, .0007));
          direction = unit_vector(transform(vec3f(.003, -1, -.001)));
        }
        float minimum = i % 5 == 0 ? float(transform(vec3f(0, 0, 2.5)).length()) : 1e-6f;
        float maximum = i % 11 == 0 ? minimum + .1f : 10000.f;
        WFGeometryProbe probe;
        probe.origin = {float(origin[0]), float(origin[1]), float(origin[2]), minimum};
        probe.direction = {float(direction[0]), float(direction[1]), float(direction[2]), maximum};
        probes.push_back(probe);
        vec3<double> o, d;
        for (int k = 0; k < 3; ++k) {
          const auto &r = q.inverse[k];
          o[k] = double(r.x) * probe.origin.x + double(r.y) * probe.origin.y +
                 double(r.z) * probe.origin.z + r.w;
          d[k] = double(r.x) * probe.direction.x + double(r.y) * probe.direction.y +
                 double(r.z) * probe.direction.z;
        }
        double a = dot(d, d), b = dot(o, d), c = dot(o, o) - 1, disc = b * b - a * c;
        double expected = INFINITY;
        if (disc >= 0) {
          double near = (-b - std::sqrt(disc)) / a, far = (-b + std::sqrt(disc)) / a;
          double candidate = near >= minimum ? near : far;
          if (candidate >= minimum && candidate <= maximum)
            expected = candidate;
        }
        distances.push_back(expected);
        vec3<double> n(0);
        if (std::isfinite(expected)) {
          auto local = unit_vector(o + expected * d);
          for (int k = 0; k < 3; ++k)
            n[k] = double((&q.inverse[0].x)[k]) * local[0] +
                   double((&q.inverse[1].x)[k]) * local[1] +
                   double((&q.inverse[2].x)[k]) * local[2];
          n = unit_vector(n) * q.measure.y;
        }
        normals.push_back(n);
      }
      auto actual = ProbeMetalGeometry(scene, probes);
      size_t classification_errors = 0;
      double distance_error = 0, normal_error = 0;
      for (size_t i = 0; i < probes.size(); ++i) {
        bool hit = actual[i].primitive == 0x80000000u;
        classification_errors += hit != std::isfinite(distances[i]);
        if (hit && std::isfinite(distances[i])) {
          distance_error = std::max(distance_error, std::abs(actual[i].point.w - distances[i]) /
                                                        std::max(1., std::abs(distances[i])));
          for (int k = 0; k < 3; ++k)
            normal_error =
                std::max(normal_error, std::abs(double((&actual[i].normal.x)[k]) - normals[i][k]));
        }
      }
      Rcpp::Rcout << "Analytic probe relative distance error: " << distance_error
                  << "; normal error: " << normal_error
                  << "; classification errors: " << classification_errors << "\n";
      expect_true(classification_errors == 0);
      expect_true(distance_error < 1e-4);
      expect_true(normal_error < 2e-4);
    }
  }
}
#endif

#if defined(NOT_CRAN) && defined(RAY_HAS_METAL_BVH)
context("Metal shared mesh intersections") {
  test_that("persistent sessions retain every indirectly referenced mesh BLAS") {
    if (!MetalBVHAvailable())
      return;
    // Three separate prototypes make two BLASes reachable only through the
    // TLAS. Binding the first BLAS is insufficient to keep the others resident.
    WavefrontScene scene;
    for (uint32_t copy = 0; copy < 3; ++copy) {
      WFTriangle face;
      face.p[0] = {-1, -1, 0, 0};
      face.p[1] = {1, -1, 0, 0};
      face.p[2] = {0, 1, 0, 0};
      face.geometric = {0, 0, 1, 0};
      for (auto &normal : face.n)
        normal = face.geometric;
      face.material = copy;
      scene.triangles.push_back(face);
      scene.meshes.push_back({copy, 1});
      WFMeshInstance placement;
      placement.triangle_start = placement.mesh = copy;
      placement.triangle_count = 1;
      placement.forward[0] = {1, 0, 0, 4.f * (int(copy) - 1)};
      placement.inverse[0] = {1, 0, 0, -placement.forward[0].w};
      placement.forward[1] = placement.inverse[1] = {0, 1, 0, 0};
      placement.forward[2] = placement.inverse[2] = {0, 0, 1, 0};
      scene.instances.push_back(placement);
      WFTexture texture;
      texture.a = {float(copy == 0), float(copy == 1), float(copy == 2), 0};
      scene.textures.push_back(texture);
      WFMaterial material;
      material.texture = copy;
      material.emissive = 1;
      material.diffuse.w = 1;
      scene.materials.push_back(material);
    }
    auto session = MakeMetalWavefront(scene, 1);
    auto other_queue = MakeMetalWavefront(scene, 1);
    WFParameters parameters;
    parameters.width = parameters.height = parameters.active = parameters.boundaries = 1;
    parameters.max_depth = 2;
    parameters.camera.orthographic = 1;
    parameters.camera.forward = {0, 0, -1, 0};
    WavefrontReport report;
    const WFPixel *pixels = nullptr;
    for (int edit = 0; edit < 9; ++edit) {
      int copy = edit % 3;
      parameters.camera.lower_left = {4.f * (copy - 1), 0, 3, 0};
      expect_true(other_queue->Render(parameters, {0}, {}, [] { return false; }, pixels, report));
      bool completed = session->Render(parameters, {0}, {}, [] { return false; }, pixels, report);
      expect_true(completed);
      if (completed) {
        const auto &actual = pixels[0].radiance;
        expect_true(std::abs(actual.x - float(copy == 0)) < 1e-6);
        expect_true(std::abs(actual.y - float(copy == 1)) < 1e-6);
        expect_true(std::abs(actual.z - float(copy == 2)) < 1e-6);
      }
    }
  }
  test_that("shared BLAS placements preserve hits, normals, UVs and containment identity") {
    if (!MetalBVHAvailable())
      return;
    WavefrontScene scene;
    WFTriangle face;
    face.p[0] = {-1, -1, 0, 0};
    face.p[1] = {1, -1, 0, 0};
    face.p[2] = {0, 1, 0, 0};
    for (auto &n : face.n)
      n = {0, 0, 1, 0};
    face.geometric = {0, 0, 1, 0};
    face.uv01 = {0, 0, 1, 0};
    face.uv2 = {.5, 1, 0, 0};
    scene.triangles.push_back(face);
    scene.meshes.push_back({0, 1});
    scene.materials.push_back(WFMaterial());
    scene.boundaries.resize(4);
    std::vector<WFGeometryProbe> probes;
    std::vector<point3f> points;
    std::vector<normal3f> normals;
    for (int copy = 0; copy < 3; ++copy) {
      Transform placement = Translate(vec3f(4 * (copy - 1), 0, 0)) * RotateY(20 * copy) *
                            Scale(copy == 1 ? -1.2 : .8, .7, 1.3);
      WFMeshInstance instance;
      instance.triangle_count = 1;
      instance.boundary = copy + 1;
      const auto &f = placement.GetMatrix(), &v = placement.GetInverseMatrix();
      for (int k = 0; k < 3; ++k) {
        instance.forward[k] = {float(f.m[k][0]), float(f.m[k][1]), float(f.m[k][2]),
                               float(f.m[k][3])};
        instance.inverse[k] = {float(v.m[k][0]), float(v.m[k][1]), float(v.m[k][2]),
                               float(v.m[k][3])};
      }
      scene.instances.push_back(instance);
      for (int j = 0; j < 51; ++j) {
        point3f p = placement(point3f((j % 7 - 3) * .07, (j % 5 - 2) * .09, 0));
        normal3f n = unit_vector(placement(normal3f(0, 0, 1)));
        vec3f direction(-n[0], -n[1], -n[2]);
        point3f origin = p - 3 * direction;
        WFGeometryProbe ray;
        ray.origin = {float(origin[0]), float(origin[1]), float(origin[2]), 0};
        ray.direction = {float(direction[0]), float(direction[1]), float(direction[2]), 4};
        ray.body = j % 3 == 0 ? 0 : (j % 3 == 1 ? copy + 1 : (copy + 1) % 3 + 1);
        probes.push_back(ray);
        points.push_back(p);
        normals.push_back(n);
      }
    }
    auto hits = ProbeMetalGeometry(scene, probes);
    double point_error = 0, normal_error = 0, distance_error = 0;
    size_t mismatches = 0;
    for (size_t i = 0; i < hits.size(); ++i) {
      const auto &hit = hits[i];
      if (i % 51 % 3 == 2) {
        mismatches += hit.primitive != UINT32_MAX;
        continue;
      }
      mismatches += hit.primitive != 0 || hit.boundary != i / 51 + 1;
      point_error =
          std::max(point_error,
                   double((point3f(hit.point.x, hit.point.y, hit.point.z) - points[i]).length()));
      normal_error = std::max(
          normal_error,
          double((normal3f(hit.normal.x, hit.normal.y, hit.normal.z) - normals[i]).length()));
      distance_error = std::max(distance_error, std::abs(double(hit.point.w) - 3));
    }
    expect_true(mismatches == 0);
    expect_true(point_error < 2e-5);
    expect_true(normal_error < 2e-5);
    expect_true(distance_error < 2e-5);
  }
}
#endif
