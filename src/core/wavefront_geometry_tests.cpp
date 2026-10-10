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
