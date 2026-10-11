#ifdef NOT_CRAN
#include "../hitables/csg.h"
#include "../volumes/boundary.h"
#include "camera.h"
#include "hlbvh.h"
#include "wavefront.h"
#include <testthat.h>

context("CSG boundary and sampling") {
  test_that("clipped and thin-shell rays find later crossings") {
    Transform identity;
    auto material = std::make_shared<diffuse_material>(std::make_shared<constant_texture>(point3f(.5)));
    auto sphere = std::make_shared<csg_sphere>(vec3f(0), 1);
    csg solid(material, sphere, &identity, &identity, false);
    random_gen rng(1);
    hit_record hit;
    expect_true(solid.hit(Ray(point3f(0, 0, 3), vec3f(0, 0, -1)), Float(2.5), 10, hit, rng));
    expect_true(std::abs(hit.t - 4) < 1e-5);
    csg shell(material, std::make_shared<csg_onion>(sphere, .0005), &identity, &identity, false);
    expect_true(shell.hit(Ray(point3f(0, 0, 1), vec3f(0, 0, 1)), 0, 1, hit, rng));
    expect_true(std::abs(hit.t - .0005) < 1e-5);
    expect_true(ValidateMediumBoundary(&shell, true));
    csg plane(material, std::make_shared<csg_plane>(vec3f(0, 1, 0), point3f(0), 1, 1), &identity, &identity, false);
    expect_true(!ValidateMediumBoundary(&plane, false));
    double pdf_sum = 0;
    for (int i = 0; i < 100; ++i) {
      vec3f direction = solid.random(point3f(0, 0, 4), rng, 0);
      pdf_sum += solid.pdf_value(point3f(0, 0, 4), direction, rng, 0);
    }
    expect_true(std::isfinite(pdf_sum));
    expect_true(pdf_sum > 100 / (4 * M_PI));
    // Diffusion walks request successive hits on an unchanged line. At grazing
    // incidence the convergence band spans more than a fixed origin nudge.
    size_t wrong_crossings = 0;
    for (int i = 0; i < 1000; ++i) {
      Ray ray(point3f(Float(i) / 1000, 0, 3), vec3f(0, 0, -1));
      Float lower = 0;
      unsigned crossings = 0;
      while (crossings < 10 && solid.hit(ray, lower, 10, hit, rng)) {
        ++crossings;
        lower = std::nextafter(hit.t, Float(INFINITY));
      }
      wrong_crossings += crossings != 2;
    }
    expect_true(wrong_crossings == 0);
  }
}
#if defined(RAY_HAS_METAL_BVH)
context("Metal CSG intersections") {
  test_that("all distance operations preserve CPU hits and normals") {
    if (!MetalBVHAvailable())
      return;
    auto material = std::make_shared<diffuse_material>(std::make_shared<constant_texture>(point3f(.5)));
    auto ball = std::make_shared<csg_sphere>(vec3f(0), .8);
    auto box = std::make_shared<csg_box>(point3f(.2, 0, 0), vec3f(1.1));
    std::vector<std::shared_ptr<ImplicitShape>> shapes = {
        ball,
        box,
        std::make_shared<csg_rounded_box>(point3f(0), vec3f(1), .15),
        std::make_shared<csg_torus>(point3f(0), .7, .2),
        std::make_shared<csg_capsule>(point3f(0, -.5, 0), point3f(.2, .5, 0), .3),
        std::make_shared<csg_cylinder>(point3f(0, -.6, 0), point3f(.2, .6, 0), .5, .1),
        std::make_shared<csg_ellipsoid>(point3f(0), vec3f(.8, .6, .4)),
        std::make_shared<csg_rounded_cone>(point3f(0, -.5, 0), point3f(0, .6, 0), .5, .2),
        std::make_shared<csg_cone>(point3f(0, -.5, 0), point3f(0, .6, 0), .6),
        std::make_shared<csg_pyramid>(point3f(0, -.4, 0), .8, 1),
        std::make_shared<csg_triangle>(point3f(-1, -.5, 0), point3f(1, -.5, 0), point3f(0, .8, 0)),
        std::make_shared<csg_plane>(vec3f(0, 0, 1), point3f(0), 1, 1),
        std::make_shared<Union>(ball, box),
        std::make_shared<Subtract>(ball, box),
        std::make_shared<Intersect>(ball, box),
        std::make_shared<Blend>(ball, box, .3f),
        std::make_shared<BlendMinus>(ball, box, .3f),
        std::make_shared<Mix>(ball, box, .4f),
        std::make_shared<csg_round>(box, .15),
        std::make_shared<csg_onion>(ball, .015),
        std::make_shared<csg_scale>(ball, .75),
        std::make_shared<csg_rotate>(box, point3f(.1), vec3f(.2, 1, .3)),
        std::make_shared<csg_translate>(ball, vec3f(.2, .1, 0)),
        std::make_shared<csg_elongate>(ball, point3f(0), vec3f(.2, 0, 0)),
        std::make_shared<csg_elongate_robust>(ball, point3f(0), vec3f(.2, 0, 0)),
        std::make_shared<csg_list>(std::vector<std::shared_ptr<ImplicitShape>>{ball, box})};
    camera cam({0, 0, 4}, {0, 0, 0}, {0, 1, 0}, 30, 1, 0, 4, 0, 1, 1);
    for (size_t k = 0; k < shapes.size(); ++k) {
      Transform transform = Translate(vec3f(.1, -.1, .2)) * RotateY(17) * Scale(-1.1, .8, 1.2);
      Transform inverse = Inverse(transform);
      auto shape = std::make_shared<csg>(material, shapes[k], &transform, &inverse, false);
      hitable_list world, lights;
      world.add(shape);
      WavefrontScene scene;
      WavefrontReport report;
      auto session = PrepareWavefront(world, lights, cam, 1, 0, scene, report);
      expect_true(bool(session));
      if (!session)
        continue;
      std::vector<WFGeometryProbe> probes;
      std::vector<hit_record> expected(169);
      std::vector<bool> hits;
      random_gen rng(1);
      for (int i = 0; i < 169; ++i) {
        point3f o = transform(point3f((i % 13 - 6) * .15, (i / 13 - 6) * .15, i % 3 == 0 ? 0 : 3));
        vec3f d = unit_vector(transform(vec3f(.03, .01, -1)));
        WFGeometryProbe p;
        p.origin = {float(o[0]), float(o[1]), float(o[2]), 0};
        p.direction = {float(d[0]), float(d[1]), float(d[2]), 10};
        probes.push_back(p);
        hits.push_back(shape->hit(Ray(o, d), 0, 10, expected[i], rng));
      }
      auto actual = ProbeMetalGeometry(scene, probes);
      size_t errors = 0;
      double distance_error = 0, normal_error = 0;
      for (size_t i = 0; i < probes.size(); ++i) {
        bool hit = actual[i].primitive == 0x80000000u;
        errors += hit != hits[i];
        if (hit && hits[i]) {
          distance_error = std::max(distance_error, std::abs(double(actual[i].point.w) - expected[i].t));
          double dotn = actual[i].normal.x * expected[i].normal[0] + actual[i].normal.y * expected[i].normal[1] +
                        actual[i].normal.z * expected[i].normal[2];
          // An unsigned triangle field has two opposite limiting normals at
          // the plane. A ray starting exactly there has no unique oriented gradient.
          if (k == 10)
            dotn = std::abs(dotn);
          normal_error = std::max(normal_error, 1 - dotn);
        }
      }
      Rcpp::Rcout << "CSG " << k << ": hit errors=" << errors << ", distance=" << distance_error
                  << ", normal=" << normal_error << "\n";
      expect_true(errors == 0);
      expect_true(distance_error < .002);
      expect_true(normal_error < .002);
    }
  }
}
#endif
#endif
