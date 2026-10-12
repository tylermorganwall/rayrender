#if defined(NOT_CRAN) && defined(RAY_HAS_METAL_BVH)
#include "../hitables/rectangle.h"
#include "../hitables/curve.h"
#include "../hitables/sphere.h"
#include "../hitables/triangle.h"
#include "../hitables/instance.h"
#include "camera.h"
#include "hlbvh.h"
#include "wavefront.h"
#include <testthat.h>

context("Metal analytic intersections") {
  test_that("mesh consistent normals preserve the native view-dependent frame") {
    if (!MetalBVHAvailable()) return;
    float vertices[] = {-2,-2,0, 2,-2,0, 0,2,0};
    float normals[] = {-.4f,0,.916515f, .6f,0,.8f, 0,.3f,.953939f};
    int indices[] = {0,1,2};
    auto material = std::make_shared<glossy>(std::make_shared<constant_texture>(point3f(.7)),
        new TrowbridgeReitzDistribution(.1,.1,nullptr,false), point3f(.06), point3f(1));
    camera cam({0,0,4},{0,0,0},{0,1,0},30,1,0,4,0,1,1);
    for (bool consistent : {false, true}) for (bool reversed : {false, true})
      for (bool instanced : {false, true}) {
      Transform transform=RotateY(23)*Scale(1.2,.7,1.5), inverse=Inverse(transform);
      TriangleMesh mesh(vertices,indices,normals,nullptr,3,3,nullptr,nullptr,material,
                        &transform,&inverse,reversed,consistent);
      auto face=std::make_shared<triangle>(&mesh,mesh.vertexIndices.data(),mesh.normalIndices.data(),
                         mesh.texIndices.data(),0,&transform,&inverse,reversed);
      hitable_list world,lights,prototype; prototype.add(face);
      Transform placement = instanced ? Translate(vec3f(.2,-.1,.3))*RotateX(17)*Scale(.8,1.4,.9)
                                     : Transform();
      std::shared_ptr<hitable> target = instanced
          ? std::static_pointer_cast<hitable>(std::make_shared<instance>(&prototype,&placement,&lights))
          : face;
      world.add(target);
      WavefrontScene scene; WavefrontReport report;
      auto session=PrepareWavefront(world,lights,cam,1,0,scene,report);
      expect_true(bool(session));
      if (!session) continue;
      std::vector<WFGeometryProbe> probes;
      std::vector<normal3f> expected;
      random_gen rng(19);
      for (float side : {-1.f,1.f}) for (float x : {-.5f,0.f,.5f})
        for (float slope : {-.8f,0.f,.8f}) {
          vec3f direction=unit_vector(placement(transform(vec3f(slope,0,-side))));
          point3f origin=placement(transform(point3f(x,0,0)))-direction*3;
          hit_record hit;
          expect_true(target->hit(Ray(origin,direction),0,10,hit,rng));
          WFGeometryProbe probe;
          probe.origin={float(origin[0]),float(origin[1]),float(origin[2]),0};
          probe.direction={float(direction[0]),float(direction[1]),float(direction[2]),10};
          probe.padding[0]=1;
          probes.push_back(probe); expected.push_back(hit.normal);
        }
      auto result=ProbeMetalGeometry(scene,probes);
      double error=0;
      for (size_t i=0;i<result.size();++i) {
        if (!result[i].normal.w) Rcpp::Rcout << "Normal probe miss at " << i << "\n";
        expect_true(result[i].normal.w!=0);
        for (int c=0;c<3;++c)
          error=std::max(error,std::abs(double((&result[i].normal.x)[c])-expected[i][c]));
      }
      Rcpp::Rcout << "Consistent normal probe: enabled=" << consistent << " reversed=" << reversed
                  << " instanced=" << instanced
                  << " max error=" << error << "\n";
      expect_true(error<3e-5);
    }
  }
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
context("Metal native curves") {
  test_that("Bezier hits preserve CPU distance, UV, normal, taper and transformations") {
    if (!MetalBVHAvailable())
      return;
    auto material = std::make_shared<hair>(point3f(.4), 1.55, .3, .4, 2);
    camera cam({0, 0, 4}, {0, 0, 0}, {0, 1, 0}, 30, 1, 0, 4, 0, 1, 1);
    random_gen rng(714);
    size_t mismatches = 0, hits = 0;
    double distance_error = 0, normal_error = 0, uv_error = 0;
    for (auto type : {CurveType::Flat, CurveType::Cylinder, CurveType::Ribbon}) {
      for (bool flipped : {false, true}) {
        for (Transform transform : {Transform(), Translate(vec3f(.3, -.2, .4)) * RotateY(31) * Scale(-1.3, .7, .9)}) {
          Transform inverse = Inverse(transform);
          for (bool curved : {false, true}) {
            point3f cp[4] = {point3f(-.4, -1, 0), point3f(.7, -.4, .3), point3f(-.5, .5, -.2), point3f(.5, 1, 0)};
            if (!curved)
              for (int k = 0; k < 4; ++k)
                cp[k] = point3f(0, -1 + 2 * k / 3.f, 0);
            vec3f normals[2] = {vec3f(.1, 0, 1), vec3f(-.1, 0, 1)};
            auto common = std::make_shared<CurveCommon>(cp, .12, .035, type, normals);
            hitable_list world, lights;
            for (int k = 0; k < 8; ++k)
              world.add(
                  std::make_shared<curve>(k / 8.f, (k + 1) / 8.f, common, material, &transform, &inverse, flipped));
            WavefrontScene scene;
            WavefrontReport report;
            auto session = PrepareWavefront(world, lights, cam, 1, 0, scene, report);
            expect_true(bool(session));
            if (!session)
              continue;
            expect_true((scene.curves.size() == 8 && scene.triangles.empty()));
            std::vector<WFGeometryProbe> probes;
            std::vector<hit_record> expected;
            std::vector<bool> found;
            for (int i = 0; i < 384; ++i) {
              // Include tilted rays, subcurve endpoints, width edges and clipping.
              Float y = -.99 + 1.98 * rng.unif_rand();
              Float u = (y + 1) * .5;
              Float center = curved ? std::pow(1 - u, 3) * cp[0][0] + 3 * (1 - u) * (1 - u) * u * cp[1][0] +
                                          3 * (1 - u) * u * u * cp[2][0] + u * u * u * cp[3][0]
                                    : 0;
              point3f o = transform(point3f(center + .16 * (rng.unif_rand() - .5), y, 3));
              vec3f d = transform(vec3f(0, 0, -1));
              WFGeometryProbe p;
              p.origin = {float(o[0]), float(o[1]), float(o[2]), i % 13 == 0 ? 3.05f : 0.f};
              p.direction = {float(d[0]), float(d[1]), float(d[2]), i % 17 == 0 ? 2.8f : 5.f};
              probes.push_back(p);
              Ray ray(point3f(p.origin.x, p.origin.y, p.origin.z), vec3f(p.direction.x, p.direction.y, p.direction.z));
              hit_record hit;
              found.push_back(world.hit(ray, p.origin.w, p.direction.w, hit, rng));
              expected.push_back(hit);
            }
            auto actual = ProbeMetalGeometry(scene, probes);
            expect_true(actual.size() == probes.size());
            for (size_t i = 0; i < actual.size(); ++i) {
              bool gpu = actual[i].normal.w != 0;
              mismatches += gpu != found[i];
              if (!(gpu && found[i]))
                continue;
              ++hits;
              distance_error = std::max(distance_error, double(std::abs(actual[i].point.w - expected[i].t)));
              vec3f n(actual[i].normal.x, actual[i].normal.y, actual[i].normal.z);
              normal_error =
                  std::max(normal_error, double((n - unit_vector(convert_to_vec3(expected[i].normal))).length()));
              uv_error = std::max(uv_error, double(std::max(std::abs(actual[i].uv.x - expected[i].u),
                                                            std::abs(actual[i].uv.y - expected[i].v))));
            }
          }
        }
      }
    }
    Rcpp::Rcout << "Curve probe hits " << hits << "; classifications " << mismatches << "; distance " << distance_error
                << "; normal " << normal_error << "; UV " << uv_error << "\n";
    expect_true(hits > 1000);
    expect_true(mismatches == 0);
    expect_true(distance_error < 2e-4);
    expect_true(normal_error < .005);
    expect_true(uv_error < .002);
  }
  test_that("direct curves keep the closest folded segment and reject degenerate curves") {
    if (!MetalBVHAvailable())
      return;
    Transform identity;
    camera cam({0, 0, 0}, {0, 0, 1}, {0, 1, 0}, 30, 1, 0, 4, 0, 1, 1);
    auto material = std::make_shared<hair>(point3f(.4), 1.55, .3, .4, 2);
    point3f cp[4] = {point3f(-1, 0, 4), point3f(2, 0, 4), point3f(-2, 0, 1), point3f(1, 0, 1)};
    auto common = std::make_shared<CurveCommon>(cp, .025, .025, CurveType::Cylinder, nullptr);
    hitable_list world, lights;
    world.add(std::make_shared<curve>(0, 1, common, material, &identity, &identity, false));
    point3f zero[4] = {point3f(0), point3f(0), point3f(0), point3f(0)};
    for (float width : {0.f, .1f}) {
      auto degenerate = std::make_shared<CurveCommon>(zero, width, width, CurveType::Flat, nullptr);
      world.add(std::make_shared<curve>(0, 1, degenerate, material, &identity, &identity, false));
    }
    WavefrontScene scene;
    WavefrontReport report;
    auto session = PrepareWavefront(world, lights, cam, 1, 0, scene, report);
    expect_true(bool(session));
    if (!session)
      return;
    std::vector<WFGeometryProbe> probes;
    for (float minimum : {0.f, 2.f, 3.f, 4.1f}) {
      WFGeometryProbe p;
      p.origin = {0, 0, 0, minimum};
      p.direction = {0, 0, 1, 10};
      probes.push_back(p);
    }
    auto actual = ProbeMetalGeometry(scene, probes);
    random_gen rng(3);
    for (size_t i = 0; i < probes.size(); ++i) {
      hit_record hit;
      bool found = world.hit(Ray(point3f(0), vec3f(0, 0, 1)), probes[i].origin.w, 10, hit, rng);
      expect_true((actual[i].normal.w != 0) == found);
      if (found)
        expect_true(std::abs(actual[i].point.w - hit.t) < 2e-5);
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
