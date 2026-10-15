#ifdef NOT_CRAN
#include "../hitables/csg.h"
#include "../volumes/boundary.h"
#include "camera.h"
#include "hlbvh.h"
#include "wavefront.h"
#include <testthat.h>

context("CSG material attribution") {
  test_that("curve pruning preserves the full signed-distance union") {
    std::vector<std::shared_ptr<ImplicitShape>> segments;
    for (int i = 0; i < 40; ++i) {
      Float t = Float(i) / 40, u = Float(i + 1) / 40;
      segments.push_back(std::make_shared<csg_rounded_cone>(
        point3f(std::cos(6*t),t,std::sin(6*t)),
        point3f(std::cos(6*u),u,std::sin(6*u)), .2-.1*t, .2-.1*u));
    }
    csg_curve_tube accelerated(segments);
    csg_list reference(segments);
    random_gen rng(121);
    double error = 0;
    for (int i = 0; i < 12000; ++i) {
      point3f p(3*rng.unif_rand()-1.5,2*rng.unif_rand()-.5,3*rng.unif_rand()-1.5);
      error = std::max(error,std::abs(double(accelerated.getDistance(p)-reference.getDistance(p))));
    }
    expect_true(error < 1e-6);
    std::vector<WFCsgInstruction> code;
    ExportWavefrontCSG(accelerated,code);
    expect_true(code[0].operation == 23);
    expect_true(code[0].padding[0] == code.size());
  }
  test_that("nested smooth materials follow analytical distance weights") {
    auto base = std::make_shared<csg_sphere>(vec3f(0), 1);
    auto red = std::make_shared<csg_material_label>(base, 1);
    auto blue = std::make_shared<csg_material_label>(base, 2);
    auto green = std::make_shared<csg_material_label>(base, 3);
    auto pair = std::make_shared<Mix>(red, blue, .25f);
    Mix nested(pair, green, .4f);
    uint32_t seed = 19;
    unsigned counts[4] = {};
    double error = 0;
    for (unsigned i = 0; i < 40000; ++i) {
      auto sample = SampleCsgSurface(nested, point3f(0,0,1), seed);
      if (sample.material < 4) ++counts[sample.material];
      error = std::max(error, std::abs(double(sample.distance)));
    }
    expect_true(error == 0);
    expect_true(counts[0] == 0);
    expect_true(std::abs(counts[1]/40000. - .45) < .015);
    expect_true(std::abs(counts[2]/40000. - .15) < .015);
    expect_true(std::abs(counts[3]/40000. - .4) < .015);
    Blend smooth(red, blue, .4f);
    smooth.blend_materials = false;
    expect_true(SampleCsgSurface(smooth,point3f(0,0,1),seed).material == 1);
    Blend zero(red,blue,0.f);
    expect_true(std::isfinite(SampleCsgSurface(zero,point3f(0,0,1),seed).distance));
    BlendMinus subtract(red,blue,.4f);
    expect_true(std::isfinite(SampleCsgSurface(subtract,point3f(0,0,1),seed).distance));
    std::vector<WFCsgInstruction> code;
    ExportWavefrontCSG(smooth,code);
    expect_true(code.size() == 5);
    expect_true((code[1].operation == 22 && code[1].padding[0] == 1));
    expect_true((code[3].operation == 22 && code[3].padding[0] == 2));
    expect_true((code[4].operation == 33 && code[4].padding[1] == 1));
  }
  test_that("material queries agree with geometry through coordinate modifiers") {
    auto sphere = std::make_shared<csg_sphere>(vec3f(0), 1);
    auto colored = std::make_shared<csg_material_label>(sphere, 7);
    std::vector<std::shared_ptr<ImplicitShape>> shapes{
      std::make_shared<csg_scale>(colored, 2),
      std::make_shared<csg_translate>(colored, vec3f(1,2,3)),
      std::make_shared<csg_rotate>(colored, point3f(1,0,0), vec3f(1,1,0)),
      std::make_shared<csg_round>(colored, .2),
      std::make_shared<csg_onion>(colored, .2),
      std::make_shared<csg_elongate>(colored, point3f(0), vec3f(.3)),
      std::make_shared<csg_elongate_robust>(colored, point3f(0), vec3f(.3))
    };
    uint32_t seed = 5;
    double error = 0;
    unsigned wrong = 0;
    for (const auto &shape : shapes) {
      for (int i = -20; i < 20; ++i) {
        point3f p(Float(i)/7, Float(i)/9, Float(i)/11);
        auto sample = SampleCsgSurface(*shape,p,seed);
        error = std::max(error,std::abs(double(sample.distance-shape->getDistance(p))));
        wrong += sample.material != 7;
      }
    }
    expect_true(error < 1e-6);
    expect_true(wrong == 0);
  }
}

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

    // A small object viewed from tens of kilometres away still has resolvable
    // local features. Every target lies strictly inside this capsule, so every
    // ray must hit regardless of the sphere-tracing implementation.
    csg distant(material, std::make_shared<csg_capsule>(point3f(-1.1, -.3, -5.6), point3f(-1.1, .4, -5.6), .3),
                &identity, &identity, false);
    size_t missed = 0;
    for (int i = 0; i < 1000; ++i) {
      Float angle = Float(i) * Float(2 * M_PI / 1000);
      point3f target(-1.1 + .24 * std::cos(angle), -.25 + .6 * (i % 17) / 16, -5.6 + .24 * std::sin(angle));
      point3f origin(35, 20, 56);
      missed += !distant.hit(Ray(origin, unit_vector(target - origin)), 0, 150, hit, rng);
    }
    expect_true(missed == 0);
  }
}
#if defined(RAY_HAS_METAL_BVH)
context("Metal CSG intersections") {
  test_that("small distant CSG solids retain local marching precision") {
    if (!MetalBVHAvailable())
      return;
    Transform identity;
    auto material = std::make_shared<diffuse_material>(std::make_shared<constant_texture>(point3f(.5)));
    auto capsule = std::make_shared<csg_capsule>(point3f(-1.1, -.3, -5.6), point3f(-1.1, .4, -5.6), .3);
    hitable_list world, lights;
    world.add(std::make_shared<csg>(material, capsule, &identity, &identity, false));
    camera cam({35, 20, 56}, {-1.1, .05, -5.6}, {0, 1, 0}, 30, 1, 0, 75, 0, 1, 1);
    WavefrontScene scene;
    WavefrontReport report;
    auto session = PrepareWavefront(world, lights, cam, 1, 0, scene, report);
    expect_true(bool(session));
    if (!session)
      return;
    std::vector<WFGeometryProbe> probes;
    for (int i = 0; i < 1000; ++i) {
      Float angle = Float(i) * Float(2 * M_PI / 1000);
      point3f target(-1.1 + .24 * std::cos(angle), -.25 + .6 * (i % 17) / 16, -5.6 + .24 * std::sin(angle));
      point3f origin(35, 20, 56);
      vec3f direction = unit_vector(target - origin);
      WFGeometryProbe probe;
      probe.origin = {35, 20, 56, 0};
      probe.direction = {direction[0], direction[1], direction[2], 150};
      probes.push_back(probe);
    }
    auto actual = ProbeMetalGeometry(scene, probes);
    size_t missed = 0;
    for (const auto &hit : actual)
      missed += hit.primitive != 0x80000000u;
    expect_true(missed == 0);
  }
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
    std::vector<std::shared_ptr<ImplicitShape>> tube_segments;
    for (int i = 0; i < 24; ++i) {
      Float a = Float(i)/24, b = Float(i+1)/24;
      tube_segments.push_back(std::make_shared<csg_rounded_cone>(
        point3f(.6*std::cos(a*4),a-.5,.6*std::sin(a*4)),
        point3f(.6*std::cos(b*4),b-.5,.6*std::sin(b*4)),.18-.1*a,.18-.1*b));
    }
    shapes.push_back(std::make_shared<csg_curve_tube>(tube_segments));
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
#ifdef NOT_CRAN
#include "../hitables/csg.h"
#include "wavefront.h"
#include <testthat.h>
context("CSG grouped curve pruning") {
  test_that("merged curve geometry preserves distance and bounds") {
    std::vector<std::shared_ptr<ImplicitShape>> curves, segments;
    for (int arm=0;arm<8;++arm) {
      std::vector<std::shared_ptr<ImplicitShape>> part;
      for (int j=0;j<24;++j) {
        Float a=j/24.f,b=(j+1)/24.f,angle=arm*.785398163f;
        point3f p(a*std::cos(angle),.2f*std::sin(a*6),a*std::sin(angle));
        point3f q(b*std::cos(angle),.2f*std::sin(b*6),b*std::sin(angle));
        part.push_back(std::make_shared<csg_rounded_cone>(p,q,.15f-.1f*a,.15f-.1f*b));
      }
      segments.insert(segments.end(),part.begin(),part.end());
      curves.push_back(std::make_shared<csg_material_label>(std::make_shared<csg_curve_tube>(part),arm+1));
    }
    csg_list original(curves);
    csg_curve_tube merged(segments);
    random_gen rng(42);
    double error=0;
    for(int i=0;i<20000;++i) {
      point3f p(4*rng.unif_rand()-2,4*rng.unif_rand()-2,4*rng.unif_rand()-2);
      error=std::max(error,std::abs(double(original.getDistance(p)-merged.getDistance(p))));
    }
    expect_true(error==0);
    aabb a,b;original.bbox(0,0,a);merged.bbox(0,0,b);
    for(int i=0;i<3;++i) {expect_true(a.min()[i]==b.min()[i]);expect_true(a.max()[i]==b.max()[i]);}
    std::vector<WFCsgInstruction> code;
    ExportWavefrontCSG(original,code);
    expect_true(code[0].operation==25);
    size_t material_start=code[0].padding[0];
    expect_true(material_start<code.size());
    expect_true(code[material_start-1].operation==26);
    expect_true(material_start-1+code[material_start-1].padding[0]==code.size());
    unsigned labels=0;for(auto &n:code) labels+=n.operation==22;
    expect_true(labels==8);

    auto nested=std::make_shared<csg_list>(curves);
    csg_list outer({nested,nested});
    code.clear();ExportWavefrontCSG(outer,code);
    unsigned dispatches=0;for(auto &n:code)dispatches+=n.operation==25;
    expect_true(dispatches==1);
    auto sphere=std::make_shared<csg_sphere>(vec3f(0),.5f);
    csg_list mixed({curves[0],sphere});
    code.clear();ExportWavefrontCSG(mixed,code);
    dispatches=0;for(auto &n:code)dispatches+=n.operation==25;
    expect_true(dispatches==0);
    Blend smooth(curves[0],curves[1],.2f);
    code.clear();ExportWavefrontCSG(smooth,code);
    dispatches=0;for(auto &n:code)dispatches+=n.operation==25;
    expect_true(dispatches==0);
    csg_list one({curves[0]});
    code.clear();ExportWavefrontCSG(one,code);
    expect_true(code[0].operation==23);
    std::shared_ptr<ImplicitShape> deep=curves[0];
    for(int i=0;i<40;++i) deep=std::make_shared<csg_list>(
      std::vector<std::shared_ptr<ImplicitShape>>{deep,curves[0]});
    code.clear();
    expect_error(ExportWavefrontCSG(*deep,code));
  }
}
#endif
