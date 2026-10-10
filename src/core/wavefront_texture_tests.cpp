#if defined(NOT_CRAN) && defined(RAY_HAS_METAL_BVH)
#include "../hitables/rectangle.h"
#include "../materials/texturegraph.h"
#include "camera.h"
#include "hlbvh.h"
#include "wavefront.h"
#include <testthat.h>

context("Metal texture evaluation") {
  test_that("graph instructions reproduce native values before shading and tonemapping") {
    if (!MetalBVHAvailable())
      return;
    auto ns = Rcpp::Environment::namespace_env("rayrender");
    Rcpp::Function coordinates = ns["texture_coordinates"], noise = ns["texture_noise"],
                   checker = ns["texture_checker"], mix = ns["texture_mix"],
                   scale = ns["texture_scale"], channel = ns["texture_channel"],
                   gradient = ns["texture_gradient"], direction = ns["texture_direction_mix"];
    for (const std::string &space : {"uv", "object", "world"}) {
      auto mapping = coordinates(space, Rcpp::Named("scale") = Rcpp::NumericVector::create(2, 3, 4),
                                 Rcpp::Named("offset") = Rcpp::NumericVector::create(-.3, .2, -.1),
                                 Rcpp::Named("rotation") = 27);
      auto pattern = noise(Rcpp::Named("coordinates") = mapping, Rcpp::Named("seed") = 81,
                           Rcpp::Named("octaves") = 5);
      auto checks = checker("coral", "steelblue", mapping);
      auto ramp = gradient(0, 1, mapping, Rcpp::Named("axis") = "x");
      auto nested = mix(checks, scale(checks, .7), mix(pattern, ramp, .3));
      auto directed =
          direction("red", "blue", Rcpp::Named("direction") = Rcpp::NumericVector::create(1, 2, -3),
                    Rcpp::Named("space") = space == "uv" ? "world" : space);
      std::vector<Rcpp::List> descriptions{Rcpp::as<Rcpp::List>(nested),
                                           Rcpp::as<Rcpp::List>(directed),
                                           Rcpp::as<Rcpp::List>(channel(nested, "luminance"))};
      for (const auto &description : descriptions) {
        TextureCache cache;
        TextureGraphBuilder builder(cache);
        auto node = builder.Build(description);
        auto texture = std::make_shared<graph_texture>(node);
        auto material = std::make_shared<diffuse_material>(texture);
        Transform identity;
        hitable_list world, lights;
        world.add(std::make_shared<xy_rect>(-2, 2, -2, 2, 0, material, nullptr, nullptr, &identity,
                                            &identity, false));
        camera cam({0, 0, 4}, {0, 0, 0}, {0, 1, 0}, 30, 1, 0, 4, 0, 1, 1);
        WavefrontScene scene;
        WavefrontReport report;
        auto session = PrepareWavefront(world, lights, cam, 1, 0, scene, report);
        expect_true(bool(session));
        if (!session)
          continue;
        std::vector<WFTextureProbe> probes;
        std::vector<point3f> expected;
        for (int i = 0; i < 97; ++i) {
          TextureEvalContext h;
          h.u = (i - 43) / 17.f;
          h.v = (i % 19 - 9) / 7.f;
          h.p = point3f(h.u, Float(.13 * i), h.v);
          h.object_p = point3f(h.v, -h.u, Float(-.07 * i));
          h.geometric_normal = unit_vector(normal3f(.3f + i * .01f, -.6f, .7f));
          h.object_normal =
              normal3f(-h.geometric_normal[2], h.geometric_normal[1], h.geometric_normal[0]);
          auto pack = [](const auto &v) {
            return WFVector{float(v[0]), float(v[1]), float(v[2]), 0};
          };
          probes.push_back({{float(h.u), float(h.v), float(scene.materials[0].texture), 0},
                            pack(h.p),
                            pack(h.object_p),
                            pack(h.geometric_normal),
                            pack(h.object_normal),
                            {}});
          expected.push_back(node->Evaluate(h));
        }
        auto actual = ProbeMetalTextures(scene, probes);
        double error = 0;
        for (size_t i = 0; i < actual.size(); ++i)
          for (int c = 0; c < 3; ++c)
            error = std::max(error, std::abs(double((&actual[i].x)[c]) - expected[i][c]));
        expect_true(error < 2e-5);
      }
    }
  }
  test_that("mapped diffuse GPU values and PDFs match the energy-conserving CPU model") {
    if (!MetalBVHAvailable())
      return;
    std::vector<WFBsdfProbe> probes;
    std::vector<WFBsdfResult> expected;
    for (double tilt : {0.0, .2, .8, 1.3})
      for (double sigma : {0.0, .7, 1.5})
        for (int j = 0; j < 13; ++j) {
          normalmap::Vector base(0, 0, 1), normal(std::sin(tilt), 0, std::cos(tilt));
          auto view = normalmap::normalized(normalmap::Vector(-.7 + .1 * j, .2, 1));
          auto incoming = normalmap::normalized(normalmap::Vector(.9 - .12 * j, -.3, 1));
          normalmap::Model model(base, normal, view);
          normalmap::DiffuseChild child(sigma);
          WFBsdfProbe p;
          double r = std::min(1.0, sigma * (2 / M_PI));
          double a = 1 / (1 + (.5 - 2 / (3 * M_PI)) * r), b = r * a;
          double avg = b * ((.5 - 2 / (3 * M_PI)) - (2. / 3 - 28 / (15 * M_PI)));
          p.material.diffuse = {float(a), float(b), float(avg), 0};
          auto pack = [](const auto &v) {
            return WFVector{float(v[0]), float(v[1]), float(v[2]), 0};
          };
          p.rho = {.2f, .5f, .9f, 0};
          p.context = pack(normal);
          p.view = pack(view);
          p.direction = pack(incoming);
          p.random = {.37f, .61f, .23f, .79f};
          WFBsdfResult result;
          for (int c = 0; c < 3; ++c)
            (&result.evaluated.x)[c] = model.eval_raw(incoming, child, (&p.rho.x)[c]) * incoming[2];
          result.evaluated.w = model.pdf(incoming);
          result.direction = pack(model.sample(p.random.x, p.random.y, p.random.z, p.random.w));
          probes.push_back(p);
          expected.push_back(result);
        }
    auto actual = ProbeMetalMaterials(probes, {});
    double value_error = 0, direction_error = 0;
    for (size_t i = 0; i < actual.size(); ++i) {
      for (int c = 0; c < 4; ++c)
        value_error = std::max(value_error, std::abs(double((&actual[i].evaluated.x)[c]) -
                                                     (&expected[i].evaluated.x)[c]));
      for (int c = 0; c < 3; ++c)
        direction_error = std::max(direction_error, std::abs(double((&actual[i].direction.x)[c]) -
                                                             (&expected[i].direction.x)[c]));
    }
    expect_true(value_error < 2e-5);
    expect_true(direction_error < 2e-5);
  }
}
#endif
