#if defined(NOT_CRAN) && defined(RAY_HAS_METAL_BVH)
#include "../hitables/rectangle.h"
#include "../materials/texturegraph.h"
#include "camera.h"
#include "hlbvh.h"
#include "wavefront.h"
#include <testthat.h>

context("Metal texture evaluation") {
  test_that("legacy noise, gradients and nested checkers match native evaluation") {
    if (!MetalBVHAvailable()) return;
    auto noise = std::make_shared<noise_texture>(1.7, point3f(.1, .3, .8),
                                                point3f(.9, .7, .2), .73, 8.3);
    std::vector<std::shared_ptr<texture>> cases{noise};
    for (bool hsv : {false, true}) {
      for (bool transpose : {false, true})
        cases.push_back(std::make_shared<gradient_texture>(point3f(.9, .1, .3),
                         point3f(.1, .2, .8), transpose, hsv));
      cases.push_back(std::make_shared<world_gradient_texture>(point3f(-.4, .2, -.7),
                        point3f(.6, 1.1, .4), point3f(.9, .1, .3), point3f(.1, .2, .8), hsv));
    }
    std::shared_ptr<texture> nested = noise;
    for (int i = 0; i < 5; ++i)
      nested = std::make_shared<checker_texture>(cases[i + 1], nested, .71 + .17 * i);
    cases.push_back(nested);
    // A graph leaf under a legacy checker exercises dispatch between the two
    // representations, including the texture-feature specialization flags.
    auto ns = Rcpp::Environment::namespace_env("rayrender");
    Rcpp::Function constant = ns["texture_constant"];
    TextureCache cache;
    TextureGraphBuilder builder(cache);
    auto graph = std::make_shared<graph_texture>(builder.Build(constant("coral")));
    cases.push_back(std::make_shared<checker_texture>(graph, nested, .93));
    Transform identity;
    camera cam({0, 0, 4}, {0, 0, 0}, {0, 1, 0}, 30, 1, 0, 4, 0, 1, 1);
    for (const auto &tex : cases) {
      hitable_list world, lights;
      world.add(std::make_shared<xy_rect>(-2, 2, -2, 2, 0,
        std::make_shared<diffuse_material>(tex), nullptr, nullptr, &identity, &identity, false));
      WavefrontScene scene;
      WavefrontReport report;
      auto session = PrepareWavefront(world, lights, cam, 1, 0, scene, report);
      expect_true(bool(session));
      if (!session) continue;
      std::vector<WFTextureProbe> probes;
      std::vector<point3f> expected;
      for (int i = 0; i < 193; ++i) {
        float u = (i % 23 - 5) / 17.f, v = (i % 31 - 8) / 19.f;
        WFVector p{(i - 96) * .173f, (i % 29 - 14) * .231f, (i % 13 - 6) * .357f, 0};
        if (i == 0) p = {}; // Checker's zero-sine convention and lattice origin.
        probes.push_back({{u, v, float(scene.materials[0].texture), 0}, p, {}, {}, {}, {}});
        expected.push_back(tex->value(u, v, point3f(p.x, p.y, p.z)));
      }
      auto actual = ProbeMetalTextures(scene, probes);
      double error = 0;
      expect_true(actual.size() == expected.size());
      for (size_t i = 0; i < actual.size(); ++i)
        for (int c = 0; c < 3; ++c) {
          double value = (&actual[i].x)[c];
          if (!std::isfinite(value)) error = INFINITY;
          error = std::max(error, std::abs(value - expected[i][c]));
        }
      expect_true(error < 5e-5);
    }
  }
  test_that("Ptex atlas lookups preserve face IDs and decode after interpolation") {
    if (!MetalBVHAvailable()) return;
    Rcpp::Function temporary = Rcpp::Environment::base_env()["tempfile"];
    Rcpp::Function remove = Rcpp::Environment::base_env()["unlink"];
    Rcpp::Function write = Rcpp::Environment::namespace_env("ptexr")["ptex_write"];
    Rcpp::Function ptex = Rcpp::Environment::namespace_env("rayrender")["texture_ptex"];
    auto filename = temporary(Rcpp::Named("fileext") = ".ptx");
    Rcpp::NumericVector image(32 * 32 * 3);
    image.attr("dim") = Rcpp::IntegerVector::create(32, 32, 3);
    for (int c = 0; c < 3; ++c)
      for (int y = 0; y < 32; ++y)
        for (int x = 0; x < 32; ++x)
          image[x + 32 * y + 1024 * c] = c == 0 ? (x + .5) / 32 :
                                         c == 1 ? (y + .5) / 32 : .375;
    write(filename, Rcpp::List::create(image, image));
    Transform identity;
    camera cam({0, 0, 4}, {0, 0, 0}, {0, 1, 0}, 30, 1, 0, 4, 0, 1, 1);
    for (bool scalar : {false, true})
      for (bool quantize : {false, true}) {
        Rcpp::List description = ptex(filename, Rcpp::Named("encoding") = "gamma 2.2",
                                      Rcpp::Named("type") = scalar ? "scalar" : "color");
        description["pbrt_encoding"] = quantize;
        TextureCache cache;
        TextureGraphBuilder builder(cache);
        auto node = builder.Build(description);
        auto material = std::make_shared<diffuse_material>(std::make_shared<graph_texture>(node));
        hitable_list world, lights;
        world.add(std::make_shared<xy_rect>(-2, 2, -2, 2, 0, material, nullptr, nullptr,
                                            &identity, &identity, false));
        WavefrontScene scene;
        WavefrontReport report;
        auto session = PrepareWavefront(world, lights, cam, 1, 0, scene, report);
        expect_true(bool(session));
        if (!session) continue;
        std::vector<WFTextureProbe> probes;
        std::vector<point3f> expected;
        for (int i = 0; i < 61; ++i) {
          TextureEvalContext context;
          context.face_index = i % 2;
          // Stay away from chart edges. An affine ramp must survive both
          // reconstruction filters; nonlinear decoding must follow filtering.
          context.u = .2f + .6f * (i % 17) / 16;
          context.v = .2f + .6f * (i % 11) / 10;
          probes.push_back({{float(context.u), float(context.v),
                             float(scene.materials[0].texture), float(context.face_index)},
                            {}, {}, {}, {}, {}});
          expected.push_back(node->Evaluate(context));
        }
        auto actual = ProbeMetalTextures(scene, probes);
        double error = 0;
        for (size_t i = 0; i < actual.size(); ++i)
          for (int c = 0; c < 3; ++c)
            error = std::max(error, std::abs(double((&actual[i].x)[c]) - expected[i][c]));
        // Quantization may straddle a single 8-bit bin after float rounding.
        expect_true(error < (quantize ? .008 : 2e-5));
      }
    remove(filename);
  }
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
