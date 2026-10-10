#if defined(NOT_CRAN) && defined(RAY_HAS_METAL_BVH)
#include "../hitables/infinite_area_light.h"
#include "../hitables/rectangle.h"
#include "camera.h"
#include "hlbvh.h"
#include "wavefront.h"
#include <cstring>
#include <testthat.h>

context("Metal image environment proposals") {
  test_that("alias frequencies, solid-angle PDFs, seams and rotation agree") {
    if (!MetalBVHAvailable())
      return;
    Transform identity;
    camera cam({0, 0, 4}, {0, 0, 0}, {0, 1, 0}, 30, 1, 0, 4, 0, 1, 1);
    for (bool black : {false, true}) {
      std::vector<Float> pixels(8 * 4 * 3, 0);
      if (!black) {
        // Bright texels at the seam and pole; bilinear spill must retain support.
        for (size_t i : {size_t(0), size_t(7), size_t(21)})
          for (int c = 0; c < 3; ++c)
            pixels[3 * i + c] = i == 21 ? 1 : 20;
      }
      auto texture = std::make_shared<latlong_image_texture>(pixels.data(), 8, 4, 3);
      auto light = std::make_shared<ImageInfiniteLight>(texture, 8, 4, 0);
      auto material =
          std::make_shared<diffuse_material>(std::make_shared<constant_texture>(point3f(.5)));
      hitable_list world, lights;
      world.add(std::make_shared<xy_rect>(-2, 2, -2, 2, 0, material, nullptr, nullptr, &identity,
                                          &identity, false));
      world.add(std::make_shared<InfiniteAreaLight>(light, 100, point3f(0), &identity, &identity));
      WavefrontScene scene;
      WavefrontReport report;
      auto session = PrepareWavefront(world, lights, cam, 1, 0, scene, report);
      expect_true(bool(session));
      if (!session)
        continue;
      auto emitter = scene.lights[0];
      expect_true((emitter.alias_offset == UINT32_MAX) == black);
      if (black) {
        auto result = ProbeMetalEnvironment(scene, emitter, {{0, 0, {}}, {UINT32_MAX, .5, {}}});
        expect_true(std::abs(result[0].w - 1 / (4 * M_PI)) < 1e-7);
        expect_true(std::abs(result[1].w - result[0].w) < 1e-7);
        continue;
      }
      double integral = 0;
      for (size_t i = 0; i < 32; ++i) {
        double area =
            2 * M_PI / 8 * (std::cos(M_PI * (i / 8) / 4) - std::cos(M_PI * (i / 8 + 1) / 4));
        integral += area * scene.texels[emitter.alias_offset + i].z;
      }
      expect_true(std::abs(integral - 1) < 1e-6);
      const size_t count = 262144;
      std::vector<WFEnvironmentProbe> random(count);
      for (size_t i = 0; i < count; ++i)
        random[i] = {uint32_t((i << 14) | ((i * 7919) & 16383)),
                     float(((i * 7919) % count + .5) / count),
                     {}};
      auto result = ProbeMetalEnvironment(scene, emitter, random);
      std::vector<size_t> frequency(32);
      double maximum_pdf_error = 0, length_error = 0;
      size_t edge_samples = 0;
      for (const auto &sample : result) {
        double u = .5 - std::atan2(-sample.x, -sample.z) / (2 * M_PI);
        double v = .5 + std::asin(std::clamp(double(sample.y), -1., 1.)) / M_PI;
        size_t cell = std::min(size_t(u * 8), size_t(7)) + 8 * std::min(size_t(v * 4), size_t(3));
        ++frequency[cell];
        // Density is discontinuous at cell edges. Host double atan2 and GPU
        // float atan2 can classify a float direction on opposite sides; bound
        // that ambiguity explicitly instead of comparing unequal neighboring PDFs.
        double edge =
            std::min(std::abs(u * 8 - std::round(u * 8)), std::abs(v * 4 - std::round(v * 4)));
        if (edge < 1e-6)
          ++edge_samples;
        else
          maximum_pdf_error =
              std::max(maximum_pdf_error,
                       std::abs(double(sample.w) - scene.texels[emitter.alias_offset + cell].z));
        length_error = std::max(length_error,
                                std::abs(double(sample.x) * sample.x + double(sample.y) * sample.y +
                                         double(sample.z) * sample.z - 1));
      }
      double maximum_mass_error = 0;
      for (size_t i = 0; i < 32; ++i) {
        double area =
            2 * M_PI / 8 * (std::cos(M_PI * (i / 8) / 4) - std::cos(M_PI * (i / 8 + 1) / 4));
        maximum_mass_error =
            std::max(maximum_mass_error, std::abs(double(frequency[i]) / count -
                                                  area * scene.texels[emitter.alias_offset + i].z));
      }
      expect_true(maximum_mass_error < .0001);
      expect_true(maximum_pdf_error < 1e-5);
      expect_true(edge_samples < 16);
      expect_true(length_error < 1e-5);
      // Rotate the emitter without rebuilding the alias table; samples and PDFs
      // must rotate together. This is the interactive environment edit path.
      emitter.forward = {0, 0, -1, 0};
      emitter.up = {1, 0, 0, 0};
      auto rotated = ProbeMetalEnvironment(scene, emitter, random);
      double error = 0;
      for (size_t i = 0; i < count; ++i)
        error = std::max(error, double(std::abs(rotated[i].x - result[i].z) +
                                       std::abs(rotated[i].z + result[i].x) +
                                       std::abs(rotated[i].w - result[i].w)));
      expect_true(error < 1e-5);
    }
  }
}
#endif
