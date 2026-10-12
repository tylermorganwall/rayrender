#ifdef NOT_CRAN
#include "pdf.h"
#include <testthat.h>

context("Glossy sampling frames") {
  test_that("cosine samples stay on the view side of rotated shading normals") {
    TrowbridgeReitzDistribution distribution(.2f, .35f, nullptr, false);
    for (normal3f normal : {normal3f(0, 0, 1), normal3f(1, 0, 0),
                            unit_vector(normal3f(.4f, -.8f, .3f))}) {
      onb frame;
      frame.build_from_w_normalized(normal);
      for (Float side : {-1.f, 1.f}) {
        vec3f view = frame.local_to_world(vec3f(.6f, 0, side * .8f));
        glossy_pdf density(normal, -view, &distribution, 0, 0);
        for (bool use_sampler : {false, true}) {
          random_gen rng(392);
          RandomSampler sampler(rng);
          int count = 0, wrong_side = 0;
          vec3f mean(0);
          for (int i = 0; i < 16384; ++i) {
            bool diffuse = false;
            vec3f sampled = use_sampler ? density.generate(&sampler, diffuse)
                                       : density.generate(rng, diffuse);
            if (!diffuse) continue;
            vec3f local = frame.world_to_local(sampled);
            wrong_side += local[2] * side <= 0;
            mean += local;
            ++count;
          }
          expect_true((count > 7000 && count < 9400));
          expect_true(wrong_side == 0);
          mean /= Float(count);
          // A cosine hemisphere has zero tangent moments and E[|cos|] = 2/3.
          expect_true((std::abs(mean[0]) < .025f && std::abs(mean[1]) < .025f));
          expect_true(std::abs(mean[2] * side - 2.f / 3) < .025f);
        }
      }
    }
  }
}
#endif
