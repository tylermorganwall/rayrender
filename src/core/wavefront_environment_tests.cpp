#if defined(NOT_CRAN) && defined(RAY_HAS_METAL_BVH)
#include "../hitables/infinite_area_light.h"
#include "../hitables/rectangle.h"
#include "camera.h"
#include "hlbvh.h"
#include "wavefront.h"
#include <cstring>
#include <testthat.h>

context("Metal image environment proposals") {
  test_that("tiny solar-cone samples retain their solar mixture density") {
    if (!MetalBVHAvailable()) return;
    WavefrontScene scene;
    WFLight light; light.flags=128;
    double radius=.53*M_PI/360;
    double elevation=28*M_PI/180, azimuth=100*M_PI/180;
    vec3f sun(-std::sin(azimuth)*std::cos(elevation),std::sin(elevation),
               std::cos(azimuth)*std::cos(elevation));
    vec3f right=unit_vector(cross(sun,vec3f(1,0,0))), up=cross(sun,right);
    light.forward={float(sun[0]),float(sun[1]),float(sun[2]),.5};
    light.right={float(right[0]),float(right[1]),float(right[2]),float(2*std::pow(std::sin(radius/2),2))};
    light.up={float(up[0]),float(up[1]),float(up[2]),0};
    const size_t count=262144;
    std::vector<WFEnvironmentProbe> samples(count);
    for (size_t i=0;i<count;++i) {
      // Concentrate on the outer 1% of the solar proposal. Cancellation in
      // dot(w,sun) used to assign sky-only density to some of these samples.
      double u=.5*(.99+.01*(i+.5)/count);
      // Keep the float variate strictly inside the solar mixture branch:
      // integer-to-float conversion can otherwise round the last value to .5.
      samples[i]={std::min(uint32_t(u*4294967296.),uint32_t(0x7fffff80)),
                  float(((i*7919)%count+.5)/count),{}};
    }
    auto result=ProbeMetalEnvironment(scene,light,samples);
    double expected=.5/(4*M_PI)+.5/(2*M_PI*light.right.w), error=0;
    size_t lost=0;
    for (auto value : result) {
      lost+=value.w<expected*.5;
      error=std::max(error,std::abs(value.w/expected-1));
    }
    Rcpp::Rcout << "Solar rim samples assigned sky-only PDF: " << lost << "/" << count << "\n";
    expect_true(lost==0);
    expect_true(error<2e-6);
    for (auto &sample : samples) sample.padding[0]=1;
    auto old=ProbeMetalEnvironment(scene,light,samples);
    size_t old_lost=0;
    for (auto value : old) old_lost+=value.w<expected*.5;
    Rcpp::Rcout << "Old dot-product solar PDF loses " << old_lost << " rim samples\n";
    expect_true(old_lost>0);
  }
  test_that("finite light BVH sampling and emitter-hit probabilities agree") {
    if (!MetalBVHAvailable()) return;
    Transform identity;
    camera cam({0,0,4},{0,0,0},{0,1,0},30,1,0,4,0,1,1);
    for (bool black : {false,true}) {
      hitable_list world,lights;
      for (int i=0;i<10;++i) {
        auto emission=std::make_shared<constant_texture>(point3f(black ? 0 : (i==5 ? 100 : 1)));
        auto material=std::make_shared<diffuse_light>(emission,1,false);
        world.add(std::make_shared<xy_rect>(2*i,2*i+1,-1,1,0,material,nullptr,nullptr,
                                            &identity,&identity,false));
      }
      auto sky=std::make_shared<ImageInfiniteLight>(std::make_shared<constant_texture>(point3f(.1)),1,1,0);
      world.add(std::make_shared<InfiniteAreaLight>(sky,100,point3f(0),&identity,&identity));
      WavefrontScene scene; WavefrontReport report;
      auto session=PrepareWavefront(world,lights,cam,1,0,scene,report);
      expect_true(bool(session)); if (!session) continue;
      expect_true(scene.light_nodes.size()==39);
      const size_t count=65536;
      for (float x : {0.f,10.f,100.f}) {
        std::vector<WFVector> inputs(count);
        for (size_t i=0;i<count;++i) inputs[i]={x,0,3,float((i+.5)/count)};
        auto result=ProbeMetalLightSelection(scene,inputs);
        std::vector<size_t> frequencies(scene.lights.size());
        double error=0;
        bool valid=true;
        for (const auto &v : result) {
          valid &= v.x>=0 && v.x<scene.lights.size() && std::isfinite(v.y) && v.y>0;
          if (valid) ++frequencies[size_t(v.x)];
          error=std::max(error,std::abs(double(v.z)-1));
        }
        expect_true(valid); expect_true(error<2e-6);
        for (size_t i=0;i<frequencies.size();++i)
          error=std::max(error,std::abs(double(frequencies[i])/count-result[i].w));
        // Stratified inverse-CDF counts differ by at most a few endpoints:
        // the tree branch and 5% uniform component each partition [0,1).
        expect_true(error<5./count);
        expect_true(std::abs(result[0].w-1./scene.lights.size())<1e-7);
      }
    }
  }
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
