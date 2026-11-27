#ifdef NOT_CRAN
#include "material.h"
#include "../hitables/rectangle.h"
#include <testthat.h>
#include <iostream>

namespace {
vec3f glossy_direction(double cosine, double phi) {
  const double radius = std::sqrt(std::max(0.0, 1-cosine*cosine));
  return vec3f(radius*std::cos(phi), radius*std::sin(phi), cosine);
}
hit_record glossy_hit(double tilt) {
  hit_record hit;
  hit.p = point3f(0);
  hit.normal = hit.geometric_normal = hit.base_shading_normal = normal3f(0, 0, 1);
  hit.bump_normal = hit.physical_shading_normal = normal3f(std::sin(tilt), 0, std::cos(tilt));
  hit.has_bump = tilt != 0;
  hit.u = hit.v = 0;
  return hit;
}
MicrofacetDistribution* glossy_distribution(bool beckmann) {
  return beckmann ? static_cast<MicrofacetDistribution*>(
                       new BeckmannDistribution(.25f, .6f, nullptr, false))
                  : new TrowbridgeReitzDistribution(.25f, .6f, nullptr, false);
}
}

context("Glossy energy and sampling") {
  test_that("a height image reaches the glossy physical frame through a real intersection") {
    auto mat = std::make_shared<glossy>(std::make_shared<constant_texture>(point3f(1)),
        glossy_distribution(false), point3f(.035), point3f(1));
    std::vector<Float> pixels(256*256);
    for (int y=0; y<256; ++y) for (int x=0; x<256; ++x)
      pixels[y*256+x] = Float(x)/255;
    auto image = std::make_shared<HeightImage>(256,256,pixels);
    auto bump = std::make_shared<bump_texture>(image,256,256,1,4,1,1,0,0);
    Transform identity;
    xy_rect plane(-2,2,-2,2,0,mat,nullptr,bump,&identity,&identity,false);
    random_gen rng(34);
    hit_record hit;
    Ray ray(point3f(3,0,1),vec3f(-3,0,-1));
    expect_true(plane.hit(ray,0,100,hit,rng));
    expect_true(hit.has_bump);
    expect_true(hit.physical_shading_normal[0] < -.7);
    expect_true(hit.physical_shading_normal[2] > .7);
    expect_true((hit.bump_normal-hit.physical_shading_normal).length() < 1e-6);
  }
  test_that("the smooth lemon agrees with the Fresnel-blend reference energy") {
    glossy material(std::make_shared<constant_texture>(point3f(1,.8,.02)),
        new TrowbridgeReitzDistribution(.050625f,.050625f,nullptr,false),
        point3f(.035),point3f(1));
    const auto hit = glossy_hit(0);
    const Ray ray(point3f(0),vec3f(0,0,-1));
    double energy[3]{};
    // At normal incidence this isotropic lobe is independent of azimuth.
    for (int i=0;i<8192;++i) {
      const auto f=material.f(ray,hit,glossy_direction((i+.5)/8192,0));
      for (int c=0;c<3;++c) energy[c]+=f[c]*(2*M_PI/8192);
    }
    const double reference[]={.967802,.780833,.0516555};
    for (int c=0;c<3;++c) expect_true(std::abs(energy[c]-reference[c])<.0002);
    std::cout << "Smooth lemon furnace RGB: " << energy[0] << ", " << energy[1]
              << ", " << energy[2] << "\n";
  }
  test_that("smooth and tilted glossy surfaces do not exceed a white furnace") {
    double maximum = 0;
    auto white = std::make_shared<constant_texture>(point3f(1));
    for (bool beckmann : {false, true})
      for (Float reflectance : {0.f, .035f, 1.f}) {
        glossy material(white, glossy_distribution(beckmann), point3f(reflectance), point3f(1));
        for (double tilt : {0., .7, 1.3})
          for (double view_x : {0., -.95, .95}) {
            auto hit = glossy_hit(tilt);
            const vec3f view(view_x, 0, std::sqrt(1-view_x*view_x));
            const Ray ray(point3f(0), -view);
            double energy = 0;
            bool finite = true;
            const int nz = 192, np = 384;
            for (int z = 0; z < nz; ++z)
              for (int p = 0; p < np; ++p) {
                const vec3f wi = glossy_direction((z+.5)/nz, 2*M_PI*(p+.5)/np);
                const double value = material.f(ray, hit, wi)[0];
                finite = finite && std::isfinite(value) && value >= 0;
                energy += value * (2*M_PI/(nz*np));
              }
            maximum = std::max(maximum, energy);
            expect_true(finite);
            expect_true(energy <= 1.003);
            expect_true(material.f(ray, hit, vec3f(0,0,-1))[0] == 0);
          }
      }
    std::cout << "Glossy maximum furnace response: " << maximum << "\n";
  }

  test_that("glossy sample frequencies match the PDF including facet paths and null mass") {
    auto white = std::make_shared<constant_texture>(point3f(1));
    for (bool beckmann : {false, true}) for (double tilt : {0., 1.2})
      for (double view_x : {-.95, .6}) {
        glossy material(white, glossy_distribution(beckmann), point3f(.035), point3f(1));
        auto hit = glossy_hit(tilt);
        const Ray ray(point3f(0), -vec3f(view_x, 0, std::sqrt(1-view_x*view_x)));
        random_gen rng(361);
        RandomSampler sampler(rng);
        scatter_record scattering;
        expect_true(material.scatter(ray, hit, scattering, rng));
        double expected[8]{}, observed[8]{}, integrated_energy = 0, sampled_energy = 0;
        const int nz = 256, np = 512, count = 131072;
        for (int z = 0; z < nz; ++z) for (int p = 0; p < np; ++p) {
          const vec3f wi = glossy_direction((z+.5)/nz, 2*M_PI*(p+.5)/np);
          const int bin = 4*(z >= nz/2) + 2*(wi[0] > 0) + (wi[1] > 0);
          expected[bin] += scattering.pdf_ptr->value(wi, rng) * (2*M_PI/(nz*np));
          integrated_energy += scattering.distribution->Evaluate(wi)[0] * (2*M_PI/(nz*np));
        }
        bool finite = true;
        for (int i = 0; i < count; ++i) {
          bool diffuse = false;
          const vec3f wi = i % 2 ? scattering.pdf_ptr->generate(&sampler, diffuse)
                                  : scattering.pdf_ptr->generate(rng, diffuse);
          if (!(wi.squared_length() > 0)) continue;
          const double pdf = scattering.pdf_ptr->value(wi, rng);
          finite = finite && std::isfinite(pdf) && pdf > 0 && wi[2] > 0;
          if (!(pdf > 0)) continue;
          const int bin = 4*(wi[2] >= .5) + 2*(wi[0] > 0) + (wi[1] > 0);
          observed[bin] += 1.0/count;
          sampled_energy += scattering.distribution->Evaluate(wi)[0] / (pdf*count);
        }
        double mass = 0, frequency_error = 0;
        for (int bin = 0; bin < 8; ++bin) {
          mass += expected[bin];
          frequency_error = std::max(frequency_error, std::abs(expected[bin]-observed[bin]));
        }
        expect_true(finite);
        expect_true(mass <= 1.003);
        expect_true(frequency_error < .006);
        expect_true(std::abs(sampled_energy-integrated_energy) < .012);
        std::cout << "Glossy PDF: Beckmann=" << beckmann << " tilt=" << tilt
                  << " view=" << view_x << " mass=" << mass << " bin error=" << frequency_error
                  << " quadrature=" << integrated_energy << " sampled=" << sampled_energy << "\n";
      }
  }
}
#endif
