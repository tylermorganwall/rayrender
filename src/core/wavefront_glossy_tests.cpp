#if defined(NOT_CRAN) && defined(RAY_HAS_METAL_BVH)
#include "../hitables/rectangle.h"
#include "camera.h"
#include "hlbvh.h"
#include "wavefront.h"
#include <testthat.h>
#include <iostream>

context("Metal glossy sampling") {
  test_that("mapped glossy GPU samples reproduce solid-angle density and furnace energy") {
    if (!MetalBVHAvailable()) return;
    Transform identity;
    camera cam({0,0,4}, {0,0,0}, {0,1,0}, 30, 1, 0, 4, 0, 1, 1);
    for (bool beckmann : {false, true}) {
      auto distribution = beckmann ? static_cast<MicrofacetDistribution*>(
          new BeckmannDistribution(.25f, .6f, nullptr, false))
          : new TrowbridgeReitzDistribution(.25f, .6f, nullptr, false);
      auto material = std::make_shared<glossy>(
          std::make_shared<constant_texture>(point3f(1)), distribution, point3f(.035), point3f(1));
      hitable_list world, lights;
      world.add(std::make_shared<xy_rect>(-1,1,-1,1,0,material,nullptr,nullptr,
                                          &identity,&identity,false));
      WavefrontScene scene;
      WavefrontReport report;
      auto session = PrepareWavefront(world,lights,cam,1,0,scene,report);
      expect_true(bool(session));
      if (!session) continue;
      for (float view_x : {-.95f, .6f}) {
        const int nz = 256, np = 256, count = nz*np;
        std::vector<WFBsdfProbe> inputs(count);
        random_gen rng(919);
        for (int z = 0; z < nz; ++z) for (int p = 0; p < np; ++p) {
          auto& input = inputs[z*np+p];
          float cosine = (z+.5f)/nz, phi = 2*M_PI*(p+.5f)/np;
          input.material = scene.materials[0];
          input.rho = {1,1,1,0};
          input.view = {view_x,0,std::sqrt(1-view_x*view_x),0};
          input.direction = {std::sqrt(1-cosine*cosine)*std::cos(phi),
                             std::sqrt(1-cosine*cosine)*std::sin(phi),cosine,0};
          input.context = {float(std::sin(1.2)),0,float(std::cos(1.2)),1};
          input.random = {float(rng.unif_rand()),float(rng.unif_rand()),
                          float(rng.unif_rand()),float(rng.unif_rand())};
        }
        const auto result = ProbeMetalMaterials(inputs, scene.material_data);
        double density[8]{}, frequency[8]{}, integral=0, estimate=0, largest_weight=0;
        bool finite=true;
        for (int i=0;i<count;++i) {
          const auto wi=inputs[i].direction, sampled=result[i].direction;
          int bin=4*(wi.z>=.5f)+2*(wi.x>0)+(wi.y>0);
          density[bin]+=result[i].evaluated.w*(2*M_PI/count);
          integral+=result[i].evaluated.x*(2*M_PI/count);
          finite=finite && std::isfinite(result[i].weight.x) && result[i].weight.x>=0;
          if (sampled.w<=0) continue; // Includes primary-facet and macro null mass.
          bin=4*(sampled.z>=.5f)+2*(sampled.x>0)+(sampled.y>0);
          frequency[bin]+=1.0/count;
          estimate+=result[i].weight.x/count;
          largest_weight=std::max(largest_weight,double(result[i].weight.x));
        }
        double error=0,mass=0;
        for (int bin=0;bin<8;++bin) {
          error=std::max(error,std::abs(density[bin]-frequency[bin]));
          mass+=density[bin];
        }
        expect_true(finite);
        expect_true(mass<=1.003);
        expect_true(integral<=1.003);
        expect_true(error<.008);
        expect_true(std::abs(estimate-integral)<.015);
        std::cout << "Metal glossy: Beckmann=" << beckmann << " view=" << view_x
                  << " quadrature=" << integral << " sampled=" << estimate
                  << " PDF bin error=" << error << " largest weight=" << largest_weight << "\n";
      }
    }
  }
}
#endif
#if defined(NOT_CRAN) && defined(RAY_HAS_METAL_BVH)
context("Metal order-16 normal mapping") {
  test_that("conditional GPU closures recover the analytic white and absorbing diffuse series") {
    if (!MetalBVHAvailable()) return;
    const int count=131072;
    for (float tilt : {.5235987756f, 1.0471975512f}) for (float rho : {1.f,.5f}) {
      std::vector<WFBsdfProbe> inputs(count);
      random_gen rng(334);
      for (int i=0;i<count;++i) {
        auto& p=inputs[i];
        p.material={};
        p.rho={rho,rho,rho,0}; p.view={0,0,1,0};
        float z=rng.unif_rand(), phi=2*M_PI*rng.unif_rand();
        p.direction={std::sqrt(1-z*z)*std::cos(phi),std::sqrt(1-z*z)*std::sin(phi),z,0};
        p.random={float(rng.unif_rand()),float(rng.unif_rand()),float(rng.unif_rand()),float(rng.unif_rand())};
        p.context={std::sin(tilt),0,std::cos(tilt),2};
      }
      const auto outputs=ProbeMetalMaterials(inputs, {});
      double sum=0,squared=0,direct=0,direct_squared=0;
      bool finite=true;
      for(const auto& out:outputs) {
        const double sample=out.direction.w>0?out.weight.x:0;
        const double light=2*M_PI*out.evaluated.x;
        sum+=sample;squared+=sample*sample;direct+=light;direct_squared+=light*light;
        finite=finite&&std::isfinite(sample)&&sample>=0&&std::isfinite(light)&&light>=0;
      }
      const double ratio=rho*(1-std::cos(tilt));
      const double expected=rho*std::cos(tilt)*(1-std::pow(ratio,8))/(1-ratio);
      const double se=std::sqrt(std::max(0.0,squared/count-(sum/count)*(sum/count))/count);
      const double direct_se=std::sqrt(std::max(0.0,direct_squared/count-(direct/count)*(direct/count))/count);
      expect_true(finite);
      expect_true(std::abs(sum/count-expected)<6*se+.001);
      expect_true(std::abs(direct/count-expected)<6*direct_se+.001);
      std::cout << "Metal order 16: tilt=" << tilt << " rho=" << rho << " reference=" << expected
                << " sample=" << sum/count << " direct=" << direct/count << "\n";
    }
  }
}
#endif
