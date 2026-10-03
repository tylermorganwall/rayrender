#ifdef NOT_CRAN
#include "rough_dielectric.h"
#include "material.h"
#include <testthat.h>

context("Rough dielectric radiance transport") {
  test_that("normal-incidence reflection and transmission have physical path weights") {
    TrowbridgeReitzDistribution dist(.3f,.3f,nullptr,false);
    point2f alpha(.1f,.3f);
    for(Float side : {-1.f,1.f}) {
      vec3f view(0,0,side);
      auto reflect=EvaluateRoughDielectric(view,view,1.5f,dist,alpha);
      auto transmit=EvaluateRoughDielectric(view,-view,1.5f,dist,alpha);
      Float ratio=side > 0 ? 1/1.5f : 1.5f;
      expect_true(reflect.pdf > 0);
      expect_true(transmit.pdf > 0);
      expect_true(std::abs(reflect.f_cos/reflect.pdf-1) < 1e-5f);
      expect_true(std::abs(transmit.f_cos/transmit.pdf-Sqr(ratio)) < 1e-5f);
      // Fresnel R=.04 and the transmission half-vector Jacobian are known
      // analytically here; this checks values independently of their ratio.
      Float D=1/(Float(M_PI)*alpha[0]*alpha[1]);
      expect_true(std::abs(reflect.pdf/(D*.04f/4)-1) < 1e-5f);
      expect_true(std::abs(transmit.pdf/(D*.96f/Sqr(1-ratio))-1) < 1e-5f);
    }
  }

  test_that("sampling matches integrated reflection and transmission densities") {
    // Integrate solid-angle PDFs independently of sampling; invalid VNDF
    // reflections are allowed to lose probability mass, not to be resampled.
    for(bool beckmann : {false,true}) for(Float side : {-1.f,1.f}) {
      std::unique_ptr<MicrofacetDistribution> dist(beckmann ?
        static_cast<MicrofacetDistribution*>(new BeckmannDistribution(.3f,.3f,nullptr,false)) :
        static_cast<MicrofacetDistribution*>(new TrowbridgeReitzDistribution(.3f,.3f,nullptr,false)));
      point2f alpha(.35f,.6f);
      random_gen rng(382); RandomSampler sampler(rng);
      vec3f view(.6f,0,side*.8f);
      micro_transmission_pdf density(normal3f(0,0,1),-view,dist.get(),1.5f,0,0);
      density.alphas=alpha;
      double integrated[2]={0,0}, observed[2]={0,0};
      const int nz=600, nphi=1200, samples=100000;
      for(int z=0;z<nz;++z) for(int p=0;p<nphi;++p) {
        Float cosine=-1+2*(z+.5f)/nz, phi=2*Float(M_PI)*(p+.5f)/nphi;
        vec3f d(std::sqrt(1-cosine*cosine)*std::cos(phi),
                std::sqrt(1-cosine*cosine)*std::sin(phi),cosine);
        integrated[cosine*side < 0] += density.value(d,rng)*4*M_PI/(nz*nphi);
      }
      bool bounce=false; Float largest_weight=0;
      for(int i=0;i<samples;++i) {
        vec3f d=i%2 ? density.generate(rng,bounce) : density.generate(&sampler,bounce);
        if (!(d.squared_length()>0)) continue;
        ++observed[d[2]*side < 0];
        Float pdf=density.value(d,rng);
        auto f=EvaluateRoughDielectric(density.wi,unit_vector(density.uvw.world_to_local(d)),
          1.5f,*dist,alpha);
        expect_true((std::isfinite(pdf) && pdf>0));
        expect_true(std::abs(pdf-density.value(d,&sampler)) < 1e-6f);
        Float weight=f.f_cos/pdf;
        if(f.transmission) weight /= Sqr(side>0 ? 1/1.5f : 1.5f);
        largest_weight=std::max(largest_weight,weight);
      }
      expect_true(largest_weight <= 1.0001f);
      expect_true(std::abs(integrated[0]-observed[0]/samples) < .008);
      expect_true(std::abs(integrated[1]-observed[1]/samples) < .008);
    }
  }

  test_that("total internal reflection, grazing, and index matching stay finite") {
    TrowbridgeReitzDistribution dist(.3f,.3f,nullptr,false);
    vec3f view(.9f,0,-std::sqrt(.19f)), reflected(-.9f,0,-std::sqrt(.19f));
    auto value=EvaluateRoughDielectric(view,reflected,1.5f,dist,point2f(.1f,.1f));
    expect_true(value.pdf>0);
    expect_true(FrDielectric(std::abs(view[2]),1.5f)==1);
    expect_true(EvaluateRoughDielectric(vec3f(1,0,0),reflected,1.5f,dist,point2f(.1f,.1f)).pdf==0);
    auto tex=std::make_shared<constant_texture>(point3f(1));
    MicrofacetTransmission mat(tex,new TrowbridgeReitzDistribution(.3f,.3f,nullptr,false),point3f(1),point3f(0));
    hit_record h; h.normal=h.geometric_normal=normal3f(0,0,1);h.has_bump=false;
    h.p=point3f(0);h.pError=vec3f(0);h.u=h.v=0;
    Ray ray(point3f(0,0,1),-vec3f(.6f,0,.8f));
    scatter_record s; random_gen rng(77); RandomSampler sampler(rng);
    expect_true(mat.scatter(ray,h,s,rng));
    expect_true((s.is_specular && s.is_transmission && s.eta==1));
    expect_true((s.specular_ray.direction()-ray.direction()).length()<1e-6f);
    expect_true(mat.scatter(ray,h,s,&sampler));
    expect_true((s.specular_ray.direction()-ray.direction()).length()<1e-6f);
  }

  test_that("material and both PDF adapters retain the radiance factor") {
    auto tex=std::make_shared<constant_texture>(point3f(1));
    MicrofacetTransmission mat(tex,new TrowbridgeReitzDistribution(.3f,.3f,nullptr,false),point3f(1.5f),point3f(0));
    hit_record hit; hit.p=point3f(0); hit.pError=vec3f(0); hit.u=hit.v=0;
    hit.normal=hit.geometric_normal=normal3f(0,0,1); hit.has_bump=false;
    random_gen rng(90); RandomSampler sampler(rng);
    for(Float side : {-1.f,1.f}) for(bool sampled : {false,true}) {
      Ray ray(point3f(0,0,side),vec3f(0,0,-side));
      scatter_record s;
      if(sampled) mat.scatter(ray,hit,s,&sampler); else mat.scatter(ray,hit,s,rng);
      vec3f outgoing(0,0,-side);
      auto value=mat.f(ray,hit,outgoing);
      Float pdf=s.pdf_ptr->value(outgoing,rng);
      expect_true(pdf>0);
      expect_true(std::abs(value[0]/pdf-Sqr(side>0 ? 1/1.5f : 1.5f))<1e-5f);
    }
  }
}
#endif
