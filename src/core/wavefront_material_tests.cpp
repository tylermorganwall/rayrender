#if defined(NOT_CRAN) && defined(RAY_HAS_METAL_BVH)
#include "../hitables/hitablelist.h"
#include "../hitables/rectangle.h"
#include "../materials/openpbr.h"
#include "../volumes/subsurface.h"
#include "camera.h"
#include "hlbvh.h"
#include "wavefront.h"
#include <testthat.h>

namespace {
float probe_error(float actual, float expected) {
  return std::isfinite(actual) && std::isfinite(expected)
             ? std::abs(actual - expected) / std::max(1.f, std::abs(expected))
             : INFINITY;
}
} // namespace
context("Metal material evaluation") {
  test_that("guided SSS uses normalized RGB proposal mixtures") {
    if (!MetalBVHAvailable()) return;
    const size_t count=65536;
    std::vector<WFSubsurfaceProbe> inputs(count);
    SubsurfaceProposal proposal;
    proposal.pole={0,.3,.95}; proposal.axis=vec3f(0,0,1); proposal.guided=true;
    for (size_t i=0;i<count;++i) {
      double cosine=2*(i+.5)/count-1;
      inputs[i].pole={0,.3,.95,0}; inputs[i].axis={0,0,1,0};
      inputs[i].wi={float(std::sqrt(1-cosine*cosine)),0,float(cosine),0};
      inputs[i].rates={1,3,10,0};
      inputs[i].random={float((i+.5)/count),float(((i*7919)%count+.5)/count),
                        float(((i*3571)%count+.5)/count),0};
    }
    auto result=ProbeMetalSubsurface(inputs);
    double error=0, integral[3]{}, moment=0;
    for (size_t i=0;i<count;++i) {
      const auto &sample=result[i].origin;
      const auto &guide=result[i].direction;
      vec3f wi(inputs[i].wi.x,inputs[i].wi.y,inputs[i].wi.z);
      vec3f direction(sample.x,sample.y,sample.z);
      double pdf=0;
      for (int c=0;c<3;++c) {
        pdf+=proposal.DirectionPdf(c,direction)/3;
        error=std::max(error,std::abs(double((&guide.x)[c])-proposal.GuidePdf(c,wi)));
        integral[c]+=(&guide.x)[c]*4*M_PI/count;
      }
      error=std::max(error,std::abs(double(sample.w)-pdf));
      error=std::max(error,std::abs(double(guide.w)-1));
      error=std::max(error,std::abs(double(direction.length())-1));
      moment+=sample.z/count;
    }
    expect_true(error<2e-5);
    for (double value : integral) expect_true(std::abs(value-1)<2e-5);
    // E[cos] of each pole component is 1/k - 1/atanh(k).
    double expected=((1/.3-1/std::atanh(.3))+(1/.95-1/std::atanh(.95)))/6;
    expect_true(std::abs(moment-expected)<.001);
  }
  test_that("GPU BSDF values and reported PDFs match the native models") {
    if (!MetalBVHAvailable())
      return;
    auto color = std::make_shared<constant_texture>(point3f(.3, .5, .7));
    auto transmission = std::make_shared<constant_texture>(point3f(.4, .3, .2));
    std::vector<std::shared_ptr<material>> models;
    for (bool beckmann : {false, true}) {
      auto distribution = [&]() -> MicrofacetDistribution * {
        return beckmann ? static_cast<MicrofacetDistribution *>(
                              new BeckmannDistribution(.2, .35, nullptr, false))
                        : static_cast<MicrofacetDistribution *>(
                              new TrowbridgeReitzDistribution(.2, .35, nullptr, false));
      };
      models.push_back(std::make_shared<MicrofacetReflection>(
          color, distribution(), point3f(.4, .7, 1.2), point3f(2, 3, 4)));
      models.push_back(std::make_shared<MicrofacetTransmission>(color, distribution(), point3f(1.4),
                                                                point3f(.1, .2, .3)));
      models.push_back(std::make_shared<glossy>(color, distribution(), point3f(.08), point3f(.7)));
    }
    models.push_back(std::make_shared<translucent_material>(color, transmission));
    models.push_back(std::make_shared<hair>(point3f(.2, .5, .9), 1.55, .3, .4, 2));
    Rcpp::Function make_pbr = Rcpp::Environment::namespace_env("rayrender")["openpbr"];
    Rcpp::List descriptors = make_pbr();
    Rcpp::List descriptor = descriptors[0];
    for (int variant = 0; variant < 4; ++variant) {
      auto parameters = Rcpp::clone(Rcpp::as<Rcpp::List>(descriptor["openpbr"]));
      if (variant == 0) {
        parameters["base_metalness"] = 1.;
        parameters["thin_film_weight"] = .8;
      }
      if (variant == 1) {
        parameters["coat_weight"] = .7;
        parameters["coat_roughness"] = .2;
      }
      if (variant == 2)
        parameters["fuzz_weight"] = .8;
      if (variant == 3) {
        parameters["geometry_thin_walled"] = true;
        parameters["transmission_weight"] = .8;
      }
      models.push_back(std::make_shared<OpenPBRMaterial>(parameters, color, nullptr, false));
    }
    for (const auto &model : models) {
      Transform identity;
      hitable_list world, lights;
      world.add(std::make_shared<xy_rect>(-1, 1, -1, 1, 0, model, nullptr, nullptr, &identity,
                                          &identity, false));
      camera cam({0, 0, 4}, {0, 0, 0}, {0, 1, 0}, 30, 1, 0, 4, 0, 1, 1);
      WavefrontScene scene;
      WavefrontReport report;
      auto session = PrepareWavefront(world, lights, cam, 1, 0, scene, report);
      expect_true(bool(session));
      if (!session)
        continue;
      std::vector<WFBsdfProbe> probes;
      for (float view_sine : {0.f, .4f, .85f})
        for (float cosine : {-.8f, -.2f, .2f, .8f})
          for (float phi : {0.f, 1.1f, 2.8f}) {
            WFBsdfProbe probe;
            probe.material = scene.materials[0];
            probe.rho = {.3, .5, .7, 0};
            probe.transmission = {.4, .3, .2, 0};
            probe.view = {view_sine, 0, std::sqrt(1 - view_sine * view_sine), 0};
            probe.direction = {std::sqrt(1 - cosine * cosine) * std::cos(phi),
                               std::sqrt(1 - cosine * cosine) * std::sin(phi), cosine, 0};
            probe.context = {.7, .37, 1.2, 0};
            probe.random = {.17f + view_sine * .7f, .23f + (cosine + 1) * .2f, .11f + phi * .2f,
                            .79f};
            probes.push_back(probe);
          }
      auto results = ProbeMetalMaterials(probes, scene.material_data);
      float eval_error = 0, pdf_error = 0, weight_error = 0;
      random_gen rng(17);
      for (size_t i = 0; i < probes.size(); ++i) {
        const auto &p = probes[i];
        const auto &r = results[i];
        vec3f view(p.view.x, p.view.y, p.view.z), wi(p.direction.x, p.direction.y, p.direction.z);
        hit_record hit;
        hit.p = point3f(0);
        hit.normal = hit.geometric_normal = hit.physical_shading_normal = normal3f(0, 0, 1);
        hit.dpdu = vec3f(1, 0, 0);
        hit.dpdv = vec3f(0, 1, 0);
        hit.has_bump = false;
        hit.u = .3;
        hit.v = (p.context.y + 1) / 2;
        hit.mat_ptr = model.get();
        Ray ray(point3f(view * .7f), -view);
        point3f value;
        float pdf;
        vec3f sampled(r.direction.x, r.direction.y, r.direction.z);
        point3f sampled_value;
        float sampled_pdf = 0;
        if (const auto *pbr = dynamic_cast<const OpenPBRMaterial *>(model.get())) {
          auto b = pbr->Prepare(ray, hit, p.context.z);
          value = b.Evaluate(wi);
          pdf = b.Pdf(wi);
          if (r.direction.w > 0) {
            sampled_value = b.Evaluate(sampled);
            sampled_pdf = b.Pdf(sampled);
          }
        } else {
          scatter_record scattering;
          model->scatter(ray, hit, scattering, rng);
          value = model->f(ray, hit, wi);
          pdf = scattering.pdf_ptr->value(wi, rng);
          // The legacy reflection class can return negative f below the
          // reflecting hemisphere. Native NEE rejects that connection before
          // using its PDF; the GPU makes the same rejection in evaluation.
          if (dynamic_cast<const MicrofacetReflection *>(model.get()) && wi[2] * view[2] <= 0) {
            value = point3f(0);
            pdf = 0;
          }
          if (r.direction.w > 0) {
            sampled_value = model->f(ray, hit, sampled);
            sampled_pdf = scattering.pdf_ptr->value(sampled, rng);
          }
        }
        const float channels[] = {r.evaluated.x, r.evaluated.y, r.evaluated.z};
        const float weights[] = {r.weight.x, r.weight.y, r.weight.z};
        for (int c = 0; c < 3; ++c)
          eval_error = std::max(eval_error, probe_error(channels[c], value[c]));
        pdf_error = std::max(pdf_error, probe_error(r.evaluated.w, pdf));
        if (r.direction.w > 0 && r.weight.w == 0) {
          pdf_error = std::max(pdf_error, probe_error(r.direction.w, sampled_pdf));
          for (int c = 0; c < 3; ++c)
            weight_error = std::max(
                weight_error,
                probe_error(weights[c], sampled_pdf > 0 ? sampled_value[c] / sampled_pdf : 0));
        }
      }
      if (eval_error > .001 || pdf_error > .001 || weight_error > .001)
        throw std::runtime_error(
            model->GetName() + " probe errors: f=" + std::to_string(eval_error) +
            " pdf=" + std::to_string(pdf_error) + " weight=" + std::to_string(weight_error));
      expect_true(eval_error < .001);
      expect_true(pdf_error < .001);
      expect_true(weight_error < .001);
    }
  }
}
#endif
