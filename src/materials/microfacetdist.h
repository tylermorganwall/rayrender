#ifndef MICROFACETDIST
#define MICROFACETDIST

#include "../math/mathinline.h"
#include "../math/rng.h"
#include "../math/vec2.h"
#include "../materials/texturegraph.h"

class MicrofacetDistribution {
public:
  virtual ~MicrofacetDistribution() {}
  std::shared_ptr<const TextureNode> roughness_graph;
  std::shared_ptr<const TextureNode> roughness_graph_v;
  bool graph_is_alpha = false;
  point2f Resolve(const hit_record& hit) const;
  virtual Float D(const vec3f &wh, point2f alphas) const = 0;
  virtual vec3f Sample_wh(const vec3f &wi, Float u1, Float u2, point2f alphas) const = 0;
  // The visible-normal sampler is conditioned on the direction toward the
  // previous vertex. Its density uses Smith G1, not the two-direction G term.
  virtual Float Lambda(const vec3f &w, point2f alphas) const = 0;
  Float VisibleNormalPdf(const vec3f &view, const vec3f &normal, point2f alphas) const {
    if (view[2] == 0 || dot(view, normal) * view[2] <= 0) return 0;
    return D(normal, alphas) * AbsDot(view, normal) /
      (AbsCosTheta(view) * (1 + Lambda(view, alphas)));
  }
  Float SmithG(const vec3f &view, const vec3f &outgoing, point2f alphas) const {
    return 1 / (1 + Lambda(view, alphas) + Lambda(outgoing, alphas));
  }
  Float Pdf(const vec3f &wo, const vec3f &wi, const vec3f &wh, point2f alphas) const {
    return D(wh, alphas) * G(wo, wi, wh) * AbsDot(wo, wh) / AbsCosTheta(wo);
  }
  virtual Float D(const vec3f &wh) const = 0;
  virtual Float Lambda(const vec3f &w) const = 0;
  virtual Float D(const vec3f &wh, Float u, Float v) const = 0;
  virtual Float Lambda(const vec3f &w, Float u, Float v) const = 0;
  Float G1(const vec3f &w, Float u, Float v) const {
    return(1 / (1 + Lambda(w, u, v)));
  }
  Float G(const vec3f &wo, const vec3f &wi, const vec3f &wh) const {
    Float NdotWh = AbsCosTheta(wh);
    Float NdotWo = AbsCosTheta(wo);
    Float NdotWi = AbsCosTheta(wi);
    Float WOdotWh = std::fabs(dot(wo, wh));
    return(std::fmin(1.f, std::fmin((2.0f * NdotWh * NdotWo / WOdotWh), (2.0f * NdotWh * NdotWi / WOdotWh))));
  }
  virtual Float GetAlpha(Float u, Float v) const = 0;
  virtual point2f GetAlphas(Float u, Float v) const = 0;
  virtual vec3f Sample_wh(const vec3f &wi, const Float u1, const Float u2) const = 0;
  virtual vec3f Sample_wh(const vec3f &wi, const Float u1, const Float u2, Float u, Float v) const = 0;
  
  Float Pdf(const vec3f &wo,const vec3f &wi, const vec3f &wh, Float u, Float v) {
    return(D(wh, u, v) * G(wo, wi, wh) * AbsDot(wo, wh) / AbsCosTheta(wo));
  }
protected:
  MicrofacetDistribution(bool sampleVisibleArea) 
      : sampleVisibleArea(sampleVisibleArea) {}
  const bool sampleVisibleArea;
};

class BeckmannDistribution : public MicrofacetDistribution {
public:
  static Float RoughnessToAlpha(Float roughness) {
    roughness = std::fmax(roughness, (Float)0.0001550155);
    Float x = std::log(roughness);
    return(1.62142f + 0.819955f * x + 0.1734f * x * x +
           0.0171201f * x * x * x + 0.000640711f * x * x * x * x );
  }
  static point2f RoughnessToAlpha(point2f roughness) {
    Float x = std::log(std::fmax(roughness.xy.x, (Float)0.0001550155));
    Float y = std::log(std::fmax(roughness.xy.y, (Float)0.0001550155));
    
    return(point2f(1.62142f + 0.819955f * x + 0.1734f * x * x +
           0.0171201f * x * x * x + 0.000640711f * x * x * x * x,
           1.62142f + 0.819955f * y + 0.1734f * y * y +
             0.0171201f * y * y * y + 0.000640711f * y * y * y * y ));
  }
  BeckmannDistribution(Float alphax_, Float alphay_, std::shared_ptr<roughness_texture> roughness,
                       bool has_roughness, bool samplevis = true)
    : MicrofacetDistribution(samplevis), roughness(roughness), has_roughness(has_roughness) {
    alphax_constant = RoughnessToAlpha(alphax_);
    alphay_constant = RoughnessToAlpha(alphay_);
    alphax_constant *= alphax_constant;
    alphay_constant *= alphay_constant;
  }
  ~BeckmannDistribution() override {}
  Float D(const vec3f &wh) const override;
  Float D(const vec3f &wh, point2f alphas) const override;
  Float Lambda(const vec3f &w, point2f alphas) const override;
  vec3f Sample_wh(const vec3f &wi, Float u1, Float u2, point2f alphas) const override;
  Float GetAlpha(Float u, Float v) const override;
  point2f GetAlphas(Float u, Float v) const override;
  vec3f Sample_wh(const vec3f &wi, const Float u1, const Float u2) const override;
  Float D(const vec3f &wh, Float u, Float v) const override;
  vec3f Sample_wh(const vec3f &wi, const Float u1, const Float u2, Float u, Float v) const override;
  Float Lambda(const vec3f &w) const override;
  Float Lambda(const vec3f &w, Float u, Float v) const override;
private:
  
  Float alphax_constant, alphay_constant;
  std::shared_ptr<roughness_texture> roughness;
  bool has_roughness;
};


class TrowbridgeReitzDistribution : public MicrofacetDistribution {
public:
  TrowbridgeReitzDistribution(const Float alphax_, const Float alphay_, std::shared_ptr<roughness_texture> roughness,
                              bool has_roughness, bool samplevis = true)
    : MicrofacetDistribution(samplevis), roughness(roughness), has_roughness(has_roughness) {
    alphax_constant = RoughnessToAlpha(alphax_);
    alphay_constant = RoughnessToAlpha(alphay_);
    alphax_constant *= alphax_constant;
    alphay_constant *= alphay_constant;
  }
  ~TrowbridgeReitzDistribution() override {}
  static Float RoughnessToAlpha(Float roughness) {
    roughness = std::fmax(roughness, (Float)0.0001550155);
    Float x = std::log(roughness);
    return(1.62142f + 0.819955f * x + 0.1734f * x * x +
           0.0171201f * x * x * x + 0.000640711f * x * x * x * x );
  }
  Float D(const vec3f &w) const override;
  Float D(const vec3f &wh, point2f alphas) const override;
  Float Lambda(const vec3f &w, point2f alphas) const override;
  vec3f Sample_wh(const vec3f &wi, Float u1, Float u2, point2f alphas) const override;
  Float GetAlpha(Float u, Float v) const override;
  point2f GetAlphas(Float u, Float v) const override;
  vec3f Sample_wh(const vec3f &wi, const Float u1, const Float u2) const override;
  
  Float D(const vec3f &wh, Float u, Float v) const override;
  vec3f Sample_wh(const vec3f &wi, const Float u1, const Float u2, Float u, Float v) const override;
  Float Lambda(const vec3f &w) const override;
  Float Lambda(const vec3f &w, Float u, Float v) const override;

private:
  Float alphax_constant, alphay_constant;
  std::shared_ptr<roughness_texture> roughness;
  bool has_roughness;
};

#endif
