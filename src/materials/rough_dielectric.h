#ifndef RAYRENDER_ROUGH_DIELECTRIC_H
#define RAYRENDER_ROUGH_DIELECTRIC_H

#include "microfacetdist.h"

struct RoughDielectricEvaluation {
  Float f_cos = 0;
  Float pdf = 0;
  bool transmission = false;
};

// Both directions point away from the interface. View points toward the
// previous path vertex; outgoing points toward the next. eta is inside/outside.
// Return f * abs(cos(outgoing)) in radiance mode and its solid-angle density.
// Keep half-vector reconstruction, Fresnel, and the Jacobian identical in the
// material and PDF paths (PBRT v4 DielectricBxDF convention).
inline RoughDielectricEvaluation EvaluateRoughDielectric(
    const vec3f &view, const vec3f &outgoing, Float eta,
    const MicrofacetDistribution &distribution, point2f alphas) {
  RoughDielectricEvaluation result;
  if (view[2] == 0 || outgoing[2] == 0 || eta == 1) return result;
  const bool reflect = view[2] * outgoing[2] > 0;
  const Float ratio = view[2] > 0 ? 1 / eta : eta; // eta_i / eta_t
  vec3f normal = reflect ? view + outgoing : ratio * view + outgoing;
  if (!(normal.squared_length() > 0)) return result;
  normal = Faceforward(unit_vector(normal), normal3f(0, 0, 1));
  const Float view_dot = dot(view, normal), out_dot = dot(outgoing, normal);
  // A sampled microfacet must face both directions on their respective sides.
  if (view_dot * view[2] <= 0 || out_dot * outgoing[2] <= 0) return result;
  const Float fresnel = FrDielectric(std::abs(view_dot), ratio);
  const Float D = distribution.D(normal, alphas);
  const Float G = distribution.SmithG(view, outgoing, alphas);
  const Float normal_pdf = distribution.VisibleNormalPdf(view, normal, alphas);
  if (reflect) {
    result.f_cos = fresnel * D * G / (4 * AbsCosTheta(view));
    result.pdf = fresnel * normal_pdf / (4 * std::abs(view_dot));
  } else {
    const Float denom = Sqr(out_dot + ratio * view_dot);
    if (!(denom > 0)) return result;
    result.transmission = true;
    result.pdf = (1 - fresnel) * normal_pdf * std::abs(out_dot) / denom;
    // Refraction changes radiance by (eta_i/eta_t)^2. This is separate from
    // the microfacet-normal to outgoing-solid-angle Jacobian in the PDF.
    result.f_cos = (1 - fresnel) * D * G * std::abs(view_dot * out_dot) *
      Sqr(ratio) / (AbsCosTheta(view) * denom);
  }
  return result;
}

#endif
