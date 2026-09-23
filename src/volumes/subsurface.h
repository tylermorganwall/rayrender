#ifndef RAYRENDER_SUBSURFACE_H
#define RAYRENDER_SUBSURFACE_H
#include "medium.h"
#include <algorithm>
#include <cmath>

// Keep a collision within the incoming side of an uncertain geometric endpoint.
// Its sampled world distance and density are unchanged.
point3f SubsurfaceCollisionPoint(const point3f &p, const Ray &ray, double t,
                                const hit_record &endpoint, bool &adjusted);

// Sampling data only: none of these quantities replace physical coefficients.
// d'Eon and Krivanek (2020), section 6.6. Original implementation.
struct SubsurfaceProposal {
  std::array<double, 3> extinction{}, pole{};
  HGPhaseFunction phase;
  vec3f wo{0, 0, 1}, axis{0, 0, 1};
  bool guided = false;
  static double Pole(double albedo);
  double GuidePdf(int c, const vec3f &wi) const;
  double DirectionPdf(int c, const vec3f &wi) const;
  vec3f SampleDirection(int hero, double strategy, double u, double v) const;
  double GuideRate(int c, const vec3f &wi) const;
  // Conditional on the sampled direction. Includes both strategy masses.
  double LogDistancePdf(int c, const vec3f &wi, double t, bool collision) const;
  double SampleDistance(int hero, const vec3f &wi, double strategy, double u) const;
  static SubsurfaceProposal Ordinary(const Medium &medium);
};

struct SubsurfaceBoundarySample {
  vec3f wi{0};
  double weight = 0, pdf = 0;
  bool transmission = false, specular = false;
};
// Neutral single-scattering GGX dielectric; alpha = roughness^2. Smith masking.
// Directions point away from the interface. eta = transmitted / incident IOR.
// Evaluate returns f * abs(cos(theta_i)) in radiance transport mode.
class SubsurfaceBoundaryBSDF {
public:
  SubsurfaceBoundaryBSDF(const vec3f &incoming, const normal3f &normal,
                         double eta, double roughness);
  double Evaluate(const vec3f &wi) const;
  double Pdf(const vec3f &wi) const;
  SubsurfaceBoundarySample Sample(double branch, double u, double v) const;
  bool IsSpecular() const { return alpha == 0 || eta == 1; }
  double Eta() const { return eta; }
private:
  onb frame;
  vec3f wo;
  double eta, alpha;
  double D(const vec3f &m) const;
  double G1(const vec3f &w) const;
  double EvaluateLocal(const vec3f &wi, bool density) const;
};
#endif
