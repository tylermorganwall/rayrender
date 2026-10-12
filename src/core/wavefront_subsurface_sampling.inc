R"MSL(
// RGB marginal of the CPU pole proposal. The physical phase is isotropic;
// neither anisotropy nor extinction is changed by this sampling strategy.
float3 subsurface_guide_pdf(float3 pole, float cosine) {
  float3 k = max(pole, float3(1e-6f));
  float3 pdf = k / (4 * pi * atanh(k) * (1 - k * cosine));
  return select(pdf, float3(1 / (4 * pi)), pole < 1e-6f);
}
float subsurface_direction_pdf(float3 pole, float3 axis, float3 wi) {
  return .5f / (4 * pi) + dot(subsurface_guide_pdf(pole, dot(axis, wi)), float3(1.f / 6));
}
float3 subsurface_direction(float3 pole, float3 axis, float branch, float u, float v,
                            float3 ordinary) {
  if (branch < .5f)
    return ordinary;
  uint channel = min(2u, uint((branch - .5f) * 6));
  float k = pole[channel];
  // Metal lacks log1p/expm1. Use log(1+k)=atanh(k)+log(1-k*k)/2
  // and short Taylor forms near zero to avoid cancellation for white media.
  float k2 = k * k;
  float log_half = k < .01f ? -.5f * k2 * (1 + .5f * k2) : .5f * log(1 - k2);
  float x = (1 - 2 * u) * atanh(k) + log_half;
  float exponential_minus_one =
      abs(x) < .01f ? x * (1 + x * (.5f + x * (1.f / 6 + x * (1.f / 24 + x / 120)))) : exp(x) - 1;
  float cosine = k < 1e-6f ? 2 * u - 1 : -exponential_minus_one / k;
  cosine = clamp(cosine, -1.f, 1.f);
  Frame frame = normal_frame(axis);
  float sine = sqrt(max(0.f, 1 - cosine * cosine)), phi = 2 * pi * v;
  return cosine * frame.z + sine * (cos(phi) * frame.x + sin(phi) * frame.y);
}
// Condition the joint angular/free-flight proposal on the sampled direction.
// Its six components are RGB x (ordinary, pole); their weights sum to one.
void subsurface_flight_mixture(float3 pole, float3 axis, float3 wi, float3 rates,
                               thread float3 &guided_rates, thread float3 &ordinary_weights,
                               thread float3 &guided_weights) {
  float cosine = dot(axis, wi);
  float3 guide = subsurface_guide_pdf(pole, cosine);
  float direction_pdf = .5f / (4 * pi) + dot(guide, float3(1.f / 6));
  ordinary_weights = float3(1 / (24 * pi * direction_pdf));
  guided_weights = guide / (6 * direction_pdf);
  guided_rates = rates * (1 - pole * cosine);
}
)MSL"
