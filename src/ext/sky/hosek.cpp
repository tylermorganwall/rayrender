/**
 Copyright (c) 2015 Eric Bruneton
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include "../../ext/sky/hosek.h"

#include <algorithm>

#include "../../ext/sky/ArHosekSkyModel.h"

namespace {

constexpr SolidAngle kSunSolidAngle = 6.87e-5 * sr;

// Copied from ArHosekSkyModel.cc.
const double originalSolarRadianceTable[] = {
  7500.0,
  12500.0,
  21127.5,
  26760.5,
  30663.7,
  27825.0,
  25503.8,
  25134.2,
  23212.1,
  21526.7,
  19870.8
};

// The Hosek sky radiance was computed uding 'originalSolarRadianceTable' as
// the solar radiance, but we want the result for SolarSpectrum(), used in all
// other models. Thus we correct the result for each wavelength by the ratio
// between the two spectra at this wavelength (including also the Sun solid
// angle to convert radiance to irradiance values).
double SolarSpectrumCorrectionFactor(Wavelength lambda) {
  double x = (lambda.to(nm) - 320.0) / 40.0;
  if (x < 0.0 || x > 10.0) {
    return 1.0;
  }
  int i = floor(x);
  double u = x - i;
  SpectralIrradiance new_irradiance = SolarSpectrum()(lambda);
  SpectralRadiance old_radiance = (originalSolarRadianceTable[i] * (1.0 - u) +
      originalSolarRadianceTable[std::min(i + 1, 10)] * u) *
      watt_per_square_meter_per_sr_per_nm;
  return (new_irradiance / (old_radiance * kSunSolidAngle))();
}

}  // anonymous namespace

Hosek::Hosek(double turbidity) : turbidity_(turbidity) {
  for (unsigned int i = 0; i < DimensionlessSpectrum::SIZE; ++i) {
    sky_model_state_[i] = NULL;
  }
}

Hosek::~Hosek() {
  for (unsigned int i = 0; i < DimensionlessSpectrum::SIZE; ++i) {
    if (sky_model_state_[i] != NULL) {
      arhosekskymodelstate_free(sky_model_state_[i]);
    }
  }
}

IrradianceSpectrum Hosek::GetSunIrradiance(Length altitude,
    Angle sun_zenith) const {
  MaybeInitSkyModelState(sun_zenith);
  IrradianceSpectrum result;
  for (unsigned int i = 0; i < result.size(); ++i) {
    Wavelength lambda = result.GetSample(i);
    if (lambda < 320.0 * nm || lambda > 720.0 * nm) {
      // arhosekskymodel_solar_radiance does not work for wavelengths outside
      // the 320-720nm range.
      result[i] = 0.0 * watt_per_square_meter_per_nm;
      continue;
    }
    double sun_and_sky_radiance = arhosekskymodel_solar_radiance(
        sky_model_state_[i], sun_zenith.to(rad), 0.0, lambda.to(nm));
    double sky_radiance = arhosekskymodel_radiance(
        sky_model_state_[i], sun_zenith.to(rad), 0.0, lambda.to(nm));
    double sun_radiance = (sun_and_sky_radiance - sky_radiance) *
        SolarSpectrumCorrectionFactor(lambda);
    result[i] =
        sun_radiance * watt_per_square_meter_per_sr_per_nm * kSunSolidAngle;
  }
  return result;
}

RadianceSpectrum Hosek::GetSkyRadiance(Length altitude, Angle sun_zenith,
    Angle view_zenith, Angle view_sun_azimuth) const {
  MaybeInitSkyModelState(sun_zenith);
  double gamma =
      GetViewSunAngle(sun_zenith, view_zenith, view_sun_azimuth).to(rad);
  RadianceSpectrum result;
  for (unsigned int i = 0; i < result.size(); ++i) {
    Wavelength lambda = result.GetSample(i);
    result[i] = arhosekskymodel_radiance(sky_model_state_[i],
        view_zenith.to(rad), gamma, lambda.to(nm)) *
        SolarSpectrumCorrectionFactor(lambda) *
        watt_per_square_meter_per_sr_per_nm;
  }
  return result;
}

void Hosek::MaybeInitSkyModelState(Angle sun_zenith) const {
  if (sky_model_state_[0] != NULL && sun_zenith == current_sun_zenith_) {
    return;
  }
  if (sky_model_state_[0] != NULL) {
    for (unsigned int i = 0; i < DimensionlessSpectrum::SIZE; ++i) {
      arhosekskymodelstate_free(sky_model_state_[i]);
    }
  }
  for (unsigned int i = 0; i < DimensionlessSpectrum::SIZE; ++i) {
    sky_model_state_[i] = arhosekskymodelstate_alloc_init(
        (pi / 2.0 - sun_zenith).to(rad), turbidity_, GroundAlbedo()[i]());
  }
  current_sun_zenith_ = sun_zenith;
}
