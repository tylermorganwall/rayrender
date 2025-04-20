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
#include "../../ext/sky/cie.h"

#include <algorithm>
#include <vector>

namespace {

struct Tables {
  std::vector<Number> x_bar_samples;
  std::vector<Number> y_bar_samples;
  std::vector<Number> z_bar_samples;
  std::vector<Number> S0_samples;
  std::vector<Number> S1_samples;
  std::vector<Number> S2_samples;

  Tables() {
    // Values from "CIE (1931) 2-deg color matching functions", see
    // "http://web.archive.org/web/20081228084047/
    //    http://www.cvrl.org/database/data/cmfs/ciexyz31.txt".
    Add(0.000000000000, 0.000000000000, 0.000000000000);
    Add(0.000129900000, 0.000003917000, 0.000606100000);
    Add(0.000232100000, 0.000006965000, 0.001086000000);
    Add(0.000414900000, 0.000012390000, 0.001946000000);
    Add(0.000741600000, 0.000022020000, 0.003486000000);
    Add(0.001368000000, 0.000039000000, 0.006450001000);
    Add(0.002236000000, 0.000064000000, 0.010549990000);
    Add(0.004243000000, 0.000120000000, 0.020050010000);
    Add(0.007650000000, 0.000217000000, 0.036210000000);
    Add(0.014310000000, 0.000396000000, 0.067850010000);
    Add(0.023190000000, 0.000640000000, 0.110200000000);
    Add(0.043510000000, 0.001210000000, 0.207400000000);
    Add(0.077630000000, 0.002180000000, 0.371300000000);
    Add(0.134380000000, 0.004000000000, 0.645600000000);
    Add(0.214770000000, 0.007300000000, 1.039050100000);
    Add(0.283900000000, 0.011600000000, 1.385600000000);
    Add(0.328500000000, 0.016840000000, 1.622960000000);
    Add(0.348280000000, 0.023000000000, 1.747060000000);
    Add(0.348060000000, 0.029800000000, 1.782600000000);
    Add(0.336200000000, 0.038000000000, 1.772110000000);
    Add(0.318700000000, 0.048000000000, 1.744100000000);
    Add(0.290800000000, 0.060000000000, 1.669200000000);
    Add(0.251100000000, 0.073900000000, 1.528100000000);
    Add(0.195360000000, 0.090980000000, 1.287640000000);
    Add(0.142100000000, 0.112600000000, 1.041900000000);
    Add(0.095640000000, 0.139020000000, 0.812950100000);
    Add(0.057950010000, 0.169300000000, 0.616200000000);
    Add(0.032010000000, 0.208020000000, 0.465180000000);
    Add(0.014700000000, 0.258600000000, 0.353300000000);
    Add(0.004900000000, 0.323000000000, 0.272000000000);
    Add(0.002400000000, 0.407300000000, 0.212300000000);
    Add(0.009300000000, 0.503000000000, 0.158200000000);
    Add(0.029100000000, 0.608200000000, 0.111700000000);
    Add(0.063270000000, 0.710000000000, 0.078249990000);
    Add(0.109600000000, 0.793200000000, 0.057250010000);
    Add(0.165500000000, 0.862000000000, 0.042160000000);
    Add(0.225749900000, 0.914850100000, 0.029840000000);
    Add(0.290400000000, 0.954000000000, 0.020300000000);
    Add(0.359700000000, 0.980300000000, 0.013400000000);
    Add(0.433449900000, 0.994950100000, 0.008749999000);
    Add(0.512050100000, 1.000000000000, 0.005749999000);
    Add(0.594500000000, 0.995000000000, 0.003900000000);
    Add(0.678400000000, 0.978600000000, 0.002749999000);
    Add(0.762100000000, 0.952000000000, 0.002100000000);
    Add(0.842500000000, 0.915400000000, 0.001800000000);
    Add(0.916300000000, 0.870000000000, 0.001650001000);
    Add(0.978600000000, 0.816300000000, 0.001400000000);
    Add(1.026300000000, 0.757000000000, 0.001100000000);
    Add(1.056700000000, 0.694900000000, 0.001000000000);
    Add(1.062200000000, 0.631000000000, 0.000800000000);
    Add(1.045600000000, 0.566800000000, 0.000600000000);
    Add(1.002600000000, 0.503000000000, 0.000340000000);
    Add(0.938400000000, 0.441200000000, 0.000240000000);
    Add(0.854449900000, 0.381000000000, 0.000190000000);
    Add(0.751400000000, 0.321000000000, 0.000100000000);
    Add(0.642400000000, 0.265000000000, 0.000049999990);
    Add(0.541900000000, 0.217000000000, 0.000030000000);
    Add(0.447900000000, 0.175000000000, 0.000020000000);
    Add(0.360800000000, 0.138200000000, 0.000010000000);
    Add(0.283500000000, 0.107000000000, 0.000000000000);
    Add(0.218700000000, 0.081600000000, 0.000000000000);
    Add(0.164900000000, 0.061000000000, 0.000000000000);
    Add(0.121200000000, 0.044580000000, 0.000000000000);
    Add(0.087400000000, 0.032000000000, 0.000000000000);
    Add(0.063600000000, 0.023200000000, 0.000000000000);
    Add(0.046770000000, 0.017000000000, 0.000000000000);
    Add(0.032900000000, 0.011920000000, 0.000000000000);
    Add(0.022700000000, 0.008210000000, 0.000000000000);
    Add(0.015840000000, 0.005723000000, 0.000000000000);
    Add(0.011359160000, 0.004102000000, 0.000000000000);
    Add(0.008110916000, 0.002929000000, 0.000000000000);
    Add(0.005790346000, 0.002091000000, 0.000000000000);
    Add(0.004109457000, 0.001484000000, 0.000000000000);
    Add(0.002899327000, 0.001047000000, 0.000000000000);
    Add(0.002049190000, 0.000740000000, 0.000000000000);
    Add(0.001439971000, 0.000520000000, 0.000000000000);
    Add(0.000999949300, 0.000361100000, 0.000000000000);
    Add(0.000690078600, 0.000249200000, 0.000000000000);
    Add(0.000476021300, 0.000171900000, 0.000000000000);
    Add(0.000332301100, 0.000120000000, 0.000000000000);
    Add(0.000234826100, 0.000084800000, 0.000000000000);
    Add(0.000166150500, 0.000060000000, 0.000000000000);
    Add(0.000117413000, 0.000042400000, 0.000000000000);
    Add(0.000083075270, 0.000030000000, 0.000000000000);
    Add(0.000058706520, 0.000021200000, 0.000000000000);
    Add(0.000041509940, 0.000014990000, 0.000000000000);
    Add(0.000029353260, 0.000010600000, 0.000000000000);
    Add(0.000020673830, 0.000007465700, 0.000000000000);
    Add(0.000014559770, 0.000005257800, 0.000000000000);
    Add(0.000010253980, 0.000003702900, 0.000000000000);
    Add(0.000007221456, 0.000002607800, 0.000000000000);
    Add(0.000005085868, 0.000001836600, 0.000000000000);
    Add(0.000003581652, 0.000001293400, 0.000000000000);
    Add(0.000002522525, 0.000000910930, 0.000000000000);
    Add(0.000001776509, 0.000000641530, 0.000000000000);
    Add(0.000001251141, 0.000000451810, 0.000000000000);
    Add(0.000000000000, 0.000000000000, 0.000000000000);

    // Values from Table T.2 in CIE TR "Colorimetry" 15:2004 Third Edition
    // (ISBN 3 901 906 33 9).
    AddS(+61.5,  38.0,  5.3);
    AddS(+68.8,  42.4,  6.1);
    AddS(+63.4,  38.5,  3.0);
    AddS(+65.8,  35.0,  1.2);
    AddS(+94.8,  43.4, -1.1);
    AddS(104.8,  46.3, -0.5);
    AddS(105.9,  43.9, -0.7);
    AddS(+96.8,  37.1, -1.2);
    AddS(113.9,  36.7, -2.6);
    AddS(125.6,  35.9, -2.9);
    AddS(125.5,  32.6, -2.8);
    AddS(121.3,  27.9, -2.6);
    AddS(121.3,  24.3, -2.6);
    AddS(113.5,  20.1, -1.8);
    AddS(113.1,  16.2, -1.5);
    AddS(110.8,  13.2, -1.3);
    AddS(106.5,   8.6, -1.2);
    AddS(108.8,   6.1, -1.0);
    AddS(105.3,   4.2, -0.5);
    AddS(104.4,   1.9, -0.3);
    AddS(100.0,   0.0,  0.0);
    AddS(+96.0,  -1.6,  0.2);
    AddS(+95.1,  -3.5,  0.5);
    AddS(+89.1,  -3.5,  2.1);
    AddS(+90.5,  -5.8,  3.2);
    AddS(+90.3,  -7.2,  4.1);
    AddS(+88.4,  -8.6,  4.7);
    AddS(+84.0,  -9.5,  5.1);
    AddS(+85.1, -10.9,  6.7);
    AddS(+81.9, -10.7,  7.3);
    AddS(+82.6, -12.0,  8.6);
    AddS(+84.9, -14.0,  9.8);
    AddS(+81.3, -13.6, 10.2);
    AddS(+71.9, -12.0,  8.3);
    AddS(+74.3, -13.3,  9.6);
    AddS(+76.4, -12.9,  8.5);
    AddS(+63.3, -10.6,  7.0);
    AddS(+71.7, -11.6,  7.6);
    AddS(+77.0, -12.2,  8.0);
    AddS(+65.2, -10.2,  6.7);
    AddS(+47.7,  -7.8,  5.2);
    AddS(+68.6, -11.2,  7.4);
    AddS(+65.0, -10.4,  6.8);
    AddS(+66.0, -10.6,  7.0);
    AddS(+61.0,  -9.7,  6.4);
    AddS(+53.3,  -8.3,  5.5);
    AddS(+58.9,  -9.3,  6.1);
    AddS(+61.9,  -9.8,  6.5);
  }

 private:
  void Add(double x_bar_sample, double y_bar_sample, double z_bar_sample) {
    x_bar_samples.push_back(x_bar_sample);
    y_bar_samples.push_back(y_bar_sample);
    z_bar_samples.push_back(z_bar_sample);
  }

  void AddS(double s0_sample, double s1_sample, double s2_sample) {
    S0_samples.push_back(s0_sample);
    S1_samples.push_back(s1_sample);
    S2_samples.push_back(s2_sample);
  }
};

const Tables tables = Tables();
const DimensionlessSpectrum x_bar(355.0 * nm, 835.0 * nm, tables.x_bar_samples);
const DimensionlessSpectrum y_bar(355.0 * nm, 835.0 * nm, tables.y_bar_samples);
const DimensionlessSpectrum z_bar(355.0 * nm, 835.0 * nm, tables.z_bar_samples);
const DimensionlessSpectrum S0(360.0 * nm, 830.0 * nm, tables.S0_samples);
const DimensionlessSpectrum S1(360.0 * nm, 830.0 * nm, tables.S1_samples);
const DimensionlessSpectrum S2(360.0 * nm, 830.0 * nm, tables.S2_samples);

}  // anonymous namespace

const DimensionlessSpectrum& cie_x_bar_function() { return x_bar; }
const DimensionlessSpectrum& cie_y_bar_function() { return y_bar; }
const DimensionlessSpectrum& cie_z_bar_function() { return z_bar; }

const DimensionlessSpectrum& S0_function() { return S0; }
const DimensionlessSpectrum& S1_function() { return S1; }
const DimensionlessSpectrum& S2_function() { return S2; }

Luminance GetLuminance(const RadianceSpectrum& spectrum) {
  return MaxLuminousEfficacy * Integral(y_bar * spectrum);
}

Color GetSrgbColor(const RadianceSpectrum& spectrum) {
  Color RGB = XYZ_to_sRGB * Color(
      MaxLuminousEfficacy * Integral(cie_x_bar_function() * spectrum),
      MaxLuminousEfficacy * Integral(cie_y_bar_function() * spectrum),
      MaxLuminousEfficacy * Integral(cie_z_bar_function() * spectrum));
  return Color(
      std::max(RGB.x, 0.0 * cd_per_square_meter),
      std::max(RGB.y, 0.0 * cd_per_square_meter),
      std::max(RGB.z, 0.0 * cd_per_square_meter));
}

Color GetSrgbColor(const RadianceSpectrum& spectrum, Wavelength min_wavelength,
    Wavelength max_wavelength, int number_of_wavelengths) {
  Color RGB = XYZ_to_sRGB * Color(
      MaxLuminousEfficacy * Integral(cie_x_bar_function() * spectrum,
          min_wavelength, max_wavelength, number_of_wavelengths),
      MaxLuminousEfficacy * Integral(cie_y_bar_function() * spectrum,
          min_wavelength, max_wavelength, number_of_wavelengths),
      MaxLuminousEfficacy * Integral(cie_z_bar_function() * spectrum,
          min_wavelength, max_wavelength, number_of_wavelengths));
  return Color(
      std::max(RGB.x, 0.0 * cd_per_square_meter),
      std::max(RGB.y, 0.0 * cd_per_square_meter),
      std::max(RGB.z, 0.0 * cd_per_square_meter));
}
