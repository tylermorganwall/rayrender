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
#ifndef PHYSICS_CIE_H_
#define PHYSICS_CIE_H_

#include "../../ext/sky/math/matrix.h"
#include "../../ext/sky/math/vector.h"
#include "../../ext/sky/units.h"

typedef dimensional::Vector3<Luminance> Color;

// The conversion factor between watts and lumens.
constexpr auto MaxLuminousEfficacy = 683.0 * lm / watt;

// The CIE color matching functions x_bar, y_bar and z_bar. See
// https://en.wikipedia.org/wiki/CIE_1931_color_space#Color_matching_functions.
const DimensionlessSpectrum& cie_x_bar_function();
const DimensionlessSpectrum& cie_y_bar_function();
const DimensionlessSpectrum& cie_z_bar_function();

// The conversion matrix from XYZ to linear sRGB color spaces.
constexpr dimensional::Matrix3<Number> XYZ_to_sRGB(
    +3.2404542, -1.5371385, -0.4985314,
    -0.9692660, +1.8760108, +0.0415560,
    +0.0556434, -0.2040259, +1.0572252);

// The CIE S0, S1 and S2 components, allowing to compute the relative spectral
// power distribution of a D65 illuminant with S = S0 + M1 S1 + M2 S2, where
// M1 and M2 depend on the CIE xy chromaticities.
// See https://en.wikipedia.org/wiki/Standard_illuminant for more details.
const DimensionlessSpectrum& S0_function();
const DimensionlessSpectrum& S1_function();
const DimensionlessSpectrum& S2_function();

// Converts the given radiance spectrum to a luminance value, using the relation
//   L = MaxLuminousEfficacy * integral of cie_y_bar_function * spectrum
Luminance GetLuminance(const RadianceSpectrum& spectrum);

// Converts the given spectrum to a RGB color, by first converting the spectrum
// to XYZ values with the x_bar, y_bar, z_bar functions, and then by converting
// from XYZ to linear sRGB.
Color GetSrgbColor(const RadianceSpectrum& spectrum);

// Same as above, but using only number_of_wavelengths wavelengths, between the
// given min and max.
Color GetSrgbColor(const RadianceSpectrum& spectrum, Wavelength min_wavelength,
    Wavelength max_wavelength, int number_of_wavelengths);

#endif  // PHYSICS_CIE_H_
