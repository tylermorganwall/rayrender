#ifndef NOISE_H
#define NOISE_H

#include "float.h"
#include "vectypes.h"

Float Noise(Float x, Float y = .5f, Float z = .5f);

Float Noise(point3f p);

vec3f DNoise(point3f p);

Float FBm(point3f p, vec3f dpdx, vec3f dpdy, Float omega, int octaves);

Float Turbulence(point3f p, vec3f dpdx, vec3f dpdy, Float omega, int octaves);


#endif  // NOISE_H