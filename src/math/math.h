#ifndef MATHH
#define MATHH

#include "float.h"
#include <span>

// Spline Interpolation Function Definitions
bool CatmullRomWeights(std::span<const Float> nodes, Float x, int *offset,
                       std::span<Float> weights);

Float CatmullRom(std::span<const Float> nodes, std::span<const Float> f, Float x);

Float InvertCatmullRom(std::span<const Float> nodes, std::span<const Float> f,
                       Float u);

Float IntegrateCatmullRom(std::span<const Float> nodes, std::span<const Float> f,
                          std::span<Float> cdf);

#endif