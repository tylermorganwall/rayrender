#include "math.h"
#include "mathinline.h"

// Spline Interpolation Function Definitions
bool CatmullRomWeights(std::span<const Float> nodes, Float x, int *offset,
                       std::span<Float> weights) {
    //CHECK_GE(weights.size(), 4);
    // Return _false_ if _x_ is out of bounds
    if (!(x >= nodes.front() && x <= nodes.back()))
        return false;

    // Search for the interval _idx_ containing _x_
    int idx = FindInterval(nodes.size(), [&](int i) { return nodes[i] <= x; });
    *offset = idx - 1;
    Float x0 = nodes[idx], x1 = nodes[idx + 1];

    // Compute the $t$ parameter and powers
    Float t = (x - x0) / (x1 - x0), t2 = t * t, t3 = t2 * t;

    // Compute initial node weights $w_1$ and $w_2$
    weights[1] = 2 * t3 - 3 * t2 + 1;
    weights[2] = -2 * t3 + 3 * t2;

    // Compute first node weight $w_0$
    if (idx > 0) {
        Float w0 = (t3 - 2 * t2 + t) * (x1 - x0) / (x1 - nodes[idx - 1]);
        weights[0] = -w0;
        weights[2] += w0;
    } else {
        Float w0 = t3 - 2 * t2 + t;
        weights[0] = 0;
        weights[1] -= w0;
        weights[2] += w0;
    }

    // Compute last node weight $w_3$
    if (idx + 2 < nodes.size()) {
        Float w3 = (t3 - t2) * (x1 - x0) / (nodes[idx + 2] - x0);
        weights[1] -= w3;
        weights[3] = w3;
    } else {
        Float w3 = t3 - t2;
        weights[1] -= w3;
        weights[2] += w3;
        weights[3] = 0;
    }

    return true;
}

Float CatmullRom(std::span<const Float> nodes, std::span<const Float> f, Float x) {
    // CHECK_EQ(nodes.size(), f.size());
    if (!(x >= nodes.front() && x <= nodes.back()))
        return 0;
    int idx = FindInterval(nodes.size(), [&](int i) { return nodes[i] <= x; });
    Float x0 = nodes[idx], x1 = nodes[idx + 1];
    Float f0 = f[idx], f1 = f[idx + 1];
    Float width = x1 - x0;
    Float d0, d1;
    if (idx > 0)
        d0 = width * (f1 - f[idx - 1]) / (x1 - nodes[idx - 1]);
    else
        d0 = f1 - f0;

    if (idx + 2 < nodes.size())
        d1 = width * (f[idx + 2] - f0) / (nodes[idx + 2] - x0);
    else
        d1 = f1 - f0;

    Float t = (x - x0) / (x1 - x0), t2 = t * t, t3 = t2 * t;
    return (2 * t3 - 3 * t2 + 1) * f0 + (-2 * t3 + 3 * t2) * f1 + (t3 - 2 * t2 + t) * d0 +
           (t3 - t2) * d1;
}

Float InvertCatmullRom(std::span<const Float> nodes, std::span<const Float> f,
                       Float u) {
    // Stop when _u_ is out of bounds
    if (!(u > f.front()))
        return nodes.front();
    else if (!(u < f.back()))
        return nodes.back();

    // Map _u_ to a spline interval by inverting _f_
    int i = FindInterval(f.size(), [&](int i) { return f[i] <= u; });

    // Look up $x_i$ and function values of spline segment _i_
    Float x0 = nodes[i], x1 = nodes[i + 1];
    Float f0 = f[i], f1 = f[i + 1];
    Float width = x1 - x0;

    // Approximate derivatives using finite differences
    Float d0 = (i > 0) ? width * (f1 - f[i - 1]) / (x1 - nodes[i - 1]) : (f1 - f0);
    Float d1 = (i + 2 < nodes.size()) ? width * (f[i + 2] - f0) / (nodes[i + 2] - x0)
                                      : (f1 - f0);

    // Invert the spline interpolant using Newton-Bisection
    auto eval = [&](Float t) -> std::pair<Float, Float> {
        // Compute powers of _t_
        Float t2 = t * t, t3 = t2 * t;

        // Set _Fhat_ using Equation (\ref{eq:cubicspline-as-basisfunctions})
        Float Fhat = (2 * t3 - 3 * t2 + 1) * f0 + (-2 * t3 + 3 * t2) * f1 +
                     (t3 - 2 * t2 + t) * d0 + (t3 - t2) * d1;

        // Set _fhat_ using Equation (\ref{eq:cubicspline-derivative})
        Float fhat = (6 * t2 - 6 * t) * f0 + (-6 * t2 + 6 * t) * f1 +
                     (3 * t2 - 4 * t + 1) * d0 + (3 * t2 - 2 * t) * d1;

        return {Fhat - u, fhat};
    };
    Float t = NewtonBisection(0, 1, eval);
    return x0 + t * width;
}

Float IntegrateCatmullRom(std::span<const Float> nodes, std::span<const Float> f,
                          std::span<Float> cdf) {
    // CHECK_EQ(nodes.size(), f.size());
    Float sum = 0;
    cdf[0] = 0;
    for (int i = 0; i < nodes.size() - 1; ++i) {
        // Look up $x_i$ and function values of spline segment _i_
        Float x0 = nodes[i], x1 = nodes[i + 1];
        Float f0 = f[i], f1 = f[i + 1];
        Float width = x1 - x0;

        // Approximate derivatives using finite differences
        Float d0 = (i > 0) ? width * (f1 - f[i - 1]) / (x1 - nodes[i - 1]) : (f1 - f0);
        Float d1 = (i + 2 < nodes.size()) ? width * (f[i + 2] - f0) / (nodes[i + 2] - x0)
                                          : (f1 - f0);

        // Keep a running sum and build a cumulative distribution function
        sum += width * ((f0 + f1) / 2 + (d0 - d1) / 12);
        cdf[i + 1] = sum;
    }
    return sum;
}