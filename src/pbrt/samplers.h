// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#ifndef PBRT_SAMPLERS_H
#define PBRT_SAMPLERS_H

#include <pbrt/pbrt.h>

#include <pbrt/base/sampler.h>
#include <pbrt/filters.h>
#include <pbrt/options.h>
#include <pbrt/util/bluenoise.h>
#include <pbrt/util/check.h>
#include <pbrt/util/error.h>
#include <pbrt/util/hash.h>
#include <pbrt/util/lowdiscrepancy.h>
#include <pbrt/util/math.h>
#include <pbrt/util/pmj02tables.h>
#include <pbrt/util/primes.h>
#include <pbrt/util/pstd.h>
#include <pbrt/util/rng.h>
#include <pbrt/util/vecmath.h>

#include <algorithm>
#include <limits>
#include <memory>
#include <string>

namespace pbrt {

// HaltonSampler Definition
class HaltonSampler {
  public:
    // HaltonSampler Public Methods
    HaltonSampler(int samplesPerPixel, Point2i fullResolution,
                  RandomizeStrategy randomize = RandomizeStrategy::PermuteDigits,
                  int seed = 0, Allocator alloc = {});

    PBRT_CPU_GPU
    static constexpr const char *Name() { return "HaltonSampler"; }
    static HaltonSampler *Create(const ParameterDictionary &parameters,
                                 Point2i fullResolution, const FileLoc *loc,
                                 Allocator alloc);

    PBRT_CPU_GPU
    int SamplesPerPixel() const { return samplesPerPixel; }

    PBRT_CPU_GPU
    RandomizeStrategy GetRandomizeStrategy() const { return randomize; }

    PBRT_CPU_GPU
    void StartPixelSample(Point2i p, int sampleIndex, int dim) {
        haltonIndex = 0;
        int sampleStride = baseScales[0] * baseScales[1];
        // Compute Halton sample index for first sample in pixel _p_
        if (sampleStride > 1) {
            Point2i pm(Mod(p[0], MaxHaltonResolution), Mod(p[1], MaxHaltonResolution));
            for (int i = 0; i < 2; ++i) {
                uint64_t dimOffset =
                    (i == 0) ? InverseRadicalInverse(pm[i], 2, baseExponents[i])
                             : InverseRadicalInverse(pm[i], 3, baseExponents[i]);
                haltonIndex +=
                    dimOffset * (sampleStride / baseScales[i]) * multInverse[i];
            }
            haltonIndex %= sampleStride;
        }

        haltonIndex += sampleIndex * sampleStride;
        dimension = std::max(2, dim);
    }

    PBRT_CPU_GPU
    Float Get1D() {
        if (dimension >= PrimeTableSize)
            dimension = 2;
        return SampleDimension(dimension++);
    }

    PBRT_CPU_GPU
    Point2f Get2D() {
        if (dimension + 1 >= PrimeTableSize)
            dimension = 2;
        int dim = dimension;
        dimension += 2;
        return {SampleDimension(dim), SampleDimension(dim + 1)};
    }

    PBRT_CPU_GPU
    Point2f GetPixel2D() {
        return {RadicalInverse(0, haltonIndex >> baseExponents[0]),
                RadicalInverse(1, haltonIndex / baseScales[1])};
    }

    Sampler Clone(Allocator alloc);
    std::string ToString() const;

  private:
    // HaltonSampler Private Methods
    static uint64_t multiplicativeInverse(int64_t a, int64_t n) {
        int64_t x, y;
        extendedGCD(a, n, &x, &y);
        return Mod(x, n);
    }

    static void extendedGCD(uint64_t a, uint64_t b, int64_t *x, int64_t *y) {
        if (b == 0) {
            *x = 1;
            *y = 0;
            return;
        }
        int64_t d = a / b, xp, yp;
        extendedGCD(b, a % b, &xp, &yp);
        *x = yp;
        *y = xp - (d * yp);
    }

    PBRT_CPU_GPU
    Float SampleDimension(int dimension) const {
        if (randomize == RandomizeStrategy::None)
            return RadicalInverse(dimension, haltonIndex);
        else if (randomize == RandomizeStrategy::PermuteDigits)
            return ScrambledRadicalInverse(dimension, haltonIndex,
                                           (*digitPermutations)[dimension]);
        else {
            DCHECK_EQ(randomize, RandomizeStrategy::Owen);
            return OwenScrambledRadicalInverse(dimension, haltonIndex,
                                               MixBits(1 + (dimension << 4)));
        }
    }

    // HaltonSampler Private Members
    int samplesPerPixel;
    RandomizeStrategy randomize;
    pstd::vector<DigitPermutation> *digitPermutations = nullptr;
    static constexpr int MaxHaltonResolution = 128;
    Point2i baseScales, baseExponents;
    int multInverse[2];
    int64_t haltonIndex = 0;
    int dimension = 0;
};

// PaddedSobolSampler Definition
class PaddedSobolSampler {
  public:
    // PaddedSobolSampler Public Methods
    PBRT_CPU_GPU
    static constexpr const char *Name() { return "PaddedSobolSampler"; }
    static PaddedSobolSampler *Create(const ParameterDictionary &parameters,
                                      const FileLoc *loc, Allocator alloc);

    PaddedSobolSampler(int samplesPerPixel, RandomizeStrategy randomizer, int seed = 0)
        : samplesPerPixel(samplesPerPixel), randomize(randomizer), seed(seed) {
        if (!IsPowerOf2(samplesPerPixel))
            Warning(
                "Sobol samplers with non power-of-two sample counts (%d) are suboptimal.",
                samplesPerPixel);
    }

    PBRT_CPU_GPU
    int SamplesPerPixel() const { return samplesPerPixel; }

    PBRT_CPU_GPU
    void StartPixelSample(Point2i p, int index, int dim) {
        pixel = p;
        sampleIndex = index;
        dimension = dim;
    }

    PBRT_CPU_GPU
    Float Get1D() {
        // Get permuted index for current pixel sample
        uint64_t hash = Hash(pixel, dimension, seed);
        int index = PermutationElement(sampleIndex, samplesPerPixel, hash);

        int dim = dimension++;
        // Return randomized 1D van der Corput sample for dimension _dim_
        return SampleDimension(0, index, hash >> 32);
    }

    PBRT_CPU_GPU
    Point2f Get2D() {
        // Get permuted index for current pixel sample
        uint64_t hash = Hash(pixel, dimension, seed);
        int index = PermutationElement(sampleIndex, samplesPerPixel, hash);

        int dim = dimension;
        dimension += 2;
        // Return randomized 2D Sobol\+$'$ sample
        return Point2f(SampleDimension(0, index, uint32_t(hash)),
                       SampleDimension(1, index, hash >> 32));
    }

    PBRT_CPU_GPU
    Point2f GetPixel2D() { return Get2D(); }

    PBRT_CPU_GPU
    RandomizeStrategy GetRandomizeStrategy() const { return randomize; }

    Sampler Clone(Allocator alloc);
    std::string ToString() const;

  private:
    // PaddedSobolSampler Private Methods
    PBRT_CPU_GPU
    Float SampleDimension(int dimension, uint32_t a, uint32_t hash) const {
        if (randomize == RandomizeStrategy::None)
            return SobolSample(a, dimension, NoRandomizer());
        else if (randomize == RandomizeStrategy::PermuteDigits)
            return SobolSample(a, dimension, BinaryPermuteScrambler(hash));
        else if (randomize == RandomizeStrategy::FastOwen)
            return SobolSample(a, dimension, FastOwenScrambler(hash));
        else
            return SobolSample(a, dimension, OwenScrambler(hash));
    }

    // PaddedSobolSampler Private Members
    int samplesPerPixel, seed;
    RandomizeStrategy randomize;
    Point2i pixel;
    int sampleIndex, dimension;
};

// ZSobolSampler Definition
class ZSobolSampler {
  public:
    // ZSobolSampler Public Methods
    ZSobolSampler(int samplesPerPixel, Point2i fullResolution,
                  RandomizeStrategy randomize, int seed = 0)
        : randomize(randomize), seed(seed) {
        if (!IsPowerOf2(samplesPerPixel))
            Warning(
                "Sobol samplers with non power-of-two sample counts (%d) are suboptimal.",
                samplesPerPixel);
        log2SamplesPerPixel = Log2Int(samplesPerPixel);
        int res = RoundUpPow2(std::max(fullResolution.x, fullResolution.y));
        int log4SamplesPerPixel = (log2SamplesPerPixel + 1) / 2;
        nBase4Digits = Log2Int(res) + log4SamplesPerPixel;
    }

    PBRT_CPU_GPU
    static constexpr const char *Name() { return "ZSobolSampler"; }

    static ZSobolSampler *Create(const ParameterDictionary &parameters,
                                 Point2i fullResolution, const FileLoc *loc,
                                 Allocator alloc);

    PBRT_CPU_GPU
    int SamplesPerPixel() const { return 1 << log2SamplesPerPixel; }

    PBRT_CPU_GPU
    void StartPixelSample(Point2i p, int index, int dim) {
        dimension = dim;
        mortonIndex = (EncodeMorton2(p.x, p.y) << log2SamplesPerPixel) | index;
    }

    PBRT_CPU_GPU
    Float Get1D() {
        uint64_t sampleIndex = GetSampleIndex();
        ++dimension;
        // Generate 1D Sobol\+$'$ sample at _sampleIndex_
        uint32_t sampleHash = Hash(dimension, seed);
        if (randomize == RandomizeStrategy::None)
            return SobolSample(sampleIndex, 0, NoRandomizer());
        else if (randomize == RandomizeStrategy::PermuteDigits)
            return SobolSample(sampleIndex, 0, BinaryPermuteScrambler(sampleHash));
        else if (randomize == RandomizeStrategy::FastOwen)
            return SobolSample(sampleIndex, 0, FastOwenScrambler(sampleHash));
        else
            return SobolSample(sampleIndex, 0, OwenScrambler(sampleHash));
    }

    PBRT_CPU_GPU
    Point2f Get2D() {
        uint64_t sampleIndex = GetSampleIndex();
        dimension += 2;
        // Generate 2D Sobol\+$'$ sample at _sampleIndex_
        uint64_t bits = Hash(dimension, seed);
        uint32_t sampleHash[2] = {uint32_t(bits), uint32_t(bits >> 32)};
        if (randomize == RandomizeStrategy::None)
            return {SobolSample(sampleIndex, 0, NoRandomizer()),
                    SobolSample(sampleIndex, 1, NoRandomizer())};
        else if (randomize == RandomizeStrategy::PermuteDigits)
            return {SobolSample(sampleIndex, 0, BinaryPermuteScrambler(sampleHash[0])),
                    SobolSample(sampleIndex, 1, BinaryPermuteScrambler(sampleHash[1]))};
        else if (randomize == RandomizeStrategy::FastOwen)
            return {SobolSample(sampleIndex, 0, FastOwenScrambler(sampleHash[0])),
                    SobolSample(sampleIndex, 1, FastOwenScrambler(sampleHash[1]))};
        else
            return {SobolSample(sampleIndex, 0, OwenScrambler(sampleHash[0])),
                    SobolSample(sampleIndex, 1, OwenScrambler(sampleHash[1]))};
    }

    PBRT_CPU_GPU
    Point2f GetPixel2D() { return Get2D(); }

    Sampler Clone(Allocator alloc);
    std::string ToString() const;

    PBRT_CPU_GPU
    uint64_t GetSampleIndex() const {
        // Define the full set of 4-way permutations in _permutations_
        static const uint8_t permutations[24][4] = {
            {0, 1, 2, 3},
            {0, 1, 3, 2},
            {0, 2, 1, 3},
            {0, 2, 3, 1},
            // Define remaining 20 4-way permutations
            {0, 3, 2, 1},
            {0, 3, 1, 2},
            {1, 0, 2, 3},
            {1, 0, 3, 2},
            {1, 2, 0, 3},
            {1, 2, 3, 0},
            {1, 3, 2, 0},
            {1, 3, 0, 2},
            {2, 1, 0, 3},
            {2, 1, 3, 0},
            {2, 0, 1, 3},
            {2, 0, 3, 1},
            {2, 3, 0, 1},
            {2, 3, 1, 0},
            {3, 1, 2, 0},
            {3, 1, 0, 2},
            {3, 2, 1, 0},
            {3, 2, 0, 1},
            {3, 0, 2, 1},
            {3, 0, 1, 2}

        };

        uint64_t sampleIndex = 0;
        // Apply random permutations to full base-4 digits
        bool pow2Samples = log2SamplesPerPixel & 1;
        int lastDigit = pow2Samples ? 1 : 0;
        for (int i = nBase4Digits - 1; i >= lastDigit; --i) {
            // Randomly permute $i$th base-4 digit in _mortonIndex_
            int digitShift = 2 * i - (pow2Samples ? 1 : 0);
            int digit = (mortonIndex >> digitShift) & 3;
            // Choose permutation _p_ to use for _digit_
            uint64_t higherDigits = mortonIndex >> (digitShift + 2);
            int p = (MixBits(higherDigits ^ (0x55555555u * dimension)) >> 24) % 24;

            digit = permutations[p][digit];
            sampleIndex |= uint64_t(digit) << digitShift;
        }

        // Handle power-of-2 (but not 4) sample count
        if (pow2Samples) {
            int digit = mortonIndex & 1;
            sampleIndex |=
                digit ^ (MixBits((mortonIndex >> 1) ^ (0x55555555u * dimension)) & 1);
        }

        return sampleIndex;
    }

  private:
    // ZSobolSampler Private Members
    RandomizeStrategy randomize;
    int seed, log2SamplesPerPixel, nBase4Digits;
    uint64_t mortonIndex;
    int dimension;
};

// PMJ02BNSampler Definition
class PMJ02BNSampler {
  public:
    // PMJ02BNSampler Public Methods
    PMJ02BNSampler(int samplesPerPixel, int seed = 0, Allocator alloc = {});

    PBRT_CPU_GPU
    static constexpr const char *Name() { return "PMJ02BNSampler"; }
    static PMJ02BNSampler *Create(const ParameterDictionary &parameters,
                                  const FileLoc *loc, Allocator alloc);

    PBRT_CPU_GPU
    int SamplesPerPixel() const { return samplesPerPixel; }

    PBRT_CPU_GPU
    void StartPixelSample(Point2i p, int index, int dim) {
        pixel = p;
        sampleIndex = index;
        dimension = std::max(2, dim);
    }

    PBRT_CPU_GPU
    Float Get1D() {
        // Find permuted sample index for 1D PMJ02BNSampler sample
        uint64_t hash = Hash(pixel, dimension, seed);
        int index = PermutationElement(sampleIndex, samplesPerPixel, hash);

        Float delta = BlueNoise(dimension, pixel);
        ++dimension;
        return std::min((index + delta) / samplesPerPixel, OneMinusEpsilon);
    }

    PBRT_CPU_GPU
    Point2f GetPixel2D() {
        int px = pixel.x % pixelTileSize, py = pixel.y % pixelTileSize;
        int offset = (px + py * pixelTileSize) * samplesPerPixel;
        return (*pixelSamples)[offset + sampleIndex];
    }

    PBRT_CPU_GPU
    Point2f Get2D() {
        // Compute index for 2D pmj02bn sample
        int index = sampleIndex;
        int pmjInstance = dimension / 2;
        if (pmjInstance >= nPMJ02bnSets) {
            // Permute index to be used for pmj02bn sample array
            uint64_t hash = Hash(pixel, dimension, seed);
            index = PermutationElement(sampleIndex, samplesPerPixel, hash);
        }

        // Return randomized pmj02bn sample for current dimension
        Point2f u = GetPMJ02BNSample(pmjInstance, index);
        // Apply Cranley-Patterson rotation to pmj02bn sample _u_
        u += Vector2f(BlueNoise(dimension, pixel), BlueNoise(dimension + 1, pixel));
        if (u.x >= 1)
            u.x -= 1;
        if (u.y >= 1)
            u.y -= 1;

        dimension += 2;
        return {std::min(u.x, OneMinusEpsilon), std::min(u.y, OneMinusEpsilon)};
    }

    Sampler Clone(Allocator alloc);
    std::string ToString() const;

  private:
    // PMJ02BNSampler Private Members
    int samplesPerPixel, seed;
    int pixelTileSize;
    pstd::vector<Point2f> *pixelSamples;
    Point2i pixel;
    int sampleIndex, dimension;
};

// IndependentSampler Definition
class IndependentSampler {
  public:
    // IndependentSampler Public Methods
    IndependentSampler(int samplesPerPixel, int seed = 0)
        : samplesPerPixel(samplesPerPixel), seed(seed) {}

    static IndependentSampler *Create(const ParameterDictionary &parameters,
                                      const FileLoc *loc, Allocator alloc);
    PBRT_CPU_GPU
    static constexpr const char *Name() { return "IndependentSampler"; }

    PBRT_CPU_GPU
    int SamplesPerPixel() const { return samplesPerPixel; }

    PBRT_CPU_GPU
    void StartPixelSample(Point2i p, int sampleIndex, int dimension) {
        rng.SetSequence(Hash(p, seed));
        rng.Advance(sampleIndex * 65536ull + dimension);
    }

    PBRT_CPU_GPU
    Float Get1D() { return rng.Uniform<Float>(); }
    PBRT_CPU_GPU
    Point2f Get2D() { return {rng.Uniform<Float>(), rng.Uniform<Float>()}; }
    PBRT_CPU_GPU
    Point2f GetPixel2D() { return Get2D(); }

    Sampler Clone(Allocator alloc);
    std::string ToString() const;

  private:
    // IndependentSampler Private Members
    int samplesPerPixel, seed;
    RNG rng;
};

// SobolSampler Definition
class SobolSampler {
  public:
    // SobolSampler Public Methods
    SobolSampler(int samplesPerPixel, Point2i fullResolution, RandomizeStrategy randomize,
                 int seed = 0)
        : samplesPerPixel(samplesPerPixel), seed(seed), randomize(randomize) {
        if (!IsPowerOf2(samplesPerPixel))
            Warning("Non power-of-two sample count %d will perform suboptimally with the "
                    "SobolSampler.",
                    samplesPerPixel);
        scale = RoundUpPow2(std::max(fullResolution.x, fullResolution.y));
    }

    PBRT_CPU_GPU
    static constexpr const char *Name() { return "SobolSampler"; }
    static SobolSampler *Create(const ParameterDictionary &parameters,
                                Point2i fullResolution, const FileLoc *loc,
                                Allocator alloc);

    PBRT_CPU_GPU
    int SamplesPerPixel() const { return samplesPerPixel; }

    PBRT_CPU_GPU
    void StartPixelSample(Point2i p, int sampleIndex, int dim) {
        pixel = p;
        dimension = std::max(2, dim);
        sobolIndex = SobolIntervalToIndex(Log2Int(scale), sampleIndex, pixel);
    }

    PBRT_CPU_GPU
    Float Get1D() {
        if (dimension >= NSobolDimensions)
            dimension = 2;
        return SampleDimension(dimension++);
    }

    PBRT_CPU_GPU
    Point2f Get2D() {
        if (dimension + 1 >= NSobolDimensions)
            dimension = 2;
        Point2f u(SampleDimension(dimension), SampleDimension(dimension + 1));
        dimension += 2;
        return u;
    }

    PBRT_CPU_GPU
    Point2f GetPixel2D() {
        Point2f u(SobolSample(sobolIndex, 0, NoRandomizer()),
                  SobolSample(sobolIndex, 1, NoRandomizer()));
        // Remap Sobol\+$'$ dimensions used for pixel samples
        for (int dim = 0; dim < 2; ++dim) {
            DCHECK_RARE(1e-7, u[dim] * scale - pixel[dim] < 0);
            DCHECK_RARE(1e-7, u[dim] * scale - pixel[dim] > 1);
            u[dim] = Clamp(u[dim] * scale - pixel[dim], 0, OneMinusEpsilon);
        }

        return u;
    }

    Sampler Clone(Allocator alloc);
    std::string ToString() const;

  private:
    // SobolSampler Private Methods
    PBRT_CPU_GPU
    Float SampleDimension(int dimension) const {
        // Return un-randomized Sobol\+$'$ sample if appropriate
        if (randomize == RandomizeStrategy::None)
            return SobolSample(sobolIndex, dimension, NoRandomizer());

        // Return randomized Sobol\+$'$ sample using _randomize_
        uint32_t hash = Hash(dimension, seed);
        if (randomize == RandomizeStrategy::PermuteDigits)
            return SobolSample(sobolIndex, dimension, BinaryPermuteScrambler(hash));
        else if (randomize == RandomizeStrategy::FastOwen)
            return SobolSample(sobolIndex, dimension, FastOwenScrambler(hash));
        else
            return SobolSample(sobolIndex, dimension, OwenScrambler(hash));
    }

    // SobolSampler Private Members
    int samplesPerPixel, scale, seed;
    RandomizeStrategy randomize;
    Point2i pixel;
    int dimension;
    int64_t sobolIndex;
};

// StratifiedSampler Definition
class StratifiedSampler {
  public:
    // StratifiedSampler Public Methods
    StratifiedSampler(int xPixelSamples, int yPixelSamples, bool jitter, int seed = 0)
        : xPixelSamples(xPixelSamples),
          yPixelSamples(yPixelSamples),
          seed(seed),
          jitter(jitter) {}

    static StratifiedSampler *Create(const ParameterDictionary &parameters,
                                     const FileLoc *loc, Allocator alloc);
    PBRT_CPU_GPU
    static constexpr const char *Name() { return "StratifiedSampler"; }

    PBRT_CPU_GPU
    int SamplesPerPixel() const { return xPixelSamples * yPixelSamples; }

    PBRT_CPU_GPU
    void StartPixelSample(Point2i p, int index, int dim) {
        pixel = p;
        sampleIndex = index;
        dimension = dim;
        rng.SetSequence(Hash(p, seed));
        rng.Advance(sampleIndex * 65536ull + dimension);
    }

    PBRT_CPU_GPU
    Float Get1D() {
        // Compute _stratum_ index for current pixel and dimension
        uint64_t hash = Hash(pixel, dimension, seed);
        int stratum = PermutationElement(sampleIndex, SamplesPerPixel(), hash);

        ++dimension;
        Float delta = jitter ? rng.Uniform<Float>() : 0.5f;
        return (stratum + delta) / SamplesPerPixel();
    }

    PBRT_CPU_GPU
    Point2f Get2D() {
        // Compute _stratum_ index for current pixel and dimension
        uint64_t hash = Hash(pixel, dimension, seed);
        int stratum = PermutationElement(sampleIndex, SamplesPerPixel(), hash);

        dimension += 2;
        int x = stratum % xPixelSamples, y = stratum / xPixelSamples;
        Float dx = jitter ? rng.Uniform<Float>() : 0.5f;
        Float dy = jitter ? rng.Uniform<Float>() : 0.5f;
        return {(x + dx) / xPixelSamples, (y + dy) / yPixelSamples};
    }

    PBRT_CPU_GPU
    Point2f GetPixel2D() { return Get2D(); }

    Sampler Clone(Allocator alloc);
    std::string ToString() const;

  private:
    // StratifiedSampler Private Members
    int xPixelSamples, yPixelSamples, seed;
    bool jitter;
    RNG rng;
    Point2i pixel;
    int sampleIndex = 0, dimension = 0;
};

// MLTSampler Definition
class MLTSampler {
  public:
    // MLTSampler Public Methods
    MLTSampler(int mutationsPerPixel, int rngSequenceIndex, Float sigma,
               Float largeStepProbability, int streamCount)
        : mutationsPerPixel(mutationsPerPixel),
          rng(MixBits(rngSequenceIndex) ^ MixBits(Options->seed)),
          sigma(sigma),
          largeStepProbability(largeStepProbability),
          streamCount(streamCount) {}

    PBRT_CPU_GPU
    void StartIteration();

    PBRT_CPU_GPU
    void Reject();

    PBRT_CPU_GPU
    void StartStream(int index);

    PBRT_CPU_GPU
    int GetNextIndex() { return streamIndex + streamCount * sampleIndex++; }

    PBRT_CPU_GPU
    int SamplesPerPixel() const { return mutationsPerPixel; }

    PBRT_CPU_GPU
    void StartPixelSample(Point2i p, int sampleIndex, int dim) {
        rng.SetSequence(Hash(p));
        rng.Advance(sampleIndex * 65536 + dim * 8192);
    }

    PBRT_CPU_GPU
    Float Get1D();

    PBRT_CPU_GPU
    Point2f Get2D();

    PBRT_CPU_GPU
    Point2f GetPixel2D();

    Sampler Clone(Allocator alloc);

    PBRT_CPU_GPU
    void Accept();

    std::string DumpState() const;

    std::string ToString() const {
        return StringPrintf(
            "[ MLTSampler rng: %s sigma: %f largeStepProbability: %f "
            "streamCount: %d X: %s currentIteration: %d largeStep: %s "
            "lastLargeStepIteration: %d streamIndex: %d sampleIndex: %d ] ",
            rng, sigma, largeStepProbability, streamCount, X, currentIteration, largeStep,
            lastLargeStepIteration, streamIndex, sampleIndex);
    }

  protected:
    // MLTSampler Private Declarations
    struct PrimarySample {
        Float value = 0;
        // PrimarySample Public Methods
        PBRT_CPU_GPU
        void Backup() {
            valueBackup = value;
            modifyBackup = lastModificationIteration;
        }
        PBRT_CPU_GPU
        void Restore() {
            value = valueBackup;
            lastModificationIteration = modifyBackup;
        }

        std::string ToString() const {
            return StringPrintf("[ PrimarySample lastModificationIteration: %d "
                                "valueBackup: %f modifyBackup: %d ]",
                                lastModificationIteration, valueBackup, modifyBackup);
        }

        // PrimarySample Public Members
        int64_t lastModificationIteration = 0;
        Float valueBackup = 0;
        int64_t modifyBackup = 0;
    };

    // MLTSampler Private Methods
    PBRT_CPU_GPU
    void EnsureReady(int index);

    // MLTSampler Private Members
    int mutationsPerPixel;
    RNG rng;
    Float sigma, largeStepProbability;
    int streamCount;
    pstd::vector<PrimarySample> X;
    int64_t currentIteration = 0;
    bool largeStep = true;
    int64_t lastLargeStepIteration = 0;
    int streamIndex, sampleIndex;
};

class DebugMLTSampler : public MLTSampler {
  public:
    static DebugMLTSampler Create(pstd::span<const std::string> state,
                                  int nSampleStreams);

    PBRT_CPU_GPU
    Float Get1D() {
        int index = GetNextIndex();
#ifdef PBRT_IS_GPU_CODE
        return 0;
#else
        DCHECK_LT(index, u.size());
        return u[index];
#endif
    }

    PBRT_CPU_GPU
    Point2f Get2D() { return {Get1D(), Get1D()}; }

    PBRT_CPU_GPU
    Point2f GetPixel2D() { return Get2D(); }

    std::string ToString() const {
        return StringPrintf("[ DebugMLTSampler %s u: %s ]",
                            ((const MLTSampler *)this)->ToString(), u);
    }

  private:
    DebugMLTSampler(int nSampleStreams) : MLTSampler(1, 0, 0.5, 0.5, nSampleStreams) {}

    std::vector<Float> u;
};

PBRT_CPU_GPU inline void Sampler::StartPixelSample(Point2i p, int sampleIndex, int dimension) {
    auto start = [&](auto ptr) {
        return ptr->StartPixelSample(p, sampleIndex, dimension);
    };
    return Dispatch(start);
}

PBRT_CPU_GPU inline int Sampler::SamplesPerPixel() const {
    auto spp = [&](auto ptr) { return ptr->SamplesPerPixel(); };
    return Dispatch(spp);
}

PBRT_CPU_GPU inline Float Sampler::Get1D() {
    auto get = [&](auto ptr) { return ptr->Get1D(); };
    return Dispatch(get);
}

PBRT_CPU_GPU inline Point2f Sampler::Get2D() {
    auto get = [&](auto ptr) { return ptr->Get2D(); };
    return Dispatch(get);
}

PBRT_CPU_GPU inline Point2f Sampler::GetPixel2D() {
    auto get = [&](auto ptr) { return ptr->GetPixel2D(); };
    return Dispatch(get);
}

// Sampler Inline Functions
template <typename S>
inline PBRT_CPU_GPU CameraSample GetCameraSample(S sampler, Point2i pPixel,
                                                 Filter filter) {
    FilterSample fs = filter.Sample(sampler.GetPixel2D());
    CameraSample cs;
    // Initialize _CameraSample_ member variables
    cs.pFilm = pPixel + fs.p + Vector2f(0.5f, 0.5f);
    cs.time = sampler.Get1D();
    cs.pLens = sampler.Get2D();
    cs.filterWeight = fs.weight;

    if (GetOptions().disablePixelJitter) {
        cs.pFilm = pPixel + Vector2f(0.5f, 0.5f);
        cs.time = 0.5f;
        cs.pLens = Point2f(0.5f, 0.5f);
        cs.filterWeight = 1;
    }
    return cs;
}

}  // namespace pbrt

#endif  // PBRT_SAMPLERS_H
