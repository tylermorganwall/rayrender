#ifndef SAMPLERH
#define SAMPLERH

#include "rng.h"
#include "vec2.h"
#include <memory>

class Sampler {
public:
  Sampler(size_t number_pixel_samples) : samplesPerPixel(number_pixel_samples) {}
  virtual ~Sampler() = 0;
  virtual void StartPixel(const vec2 &p);
  virtual Float Get1D() = 0;
  virtual vec2 Get2D() = 0;
  virtual int RoundCount(int n) const {
    return n;
  }
  const Float *Get1DArray(int n);
  const vec2 *Get2DArray(int n);
  void Request1DArray(int n);
  void Request2DArray(int n);
  
  virtual bool StartNextSample();
  virtual std::unique_ptr<Sampler> Clone(int seed) = 0;
  virtual bool SetSampleNumber(size_t sampleNum);
  
  const size_t samplesPerPixel;
  
protected:
  //Add int vec2 class
  vec2 currentPixel;
  size_t currentPixelSampleIndex;
  std::vector<int> samples1DArraySizes, samples2DArraySizes;
  std::vector<std::vector<Float>> sampleArray1D;
  std::vector<std::vector<vec2>> sampleArray2D;
  
private: 
  size_t array1DOffset, array2DOffset;
};

class PixelSampler : public Sampler {
public:
  PixelSampler(size_t samplesPerPixel, size_t nSampledDimensions, random_gen& rng)
    : Sampler(samplesPerPixel), current1DDimension(0), current2DDimension(0), rng(rng) {
    for (size_t i = 0; i < nSampledDimensions; ++i) {
      samples1D.push_back(std::vector<Float>(samplesPerPixel));
      samples2D.push_back(std::vector<vec2>(samplesPerPixel));
    }
  }
  bool StartNextSample();
  bool SetSampleNumber(size_t);
  virtual Float Get1D();
  virtual vec2 Get2D();
  
protected:
  std::vector<std::vector<Float> > samples1D;
  std::vector<std::vector<vec2> > samples2D;
  size_t current1DDimension;
  size_t current2DDimension;
  random_gen rng;
};

class StratifiedSampler : public PixelSampler {
public:
  StratifiedSampler(int xPixelSamples, int yPixelSamples, 
                  bool jitterSamples, size_t nSampledDimensions, random_gen& rng)
  : PixelSampler(xPixelSamples * yPixelSamples, nSampledDimensions, rng),
    xPixelSamples(xPixelSamples), yPixelSamples(yPixelSamples), jitterSamples(jitterSamples) { }
  void StartPixel(const vec2 &p);
  std::unique_ptr<Sampler> Clone(int seed);

private:
  const int xPixelSamples, yPixelSamples;
  const bool jitterSamples;
};

class RandomSampler : public PixelSampler {
public:
  RandomSampler(random_gen& rng) : PixelSampler(1 * 1, 0, rng) {}
  void StartPixel(const vec2 &p) {};
  Float Get1D() {
    return(rng.unif_rand());
  }
  vec2 Get2D() {
    return(vec2(rng.unif_rand(), rng.unif_rand()));
  }
  bool StartNextSample() {
    return(true);
  }
  bool SetSampleNumber(size_t) {
    return(true);
  }
  std::unique_ptr<Sampler> Clone(int seed);
};

#endif
