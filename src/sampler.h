#ifndef SAMPLERH
#define SAMPLERH

#include "rng.h"
#include "vec2.h"

class Sampler {
public:
  Sampler(size_t number_pixel_samples) : samplesPerPixel(number_pixel_samples) {}
  virtual ~Sampler() = 0;
  // void markOccupiedStrata1(vec2 pt, int NN) {
  //   int shape = 0;
  //   int xdivs = NN;
  //   int ydivs = 1;
  //   while(xdivs != 0) {
  //     int xstratum = (xdivs * pt.x());
  //     int ystratum = (ydivs * pt.y());
  //   }
  // }
  // std::vector<bool> occupiedStrata;
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
  PixelSampler(size_t samplesPerPixel, size_t nSampledDimensions, unsigned int seed)
    : Sampler(samplesPerPixel) {
    for (size_t i = 0; i < nSampledDimensions; ++i) {
      samples1D.push_back(std::vector<Float>(samplesPerPixel));
      samples2D.push_back(std::vector<vec2>(samplesPerPixel));
    }
    random_gen rng2(seed);
    rng = rng2;
    current1DDimension = 0;
    current2DDimension = 0;
  }
  bool StartNextSample();
  bool SetSampleNumber(size_t);
  Float Get1D();
  vec2 Get2D();
  
protected:
  std::vector<std::vector<Float>> samples1D;
  std::vector<std::vector<vec2>> samples2D;
  int current1DDimension;
  int current2DDimension;
  random_gen rng;
};

class StratifiedSampler : public PixelSampler {
public:
  StratifiedSampler(int xPixelSamples, int yPixelSamples,
                  bool jitterSamples, size_t nSampledDimensions, unsigned int seed)
  : PixelSampler(xPixelSamples * yPixelSamples, nSampledDimensions, seed),
    xPixelSamples(xPixelSamples), yPixelSamples(yPixelSamples), jitterSamples(jitterSamples) { }
  void StartPixel(const vec2 &p);
  std::unique_ptr<Sampler> Clone(int seed);

private:
  const int xPixelSamples, yPixelSamples;
  const bool jitterSamples;
};

#endif
