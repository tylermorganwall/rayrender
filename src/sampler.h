#ifndef SAMPLERH
#define SAMPLERH

#include "rng.h"
#include "vec2.h"
#include <memory>
#include "single_sample.h"

class Sampler {
public:
  Sampler(size_t number_pixel_samples) : samplesPerPixel(number_pixel_samples) {}
  virtual ~Sampler() = 0;
  virtual void StartPixel(const unsigned int i, const unsigned int j);
  virtual Float Get1D() = 0;
  virtual vec2f Get2D() = 0;
  virtual int RoundCount(int n) const {
    return n;
  }
  const Float *Get1DArray(int n);
  const vec2f *Get2DArray(int n);
  void Request1DArray(int n);
  void Request2DArray(int n);
  
  virtual bool StartNextSample();
  virtual std::unique_ptr<Sampler> Clone(int seed) = 0;
  virtual bool SetSampleNumber(size_t sampleNum);
  
  const size_t samplesPerPixel;
  
protected:
  //Add int vec2f class
  unsigned int currentPixelx, currentPixely;
  size_t currentPixelSampleIndex;
  std::vector<int> samples1DArraySizes, samples2DArraySizes;
  std::vector<std::vector<Float>> sampleArray1D;
  std::vector<std::vector<vec2f>> sampleArray2D;
  
private: 
  size_t array1DOffset, array2DOffset;
};

class PixelSampler : public Sampler {
public:
  PixelSampler(size_t samplesPerPixel, size_t nSampledDimensions, random_gen& rng)
    : Sampler(samplesPerPixel), current1DDimension(0), current2DDimension(0), rng(rng) {
    for (size_t i = 0; i < nSampledDimensions; ++i) {
      samples1D.push_back(std::vector<Float>(samplesPerPixel));
      samples2D.push_back(std::vector<vec2f>(samplesPerPixel));
    }
  }
  bool StartNextSample();
  bool SetSampleNumber(size_t);
  virtual Float Get1D();
  virtual vec2f Get2D();
  
protected:
  std::vector<std::vector<Float> > samples1D;
  std::vector<std::vector<vec2f> > samples2D;
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
  void StartPixel( unsigned int i,  unsigned int j);
  std::unique_ptr<Sampler> Clone(int seed);

private:
  const int xPixelSamples, yPixelSamples;
  const bool jitterSamples;
};

class RandomSampler : public PixelSampler {
public:
  RandomSampler(random_gen& rng) : PixelSampler(1 * 1, 0, rng) {}
  void StartPixel( unsigned int i,  unsigned int j) {};
  Float Get1D();
  vec2f Get2D();
  bool StartNextSample();
  bool SetSampleNumber(size_t);
  std::unique_ptr<Sampler> Clone(int seed);
};

class SobolSampler : public PixelSampler {
  public:
    SobolSampler(unsigned int xPixelSamples, unsigned int yPixelSamples,
                 unsigned int maxSamples,
                 random_gen& rng);
    void StartPixel( unsigned int i,  unsigned int j);
    Float Get1D();
    vec2f Get2D();
    bool StartNextSample();
    bool SetSampleNumber(size_t);
    std::unique_ptr<Sampler> Clone(int seed);
  private:
    const unsigned int xPixelSamples, yPixelSamples;
    unsigned int current1Dsample, current2Dsample;
    unsigned int pixelseed;
};

class SobolBlueNoiseSampler : public PixelSampler {
  public:
    SobolBlueNoiseSampler(random_gen& rng);
    void StartPixel( unsigned int i,  unsigned int j);
    Float Get1D();
    vec2f Get2D();
    bool StartNextSample();
    bool SetSampleNumber(size_t);
    std::unique_ptr<Sampler> Clone(int seed);
    private:
      unsigned int current1Dsample, current2Dsample;

};

#endif
