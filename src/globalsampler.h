#ifndef GLOBALSAMPLERH
#define GLOBALSAMPLERH

#include "sampler.h"

class GlobalSampler : public Sampler {
  public:
    bool StartNextSample();
  void StartPixel(const vec2f &);
  bool SetSampleNumber(int64_t sampleNum);
  Float Get1D();
  vec2f Get2D();
  GlobalSampler(int64_t samplesPerPixel) : Sampler(samplesPerPixel) { }
  virtual int64_t GetIndexForSample(int64_t sampleNum) const = 0;
  virtual Float SampleDimension(int64_t index, int dimension) const = 0;
  
  private:
    int dimension;
  int64_t intervalSampleIndex;
  static const int arrayStartDim = 5;
  int arrayEndDim;
};

void GlobalSampler::StartPixel(const vec2f &p) {
  Sampler::StartPixel(p);
  dimension = 0;
  intervalSampleIndex = GetIndexForSample(0);
  // Compute _arrayEndDim_ for dimensions used for array samples
  arrayEndDim = arrayStartDim + sampleArray1D.size() + 2 * sampleArray2D.size();
  
  // Compute 1D array samples for _GlobalSampler_
  for (size_t i = 0; i < samples1DArraySizes.size(); ++i) {
    int nSamples = samples1DArraySizes[i] * samplesPerPixel;
    for (int j = 0; j < nSamples; ++j) {
      int64_t index = GetIndexForSample(j);
      sampleArray1D[i][j] = SampleDimension(index, arrayStartDim + i);
    }
  }
  
  // Compute 2D array samples for _GlobalSampler_
  int dim = arrayStartDim + samples1DArraySizes.size();
  for (size_t i = 0; i < samples2DArraySizes.size(); ++i) {
    int nSamples = samples2DArraySizes[i] * samplesPerPixel;
    for (int j = 0; j < nSamples; ++j) {
      int64_t idx = GetIndexForSample(j);
      sampleArray2D[i][j].e[0] = SampleDimension(idx, dim);
      sampleArray2D[i][j].e[1] = SampleDimension(idx, dim + 1);
    }
    dim += 2;
  }
}

bool GlobalSampler::StartNextSample() {
  dimension = 0;
  intervalSampleIndex = GetIndexForSample(currentPixelSampleIndex + 1);
  return(Sampler::StartNextSample());
}

bool GlobalSampler::SetSampleNumber(int64_t sampleNum) {
  dimension = 0;
  intervalSampleIndex = GetIndexForSample(sampleNum);
  return(Sampler::SetSampleNumber(sampleNum));
}

Float GlobalSampler::Get1D() {
  if (dimension >= arrayStartDim && dimension < arrayEndDim) {
    dimension = arrayEndDim;
  }
  return(SampleDimension(intervalSampleIndex, dimension++));
}

vec2f GlobalSampler::Get2D() {
  if (dimension + 1 >= arrayStartDim && dimension < arrayEndDim) {
    dimension = arrayEndDim;
  }
  vec2f p(SampleDimension(intervalSampleIndex, dimension),
         SampleDimension(intervalSampleIndex, dimension + 1));
  dimension += 2;
  return(p);
}

#endif
