#ifndef SAMPLERH
#define SAMPLERH

#include "rng.h"
#include "vec2.h"

struct CameraSample {
  vec2 pFilm;
  vec2 pLens;
  Float time;
};

template <typename T>
void Shuffle(T *samp, int count, int nDimensions, random_gen &rng) {
  for (int i = 0; i < count; ++i) {
    int other = i + rng.UniformUInt32(count - i);
    for (int j = 0; j < nDimensions; ++j)
      std::swap(samp[nDimensions * i + j],
                samp[nDimensions * other + j]);
  }
}

class Sampler {
public:
  Sampler(int number_pixel_samples) : samplesPerPixel(number_pixel_samples) {}
  
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
  CameraSample GetCameraSample(const vec2 &pRaster) ;
  virtual int RoundCount(int n) const {
    return n;
  }
  const Float *Get1DArray(int n);
  const vec2 *Get2DArray(int n);
  void Request1DArray(int n);
  void Request2DArray(int n);
  
  virtual bool StartNextSample();
  virtual std::unique_ptr<Sampler> Clone(int seed) = 0;
  virtual bool SetSampleNumber(int64_t sampleNum);
  
  const int64_t samplesPerPixel;
  
protected:
  vec2 currentPixel;
  int64_t currentPixelSampleIndex;
  std::vector<int> samples1DArraySizes, samples2DArraySizes;
  std::vector<std::vector<Float>> sampleArray1D;
  std::vector<std::vector<vec2>> sampleArray2D;
  
private: 
  size_t array1DOffset, array2DOffset;
};

void Sampler::StartPixel(const vec2 &p) {
  currentPixel = p;
  currentPixelSampleIndex = 0;
  array1DOffset = array2DOffset = 0;
}

bool Sampler::StartNextSample() {
  array1DOffset = array2DOffset = 0;
  return(++currentPixelSampleIndex < samplesPerPixel);
}

bool Sampler::SetSampleNumber(int64_t sampleNum) {
  currentPixelSampleIndex = sampleNum;
  return(currentPixelSampleIndex < samplesPerPixel);
}

void Sampler::Request1DArray(int n) {
  samples1DArraySizes.push_back(n);
  sampleArray1D.push_back(std::vector<Float>(n * samplesPerPixel));
}

void Sampler::Request2DArray(int n) {
  samples2DArraySizes.push_back(n);
  sampleArray2D.push_back(std::vector<vec2>(n * samplesPerPixel));
}

const Float *Sampler::Get1DArray(int n) {
  if (array1DOffset == sampleArray1D.size()) {
    return(nullptr);
  }
  return(&sampleArray1D[array1DOffset++][currentPixelSampleIndex * n]);
}

const vec2 *Sampler::Get2DArray(int n) {
  if (array2DOffset == sampleArray2D.size()) {
    return(nullptr);
  }
  return(&sampleArray2D[array2DOffset++][currentPixelSampleIndex * n]);
}

CameraSample Sampler::GetCameraSample(const vec2 &pRaster) {
  CameraSample cs;
  cs.pFilm = (vec2)pRaster + Get2D();
  cs.time = Get1D();
  cs.pLens = Get2D();
  return cs;
}


class PixelSampler : public Sampler {
public:
  PixelSampler(int64_t samplesPerPixel,int nSampledDimensions)
    : Sampler(samplesPerPixel) {
    for (int i = 0; i < nSampledDimensions; ++i) {
      samples1D.push_back(std::vector<Float>(samplesPerPixel));
      samples2D.push_back(std::vector<vec2>(samplesPerPixel));
    }
  }
  bool StartNextSample();
  bool SetSampleNumber(int64_t);
  Float Get1D();
  vec2 Get2D();
  
protected:
  std::vector<std::vector<Float>> samples1D;
  std::vector<std::vector<vec2>> samples2D;
  int current1DDimension = 0, current2DDimension = 0;
  random_gen rng;
};

bool PixelSampler::StartNextSample() {
  current1DDimension = current2DDimension = 0;
  return(Sampler::StartNextSample());
}

bool PixelSampler::SetSampleNumber(int64_t sampleNum) {
  current1DDimension = current2DDimension = 0;
  return(Sampler::SetSampleNumber(sampleNum));
}

Float PixelSampler::Get1D() {
  if (current1DDimension < samples1D.size()) {
    return(samples1D[current1DDimension++][currentPixelSampleIndex]);
  } else {
    return(rng.unif_rand());
  }
}

vec2 PixelSampler::Get2D() {
  if (current2DDimension < samples2D.size()) {
    return(samples2D[current2DDimension++][currentPixelSampleIndex]);
  } else {
    return(vec2(rng.unif_rand(), rng.unif_rand()));
  }
}

class GlobalSampler : public Sampler {
public:
  bool StartNextSample();
  void StartPixel(const vec2 &);
  bool SetSampleNumber(int64_t sampleNum);
  Float Get1D();
  vec2 Get2D();
  GlobalSampler(int64_t samplesPerPixel) : Sampler(samplesPerPixel) { }
  virtual int64_t GetIndexForSample(int64_t sampleNum) const = 0;
  virtual Float SampleDimension(int64_t index, int dimension) const = 0;
  
private:
  int dimension;
  int64_t intervalSampleIndex;
  static const int arrayStartDim = 5;
  int arrayEndDim;
};

void GlobalSampler::StartPixel(const vec2 &p) {
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

vec2 GlobalSampler::Get2D() {
  if (dimension + 1 >= arrayStartDim && dimension < arrayEndDim) {
    dimension = arrayEndDim;
  }
  vec2 p(SampleDimension(intervalSampleIndex, dimension),
         SampleDimension(intervalSampleIndex, dimension + 1));
  dimension += 2;
  return(p);
}

class StratifiedSampler : public PixelSampler {
public:
  StratifiedSampler(int xPixelSamples, int yPixelSamples,
                  bool jitterSamples, int nSampledDimensions)
  : PixelSampler(xPixelSamples * yPixelSamples, nSampledDimensions),
    xPixelSamples(xPixelSamples), yPixelSamples(yPixelSamples), jitterSamples(jitterSamples) { }
  void StartPixel(const vec2 &);
  std::unique_ptr<Sampler> Clone(int seed);
  
private:
  const int xPixelSamples, yPixelSamples;
  const bool jitterSamples;
  
};

std::unique_ptr<Sampler> StratifiedSampler::Clone(int seed) {
  StratifiedSampler *ss = new StratifiedSampler(*this);
  ss->rng.SetSequence(seed);
  return(std::unique_ptr<Sampler>(ss));
}

void LatinHypercube(Float *samples, int nSamples, int nDim, random_gen &rng) {
Float invNSamples = (Float)1 / nSamples;
  for (int i = 0; i < nSamples; ++i) {
    for (int j = 0; j < nDim; ++j) {
      Float sj = (i + (rng.unif_rand())) * invNSamples;
      samples[nDim * i + j] = std::min(sj, OneMinusEpsilon);
    }
  }
  for (int i = 0; i < nDim; ++i) {
    for (int j = 0; j < nSamples; ++j) {
      int other = j + rng.UniformUInt32(nSamples - j);
      std::swap(samples[nDim * j + i], samples[nDim * other + i]);
    }
  }
}

void StratifiedSample1D(Float *samp, int nSamples, random_gen &rng,
                        bool jitter) {
  Float invNSamples = (Float)1 / nSamples;
  for (int i = 0; i < nSamples; ++i) {
    Float delta = jitter ? rng.unif_rand() : 0.5f;
    samp[i] = std::min((i + delta) * invNSamples, OneMinusEpsilon);
  }
}

void StratifiedSample2D(vec2 *samp, int nx, int ny, random_gen &rng,
                        bool jitter) {
  Float dx = (Float)1 / nx, dy = (Float)1 / ny;
  for (int y = 0; y < ny; ++y) {
    for (int x = 0; x < nx; ++x) {
      Float jx = jitter ? rng.unif_rand() : 0.5f;
      Float jy = jitter ? rng.unif_rand() : 0.5f;
      samp->e[0] = std::min((x + jx) * dx, OneMinusEpsilon);
      samp->e[1] = std::min((y + jy) * dy, OneMinusEpsilon);
      ++samp;
    }
  }
}

void StratifiedSampler::StartPixel(const vec2 &p) {
  random_gen rng;
  // Generate single stratified samples for the pixel
  for (size_t i = 0; i < samples1D.size(); ++i) {
    StratifiedSample1D(&samples1D[i][0], xPixelSamples * yPixelSamples, rng, jitterSamples);
    Shuffle(&samples1D[i][0], xPixelSamples * yPixelSamples, 1, rng);
  }
  for (size_t i = 0; i < samples2D.size(); ++i) {
    StratifiedSample2D(&samples2D[i][0], xPixelSamples, yPixelSamples, rng, jitterSamples);
    Shuffle(&samples2D[i][0], xPixelSamples * yPixelSamples, 1, rng);
  }
  
  // Generate arrays of stratified samples for the pixel
  for (size_t i = 0; i < samples1DArraySizes.size(); ++i) {
    for (int64_t j = 0; j < samplesPerPixel; ++j) {
      int count = samples1DArraySizes[i];
      StratifiedSample1D(&sampleArray1D[i][j * count], count, rng, jitterSamples);
      Shuffle(&sampleArray1D[i][j * count], count, 1, rng);
    }
  }
  for (size_t i = 0; i < samples2DArraySizes.size(); ++i) {
    for (int64_t j = 0; j < samplesPerPixel; ++j) {
      int count = samples2DArraySizes[i];
      LatinHypercube(&sampleArray2D[i][j * count].e[0], count, 2, rng);
    }
    PixelSampler::StartPixel(p);
  }
}

#endif
