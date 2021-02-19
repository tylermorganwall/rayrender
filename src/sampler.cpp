#include "sampler.h"
#include "samplerBlueNoise.h"

template <typename T>
void Shuffle(T *samp, int count, int nDimensions, random_gen &rng) {
  for (int i = 0; i < count; ++i) {
    int other = i + rng.UniformUInt32(count - i);
    for (int j = 0; j < nDimensions; ++j)
      std::swap(samp[nDimensions * i + j],
                samp[nDimensions * other + j]);
  }
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


Sampler::~Sampler() {};

void Sampler::StartPixel(unsigned int i, unsigned int j) {
  currentPixelx = i;
  currentPixely = j;
  
  currentPixelSampleIndex = 0;
  array1DOffset = array2DOffset = 0;
}

bool Sampler::SetSampleNumber(size_t sampleNum) {
  array1DOffset = array2DOffset = 0;
  currentPixelSampleIndex = sampleNum;
  return(currentPixelSampleIndex < samplesPerPixel);
}

bool Sampler::StartNextSample() {
  array1DOffset = array2DOffset = 0;
  return(++currentPixelSampleIndex < samplesPerPixel);
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



bool PixelSampler::StartNextSample() {
  current1DDimension = current2DDimension = 0;
  return(Sampler::StartNextSample());
}

bool PixelSampler::SetSampleNumber(size_t sampleNum) {
  current1DDimension = current2DDimension = 0;
  return(Sampler::SetSampleNumber(sampleNum));
}

//This both fetches and increments `current1DDimension`
Float PixelSampler::Get1D() {
  if (current1DDimension < samples1D.size() && 
      currentPixelSampleIndex < samples1D[current1DDimension].size()) {
    return(samples1D[current1DDimension++][currentPixelSampleIndex]);
  } else {
    return(rng.unif_rand());
  }
}

//This both fetches and increments `current2DDimension`
vec2 PixelSampler::Get2D() {
  if (current2DDimension < samples2D.size() && 
      currentPixelSampleIndex < samples2D[current2DDimension].size()) {
    return(samples2D[current2DDimension++][currentPixelSampleIndex]);
  } else {
    return(vec2(rng.unif_rand(), rng.unif_rand()));
  }
}


std::unique_ptr<Sampler> StratifiedSampler::Clone(int seed) {
  StratifiedSampler *ss = new StratifiedSampler(*this);
  ss->rng.SetSequence(seed);
  return(std::unique_ptr<Sampler>(ss));
}

void StratifiedSampler::StartPixel(unsigned int i, unsigned int j) {
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
    for (size_t j = 0; j < samplesPerPixel; ++j) {
      int count = samples1DArraySizes[i];
      StratifiedSample1D(&sampleArray1D[i][j * count], count, rng, jitterSamples);
      Shuffle(&sampleArray1D[i][j * count], count, 1, rng);
    }
  }
  for (size_t i = 0; i < samples2DArraySizes.size(); ++i) {
    for (size_t j = 0; j < samplesPerPixel; ++j) {
      int count = samples2DArraySizes[i];
      LatinHypercube(&sampleArray2D[i][j * count].e[0], count, 2, rng);
    }
  }
  PixelSampler::StartPixel(i,j);
}

std::unique_ptr<Sampler> RandomSampler::Clone(int seed) {
  RandomSampler *ss = new RandomSampler(*this);
  ss->rng.SetSequence(seed);
  return(std::unique_ptr<Sampler>(ss));
}


Float RandomSampler::Get1D() {
  return(rng.unif_rand());
}
vec2 RandomSampler::Get2D() {
  return(vec2(rng.unif_rand(), rng.unif_rand()));
}
bool RandomSampler::StartNextSample() {
  return(true);
}

bool RandomSampler::SetSampleNumber(size_t) {
  return(true);
}


static inline float sobol_calc_single(unsigned long long  i, unsigned int dim, unsigned int scramble) {
  return(spacefillr::sobol_owen_single(i, dim, scramble));
} 

static inline vec2 sobol_calc_double(unsigned long long  i, unsigned int dim, unsigned int scramble) {
  return(vec2(spacefillr::sobol_owen_single(i, dim, scramble), 
              spacefillr::sobol_owen_single(i, dim+1, scramble)));
}

static inline double sobol_calc_single_bluenoise(int  x, int  y, int i, int dim) {
  return(spacefillr::samplerBlueNoise(x,y, i, dim)); 
}

static inline vec2 sobol_calc_double_bluenoise(int  x, int  y, int i, int dim) {
  return(vec2(spacefillr::samplerBlueNoise(x,y, i, dim), 
              spacefillr::samplerBlueNoise(x,y, i, dim+1)));
}

SobolSampler::SobolSampler(unsigned int xPixelSamples, unsigned int yPixelSamples,
                           unsigned int maxSamples,
             random_gen& rng) : PixelSampler(1000000000,0,rng),
             xPixelSamples(xPixelSamples), yPixelSamples(yPixelSamples),
             current1Dsample(0), current2Dsample(0)  {
  pixelseed = rng.unif_rand()*std::numeric_limits<unsigned int>::max();
}

void SobolSampler::StartPixel( unsigned int i,  unsigned int j) { 
  currentPixelx = i;
  currentPixely = j;
  current2DDimension = 0;
  current1DDimension = 0;
}

#include "RcppThread.h"

Float SobolSampler::Get1D() {
  double temp = sobol_calc_single(current1Dsample,
                                  2,
                                  pixelseed + current1DDimension
                                  );
  current1DDimension += 1;
  return(temp);
}


vec2 SobolSampler::Get2D() {
  vec2 temp = sobol_calc_double(current2Dsample,
                                0,
                                pixelseed + current2DDimension
                                );
  current2DDimension += 1;
  return(temp);
}


bool SobolSampler::StartNextSample() {
  current1Dsample++;
  current2Dsample++;
  current1DDimension = current2DDimension = 0;
  return(true);
}
bool SobolSampler::SetSampleNumber(size_t) {
  return(true);
}

std::unique_ptr<Sampler> SobolSampler::Clone(int seed) {
  SobolSampler *ss = new SobolSampler(*this);
  ss->rng.SetSequence(seed);
  return(std::unique_ptr<Sampler>(ss));
}


Float SobolBlueNoiseSampler::Get1D() {
  double temp = sobol_calc_single_bluenoise(currentPixelx,
                                            currentPixely,
                                            current1Dsample,
                                            current1DDimension);
  
  current1DDimension += 1;
  return(temp);
}

vec2 SobolBlueNoiseSampler::Get2D() {
  vec2 temp = sobol_calc_double_bluenoise(currentPixelx,
                                          currentPixely,
                                          current2Dsample,
                                          current2DDimension);
  current2DDimension += 2;
  return(temp);
}

SobolBlueNoiseSampler::SobolBlueNoiseSampler(random_gen& rng)
  : PixelSampler(1000000000,0,rng),  current1Dsample(0), current2Dsample(0)  {}


void SobolBlueNoiseSampler::StartPixel( unsigned int i,  unsigned int j) {
  currentPixelx = i;
  currentPixely = j;
  current2DDimension = 0;
  current1DDimension = 0;
}

bool SobolBlueNoiseSampler::StartNextSample() {
  current1Dsample++;
  current2Dsample++;
  current1DDimension = current2DDimension = 0;
  return(true);
}
bool SobolBlueNoiseSampler::SetSampleNumber(size_t) {
  return(true);
}

std::unique_ptr<Sampler> SobolBlueNoiseSampler::Clone(int seed) {
  SobolBlueNoiseSampler *ss = new SobolBlueNoiseSampler(*this);
  ss->rng.SetSequence(seed);
  return(std::unique_ptr<Sampler>(ss));
}
