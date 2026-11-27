#include "../math/distributions.h"
#include <stdexcept>

Distribution1D::Distribution1D(const Float *f, int n) : func(f, f + n), cdf(n + 1) {
  // Compute integral of step function at $x_i$
  cdf[0] = 0;
  for (int i = 1; i < n + 1; ++i) {
    cdf[i] = cdf[i - 1] + func[i - 1] / n;
  }
  
  // Transform step function integral into CDF
  funcInt = cdf[n];
  if (funcInt == 0) {
    for (int i = 1; i < n + 1; ++i) {
      cdf[i] = Float(i) / Float(n);
    }
  } else {
    for (int i = 1; i < n + 1; ++i) {
      cdf[i] /= funcInt;
    }
  }
}

int Distribution1D::Count() const { 
  return((int)func.size()); 
}

Float Distribution1D::SampleContinuous(Float u, Float *pdf, int *off) const {
  // Find surrounding CDF segments and _offset_
  int offset = FindInterval((int)cdf.size(),
                            [&](int index) { return cdf[index] <= u;});
  if (off) {
    *off = offset;
  }
  // Compute offset along CDF segment
  Float du = u - cdf[offset];
  if ((cdf[offset + 1] - cdf[offset]) > 0) {
    du /= (cdf[offset + 1] - cdf[offset]);
  }
  // Compute PDF for sampled offset
  if (pdf) {
    *pdf = (funcInt > 0) ? func[offset] / funcInt : 0;
  }
  // Return $x\in{}[0,1)$ corresponding to sample
  return (offset + du) / Count();
}

int Distribution1D::SampleDiscrete(Float u, Float *pdf,
                   Float *uRemapped) const {
  // Find surrounding CDF segments and _offset_
  int offset = FindInterval((int)cdf.size(),
                            [&](int index) { return cdf[index] <= u; });
  if (pdf) {
    *pdf = (funcInt > 0) ? func[offset] / (funcInt * Count()) : 0;
  }
  if (uRemapped) {
    *uRemapped = (u - cdf[offset]) / (cdf[offset + 1] - cdf[offset]);
  }
  return offset;
}

Float Distribution1D::DiscretePDF(int index) const {
  return(func[index] / (funcInt * Count()));
}

size_t Distribution1D::GetSize()  {
  return(sizeof(*this) + sizeof(Float)*func.size() + sizeof(Float)*cdf.size());
}



vec2f Distribution2D::SampleContinuous(const vec2f &u, Float *pdf) const {
  Float pdfs[2];
  int v;
  Float d1 = pMarginal->SampleContinuous(u[1], &pdfs[1], &v);
  Float d0 = pConditionalV[v]->SampleContinuous(u[0], &pdfs[0]);
  *pdf = pdfs[0] * pdfs[1];
  return(vec2f(d0, d1));
}

Float Distribution2D::Pdf(const vec2f &p) const {
  int iu = clamp(int(p[0] * pConditionalV[0]->Count()), 0,
                 pConditionalV[0]->Count() - 1);
  int iv = clamp(int(p[1] * pMarginal->Count()), 0, pMarginal->Count() - 1);
  return(pConditionalV[iv]->func[iu] / pMarginal->funcInt);
}

Distribution2D::Distribution2D(const Float *func, int nu, int nv) {
  pConditionalV.reserve(nv);
  for (int v = 0; v < nv; ++v) {
    // Compute conditional sampling distribution for $\tilde{v}$
    pConditionalV.emplace_back(std::unique_ptr<Distribution1D>(new Distribution1D(&func[v * nu], nu)));
  }
  // Compute marginal sampling distribution $p[\tilde{v}]$
  std::vector<Float> marginalFunc;
  marginalFunc.reserve(nv);
  for (int v = 0; v < nv; ++v)
    marginalFunc.push_back(pConditionalV[v]->funcInt);
  pMarginal.reset(new Distribution1D(&marginalFunc[0], nv));
}

size_t Distribution2D::GetSize()  {
  size_t total_size = 0;
  for(size_t i = 0; i < pConditionalV.size(); i++) {
    total_size += pConditionalV[i]->GetSize();
  }
  total_size += pMarginal->GetSize();
  return(sizeof(*this) + total_size);
}

AliasDistribution2D::AliasDistribution2D(const Float *data, int w, int h)
    : width(w), height(h) {
  if (width <= 0 || height <= 0)
    throw std::invalid_argument("Alias distribution dimensions must be positive");
  columns.resize(size_t(width) * height);
  rows.resize(height);
  std::vector<double> row_weights(height), weights(width);
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      double value = data[size_t(y) * width + x];
      if (!(value >= 0) || !std::isfinite(value))
        throw std::invalid_argument("Alias distribution weights must be finite and nonnegative");
      weights[x] = value;
      row_weights[y] += value;
    }
    Build(weights, columns.data() + size_t(y) * width);
  }
  Build(row_weights, rows.data());
}

void AliasDistribution2D::Build(const std::vector<double> &weights, Entry *entries) {
  const size_t count = weights.size();
  double total = 0;
  for (double value : weights)
    total += value;
  if (!std::isfinite(total))
    throw std::invalid_argument("Alias distribution weight sum must be finite");
  std::vector<double> scaled(count);
  std::vector<uint32_t> small, large;
  small.reserve(count);
  large.reserve(count);
  for (size_t i = 0; i < count; ++i) {
    // Empty rows are unreachable unless the entire map is black. A uniform
    // fallback also gives this distribution a defined standalone black-map PDF.
    scaled[i] = total > 0 ? weights[i] / total * count : 1;
    entries[i].alias = i;
    (scaled[i] < 1 ? small : large).push_back(i);
  }
  while (!small.empty() && !large.empty()) {
    uint32_t a = small.back(), b = large.back();
    small.pop_back();
    large.pop_back();
    entries[a].cutoff = std::clamp(Float(scaled[a]), Float(0), Float(1));
    entries[a].alias = b;
    scaled[b] -= 1 - scaled[a];
    (scaled[b] < 1 ? small : large).push_back(b);
  }
  // Report the actual density of the rounded cutoffs, not the ideal input
  // weights. Each column has probability 1/count; density multiplies by count.
  std::fill(scaled.begin(), scaled.end(), 0);
  for (size_t i = 0; i < count; ++i) {
    scaled[i] += entries[i].cutoff;
    scaled[entries[i].alias] += 1 - double(entries[i].cutoff);
  }
  for (size_t i = 0; i < count; ++i)
    entries[i].density = scaled[i];
}

int AliasDistribution2D::Sample(const Entry *entries, int count, Float u, Float &residual) {
  // Double multiplication retains the available input bits for non-power-of-two
  // dimensions. Inputs at 1 are clamped before indexing, including rounded RNGs.
  double target = std::clamp(double(u), 0.0, double(std::nextafter(Float(1), Float(0)))) * count;
  int column = int(target);
  double coin = target - column;
  const auto &entry = entries[column];
  if (coin < entry.cutoff) {
    residual = coin / entry.cutoff;
  } else {
    residual = (coin - entry.cutoff) / (1 - double(entry.cutoff));
    column = entry.alias;
  }
  residual = std::min(residual, std::nextafter(Float(1), Float(0)));
  return column;
}

vec2f AliasDistribution2D::SampleContinuous(const vec2f &u, Float *pdf) const {
  Float x_residual, y_residual;
  int y = Sample(rows.data(), height, u[1], y_residual);
  const Entry *row = columns.data() + size_t(y) * width;
  int x = Sample(row, width, u[0], x_residual);
  // Do not round the returned UV into a neighboring cell (especially a zero
  // weight cell). This matters at exactly representable alias branch endpoints.
  Float uu = (x + double(x_residual)) / width;
  Float vv = (y + double(y_residual)) / height;
  if (double(uu) * width < x)
    uu = std::nextafter(uu, Float(1));
  else if (double(uu) * width >= x + 1)
    uu = std::nextafter(uu, Float(0));
  if (double(vv) * height < y)
    vv = std::nextafter(vv, Float(1));
  else if (double(vv) * height >= y + 1)
    vv = std::nextafter(vv, Float(0));
  *pdf = row[x].density * rows[y].density;
  return vec2f(uu, vv);
}

Float AliasDistribution2D::Pdf(const vec2f &uv) const {
  int x = clamp(int(double(uv[0]) * width), 0, width - 1);
  int y = clamp(int(double(uv[1]) * height), 0, height - 1);
  return columns[size_t(y) * width + x].density * rows[y].density;
}

size_t AliasDistribution2D::GetSize() const {
  return sizeof(*this) + (columns.capacity() + rows.capacity()) * sizeof(Entry);
}
