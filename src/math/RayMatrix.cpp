#include "../math/RayMatrix.h"
#include "../math/assert.h"

RayMatrix::RayMatrix() : nrow(0), ncol(0), channels(0) {
  data.resize(0);
};

RayMatrix::RayMatrix(unsigned int _rows, unsigned int _cols, unsigned int _channels, 
float start_value) : nrow(_rows), ncol(_cols), channels(_channels) {
  data.resize(channels*nrow*ncol);
  std::fill(data.begin(), data.end(), start_value);
};

Rcpp::NumericMatrix RayMatrix::ConvertRcpp() const  {
  Rcpp::NumericMatrix out(channels * nrow,ncol);
  for(unsigned int i = 0; i < nrow; i++) {
    for(unsigned int j = 0; j < ncol; j++) {
      for(unsigned int k = 0; k < channels; k++) {
        out(k + channels * i,j) = (double)data[k + channels * i + channels * nrow * j];
      }
    }
  }
  return(out);
}

Rcpp::NumericMatrix RayMatrix::ConvertRcpp(unsigned int channel) const {
  Rcpp::NumericMatrix out(nrow,ncol);
  ASSERT(channel < channels && channel >= 0);
  for(unsigned int i = 0; i < nrow; i++) {
    for(unsigned int j = 0; j < ncol; j++) {
      out(i,j) = (double)data[channel + channels * i + channels * nrow * j];
    }
  }
  return(out);
}