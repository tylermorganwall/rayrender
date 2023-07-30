#include "RayMatrix.h"

RayMatrix::RayMatrix() : ncol(0), nrow(0) {
  data.resize(0);
};

RayMatrix::RayMatrix(unsigned int _rows, unsigned int _cols, float start_value) : nrow(_rows), ncol(_cols) {
  data.resize(nrow*ncol);
  std::fill(data.begin(), data.end(), start_value);
};

Rcpp::NumericMatrix RayMatrix::ConvertRcpp() {
  Rcpp::NumericMatrix out(nrow,ncol);
  for(unsigned int i = 0; i < nrow; i++) {
    for(unsigned int j = 0; j < ncol; j++) {
      out(i,j) = (double)data[i + nrow * j];
    }
  }
  return(out);
}