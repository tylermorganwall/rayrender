#ifndef RAYMATRIXH
#define RAYMATRIXH

#include "Rcpp.h"

//Row major
class RayMatrix {
public: 
  RayMatrix();
  RayMatrix(unsigned int _rows, unsigned int _cols, unsigned int channels, float start_value = 0);
  inline float& operator()(unsigned int x, unsigned int y, unsigned int channel) { 
    return data[channel + channels * x + channels * nrow * y]; 
  }
  inline void add_one(unsigned int x, unsigned int y, unsigned int channel) { 
    data[channel + channels * x + channels * nrow * y] += 1;
  }
  inline void print() {
    for(unsigned int k = 0; k < channels; k++) {
      Rcpp::Rcout << "Channel: " << k << "\n\n";
      for(unsigned int i = 0; i < nrow; i++) {
        for(unsigned int j = 0; j < ncol; j++) {
          Rcpp::Rcout << (*this)(i,j,k) << " ";
        }
        Rcpp::Rcout << "\n";
      }
    }
  }
  unsigned int rows() {return(nrow);}
  unsigned int cols() {return(ncol);}
  unsigned int nrows() {return(nrow);}
  unsigned int ncols() {return(ncol);}
  unsigned int length() {return(ncol*nrow);}
  unsigned int size() {return(ncol*nrow);}
  void reset() {
    std::fill(data.begin(), data.end(), 0);
  }

  float* begin() {return(&data[0]);}
  float* end() {return(&data[0]+data.size());}
  
  Rcpp::NumericMatrix ConvertRcpp() const;
  Rcpp::NumericMatrix ConvertRcpp(unsigned int channel) const;

  std::vector<float> data;
  unsigned int nrow;
  unsigned int ncol;
  unsigned int channels;
};




#endif