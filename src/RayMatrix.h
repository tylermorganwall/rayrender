#ifndef RAYMATRIXH
#define RAYMATRIXH

#include "Rcpp.h"

class RayMatrix {
public: 
  RayMatrix();
  RayMatrix(unsigned int _rows, unsigned int _cols, float start_value = 0);
  inline float& operator()(unsigned int x, unsigned int y) { 
    return data[x + nrow * y]; 
  }
  inline void add_one(unsigned int x, unsigned int y) { 
    data[x + nrow * y] += 1;
  }
  unsigned int rows() {return(nrow);}
  unsigned int cols() {return(ncol);}
  unsigned int nrows() {return(nrow);}
  unsigned int ncols() {return(ncol);}
  unsigned int length() {return(ncol*nrow);}
  unsigned int size() {return(ncol*nrow);}
  
  float* begin() {return(&data[0]);}
  float* end() {return(&data[0]+data.size());}
  
  Rcpp::NumericMatrix ConvertRcpp();
  
  std::vector<float> data;
  unsigned int nrow;
  unsigned int ncol;
};




#endif