#ifndef TONEMAPH
#define TONEMAPH

#include "mathinline.h"
#include "Rcpp.h"

static const Float A = 0.15;
static const Float B = 0.50;
static const Float C = 0.10;
static const Float D = 0.20;
static const Float E = 0.02;
static const Float F = 0.30;
static const Float W = 11.2;

// [[Rcpp::export]]
Rcpp::List tonemap_image(Rcpp::NumericMatrix routput, Rcpp::NumericMatrix goutput, Rcpp::NumericMatrix boutput, 
                         int toneval);

#endif
