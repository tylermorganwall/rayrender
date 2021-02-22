#ifndef TONEMAPH
#define TONEMAPH

#include "mathinline.h"
#include "Rcpp.h"

static inline Float reinhard(Float color, Float sum);

static Float A = 0.15;
static Float B = 0.50;
static Float C = 0.10;
static Float D = 0.20;
static Float E = 0.02;
static Float F = 0.30;
static Float W = 11.2;

static Float uncharted(Float x);
static Float hable(Float color);
static Float hbd(Float color);

// [[Rcpp::export]]
Rcpp::List tonemap_image(int nx, int ny, 
                         Rcpp::NumericMatrix routput, Rcpp::NumericMatrix goutput, Rcpp::NumericMatrix boutput, 
                         int toneval);

#endif
