#ifndef DSPLINE_H
#define DSPLINE_H

#include <Rcpp.h>
using namespace Rcpp; 

/******************************************************************************/
// Simple utilities

int fact(int k);

/******************************************************************************/
// Divided differences, discrete derivatives, and discrete integrals 

void rcpp_dot_divided_diff(NumericVector f, NumericVector z);
double rcpp_divided_diff(NumericVector f, NumericVector z);
NumericVector rcpp_discrete_deriv(NumericVector f, int k, NumericVector xd, NumericVector x);
NumericVector rcpp_discrete_integ(NumericVector f, int k, NumericVector xd, NumericVector x);

/******************************************************************************/
// Bases functions

double rcpp_hj_fun(int k, NumericVector xd, int j, double x);

#endif
