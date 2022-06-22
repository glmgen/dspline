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
// Construction of discrete derivative and discrete spline basis matrices

List rcpp_d_mat(int k, NumericVector xd, bool tf_weighting, IntegerVector row_idx, bool ext);
List rcpp_h_mat(int k, NumericVector xd, bool di_weighting, IntegerVector col_idx);
List rcpp_n_mat(int k, NumericVector xd, bool normalized, IntegerVector knot_idx);
NumericMatrix rcpp_hx_mat(int k, NumericVector xd, NumericVector x, bool di_weighting, IntegerVector col_idx);
List rcpp_nx_mat(int k, NumericVector xd, NumericVector x, bool normalized, IntegerVector knot_idx);

/******************************************************************************/
// Basis and weight functions 

double rcpp_hj_fun(int k, NumericVector xd, int j, double x);
double bij(int k, NumericVector xd, int i, int j);
double dij(int k, NumericVector xd, int i, int j);

/******************************************************************************/

// Sparse matrix support?
// Consider just copying this header file, albeit for read-only access:
// https://github.com/zdebruine/RcppSparse/blob/main/inst/include/RcppSparse.h

#endif
