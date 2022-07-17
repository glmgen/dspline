#ifndef DSPLINE_H
#define DSPLINE_H

#include <Rcpp.h>
using namespace Rcpp; 

/******************************************************************************/
// Simple utilities

int fact(int k);
void CumSum(NumericVector v, int i);
void RevCumSum(NumericVector v, int i);
void Diff(NumericVector v, int i);
void RevDiff(NumericVector v, int i);
void GapWeight(NumericVector v, int i, NumericVector xd);
void InvGapWeight(NumericVector v, int i, NumericVector xd);

/******************************************************************************/
// Divided differences, discrete derivatives, and discrete integrals 
void dot_divided_diff(NumericVector f, NumericVector z, int n);
double divided_diff(NumericVector f, NumericVector z, int n);
void rcpp_dot_divided_diff(NumericVector f, NumericVector z);
double rcpp_divided_diff(NumericVector f, NumericVector z);
NumericVector rcpp_discrete_deriv(NumericVector f, int k, NumericVector xd, NumericVector x);
NumericVector rcpp_discrete_integ(NumericVector f, int k, NumericVector xd, NumericVector x);

/******************************************************************************/
// Multiplication by discrete derivative and discrete spline basis matrices

void rcpp_dot_b_mat_mult(NumericVector v, int k, NumericVector xd, bool tf_weighting, bool transpose, bool inverse);
void rcpp_dot_h_mat_mult(NumericVector v, int k, NumericVector xd, bool di_weighting, bool transpose, bool inverse);
NumericVector rcpp_d_mat_mult(NumericVector v, int k, NumericVector xd, bool tf_weighting, bool transpose);
NumericVector rcpp_b_mat_mult(NumericVector v, int k, NumericVector xd, bool tf_weighting, bool transpose, bool inverse);
NumericVector rcpp_h_mat_mult(NumericVector v, int k, NumericVector xd, bool di_weighting, bool transpose, bool inverse);

/******************************************************************************/
// Construction of discrete derivative and discrete spline basis matrices

List rcpp_b_mat(int k, NumericVector xd, bool tf_weighting, IntegerVector row_idx, bool d_only);
List rcpp_h_mat(int k, NumericVector xd, bool di_weighting, IntegerVector col_idx);
List rcpp_n_mat(int k, NumericVector xd, bool normalized, IntegerVector knot_idx);

/******************************************************************************/
// Evaluation of falling factorial at arbitrary query points

List rcpp_h_eval(int k, NumericVector xd, NumericVector x, IntegerVector col_idx);

/******************************************************************************/
// Basis and weight functions 

double newton_poly(NumericVector z, double x, int n);
double dij(int k, NumericVector xd, int i, int j);
double bij(int k, NumericVector xd, int i, int j);
double hxj(int k, NumericVector xd, double x, int j);

/******************************************************************************/
// Workhorse for computing discrete B-spline evaluations

void nj(int k, NumericVector xd, IntegerVector knot_idx, int j_col, IntegerVector i_vec, IntegerVector j_vec, NumericVector x_vec, int& l);

/******************************************************************************/
// Interpolation within polynomial and discrete spline spaces

NumericVector rcpp_newton_interp(NumericVector v, NumericVector z, NumericVector x);
NumericVector rcpp_dspline_interp(NumericVector v, int k, NumericVector xd, NumericVector x, bool implicit);

/******************************************************************************/

// Sparse matrix support?
// Consider just copying this header file, albeit for read-only access:
// https://github.com/zdebruine/RcppSparse/blob/main/inst/include/RcppSparse.h

#endif
