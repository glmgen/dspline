/******************************************************************************/
// Multiplication by discrete derivative and discrete spline basis matrices

#include <Rcpp.h>
#include "dspline.h"

// [[Rcpp::interfaces(r, cpp)]]

using namespace Rcpp;

// [[Rcpp::export]]
void rcpp_dot_b_mat_mult(NumericVector v, int k, NumericVector xd, bool tf_weighting, bool transpose, bool inverse) {
  if (!transpose && !inverse) {
    for (int i = 1; i <= k; i++) {
      Diff(v, i);
      InvGapWeight(v, i, xd);
    }
    if (tf_weighting) GapWeight(v, k, xd);
  }
  else if (transpose && !inverse) {
    if (tf_weighting) GapWeight(v, k, xd);
    for (int i = k; i >= 1; i--) {
      InvGapWeight(v, i, xd);
      RevDiff(v, i-1);
    }
  }
  else if (!transpose && inverse) {
    if (tf_weighting) InvGapWeight(v, k, xd);
    for (int i = k; i >= 1; i--) {
      GapWeight(v, i, xd);
      CumSum(v, i-1);
    }
  }
  else {
    for (int i = 1; i <= k; i++) {
      RevCumSum(v, i-1);
      GapWeight(v, i, xd);
    }
    if (tf_weighting) InvGapWeight(v, k, xd);
  }
}

// [[Rcpp::export]]
void rcpp_dot_h_mat_mult(NumericVector v, int k, NumericVector xd, bool di_weighting, bool transpose, bool inverse) {
  if (!transpose && !inverse) {
    if (di_weighting) GapWeight(v, k+1, xd);
    for (int i = k; i >= 0; i--) {
      CumSum(v, i);
      if (i != 0) GapWeight(v, i, xd);
    }
  }
  else if (transpose && !inverse) {
    for (int i = 0; i <= k; i++) {
      if (i != 0) GapWeight(v, i, xd);
      RevCumSum(v, i);
    }
    if (di_weighting) GapWeight(v, k+1, xd);
  }
  else if (!transpose && inverse) {
    for (int i = 0; i <= k; i++) {
      if (i != 0) InvGapWeight(v, i, xd);
      Diff(v, i+1);
    }
    if (di_weighting) InvGapWeight(v, k+1, xd);
  }
  else {
    if (di_weighting) InvGapWeight(v, k+1, xd);
    for (int i = k; i >= 0; i--) {
      RevDiff(v, i);
      if (i != 0) InvGapWeight(v, i, xd);
    }
  }
}

// [[Rcpp::export]]
NumericVector rcpp_d_mat_mult(NumericVector v, int k, NumericVector xd, bool tf_weighting, bool transpose) {
  if (!transpose) {
    NumericVector v_copy = clone(v);
    rcpp_dot_b_mat_mult(v_copy, k, xd, tf_weighting, transpose, 0);
    return v_copy[Range(k, v_copy.size()-1)];
  }
  else {
    int n = v.size() + k;
    NumericVector v_copy (n);
    for (int i = 0; i < k; i++) v_copy[i] = 0;
    for (int i = k; i < n; i++) v_copy[i] = v[i-k];
    rcpp_dot_b_mat_mult(v_copy, k, xd, tf_weighting, transpose, 0);
    return v_copy;
  }
}

// [[Rcpp::export]]
NumericVector rcpp_b_mat_mult(NumericVector v, int k, NumericVector xd, bool tf_weighting, bool transpose, bool inverse) {
  NumericVector v_copy = clone(v);
  rcpp_dot_b_mat_mult(v_copy, k, xd, tf_weighting, transpose, inverse);
  return v_copy;
}

// [[Rcpp::export]]
NumericVector rcpp_h_mat_mult(NumericVector v, int k, NumericVector xd, bool di_weighting, bool transpose, bool inverse) {
  NumericVector v_copy = clone(v);
  rcpp_dot_h_mat_mult(v_copy, k, xd, di_weighting, transpose, inverse);
  return v_copy;
}
