/******************************************************************************/
// Construction of discrete derivative and discrete spline basis matrices

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Sparse>
#include "dspline.h"

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::depends(RcppEigen)]]

typedef Eigen::Triplet<double> T;

using namespace Rcpp;
using Eigen::SparseMatrix;
using Eigen::VectorXd;

#define DOUBLE_INF std::numeric_limits<double>::infinity()

// Concatenate two Vectors
template<typename V>
V concat(V x, V y) {
  int n_x = x.size();
  int n_y = y.size();
  V xy(n_x+n_y);
  xy[Rcpp::Range(0, n_x-1)] = x;
  xy[Rcpp::Range(n_x, n_x+n_y-1)] = y;
  return xy;
}


// [[Rcpp::export]]
Eigen::SparseMatrix<double> rcpp_b_mat(int k, NumericVector xd, bool tf_weighting, IntegerVector row_idx, bool d_only) {
  // Compute number of nonzero elements for D. We do so by computing nonzeros
  // per row (discrete derivative vector). Start by guessing k+1 nonzeros per
  // row, then if needed, adjust for rows that give lower order derivatives
  int  n_row = row_idx.size(), N = n_row * (k+1);
  if (!d_only) {
    for (int i = 0; i < n_row; i++) {
      if (row_idx[i] < k) N += (row_idx[i]+1 - k-1);
    }
  }

  // Now construct D in sparse triplet format
  std::vector<T> b_list;
  b_list.reserve(N);
  int j_start, j_end;
  double x;
  for (int i = 0; i < n_row; i++) {
    j_start = std::max(row_idx[i] - (!d_only * k), 0);
    j_end = row_idx[i] + (d_only * k);
    for (int j = j_start; j <= j_end; j++) {
      if (!d_only) x = bij(k, xd, row_idx[i], j);
      else x = dij(k, xd, row_idx[i], j);

      // Trend filtering weighting (only applies to k >= 0)
      if (tf_weighting && k > 0) {
        if (!d_only && row_idx[i] >= k) {
          x *= (xd[row_idx[i]] - xd[row_idx[i]-k]) / k;
        }
        else if (d_only) {
          x *= (xd[row_idx[i]+k] - xd[row_idx[i]]) / k;
        }
      }
      b_list.push_back(T(i, j, x));
    }
  }
  SparseMatrix<double> b_mat(n_row, xd.size());
  b_mat.setFromTriplets(b_list.begin(), b_list.end());
  return b_mat;
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> rcpp_h_mat(int k, NumericVector xd, bool di_weighting, IntegerVector col_idx) {
  // Compute number of nonzero elements for H. We do so by computing nonzeros
  // per column (falling factorial basis vector)
  int n_row = xd.size(), n_col = col_idx.size(), N = n_col * n_row;
  for (int j = 0; j < n_col; j++) N -= col_idx[j];

  // Now construct H in sparse triplet format
  std::vector<T> h_list;
  h_list.reserve(N);
  double x;
  for (int j = 0; j < n_col; j++) {
    for (int i = col_idx[j]; i < n_row; i++) {
      x = hxj(k, xd, xd[i], col_idx[j]);

      // Discrete integration weighting
      if (di_weighting && col_idx[j] >= k+1) {
        x *= (xd[col_idx[j]] - xd[col_idx[j]-k-1]) / (k+1);
      }
      h_list.push_back(T(i, j, x));
    }
  }

  SparseMatrix<double> h_mat(n_row, n_col);
  h_mat.setFromTriplets(h_list.begin(), h_list.end());
  return h_mat;
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> rcpp_n_mat(int k, NumericVector xd, bool normalized, IntegerVector knot_idx) {
  // Add k+1 boundary knots (cf. Tibshirani, 2020, Section 8.2)
  int n = xd.size();
  int n_knots = knot_idx.size();
  double max_diff = Rcpp::max(xd[Rcpp::Range(1, n-1)]-xd[Rcpp::Range(0, n-2)]);
  NumericVector extra_knots = Rcpp::cumsum(Rcpp::rep(max_diff, k+1));
  extra_knots = extra_knots + Rcpp::max(xd);
  NumericVector ext_xd = concat<NumericVector>(xd, extra_knots);
  IntegerVector ext_knot_idx = concat<IntegerVector>(knot_idx, Rcpp::Range(n-1, n+k-1));

  // Compute number of nonzero elements for N. We do so by computing nonzeros
  // per column (discrete B-spline basis vector)
  int n_row = n, n_col = n_knots+k+1, N = 0;
  for (int j = 0; j < n_col; j++) {
    if (j < k+1) N += ext_knot_idx[j] + 1; // First k+1 functions
    else if (j < n_col-k-1) N += ext_knot_idx[j] - ext_knot_idx[j-k-1] + 1;
    else N += n_row - ext_knot_idx[j-k-1]; // Last k+1 functions
  }

  // Now construct N in sparse triplet format
  std::vector<T> n_list;
  n_list.reserve(N);
  VectorXd max_vals(n_col);
  for (int j = 0; j < n_col; j++) {
    // Solve a linear system to get the appropriate basis vector N_j, store it
    // in the triplet list
    nj(k, ext_xd, ext_knot_idx, j, n_list, max_vals);
  }

  SparseMatrix<double> n_mat(n_row, n_col);
  n_mat.setFromTriplets(n_list.begin(), n_list.end());
  // Rescale columns to have largest value 1 if necessary
  if (normalized) {
    // Canonical matrix-vector broadcasting in Eigen is done via arrays (see,
    // e.g., https://eigen.tuxfamily.org/dox/group__TutorialArrayClass.html).
    // However, sparse matrices cannot be converted to arrays, so here we
    // perform column scaling via right-multiplication by a diagonal matrix.
    n_mat = n_mat * max_vals.cwiseInverse().asDiagonal();
  }
  return n_mat;
}

/******************************************************************************/
// Evaluation of falling factorial basis at arbitrary query points

// [[Rcpp::export]]
Eigen::SparseMatrix<double> rcpp_h_eval(int k, NumericVector xd, NumericVector x, IntegerVector col_idx) {
  // Compute number of nonzero elements for H. We do so by computing nonzeros
  // per column (falling factorial basis vector). Along the way we compute the
  // index of the smallest nonzero evaluation per column
  int n_row = x.size(), n_col = col_idx.size(), N = n_row * n_col;
  IntegerVector I (n_col);
  for (int j = 0; j < n_col; j++) {
    if (col_idx[j] < k+1) {
      // Evaluate all entries for polynomial basis functions
      I[j] = 0;
    } else {
      // Find the index of smallest nonzero query evaluation for h_j
      I[j] = std::upper_bound(x.begin(), x.end(), xd[col_idx[j]-1]) - x.begin();
    }
    N -= I[j];
  }

  // Now construct H in sparse triplet format
  std::vector<T> h_list;
  h_list.reserve(N);
  for (int j = 0; j < n_col; j++) {
    for (int i = I[j]; i < n_row; i++) {
      h_list.push_back(T(i, j, hxj(k, xd, x[i], col_idx[j])));
    }
  }

  SparseMatrix<double> h_mat(n_row, n_col);
  h_mat.setFromTriplets(h_list.begin(), h_list.end());
  return h_mat;
}

/******************************************************************************/
// Evaluation of discrete spline basis at arbitrary query points

// [[Rcpp::export]]
Eigen::SparseMatrix<double> rcpp_n_eval(int k, NumericVector xd, NumericVector x, bool normalized, IntegerVector knot_idx) {
  Eigen::SparseMatrix<double> n_mat = rcpp_n_mat(k, xd, normalized, knot_idx);
  return rcpp_n_eval_precomputed(k, xd, x, knot_idx, n_mat);
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> rcpp_n_eval_precomputed(int k, NumericVector xd, NumericVector x, IntegerVector knot_idx, Eigen::SparseMatrix<double> n_mat) {
  int n_row = x.size();
  int n_col = n_mat.cols();
  std::vector<T> n_list;
  // n_list.reserve(N);
  for (int j = 0; j < n_col; j++) {
    double lower = (j >= k+1) ? xd[knot_idx[j-k-1]] : -DOUBLE_INF;
    double upper = (j <= n_col-k-2) ? xd[knot_idx[j]] : DOUBLE_INF;
    LogicalVector mask = ((lower <= x) & (x <= upper));
    // Compiler doesn't understand that a Rcpp::Range can be coerced to
    // IntegerVector and refuses to allow subsetting
    IntegerVector tmp = Rcpp::Range(0, mask.size()-1);
    IntegerVector I = tmp[mask];
    Eigen::VectorXd vx = n_mat.col(j);
    NumericVector vals = rcpp_dspline_interp(Rcpp::wrap(vx), k, xd, x[I], true);
    for (int i = 0; i < I.size(); i++) {
      n_list.push_back(T(I[i], j, vals[i]));
    }
  }

  SparseMatrix<double> n_evals(n_row, n_col);
  n_evals.setFromTriplets(n_list.begin(), n_list.end());
  return n_evals;
}
