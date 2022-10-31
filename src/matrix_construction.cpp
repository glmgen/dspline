/******************************************************************************/
// Construction of discrete derivative and discrete spline basis matrices

#include <Rcpp.h>
#include "dspline.h"
using namespace Rcpp;

// [[Rcpp::export]]
List rcpp_b_mat(int k, NumericVector xd, bool tf_weighting, IntegerVector row_idx, bool d_only) {
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
  IntegerVector i_vec (N);
  IntegerVector j_vec (N);
  NumericVector x_vec (N);
  int l = 0, j_start, j_end;
  for (int i = 0; i < n_row; i++) {
    j_start = std::max(row_idx[i] - (!d_only * k), 0);
    j_end = row_idx[i] + (d_only * k);
    for (int j = j_start; j <= j_end; j++) {
      i_vec[l] = i;
      j_vec[l] = j;
      if (!d_only) x_vec[l] = bij(k, xd, row_idx[i], j);
      else x_vec[l] = dij(k, xd, row_idx[i], j);

      // Trend filtering weighting (only applies to k >= 0)
      if (tf_weighting && k > 0) {
        if (!d_only && row_idx[i] >= k) {
          x_vec[l] *= (xd[row_idx[i]] - xd[row_idx[i]-k]) / k;
        }
        else if (d_only) {
          x_vec[l] *= (xd[row_idx[i]+k] - xd[row_idx[i]]) / k;
        }
      }
      l++;
    }
  }

  return List::create(Named("i") = i_vec, Named("j") = j_vec, Named("x") = x_vec);
}

// [[Rcpp::export]]
List rcpp_h_mat(int k, NumericVector xd, bool di_weighting, IntegerVector col_idx) {
  // Compute number of nonzero elements for H. We do so by computing nonzeros
  // per column (falling factorial basis vector)
  int n_row = xd.size(), n_col = col_idx.size(), N = n_col * n_row;
  for (int j = 0; j < n_col; j++) N -= col_idx[j];

  // Now construct H in sparse triplet format
  IntegerVector i_vec (N);
  IntegerVector j_vec (N);
  NumericVector x_vec (N);
  int l = 0;
  for (int j = 0; j < n_col; j++) {
    for (int i = col_idx[j]; i < n_row; i++) {
      i_vec[l] = i;
      j_vec[l] = j;
      x_vec[l] = hxj(k, xd, xd[i], col_idx[j]);

      // Discrete integration weighting
      if (di_weighting && col_idx[j] >= k+1) {
        x_vec[l] *= (xd[col_idx[j]] - xd[col_idx[j]-k-1]) / (k+1);
      }
      l++;
    }
  }

  return List::create(Named("i") = i_vec, Named("j") = j_vec, Named("x") = x_vec);
}

// [[Rcpp::export]]
List rcpp_n_mat(int k, NumericVector xd, bool normalized, IntegerVector knot_idx) {
  // Compute number of nonzero elements for N. We do so by computing nonzeros
  // per column (discrete B-spline basis vector)
  int n_row = xd.size() - k - 1, n_col = knot_idx.size(), N = 0;
  for (int j = 0; j < n_col; j++) {
    if (j < k+1) N += knot_idx[j] + 1; // First k+1 functions
    else if (j < n_col-k-1) N += knot_idx[j] - knot_idx[j-k-1] + 1;
    else N += n_row - knot_idx[j-k-1]; // Last k+1 functions
  }

  // Now construct N in sparse triplet format
  IntegerVector i_vec (N);
  IntegerVector j_vec (N);
  NumericVector x_vec (N);
  int l = 0, l_last;
  for (int j = 0; j < n_col; j++) {
    l_last = l; // Save this in case we need it later for the normalization

    // Solve a linear system to get the appropriate basis vector N_j, store it
    // in the triplet list, and update the counter l by reference
    nj(k, xd, knot_idx, j, i_vec, j_vec, x_vec, l);

    if (normalized) {
      double max = *std::max_element(x_vec.begin() + l_last, x_vec.begin() + l);
      for (int p = l_last; p < l; p++) x_vec[p] = x_vec[p] / max;
    }
  }

  return List::create(Named("i") = i_vec, Named("j") = j_vec, Named("x") = x_vec);
}

/******************************************************************************/
// Evaluation of falling factorial basis at arbitrary query points

// [[Rcpp::export]]
List rcpp_h_eval(int k, NumericVector xd, NumericVector x, IntegerVector col_idx) {
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
  IntegerVector i_vec (N);
  IntegerVector j_vec (N);
  NumericVector x_vec (N);
  int l = 0;
  for (int j = 0; j < n_col; j++) {
    for (int i = I[j]; i < n_row; i++) {
      i_vec[l] = i;
      j_vec[l] = j;
      x_vec[l] = hxj(k, xd, x[i], col_idx[j]);
      l++;
    }
  }

  return List::create(Named("i") = i_vec, Named("j") = j_vec, Named("x") = x_vec);
}
