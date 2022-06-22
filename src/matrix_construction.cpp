/******************************************************************************/
// Construction of discrete derivative and discrete spline basis matrices

#include <Rcpp.h>
#include "dspline.h"
using namespace Rcpp; 

// [[Rcpp::export]]
List rcpp_d_mat(int k, NumericVector xd, bool tf_weighting, IntegerVector row_idx, bool ext) {
	// Compute number of nonzero elements for D. We do so by computing nonzeros per
	// row. Start by guessing k+1 nonzeros per row, then if needed, adjust for rows 
	// that give lower order derivatives
	int	n_row = row_idx.size(), N = n_row * (k+1);
	if (ext) {
		for (int i = 0; i < n_row; i++) { 
			if (row_idx[i] < k) N += (row_idx[i]+1 - k-1);
		}
	}

	// Now construct D in sparse triplet format
	IntegerVector ivec (N);
	IntegerVector jvec (N);
	NumericVector xvec (N);
	int l = 0, j_start, j_end;
	for (int i = 0; i < n_row; i++) {
		j_start = std::max(row_idx[i] - ext * k, 0);
		j_end = row_idx[i] + (!ext * k);
		for (int j = j_start; j <= j_end; j++) {
			ivec[l] = i;
			jvec[l] = j;
			if (ext) xvec[l] = bij(k, xd, row_idx[i], j);
			else xvec[l] = dij(k, xd, row_idx[i], j);

			// Trend filtering weighting (only applies to k >= 0)
			if (tf_weighting && k > 0) {
				if (ext && row_idx[i] >= k) {
					xvec[l] *= (xd[row_idx[i]] - xd[row_idx[i]-k]) / k; 
				}
				else if (!ext) {
					xvec[l] *= (xd[row_idx[i]+k] - xd[row_idx[i]]) / k; 
				}
			}
			l++;
		}
	}

	return List::create(Named("i") = ivec, Named("j") = jvec, Named("x") = xvec);
}

// [[Rcpp::export]]
List rcpp_h_mat(int k, NumericVector xd, bool di_weighting, IntegerVector col_idx) {
	// Compute number of nonzero elements for H. We do so by computing nonzeros 
	// per column
	int n_row = xd.size(), n_col = col_idx.size(), N = n_col * n_row;
	for (int j = 0; j < n_col; j++) N -= col_idx[j];

	// Now construct H in sparse triplet format
	IntegerVector ivec (N);
	IntegerVector jvec (N);
	NumericVector xvec (N);
	int l = 0;
	for (int j = 0; j < n_col; j++) {
		for (int i = col_idx[j]; i < n_row; i++) {
			ivec[l] = i;
			jvec[l] = j;
			xvec[l] = rcpp_hj_fun(k, xd, col_idx[j], xd[i]);

			// Discrete integration weighting
			if (di_weighting && col_idx[j] >= k+1) {
				xvec[l] *= (xd[col_idx[j]] - xd[col_idx[j]-k-1]) / (k+1);
			}
			l++;
		}
	}

	return List::create(Named("i") = ivec, Named("j") = jvec, Named("x") = xvec);
}

// TODO: construct this in sparse triplet format? Probably worth sorting x ahead of time ...
// [[Rcpp::export]]
NumericMatrix rcpp_hx_mat(int k, NumericVector xd, NumericVector x, bool di_weighting, IntegerVector col_idx) {
	int n_row = x.size(), n_col = col_idx.size();
	NumericMatrix H (n_row, n_col); 
	for (int i = 0; i < n_row; i++) {
		for (int j = 0; j < n_col; j++) {
			H(i, j) = rcpp_hj_fun(k, xd, col_idx[j], x[i]);

			// Discrete integration weighting
			if (di_weighting && col_idx[j] >= k+1) {
				H(i, j) *= (xd[col_idx[j]] - xd[col_idx[j]-k-1]) / (k+1);
			}
		}
	}

	return H;
}
