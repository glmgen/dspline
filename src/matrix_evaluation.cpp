/******************************************************************************/
// Evaluation of falling factorial and discrete B-spline bases

#include <Rcpp.h>
#include "dspline.h"
using namespace Rcpp; 

// [[Rcpp::export]]
List rcpp_hx_mat(int k, NumericVector xd, NumericVector x, bool di_weighting, IntegerVector col_idx) {
	int n_row = x.size(), n_col = col_idx.size();
	NumericMatrix H (n_row, n_col); 
	for (int i = 0; i < n_row; i++) {
		for (int j = 0; j < n_col; j++) {
			H(i, j) = hxj(k, xd, x[i], col_idx[j]);

			// Discrete integration weighting
			if (di_weighting && col_idx[j] >= k+1) {
				H(i, j) *= (xd[col_idx[j]] - xd[col_idx[j]-k-1]) / (k+1);
			}
		}
	}

	//return H;

	int N = 0;
	IntegerVector i_vec (N);
	IntegerVector j_vec (N);
	NumericVector x_vec (N);

	return List::create(Named("i") = i_vec, Named("j") = j_vec, Named("x") = x_vec);
}

// [[Rcpp::export]]
List rcpp_nx_mat(int k, NumericVector xd, NumericVector x, bool normalized, IntegerVector knot_idx) {
	int N = 0;
	IntegerVector i_vec (N);
	IntegerVector j_vec (N);
	NumericVector x_vec (N);

	return List::create(Named("i") = i_vec, Named("j") = j_vec, Named("x") = x_vec);
}

