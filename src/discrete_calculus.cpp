/******************************************************************************/
// Divided differences, discrete derivatives, and discrete integrals 

#include <Rcpp.h>
#include "dspline.h"
using namespace Rcpp; 

// [[Rcpp::export]]
void rcpp_dot_divided_diff(NumericVector f, NumericVector z) {
	int n = f.size();
  for (int i = 0; i < n; i++) {
		for (int j = n-1; j >= i+1; j--) {
			f[j] = (f[j] - f[j-1]) / (z[j] - z[j-i-1]);
		}
	}
	return;
}

// [[Rcpp::export]]
double rcpp_divided_diff(NumericVector f, NumericVector z) {
	NumericVector f_copy = clone(f);
	rcpp_dot_divided_diff(f_copy, z);
	return f_copy[f_copy.size()-1];
}

// [[Rcpp::export]]
NumericVector rcpp_discrete_deriv(NumericVector f, int k, NumericVector xd, NumericVector x) {   
	// Take care of trivial case first
	int nx = x.size();
	if (k == 0) return f[Range(f.size()-nx, f.size()-1)];

	NumericVector a (nx);
	for (int l = 0; l < nx; l++) {
		// Find the largest i such that x_i < x
		int i = std::lower_bound(xd.begin(), xd.end(), x[l]) - xd.begin() - 1;

		// If x <= x_1, then return f(x)
		if (i < 0) a[l] = f[f.size()-nx+l];
	
		// Else assemble centers, evaluations, compute scaled divided difference 
		else {
			int j = std::min(k-1, i);
			NumericVector xd_sub = xd[Range(i-j, i)];
			NumericVector f_sub = f[Range(i-j, i)];
			xd_sub.push_back(x[l]);
			f_sub.push_back(f[f.size()-nx+l]);
			a[l] = rcpp_divided_diff(f_sub, xd_sub) * fact(j+1);
		}
	}
	return a;
}

// [[Rcpp::export]]
NumericVector rcpp_discrete_integ(NumericVector f, int k, NumericVector xd, NumericVector x) {   
	// Take care of trivial case first
	int nx = x.size();
	if (k == 0) return f[Range(f.size()-nx, f.size()-1)];
	
	NumericVector a (nx);
	for (int l = 0; l < nx; l++) {
		// Find the largest i such that x_i < x
		int i = std::lower_bound(xd.begin(), xd.end(), x[l]) - xd.begin() - 1;

		// If x <= x_1, then return f(x)
		if (i < 0) a[l] = f[f.size()-nx+l];
	
		// Else use linear combination form in (48) of Tibshirani (2020)
		else {
			double b, sum = 0;
			// Summands for x_1, ..., x_i
			for (int j = 0; j <= i; j++) {
				b = rcpp_hj_fun(k-1, xd, j, x[l]);
				if (j >= k) b *= (xd[j] - xd[j-k]) / k;
				b *= f[j];
				sum += b;
			}
			// Summand for x
			b = rcpp_hj_fun(k-1, xd, i+1, x[l]);
			if (i+1 >= k) b *= (x[l] - xd[i+1-k]) / k;
			b *= f[f.size()-nx+l];
			sum += b;
			a[l] = sum;
		}
	}
	return a;
}
