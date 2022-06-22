/******************************************************************************/
// Basis and weight functions

#include <Rcpp.h>
#include "dspline.h"
using namespace Rcpp; 

// [[Rcpp::export]]
double rcpp_hj_fun(int k, NumericVector xd, int j, double x) {
	double h = 1; 
	if (j <= k) {
		for (int l = 0; l < j; l++) {
			h *= (x - xd[l]) / (l+1);
		}
	}
	if (j > k) {
		if (x <= xd[j-1]) h = 0;
		else {
			for (int l = 0; l < k; l++) {
				h *= (x - xd[j-k+l]) / (l+1);
			}
		}
	}
	return h;
}

// Explicit representation in (40) of Tibshirani (2020)
double bij(int k, NumericVector xd, int i, int j) {		
	if (i < k) {
		double w = fact(i);
		for (int l = 0; l < i+1; l++) {
			if (l != j) {
				w *= 1 / (xd[j] - xd[l]);
			}
		}
		return w;
	}
	else return dij(k, xd, i-k, j);
}

// Explicit representation in (40) of Tibshirani (2020)
double dij(int k, NumericVector xd, int i, int j) {
	double w = fact(k);
	for (int l = i; l < i+k+1; l++) {
		if (l != j) {
			w *= 1 / (xd[j] - xd[l]);
		}
	}
	return w;
}
