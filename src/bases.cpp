/******************************************************************************/
// Bases functions

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
