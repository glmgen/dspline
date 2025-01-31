/******************************************************************************/
// Interpolation within polynomial and discrete spline spaces

#include <Rcpp.h>
#include "dspline.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rcpp_newton_interp(NumericVector v, NumericVector z, NumericVector x) {
  NumericVector p (x.size());
  for (int i = 0; i < x.size(); i++) {
    p[i] = 0;
    for (int l = 0; l < z.size(); l++) {
      p[i] += divided_diff(v, z, l+1) * newton_poly(z, x[i], l);
    }
  }
  return p;
}

// [[Rcpp::export]]
NumericVector rcpp_dspline_interp(NumericVector v, int k, NumericVector xd, NumericVector x, bool implicit) {
  int nx = x.size();
  NumericVector f (nx);

  // Falling factorial interpolation
  if (!implicit) {
    // Locate each x point among knot points only
    IntegerVector J (nx);
    for (int l = 0 ; l < nx; l++) {
      // Find the smallest i such that x_i \geq x
      J[l] = std::lower_bound(xd.begin() + k + 1, xd.end(), x[l]) - xd.begin();
      J[l] = std::min(J[l], (int) xd.size()-1);
    }

    // Compute basis coefficients a = B^{k+1} Z^{k+1} v
    NumericVector a = rcpp_b_mat_mult(v, k+1, xd, 1, 0, 0);

    // Multiply by falling factorial basis H^{k+1}_x at x
    for (int i = 0; i < nx; i++) {
      f[i] = 0;
      // Pure polynomials
      for (int j = 0; j < k+1; j++) {
        f[i] += hxj(k, xd, x[i], j) * a[j];
      }
      // Piecewise polynomials
      for (int j = k+1; j <= J[i]; j++) {
        if (a[j] != 0) {
          f[i] += hxj(k, xd, x[i], j) * a[j];
        }
      }
    }
  }

  // Implicit form interpolation
  else {
    // Locate each x point among all design points
    IntegerVector J (nx);
    for (int l = 0 ; l < nx; l++) {
      // Find the smallest i such that x_i \geq x
      J[l] = std::lower_bound(xd.begin(), xd.end(), x[l]) - xd.begin();
      J[l] = std::min(J[l], (int) xd.size()-1);
    }

    double s, w;
    int offset;

    for (int i = 0; i < nx; i++) {
      // Trivial case where x is one of the design points
      if (xd[J[i]] == x[i]) {
        f[i] = v[J[i]];
        continue;
      }

      // If J[i] < k, solve for f(x) in f[xd[0], ..., xd[k], x] = 0
      // Else, solve for f(x) in f[xd[J[i]-k], ..., xd[J[i]], x] = 0
      if (J[i] < k) offset = 0;
      else offset = J[i]-k;

      s = 0;
      for (int j = 0; j < k+1; j++) {
        w = 1;
        for (int l = 0; l < k+1; l++) {
          if (l != j) {
            w /= xd[offset + j] - xd[offset + l];
          }
        }
        w /= xd[offset + j] - x[i];
        s -= w * v[offset + j];
      }

      w = 1;
      for (int l = 0; l < k+1; l++) {
        w /= x[i] - xd[offset + l];
      }
      f[i] = s / w;
    }
  }

  return f;
}
