/******************************************************************************/
// Basis and weight functions

#include <Rcpp.h>
#include "dspline.h"
using namespace Rcpp;

double newton_poly(NumericVector z, double x, int n) {
  double p = 1;
  for (int i = 0; i < n; i++) {
    p *= (x - z[i]);
  }
  return p;
}

double dij(int k, NumericVector xd, int i, int j) {
  if (i <= j && j <= i+k) {
    double w = fact(k);
    for (int l = i; l < i+k+1; l++) {
      if (l != j) {
        w /= xd[j] - xd[l];
      }
    }
    return w;
  }
  else return 0;
}

double bij(int k, NumericVector xd, int i, int j) {
  if (i < k) {
    if (j <= i) {
      double w = fact(i);
      for (int l = 0; l < i+1; l++) {
        if (l != j) {
          w /= xd[j] - xd[l];
        }
      }
      return w;
    }
    else return 0;
  }
  else return dij(k, xd, i-k, j);
}

double hxj(int k, NumericVector xd, double x, int j) {
  if (j <= k) {
    double h = 1;
    for (int l = 0; l < j; l++) {
      h *= (x - xd[l]) / (l+1);
    }
    return h;
  }
  else {
    if (x <= xd[j-1]) return 0;
    else {
      double h = 1;
      for (int l = 0; l < k; l++) {
        h *= (x - xd[j-k+l]) / (l+1);
      }
      return h;
    }
  }
}

/******************************************************************************/
// Workhorse for computing discrete B-spline evaluations

#include <RcppEigen.h>
#include <Eigen/Sparse>
typedef Eigen::Triplet<double> T;
typedef Eigen::COLAMDOrdering<int> Ord;
using Eigen::SparseMatrix;
using Eigen::SparseQR;
using Eigen::VectorXd;

void nj(int k, NumericVector xd, IntegerVector knot_idx, int j_col, std::vector<T>& n_list, Eigen::VectorXd& max_vals) {
  int n_row = xd.size() - k - 1;
  int t = knot_idx[j_col]; // Last knot
  VectorXd x;
  max_vals[j_col] = 1; // First element added will be 1.

  // First k+1 basis vectors
  if (j_col < k+1) {
    // Add entry where the basis vector equals 1
    n_list.push_back(T(j_col, j_col, 1));

    // Add remaining entries only if t-k >= j_col+1
    if (t-k >= j_col+1) {
      // Compute row indices i such that i+k is NOT one of first j_col-1 knots
      IntegerVector I0 = Range(k, t-1);
      IntegerVector I (t-k-j_col);
      std::set_difference(I0.begin(), I0.end(), knot_idx.begin(), knot_idx.begin()+j_col, I.begin());
      for (int p = 0; p < t-k-j_col; p++) I[p] -= k;

      // Compute column indices that correspond to unknown indices
      IntegerVector J = Range(j_col+1, t-k);

      // Step 1a: form appropriate matrix D
      std::vector<T> d_list;
      int q_start = 0, q_end;
      d_list.reserve((t-k-j_col) * (k+2));
      for (int p = 0; p < t-k-j_col; p++) {
        // We know D is banded, find relevant start/end points
        while(J[q_start] < I[p]) q_start++;
        q_end = std::min(q_start+k+1, t-k-j_col-1);
        for (int q = q_start; q <= q_end; q++) {
          d_list.push_back(T(p, q, dij(k+1, xd, I[p], J[q])));
        }
      }
      SparseMatrix<double> D(t-k-j_col, t-k-j_col);
      D.setFromTriplets(d_list.begin(), d_list.end());

      // Step 1b: form appropriate vector b
      VectorXd b (t-k-j_col);
      for (int p = 0; p < t-k-j_col; p++) {
        b[p] = -1.0 * dij(k+1, xd, I[p], j_col);
      }

      // Step 2: solve linear system Dx = b
      SparseQR<SparseMatrix<double>, Ord> qr;
      qr.compute(D);
      x = qr.solve(b);

      // Step 3: record vector x in the list
      for (int q = 0; q < t-k-j_col; q++) {
        if (J[q] < n_row) {
          n_list.push_back(T(J[q], j_col, x[q]));
          max_vals[j_col] = (x[q] > max_vals[j_col]) ? x[q] : max_vals[j_col];
        }
      }
    }
  }

  // Last r-k-1 basis vectors
  else {
    // First and second knots
    int t1 = knot_idx[j_col-k-1];
    int t2 = knot_idx[j_col-k];

    // Add entry where basis vector equals 1
    n_list.push_back(T(t2, j_col, 1));

    // Add more entries between first and second knot
    if (t2-1 >= t1+1) {
      // Compute row indices and columns indices for the appropriate divided
      // differences and unknown indices between first and second knot
      IntegerVector I = Range(t1-k+1, t2-k-1);
      IntegerVector J = Range(t1+1, t2-1);

      // Step 1a: form appropriate matrix D
      std::vector<T> d_list;
      int q_start = 0, q_end;
      d_list.reserve((t2-t1-1) * (k+2));
      for (int p = 0; p < t2-t1-1; p++) {
        q_start = std::max(p-k, 0);
        q_end = std::min(p+1, t2-t1-2);
        for (int q = q_start; q <= q_end; q++) {
          d_list.push_back(T(p, q, dij(k+1, xd, I[p], J[q])));
        }
      }
      SparseMatrix<double> D(t2-t1-1, t2-t1-1);
      D.setFromTriplets(d_list.begin(), d_list.end());

      // Step 1b: form appropriate vector b
      VectorXd b (t2-t1-1);
      for (int p = 0; p < t2-t1-1; p++) {
        b[p] = -1.0 * dij(k+1, xd, I[p], t2);
      }

      // Step 2: solve linear system Dx = b
      SparseQR<SparseMatrix<double>, Ord> qr;
      qr.compute(D);
      x = qr.solve(b);

      // Step 3: record vector x in the list
      for (int q = 0; q < t2-t1-1; q++) {
        if (J[q] < n_row) {
          n_list.push_back(T(J[q], j_col, x[q]));
          max_vals[j_col] = (x[q] > max_vals[j_col]) ? x[q] : max_vals[j_col];
        }
      }
    }

    // Add more entries between second and last knot
    if (t-k >= t2+1) {
      // Compute row indices i such that i+k is NOT one of the local knots
      IntegerVector I0 = Range(t2+1, t-1);
      IntegerVector I (t-k-t2);
      std::set_difference(I0.begin(), I0.end(), knot_idx.begin()+j_col-k+1, knot_idx.begin()+j_col, I.begin());
      for (int p = 0; p < t-k-t2; p++) I[p] -= k;

      // Compute column indices that correspond to the unknown indices between
      // second and last knot
      IntegerVector J = Range(t2+1, t-k);

      // Step 1a: form appropriate matrix D
      std::vector<T> d_list;
      int q_start = 0, q_end;
      d_list.reserve((t-k-t2) * (k+2));
      for (int p = 0; p < t-k-t2; p++) {
        // We know D is banded, find relevant start/end points
        while(J[q_start] < I[p]) q_start++;
        q_end = std::min(q_start+k+1, t-k-t2-1);
        for (int q = q_start; q <= q_end; q++) {
          d_list.push_back(T(p, q, dij(k+1, xd, I[p], J[q])));
        }
      }
      SparseMatrix<double> D(t-k-t2, t-k-t2);
      D.setFromTriplets(d_list.begin(), d_list.end());

      // Step 1b: form appropriate vector b
      VectorXd b (t-k-t2);
      for (int p = 0; p < t-k-t2; p++) {
        b[p] = -1.0 * dij(k+1, xd, I[p], t2);

        // Contributions between first and second knot
        if (t2-t1-1 > 0) {
          for (int q = 0; q < t2-t1-1; q++) {
            b[p] += -x[q] * dij(k+1, xd, I[p], t1+1+q);
          }
        }
      }

      // Step 2: solve linear system Dx = b
      SparseQR<SparseMatrix<double>, Ord> qr;
      qr.compute(D);
      x = qr.solve(b);

      // Step 3: record vector x in the list
      for (int q = 0; q < t-k-t2; q++) {
        if (J[q] < n_row) {
          n_list.push_back(T(J[q], j_col, x[q]));
          max_vals[j_col] = (x[q] > max_vals[j_col]) ? x[q] : max_vals[j_col];
        }
      }
    }
  }

  return;
}
