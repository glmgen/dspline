/******************************************************************************/
// Simple utilities

#include <Rcpp.h>
using namespace Rcpp;

// Factorial
int fact(int k) {
  int a = 1;
  for (int i = 1; i <= k; i++) a *= i;
  return a;
}

// Overwrite v with cumulative sums, starting at i
void CumSum(NumericVector v, int i) {
  for (int j = i+1; j < v.size(); j++) {
    v[j] += v[j-1];
  }
}

// Overwrite v with reverse cumulative sums, ending at i
void RevCumSum(NumericVector v, int i) {
  for (int j = v.size()-2; j >= i; j--) {
    v[j] += v[j+1];
  }
}

// Overwrite v with pairwise differences, starting at i
void Diff(NumericVector v, int i) {
  for (int j = v.size()-1; j >= i; j--) {
    v[j] -= v[j-1];
  }
}

// Overwrite v with reverse pairwise differences, ending at i
void RevDiff(NumericVector v, int i) {
  for (int j = i; j < v.size()-1; j++) {
    v[j] -= v[j+1];
  }
}

// Apply gap weighting to v (wrt xd), starting at i
void GapWeight(NumericVector v, int i, NumericVector xd) {
  for (int j = i; j < v.size(); j++) {
    v[j] *= (xd[j] - xd[j-i]) / i;
  }
}

// Apply invesre gap weighting to v (wrt xd), starting at i
void InvGapWeight(NumericVector v, int i, NumericVector xd) {
  for (int j = i; j < v.size(); j++) {
    v[j] /= (xd[j] - xd[j-i]) / i;
  }
}
