/******************************************************************************/
// Simple utilities

#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Sparse>
using namespace Rcpp;

// Factorial
int fact(int k) {
  int a = 1;
  for (int i = 1; i <= k; i++) a *= i;
  return a;
}

// Overwrite v with cumulative sums, starting at i
template<typename T>
void CumSum(T& v, int i) {
  for (int j = i+1; j < v.size(); j++) {
    v[j] += v[j-1];
  }
}
template void CumSum<NumericVector>(NumericVector& v, int i);
template void CumSum<Eigen::VectorXd>(Eigen::VectorXd& v, int i);

// Overwrite v with reverse cumulative sums, ending at i
template<typename T>
void RevCumSum(T& v, int i) {
  for (int j = v.size()-2; j >= i; j--) {
    v[j] += v[j+1];
  }
}
template void RevCumSum<NumericVector>(NumericVector& v, int i);
template void RevCumSum<Eigen::VectorXd>(Eigen::VectorXd& v, int i);

// Overwrite v with pairwise differences, starting at i
template<typename T>
void Diff(T& v, int i) {
  for (int j = v.size()-1; j >= i; j--) {
    v[j] -= v[j-1];
  }
}
template void Diff<NumericVector>(NumericVector& v, int i);
template void Diff<Eigen::VectorXd>(Eigen::VectorXd& v, int i);

// Overwrite v with reverse pairwise differences, ending at i
template<typename T>
void RevDiff(T& v, int i) {
  for (int j = i; j < v.size()-1; j++) {
    v[j] -= v[j+1];
  }
}
template void RevDiff<NumericVector>(NumericVector& v, int i);
template void RevDiff<Eigen::VectorXd>(Eigen::VectorXd& v, int i);

// Apply gap weighting to v (wrt xd), starting at i
template<typename T, typename U>
void GapWeight(T& v, int i, U& xd) {
  for (int j = i; j < v.size(); j++) {
    v[j] *= (xd[j] - xd[j-i]) / i;
  }
}
template void GapWeight<NumericVector,NumericVector>(NumericVector& v, int i, NumericVector& xd);
template void GapWeight<NumericVector,Eigen::VectorXd>(NumericVector& v, int i, Eigen::VectorXd& xd);
template void GapWeight<Eigen::VectorXd,NumericVector>(Eigen::VectorXd& v, int i, NumericVector& xd);
template void GapWeight<Eigen::VectorXd,Eigen::VectorXd>(Eigen::VectorXd& v, int i, Eigen::VectorXd& xd);

// Apply invesre gap weighting to v (wrt xd), starting at i
template<typename T, typename U>
void InvGapWeight(T& v, int i, U& xd) {
  for (int j = i; j < v.size(); j++) {
    v[j] /= (xd[j] - xd[j-i]) / i;
  }
}
template void InvGapWeight<NumericVector,NumericVector>(NumericVector& v, int i, NumericVector& xd);
template void InvGapWeight<NumericVector,Eigen::VectorXd>(NumericVector& v, int i, Eigen::VectorXd& xd);
template void InvGapWeight<Eigen::VectorXd,NumericVector>(Eigen::VectorXd& v, int i, NumericVector& xd);
template void InvGapWeight<Eigen::VectorXd,Eigen::VectorXd>(Eigen::VectorXd& v, int i, Eigen::VectorXd& xd);
