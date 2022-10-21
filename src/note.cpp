/******************************************************************************/
// Pass by reference note

// It looks like NumericVector objects (and other vector/matrix objects) are
// always be passed reference. To check this, uncomment the functions below,
// recompile, and then run the following in R. >>> a = rnorm(5); foo(a); a

// #include <Rcpp.h>
// using namespace Rcpp;

// void bar(NumericVector a) {
//   a[0] = 5;
// }

// // [[Rcpp::export]]
// void foo(NumericVector a) {
//   bar(a);
// }
