// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rcpp_hj_fun
double rcpp_hj_fun(int k, NumericVector xd, int j, double x);
RcppExport SEXP _dspline_rcpp_hj_fun(SEXP kSEXP, SEXP xdSEXP, SEXP jSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xd(xdSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_hj_fun(k, xd, j, x));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_dot_divided_diff
void rcpp_dot_divided_diff(NumericVector f, NumericVector z);
RcppExport SEXP _dspline_rcpp_dot_divided_diff(SEXP fSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    rcpp_dot_divided_diff(f, z);
    return R_NilValue;
END_RCPP
}
// rcpp_divided_diff
double rcpp_divided_diff(NumericVector f, NumericVector z);
RcppExport SEXP _dspline_rcpp_divided_diff(SEXP fSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_divided_diff(f, z));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_discrete_deriv
NumericVector rcpp_discrete_deriv(NumericVector f, int k, NumericVector xd, NumericVector x);
RcppExport SEXP _dspline_rcpp_discrete_deriv(SEXP fSEXP, SEXP kSEXP, SEXP xdSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xd(xdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_discrete_deriv(f, k, xd, x));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_discrete_integ
NumericVector rcpp_discrete_integ(NumericVector f, int k, NumericVector xd, NumericVector x);
RcppExport SEXP _dspline_rcpp_discrete_integ(SEXP fSEXP, SEXP kSEXP, SEXP xdSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xd(xdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_discrete_integ(f, k, xd, x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dspline_rcpp_hj_fun", (DL_FUNC) &_dspline_rcpp_hj_fun, 4},
    {"_dspline_rcpp_dot_divided_diff", (DL_FUNC) &_dspline_rcpp_dot_divided_diff, 2},
    {"_dspline_rcpp_divided_diff", (DL_FUNC) &_dspline_rcpp_divided_diff, 2},
    {"_dspline_rcpp_discrete_deriv", (DL_FUNC) &_dspline_rcpp_discrete_deriv, 4},
    {"_dspline_rcpp_discrete_integ", (DL_FUNC) &_dspline_rcpp_discrete_integ, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_dspline(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}